// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUChainTrackingCompression.cxx
/// \author David Rohr

#include "GPUChainTracking.h"
#include "GPULogging.h"
#include "GPUO2DataTypes.h"
#include "GPUTrackingInputProvider.h"

#ifdef GPUCA_HAVE_O2HEADERS
#include "GPUTPCCFChainContext.h"
#include "TPCClusterDecompressor.h"
#endif
#include "utils/strtag.h"

using namespace GPUCA_NAMESPACE::gpu;
using namespace o2::tpc;

int GPUChainTracking::RunTPCCompression()
{
#ifdef GPUCA_HAVE_O2HEADERS
  mRec->PushNonPersistentMemory(qStr2Tag("TPCCOMPR"));
  RecoStep myStep = RecoStep::TPCCompression;
  bool doGPU = GetRecoStepsGPU() & RecoStep::TPCCompression;
  GPUTPCCompression& Compressor = processors()->tpcCompressor;
  GPUTPCCompression& CompressorShadow = doGPU ? processorsShadow()->tpcCompressor : Compressor;
  const auto& threadContext = GetThreadContext();
  if (mPipelineFinalizationCtx && GetProcessingSettings().doublePipelineClusterizer) {
    RecordMarker(&mEvents->single, 0);
  }

  if (ProcessingSettings().tpcCompressionGatherMode == 3) {
    mRec->AllocateVolatileDeviceMemory(0); // make future device memory allocation volatile
  }
  SetupGPUProcessor(&Compressor, true);
  new (Compressor.mMemory) GPUTPCCompression::memory;
  WriteToConstantMemory(myStep, (char*)&processors()->tpcCompressor - (char*)processors(), &CompressorShadow, sizeof(CompressorShadow), 0);
  TransferMemoryResourcesToGPU(myStep, &Compressor, 0);
  runKernel<GPUMemClean16>(GetGridAutoStep(0, RecoStep::TPCCompression), krnlRunRangeNone, krnlEventNone, CompressorShadow.mClusterStatus, Compressor.mMaxClusters * sizeof(CompressorShadow.mClusterStatus[0]));
  runKernel<GPUTPCCompressionKernels, GPUTPCCompressionKernels::step0attached>(GetGridAuto(0), krnlRunRangeNone, krnlEventNone);
  runKernel<GPUTPCCompressionKernels, GPUTPCCompressionKernels::step1unattached>(GetGridAuto(0), krnlRunRangeNone, krnlEventNone);
  TransferMemoryResourcesToHost(myStep, &Compressor, 0);
#ifdef GPUCA_TPC_GEOMETRY_O2
  if (mPipelineFinalizationCtx && GetProcessingSettings().doublePipelineClusterizer) {
    SynchronizeEventAndRelease(&mEvents->single);
    ((GPUChainTracking*)GetNextChainInQueue())->RunTPCClusterizer_prepare(false);
    ((GPUChainTracking*)GetNextChainInQueue())->mCFContext->ptrClusterNativeSave = processorsShadow()->ioPtrs.clustersNative;
  }
#endif
  SynchronizeStream(0);
  o2::tpc::CompressedClusters* O = Compressor.mOutput;
  memset((void*)O, 0, sizeof(*O));
  O->nTracks = Compressor.mMemory->nStoredTracks;
  O->nAttachedClusters = Compressor.mMemory->nStoredAttachedClusters;
  O->nUnattachedClusters = Compressor.mMemory->nStoredUnattachedClusters;
  O->nAttachedClustersReduced = O->nAttachedClusters - O->nTracks;
  O->nSliceRows = NSLICES * GPUCA_ROW_COUNT;
  O->nComppressionModes = param().rec.tpc.compressionTypeMask;
  O->solenoidBz = param().bzkG;
  O->maxTimeBin = param().par.continuousMaxTimeBin;
  size_t outputSize = AllocateRegisteredMemory(Compressor.mMemoryResOutputHost, mSubOutputControls[GPUTrackingOutputs::getIndex(&GPUTrackingOutputs::compressedClusters)]);
  Compressor.mOutputFlat->set(outputSize, *Compressor.mOutput);
  char* hostFlatPtr = (char*)Compressor.mOutput->qTotU; // First array as allocated in GPUTPCCompression::SetPointersCompressedClusters
  size_t copySize = 0;
  if (ProcessingSettings().tpcCompressionGatherMode == 3) {
    CompressorShadow.mOutputA = Compressor.mOutput;
    copySize = AllocateRegisteredMemory(Compressor.mMemoryResOutputGPU); // We overwrite Compressor.mOutput with the allocated output pointers on the GPU
  }
  const o2::tpc::CompressedClustersPtrs* P = nullptr;
  HighResTimer* gatherTimer = nullptr;
  int outputStream = 0;
  if (ProcessingSettings().doublePipeline) {
    SynchronizeStream(OutputStream()); // Synchronize output copies running in parallel from memory that might be released, only the following async copy from stacked memory is safe after the chain finishes.
    outputStream = OutputStream();
  }
  if (ProcessingSettings().tpcCompressionGatherMode >= 2) {
    if (ProcessingSettings().tpcCompressionGatherMode == 2) {
      void* devicePtr = mRec->getGPUPointer(Compressor.mOutputFlat);
      if (devicePtr != Compressor.mOutputFlat) {
        CompressedClustersPtrs& ptrs = *Compressor.mOutput; // We need to update the ptrs with the gpu-mapped version of the host address space
        for (unsigned int i = 0; i < sizeof(ptrs) / sizeof(void*); i++) {
          reinterpret_cast<char**>(&ptrs)[i] = reinterpret_cast<char**>(&ptrs)[i] + (reinterpret_cast<char*>(devicePtr) - reinterpret_cast<char*>(Compressor.mOutputFlat));
        }
      }
    }
    TransferMemoryResourcesToGPU(myStep, &Compressor, outputStream);
    constexpr unsigned int nBlocksDefault = 2;
    constexpr unsigned int nBlocksMulti = 1 + 2 * 200;
    switch (ProcessingSettings().tpcCompressionGatherModeKernel) {
      case 0:
        runKernel<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::unbuffered>(GetGridBlkStep(nBlocksDefault, outputStream, RecoStep::TPCCompression), krnlRunRangeNone, krnlEventNone);
        getKernelTimer<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::unbuffered>(RecoStep::TPCCompression, 0, outputSize);
        break;
      case 1:
        runKernel<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::buffered32>(GetGridBlkStep(nBlocksDefault, outputStream, RecoStep::TPCCompression), krnlRunRangeNone, krnlEventNone);
        getKernelTimer<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::buffered32>(RecoStep::TPCCompression, 0, outputSize);
        break;
      case 2:
        runKernel<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::buffered64>(GetGridBlkStep(nBlocksDefault, outputStream, RecoStep::TPCCompression), krnlRunRangeNone, krnlEventNone);
        getKernelTimer<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::buffered64>(RecoStep::TPCCompression, 0, outputSize);
        break;
      case 3:
        runKernel<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::buffered128>(GetGridBlkStep(nBlocksDefault, outputStream, RecoStep::TPCCompression), krnlRunRangeNone, krnlEventNone);
        getKernelTimer<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::buffered128>(RecoStep::TPCCompression, 0, outputSize);
        break;
      case 4:

        static_assert((nBlocksMulti & 1) && nBlocksMulti >= 3);
        runKernel<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::multiBlock>(GetGridBlkStep(nBlocksMulti, outputStream, RecoStep::TPCCompression), krnlRunRangeNone, krnlEventNone);
        getKernelTimer<GPUTPCCompressionGatherKernels, GPUTPCCompressionGatherKernels::multiBlock>(RecoStep::TPCCompression, 0, outputSize);
        break;
      default:
        GPUError("Invalid compression kernel selected.");
        return 1;
    }
    if (ProcessingSettings().tpcCompressionGatherMode == 3) {
      RecordMarker(&mEvents->stream[outputStream], outputStream);
      char* deviceFlatPts = (char*)Compressor.mOutput->qTotU;
      if (GetProcessingSettings().doublePipeline) {
        const size_t blockSize = CAMath::nextMultipleOf<1024>(copySize / 30);
        const unsigned int n = (copySize + blockSize - 1) / blockSize;
        for (unsigned int i = 0; i < n; i++) {
          GPUMemCpy(myStep, hostFlatPtr + i * blockSize, deviceFlatPts + i * blockSize, CAMath::Min(blockSize, copySize - i * blockSize), outputStream, false);
        }
      } else {
        GPUMemCpy(myStep, hostFlatPtr, deviceFlatPts, copySize, outputStream, false);
      }
    }
  } else {
    char direction = 0;
    if (ProcessingSettings().tpcCompressionGatherMode == 0) {
      P = &CompressorShadow.mPtrs;
    } else if (ProcessingSettings().tpcCompressionGatherMode == 1) {
      P = &Compressor.mPtrs;
      direction = -1;
      gatherTimer = &getTimer<GPUTPCCompressionKernels>("GPUTPCCompression_GatherOnCPU", 0);
      gatherTimer->Start();
    }
    GPUMemCpyAlways(myStep, O->nSliceRowClusters, P->nSliceRowClusters, NSLICES * GPUCA_ROW_COUNT * sizeof(O->nSliceRowClusters[0]), outputStream, direction);
    GPUMemCpyAlways(myStep, O->nTrackClusters, P->nTrackClusters, O->nTracks * sizeof(O->nTrackClusters[0]), outputStream, direction);
    SynchronizeStream(outputStream);
    unsigned int offset = 0;
    for (unsigned int i = 0; i < NSLICES; i++) {
      for (unsigned int j = 0; j < GPUCA_ROW_COUNT; j++) {
        unsigned int srcOffset = mIOPtrs.clustersNative->clusterOffset[i][j] * Compressor.mMaxClusterFactorBase1024 / 1024;
        GPUMemCpyAlways(myStep, O->qTotU + offset, P->qTotU + srcOffset, O->nSliceRowClusters[i * GPUCA_ROW_COUNT + j] * sizeof(O->qTotU[0]), outputStream, direction);
        GPUMemCpyAlways(myStep, O->qMaxU + offset, P->qMaxU + srcOffset, O->nSliceRowClusters[i * GPUCA_ROW_COUNT + j] * sizeof(O->qMaxU[0]), outputStream, direction);
        GPUMemCpyAlways(myStep, O->flagsU + offset, P->flagsU + srcOffset, O->nSliceRowClusters[i * GPUCA_ROW_COUNT + j] * sizeof(O->flagsU[0]), outputStream, direction);
        GPUMemCpyAlways(myStep, O->padDiffU + offset, P->padDiffU + srcOffset, O->nSliceRowClusters[i * GPUCA_ROW_COUNT + j] * sizeof(O->padDiffU[0]), outputStream, direction);
        GPUMemCpyAlways(myStep, O->timeDiffU + offset, P->timeDiffU + srcOffset, O->nSliceRowClusters[i * GPUCA_ROW_COUNT + j] * sizeof(O->timeDiffU[0]), outputStream, direction);
        GPUMemCpyAlways(myStep, O->sigmaPadU + offset, P->sigmaPadU + srcOffset, O->nSliceRowClusters[i * GPUCA_ROW_COUNT + j] * sizeof(O->sigmaPadU[0]), outputStream, direction);
        GPUMemCpyAlways(myStep, O->sigmaTimeU + offset, P->sigmaTimeU + srcOffset, O->nSliceRowClusters[i * GPUCA_ROW_COUNT + j] * sizeof(O->sigmaTimeU[0]), outputStream, direction);
        offset += O->nSliceRowClusters[i * GPUCA_ROW_COUNT + j];
      }
    }
    offset = 0;
    for (unsigned int i = 0; i < O->nTracks; i++) {
      GPUMemCpyAlways(myStep, O->qTotA + offset, P->qTotA + Compressor.mAttachedClusterFirstIndex[i], O->nTrackClusters[i] * sizeof(O->qTotA[0]), outputStream, direction);
      GPUMemCpyAlways(myStep, O->qMaxA + offset, P->qMaxA + Compressor.mAttachedClusterFirstIndex[i], O->nTrackClusters[i] * sizeof(O->qMaxA[0]), outputStream, direction);
      GPUMemCpyAlways(myStep, O->flagsA + offset, P->flagsA + Compressor.mAttachedClusterFirstIndex[i], O->nTrackClusters[i] * sizeof(O->flagsA[0]), outputStream, direction);
      GPUMemCpyAlways(myStep, O->sigmaPadA + offset, P->sigmaPadA + Compressor.mAttachedClusterFirstIndex[i], O->nTrackClusters[i] * sizeof(O->sigmaPadA[0]), outputStream, direction);
      GPUMemCpyAlways(myStep, O->sigmaTimeA + offset, P->sigmaTimeA + Compressor.mAttachedClusterFirstIndex[i], O->nTrackClusters[i] * sizeof(O->sigmaTimeA[0]), outputStream, direction);

      // First index stored with track
      GPUMemCpyAlways(myStep, O->rowDiffA + offset - i, P->rowDiffA + Compressor.mAttachedClusterFirstIndex[i] + 1, (O->nTrackClusters[i] - 1) * sizeof(O->rowDiffA[0]), outputStream, direction);
      GPUMemCpyAlways(myStep, O->sliceLegDiffA + offset - i, P->sliceLegDiffA + Compressor.mAttachedClusterFirstIndex[i] + 1, (O->nTrackClusters[i] - 1) * sizeof(O->sliceLegDiffA[0]), outputStream, direction);
      GPUMemCpyAlways(myStep, O->padResA + offset - i, P->padResA + Compressor.mAttachedClusterFirstIndex[i] + 1, (O->nTrackClusters[i] - 1) * sizeof(O->padResA[0]), outputStream, direction);
      GPUMemCpyAlways(myStep, O->timeResA + offset - i, P->timeResA + Compressor.mAttachedClusterFirstIndex[i] + 1, (O->nTrackClusters[i] - 1) * sizeof(O->timeResA[0]), outputStream, direction);
      offset += O->nTrackClusters[i];
    }
    GPUMemCpyAlways(myStep, O->qPtA, P->qPtA, O->nTracks * sizeof(O->qPtA[0]), outputStream, direction);
    GPUMemCpyAlways(myStep, O->rowA, P->rowA, O->nTracks * sizeof(O->rowA[0]), outputStream, direction);
    GPUMemCpyAlways(myStep, O->sliceA, P->sliceA, O->nTracks * sizeof(O->sliceA[0]), outputStream, direction);
    GPUMemCpyAlways(myStep, O->timeA, P->timeA, O->nTracks * sizeof(O->timeA[0]), outputStream, direction);
    GPUMemCpyAlways(myStep, O->padA, P->padA, O->nTracks * sizeof(O->padA[0]), outputStream, direction);
  }
  if (ProcessingSettings().tpcCompressionGatherMode == 1) {
    gatherTimer->Stop();
  }
  mIOPtrs.tpcCompressedClusters = Compressor.mOutputFlat;
  if (ProcessingSettings().tpcCompressionGatherMode == 3) {
    SynchronizeEventAndRelease(&mEvents->stream[outputStream]);
    mRec->ReturnVolatileDeviceMemory();
  }

  if (mPipelineFinalizationCtx == nullptr) {
    SynchronizeStream(outputStream);
  } else {
    ((GPUChainTracking*)GetNextChainInQueue())->mRec->BlockStackedMemory(mRec);
  }
  mRec->PopNonPersistentMemory(RecoStep::TPCCompression, qStr2Tag("TPCCOMPR"));
#endif
  return 0;
}

int GPUChainTracking::RunTPCDecompression()
{
#ifdef GPUCA_HAVE_O2HEADERS
  if (GetProcessingSettings().tpcUseOldCPUDecoding) {
    const auto& threadContext = GetThreadContext();
    TPCClusterDecompressor decomp;
    auto allocator = [this](size_t size) {
      this->mInputsHost->mNClusterNative = this->mInputsShadow->mNClusterNative = size;
      this->AllocateRegisteredMemory(this->mInputsHost->mResourceClusterNativeOutput, this->mSubOutputControls[GPUTrackingOutputs::getIndex(&GPUTrackingOutputs::clustersNative)]);
      return this->mInputsHost->mPclusterNativeOutput;
    };
    auto& gatherTimer = getTimer<TPCClusterDecompressor>("TPCDecompression", 0);
    gatherTimer.Start();
    if (decomp.decompress(mIOPtrs.tpcCompressedClusters, *mClusterNativeAccess, allocator, param())) {
      GPUError("Error decompressing clusters");
      return 1;
    }
    gatherTimer.Stop();
    mIOPtrs.clustersNative = mClusterNativeAccess.get();
    if (mRec->IsGPU()) {
      AllocateRegisteredMemory(mInputsHost->mResourceClusterNativeBuffer);
      processorsShadow()->ioPtrs.clustersNative = mInputsShadow->mPclusterNativeAccess;
      WriteToConstantMemory(RecoStep::TPCDecompression, (char*)&processors()->ioPtrs - (char*)processors(), &processorsShadow()->ioPtrs, sizeof(processorsShadow()->ioPtrs), 0);
      *mInputsHost->mPclusterNativeAccess = *mIOPtrs.clustersNative;
      mInputsHost->mPclusterNativeAccess->clustersLinear = mInputsShadow->mPclusterNativeBuffer;
      mInputsHost->mPclusterNativeAccess->setOffsetPtrs();
      GPUMemCpy(RecoStep::TPCDecompression, mInputsShadow->mPclusterNativeBuffer, mIOPtrs.clustersNative->clustersLinear, sizeof(mIOPtrs.clustersNative->clustersLinear[0]) * mIOPtrs.clustersNative->nClustersTotal, 0, true);
      TransferMemoryResourceLinkToGPU(RecoStep::TPCDecompression, mInputsHost->mResourceClusterNativeAccess, 0);
      SynchronizeStream(0);
    }
  } else {
    mRec->PushNonPersistentMemory(qStr2Tag("TPCDCMPR"));
    RecoStep myStep = RecoStep::TPCDecompression;
    bool doGPU = GetRecoStepsGPU() & RecoStep::TPCDecompression;
    GPUTPCDecompression& Decompressor = processors()->tpcDecompressor;
    GPUTPCDecompression& DecompressorShadow = doGPU ? processorsShadow()->tpcDecompressor : Decompressor;
    const auto& threadContext = GetThreadContext();
    CompressedClusters cmprClsHost = *mIOPtrs.tpcCompressedClusters;
    CompressedClusters& inputGPU = Decompressor.mInputGPU;
    CompressedClusters& inputGPUShadow = DecompressorShadow.mInputGPU;

    inputGPU.nAttachedClusters = cmprClsHost.nAttachedClusters;
    inputGPU.nUnattachedClusters = cmprClsHost.nUnattachedClusters;
    inputGPU.nTracks = cmprClsHost.nTracks;
    inputGPU.nAttachedClustersReduced = inputGPU.nAttachedClusters - inputGPU.nTracks;
    inputGPU.nSliceRows = NSLICES * GPUCA_ROW_COUNT;
    inputGPU.nComppressionModes = param().rec.tpc.compressionTypeMask;
    inputGPU.solenoidBz = param().bzkG;
    inputGPU.maxTimeBin = param().par.continuousMaxTimeBin;
    SetupGPUProcessor(&Decompressor, true);
    WriteToConstantMemory(myStep, (char*)&processors()->tpcDecompressor - (char*)processors(), &DecompressorShadow, sizeof(DecompressorShadow), 0);
    int inputStream = 0;
    bool toGPU = true;
    GPUMemCpy(myStep, inputGPUShadow.nTrackClusters, cmprClsHost.nTrackClusters, cmprClsHost.nTracks * sizeof(cmprClsHost.nTrackClusters[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.qTotU, cmprClsHost.qTotU, cmprClsHost.nUnattachedClusters * sizeof(cmprClsHost.qTotU[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.qMaxU, cmprClsHost.qMaxU, cmprClsHost.nUnattachedClusters * sizeof(cmprClsHost.qMaxU[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.flagsU, cmprClsHost.flagsU, cmprClsHost.nUnattachedClusters * sizeof(cmprClsHost.flagsU[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.padDiffU, cmprClsHost.padDiffU, cmprClsHost.nUnattachedClusters * sizeof(cmprClsHost.padDiffU[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.timeDiffU, cmprClsHost.timeDiffU, cmprClsHost.nUnattachedClusters * sizeof(cmprClsHost.timeDiffU[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.sigmaPadU, cmprClsHost.sigmaPadU, cmprClsHost.nUnattachedClusters * sizeof(cmprClsHost.sigmaPadU[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.sigmaTimeU, cmprClsHost.sigmaTimeU, cmprClsHost.nUnattachedClusters * sizeof(cmprClsHost.sigmaTimeU[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.nSliceRowClusters, cmprClsHost.nSliceRowClusters, NSLICES * GPUCA_ROW_COUNT * sizeof(cmprClsHost.nSliceRowClusters[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.qTotA, cmprClsHost.qTotA, cmprClsHost.nAttachedClusters * sizeof(cmprClsHost.qTotA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.qMaxA, cmprClsHost.qMaxA, cmprClsHost.nAttachedClusters * sizeof(cmprClsHost.qMaxA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.flagsA, cmprClsHost.flagsA, cmprClsHost.nAttachedClusters * sizeof(cmprClsHost.flagsA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.rowDiffA, cmprClsHost.rowDiffA, (cmprClsHost.nAttachedClusters - cmprClsHost.nTracks) * sizeof(cmprClsHost.rowDiffA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.sliceLegDiffA, cmprClsHost.sliceLegDiffA, (cmprClsHost.nAttachedClusters - cmprClsHost.nTracks) * sizeof(cmprClsHost.sliceLegDiffA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.padResA, cmprClsHost.padResA, (cmprClsHost.nAttachedClusters - cmprClsHost.nTracks) * sizeof(cmprClsHost.padResA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.timeResA, cmprClsHost.timeResA, (cmprClsHost.nAttachedClusters - cmprClsHost.nTracks) * sizeof(cmprClsHost.timeResA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.sigmaPadA, cmprClsHost.sigmaPadA, cmprClsHost.nAttachedClusters * sizeof(cmprClsHost.sigmaPadA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.sigmaTimeA, cmprClsHost.sigmaTimeA, cmprClsHost.nAttachedClusters * sizeof(cmprClsHost.sigmaTimeA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.qPtA, cmprClsHost.qPtA, cmprClsHost.nTracks * sizeof(cmprClsHost.qPtA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.rowA, cmprClsHost.rowA, cmprClsHost.nTracks * sizeof(cmprClsHost.rowA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.sliceA, cmprClsHost.sliceA, cmprClsHost.nTracks * sizeof(cmprClsHost.sliceA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.timeA, cmprClsHost.timeA, cmprClsHost.nTracks * sizeof(cmprClsHost.timeA[0]), inputStream, toGPU);
    GPUMemCpy(myStep, inputGPUShadow.padA, cmprClsHost.padA, cmprClsHost.nTracks * sizeof(cmprClsHost.padA[0]), inputStream, toGPU);
    inputGPU.nTrackClusters = cmprClsHost.nTrackClusters;
    inputGPU.qTotU = cmprClsHost.qTotU;
    inputGPU.qMaxU = cmprClsHost.qMaxU;
    inputGPU.flagsU = cmprClsHost.flagsU;
    inputGPU.padDiffU = cmprClsHost.padDiffU;
    inputGPU.timeDiffU = cmprClsHost.timeDiffU;
    inputGPU.sigmaPadU = cmprClsHost.sigmaPadU;
    inputGPU.sigmaTimeU = cmprClsHost.sigmaTimeU;
    inputGPU.nSliceRowClusters = cmprClsHost.nSliceRowClusters;
    inputGPU.qTotA = cmprClsHost.qTotA;
    inputGPU.qMaxA = cmprClsHost.qMaxA;
    inputGPU.flagsA = cmprClsHost.flagsA;
    inputGPU.rowDiffA = cmprClsHost.rowDiffA;
    inputGPU.sliceLegDiffA = cmprClsHost.sliceLegDiffA;
    inputGPU.padResA = cmprClsHost.padResA;
    inputGPU.timeResA = cmprClsHost.timeResA;
    inputGPU.sigmaPadA = cmprClsHost.sigmaPadA;
    inputGPU.sigmaTimeA = cmprClsHost.sigmaTimeA;
    inputGPU.qPtA = cmprClsHost.qPtA;
    inputGPU.rowA = cmprClsHost.rowA;
    inputGPU.sliceA = cmprClsHost.sliceA;
    inputGPU.timeA = cmprClsHost.timeA;
    inputGPU.padA = cmprClsHost.padA;

    runKernel<GPUMemClean16>(GetGridAutoStep(inputStream, RecoStep::TPCDecompression), krnlRunRangeNone, krnlEventNone, DecompressorShadow.mNativeClustersIndex, NSLICES * GPUCA_ROW_COUNT * sizeof(DecompressorShadow.mNativeClustersIndex[0]));
    runKernel<GPUTPCDecompressionKernels, GPUTPCDecompressionKernels::step0attached>(GetGridAuto(inputStream), krnlRunRangeNone, krnlEventNone);

    mInputsHost->mNClusterNative = mInputsShadow->mNClusterNative = cmprClsHost.nAttachedClusters + cmprClsHost.nUnattachedClusters;
    AllocateRegisteredMemory(mInputsHost->mResourceClusterNativeOutput, mSubOutputControls[GPUTrackingOutputs::getIndex(&GPUTrackingOutputs::clustersNative)]);
    AllocateRegisteredMemory(mInputsHost->mResourceClusterNativeBuffer);
    DecompressorShadow.mNativeClustersBuffer = mInputsShadow->mPclusterNativeBuffer;
    Decompressor.mNativeClustersBuffer = mInputsHost->mPclusterNativeOutput;
    WriteToConstantMemory(myStep, (char*)&processors()->tpcDecompressor - (char*)processors(), &DecompressorShadow, sizeof(DecompressorShadow), inputStream);
    TransferMemoryResourceLinkToHost(RecoStep::TPCDecompression, Decompressor.mResourceTmpIndexes, inputStream);
    SynchronizeStream(inputStream);

    unsigned int offset = 0;
    unsigned int decodedAttachedClusters = 0;
    for (unsigned int i = 0; i < NSLICES; i++) {
      for (unsigned int j = 0; j < GPUCA_ROW_COUNT; j++) {
        unsigned int linearIndex = i * GPUCA_ROW_COUNT + j;
        unsigned int unattachedOffset = (linearIndex >= cmprClsHost.nSliceRows) ? 0 : cmprClsHost.nSliceRowClusters[linearIndex];
        (mClusterNativeAccess->nClusters)[i][j] = Decompressor.mNativeClustersIndex[linearIndex] + unattachedOffset;
        Decompressor.mUnattachedClustersOffsets[linearIndex] = offset;
        offset += unattachedOffset;
        decodedAttachedClusters += Decompressor.mNativeClustersIndex[linearIndex];
      }
    }
    TransferMemoryResourceLinkToGPU(RecoStep::TPCDecompression, Decompressor.mResourceTmpClustersOffsets, inputStream);
    if (decodedAttachedClusters != cmprClsHost.nAttachedClusters) {
      GPUWarning("%u / %u clusters failed track model decoding (%f %%)", cmprClsHost.nAttachedClusters - decodedAttachedClusters, cmprClsHost.nAttachedClusters, 100.f * (float)(cmprClsHost.nAttachedClusters - decodedAttachedClusters) / (float)cmprClsHost.nAttachedClusters);
    }

    if (doGPU) {
      mClusterNativeAccess->clustersLinear = mInputsShadow->mPclusterNativeBuffer;
      mClusterNativeAccess->setOffsetPtrs();
      *mInputsHost->mPclusterNativeAccess = *mClusterNativeAccess;
      processorsShadow()->ioPtrs.clustersNative = mInputsShadow->mPclusterNativeAccess;
      WriteToConstantMemory(RecoStep::TPCDecompression, (char*)&processors()->ioPtrs - (char*)processors(), &processorsShadow()->ioPtrs, sizeof(processorsShadow()->ioPtrs), inputStream);
      TransferMemoryResourceLinkToGPU(RecoStep::TPCDecompression, mInputsHost->mResourceClusterNativeAccess, inputStream);
    }
    mIOPtrs.clustersNative = mClusterNativeAccess.get();
    mClusterNativeAccess->clustersLinear = mInputsHost->mPclusterNativeOutput;
    mClusterNativeAccess->setOffsetPtrs();

    runKernel<GPUTPCDecompressionKernels, GPUTPCDecompressionKernels::step1unattached>(GetGridAuto(inputStream), krnlRunRangeNone, krnlEventNone);
    if (GetProcessingSettings().deterministicGPUReconstruction || GetProcessingSettings().debugLevel >= 4) {
      runKernel<GPUTPCDecompressionUtilKernels, GPUTPCDecompressionUtilKernels::sortPerSectorRow>(GetGridAutoStep(inputStream, RecoStep::TPCDecompression), krnlRunRangeNone, krnlEventNone);
    }
    SynchronizeStream(inputStream);
    GPUMemCpy(RecoStep::TPCDecompression, (void*)mInputsHost->mPclusterNativeOutput, (void*)mInputsShadow->mPclusterNativeBuffer, sizeof(mInputsShadow->mPclusterNativeBuffer[0]) * mIOPtrs.clustersNative->nClustersTotal, inputStream, false);
    mRec->PopNonPersistentMemory(RecoStep::TPCDecompression, qStr2Tag("TPCDCMPR"));
  }
#endif
  return 0;
}
