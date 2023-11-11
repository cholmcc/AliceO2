// Copyright 2023-2099 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   AODMcProducerHelpers.h
/// @author Christian Holm Christensen <cholm@nbi.dk>
/// common helpers for AOD MC producers
#include "AODProducerWorkflow/AODMcProducerHelpers.h"
#include <SimulationDataFormat/MCUtils.h>
#include <Framework/AnalysisDataModel.h>
#include <MathUtils/Utils.h>
#include <algorithm>
namespace o2::aodmchelpers
{
//==================================================================
bool hasKeys(o2::dataformats::MCEventHeader const& header,
             const std::vector<std::string>& keys,
             bool anyNotAll)
{
  auto check = [&header](const std::string& key) { return header.hasInfo(key); };
  return (anyNotAll ? std::any_of(keys.cbegin(), keys.cend(), check) : std::all_of(keys.cbegin(), keys.cend(), check));
}
//====================================================================
short updateMCCollisions(const CollisionCursor& cursor,
                         int bcId,
                         float time,
                         o2::dataformats::MCEventHeader const& header,
                         short generatorId,
                         int sourceId,
                         unsigned int mask)
{
  using Key = o2::dataformats::MCInfoKeys;
  using GenProp = o2::mcgenid::GeneratorProperty;
  using namespace o2::math_utils;

  int subGenId = getEventInfo(header, GenProp::SUBGENERATORID, -1);
  int genId = getEventInfo(header, GenProp::GENERATORID,
                           int(generatorId));
  float weight = getEventInfo(header, Key::weight, 1.f);

  auto encodedGeneratorId = o2::mcgenid::getEncodedGenId(genId,
                                                         sourceId,
                                                         subGenId);
  cursor(0,
         bcId,
         encodedGeneratorId,
         truncateFloatFraction(header.GetX(), mask),
         truncateFloatFraction(header.GetY(), mask),
         truncateFloatFraction(header.GetZ(), mask),
         truncateFloatFraction(time, mask),
         truncateFloatFraction(weight, mask),
         header.GetB());
  return encodedGeneratorId;
}
//--------------------------------------------------------------------
bool updateHepMCXSection(const XSectionCursor& cursor,
                         int collisionID,
                         short generatorID,
                         o2::dataformats::MCEventHeader const& header,
                         HepMCUpdate when)
{
  using Key = o2::dataformats::MCInfoKeys;

  if (when == HepMCUpdate::never or
      (when != HepMCUpdate::always and
       not hasKeys(header, {Key::acceptedEvents, Key::attemptedEvents, Key::xSection, Key::xSectionError}, when == HepMCUpdate::anyKey))) {
    return false;
  }
  
  cursor(0,
         collisionID,
         generatorID,
         getEventInfo(header, Key::acceptedEvents, 0),
         getEventInfo(header, Key::attemptedEvents, 0),
         getEventInfo(header, Key::xSection, 0.f),
         getEventInfo(header, Key::xSectionError, 0.f),
         getEventInfo(header, "ptHard", 1.f),
         getEventInfo(header, "MPI", -1),
         getEventInfo(header, "processId", -1));
  return true;
}
//--------------------------------------------------------------------
bool updateHepMCPdfInfo(const PdfInfoCursor& cursor,
                        int collisionID,
                        short generatorID,
                        o2::dataformats::MCEventHeader const& header,
                        HepMCUpdate when)
{
  using Key = o2::dataformats::MCInfoKeys;

  if (when == HepMCUpdate::never or
      (when != HepMCUpdate::always and not hasKeys(header, {Key::pdfParton1Id, Key::pdfParton2Id, Key::pdfCode1, Key::pdfCode2, Key::pdfX1, Key::pdfX2, Key::pdfScale, Key::pdfXF1, Key::pdfXF2}, when == HepMCUpdate::anyKey))) {
    return false;
  }
  
  cursor(0,
         collisionID,
         generatorID,
         getEventInfo(header, Key::pdfParton1Id, int(0)),
         getEventInfo(header, Key::pdfParton2Id, int(0)),
         getEventInfo(header, Key::pdfCode1, int(0)),
         getEventInfo(header, Key::pdfCode2, int(0)),
         getEventInfo(header, Key::pdfX1, 0.f),
         getEventInfo(header, Key::pdfX2, 0.f),
         getEventInfo(header, Key::pdfScale, 1.f),
         getEventInfo(header, Key::pdfXF1, 0.f),
         getEventInfo(header, Key::pdfXF2, 0.f));

  return true;
}
//--------------------------------------------------------------------
bool updateHepMCHeavyIon(const HeavyIonCursor& cursor,
                         int collisionID,
                         short generatorID,
                         o2::dataformats::MCEventHeader const& header,
                         HepMCUpdate when)
{
  using Key = dataformats::MCInfoKeys;

  if (when == HepMCUpdate::never or
      (when != HepMCUpdate::always and not hasKeys(header, {Key::nCollHard, Key::nPartProjectile, Key::nPartTarget, Key::nColl, Key::nCollNNWounded, Key::nCollNWoundedN, Key::nCollNWoundedNwounded, Key::nSpecProjectileNeutron, Key::nSpecTargetNeutron, Key::nSpecProjectileProton, Key::nSpecTargetProton, Key::planeAngle, "eccentricity", Key::sigmaInelNN, Key::centrality}, when == HepMCUpdate::anyKey))) {
    return false;
  }
  
  int specNeutrons = (getEventInfo(header, Key::nSpecProjectileNeutron, -1) +
                      getEventInfo(header, Key::nSpecTargetNeutron, -1));
  int specProtons = (getEventInfo(header, Key::nSpecProjectileProton, -1) +
                     getEventInfo(header, Key::nSpecTargetProton, -1));
  cursor(0,
         collisionID,
         generatorID,
         getEventInfo(header, Key::nCollHard, -1),
         getEventInfo(header, Key::nPartProjectile, -1),
         getEventInfo(header, Key::nPartTarget, -1),
         getEventInfo(header, Key::nColl, -1),
         getEventInfo(header, Key::nCollNNWounded, -1),
         getEventInfo(header, Key::nCollNWoundedN, -1),
         getEventInfo(header, Key::nCollNWoundedNwounded, -1),
         specNeutrons,
         specProtons,
         header.GetB(),
         getEventInfo(header, Key::planeAngle, header.GetRotZ()),
         getEventInfo(header, "eccentricity", 0),
         getEventInfo(header, Key::sigmaInelNN, 0.),
         getEventInfo(header, Key::centrality, -1));

  return true;
}
//--------------------------------------------------------------------
void updateParticle(const ParticleCursor& cursor,
                    const TrackToIndex& toStore,
                    int collisionID,
                    o2::MCTrack const& track,
                    std::vector<MCTrack> const& tracks,
                    uint8_t flags,
                    uint32_t weightMask,
                    uint32_t momentumMask,
                    uint32_t positionMask)
{
  using o2::mcutils::MCTrackNavigator;
  using namespace o2::aod::mcparticle::enums;
  using namespace o2::math_utils;
  using namespace o2::mcgenstatus;

  auto mapping = [&toStore](int trackNo) {
    auto iter = toStore.find(trackNo);
    if (iter == toStore.end()) {
      return -1;
    }
    return iter->second;
  };

  auto statusCode = track.getStatusCode().fullEncoding;
  auto hepmc = getHepMCStatusCode(track.getStatusCode());
  if (not track.isPrimary()) {
    flags = ProducedByTransport;
    statusCode = track.getProcess();
  }
  if (MCTrackNavigator::isPhysicalPrimary(track, tracks)) {
    flags = PhysicalPrimary;
  }
  
  int daughters[2] = {-1, -1};
  std::vector<int> mothers;
  int id;
  if ((id = mapping(track.getMotherTrackId())) >= 0) {
    mothers.push_back(id);
  }
  if ((id = mapping(track.getSecondMotherTrackId())) >= 0) {
    mothers.push_back(id);
  }
  if ((id = mapping(track.getFirstDaughterTrackId())) >= 0) {
    daughters[0] = id;
  }
  if ((id = mapping(track.getFirstDaughterTrackId())) >= 0) {
    daughters[1] = id;
  }
  else {
    daughters[1] = daughters[0];
  }
  if (daughters[0] < 0 and daughters[1] >= 0) {
    LOG(error) << "Problematic daughters: " << daughters[0] << " and "
               << daughters[1];
    daughters[0] = daughters[1];
  }
  if (daughters[0] > daughters[1]) {
    std::swap(daughters[0], daughters[1]);
  }
  
  float weight = track.getWeight();
  float pX = float(track.Px());
  float pY = float(track.Py());
  float pZ = float(track.Pz());
  float energy = float(track.GetEnergy());
  float vX = float(track.Vx());
  float vY = float(track.Vy());
  float vZ = float(track.Vz());
  float time = float(track.T());

  cursor(0,
         collisionID,
         track.GetPdgCode(),
         statusCode,
         flags,
         mothers,
         daughters,
         truncateFloatFraction(weight, weightMask),
         truncateFloatFraction(pX, momentumMask),
         truncateFloatFraction(pY, momentumMask),
         truncateFloatFraction(pZ, momentumMask),
         truncateFloatFraction(energy, momentumMask),
         truncateFloatFraction(vX, positionMask),
         truncateFloatFraction(vY, positionMask),
         truncateFloatFraction(vZ, positionMask),
         truncateFloatFraction(time, positionMask));
}
//--------------------------------------------------------------------
uint32_t updateParticles(const ParticleCursor& cursor,
                         int collisionID,
                         std::vector<MCTrack> const& tracks,
                         uint32_t offset,
                         bool filter,
                         bool background,
                         const TrackToIndex& preselect,
                         uint32_t weightMask,
                         uint32_t momentumMask,
                         uint32_t positionMask)
{
  using o2::mcutils::MCTrackNavigator;
  using namespace o2::aod::mcparticle::enums;
  using namespace o2::mcgenstatus;

  // First loop over particles to find out which to store
  TrackToIndex toStore(preselect.begin(), preselect.end());

  // Mapping lambda.  This maps the track number to the index into
  // the table exported.
  auto mapping = [&toStore](int trackNo) {
    auto iter = toStore.find(trackNo);
    if (iter == toStore.end()) {
      return -1;
    }
    return iter->second;
  };

  LOG(debug) << "Got a total of " << tracks.size();
  for (int trackNo = tracks.size() - 1; trackNo >= 0; trackNo--) {
    auto& track = tracks[trackNo];
    auto hepmc = getHepMCStatusCode(track.getStatusCode());
    if (filter) {
      if (mapping(trackNo) < 0 and
          not track.isPrimary() and
          not MCTrackNavigator::isPhysicalPrimary(track, tracks) and
          not MCTrackNavigator::isKeepPhysics(track, tracks)) {
        LOG(debug) << "Skipping track " << trackNo << " " << hepmc;
        continue;
      }
    }

    // Store this particle.  We mark that putting a 1 in the
    // `toStore` mapping. This will later on be updated with the
    // actual index into the table
    toStore[trackNo] = 1;

    // If we're filtering tracks, then also mark mothers and
    // daughters(?) to be stored.
    if (filter) {
      int id;
      if ((id = track.getMotherTrackId()) >= 0) {
        toStore[id] = 1;
      }
      if ((id = track.getSecondMotherTrackId()) >= 0) {
        toStore[id] = 1;
      }
      if ((id = track.getFirstDaughterTrackId()) >= 0) {
        toStore[id] = 1;
      }
      if ((id = track.getLastDaughterTrackId()) >= 0) {
        toStore[id] = 1;
      }
    }
  }

  // Second loop to set indexes.  This is needed to be done before
  // we actually update the table, because a particle may point to a
  // later particle.
  LOG(debug) << "Will store " << toStore.size() << " particles";
  size_t index = 0;
  for (size_t trackNo = 0U; trackNo < tracks.size(); trackNo++) {
    auto storeIt = mapping(trackNo);
    if (storeIt < 0) {
      continue;
    }

    toStore[trackNo] = offset + index;
    index++;
  }

  // Make sure we have the right number
  assert(index == toStore.size());

  LOG(debug) << "Starting index " << offset << ", last index" << index;
  // Third loop to actually store the particles in the order given
  for (size_t trackNo = 0U; trackNo < tracks.size(); trackNo++) {
    auto storeIt = mapping(trackNo);
    if (storeIt < 0) {
      continue;
    }
    
    auto& track = tracks[trackNo];
    auto hepmc = getHepMCStatusCode(track.getStatusCode());
    uint8_t flags = (background ? FromBackgroundEvent : 0);
    updateParticle(cursor,
                   toStore,
                   collisionID,
                   track,
                   tracks,
                   flags,
                   weightMask,
                   momentumMask,
                   positionMask);
  }
  LOG(debug) << "Return new offset " << offset + index;
  return offset + index;
}
} // namespace o2::aodmchelpers

// Local Variables:
//   mode: C++
// End:
