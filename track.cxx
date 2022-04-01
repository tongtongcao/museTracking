#include <algorithm>
#include "track.h"

track::track(){
  fNHits = kDefault;
  fChi2 = kDefault;
  fNDF = kDefault;
  fChi2PerNDF = kDefault;

  for(Int_t i=0; i<kNTrackers; i++){
    fIsTrackerForTrack[i] = kFALSE;
    fTrackersPathMomentum[i] = kDefault;
    fTrackersPathTime[i] = kDefault;
    
    fTrackersMX[i] = kDefault;
    fTrackersMY[i] = kDefault;
    fTrackersMXErr[i] = kDefault;
    fTrackersMYErr[i] = kDefault;

    fTrackersSX[i] = kDefault;
    fTrackersSY[i] = kDefault;
    fTrackersSTx[i] = kDefault;
    fTrackersSTy[i] = kDefault;
  }
}

ClassImp(track)

