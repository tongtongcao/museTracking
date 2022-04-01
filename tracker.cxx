#include "tracker.h"
#include <vector>
using namespace std;

tracker::tracker(planeId id):plane(id, kPlanePosition[id]){
  fNHits=0;
  fHits = new TClonesArray("planeHit");
}

tracker::~tracker(){
  fHits->SetOwner(kTRUE);
  fHits->Clear("C");
  delete fHits;
}

void tracker::init(UInt_t nHitsX, UInt_t nHitsY, Double_t posX[], Double_t posY[], Double_t posXErr[], Double_t posYErr[]){
  fNHits=0;
  if(fHits->GetEntries()!=0) fHits->Clear();
  for(UInt_t i=0;i<nHitsX;i++){
    for(UInt_t j=0;j<nHitsY;j++){
      if(posXErr[i] == 0 || posYErr[j] ==0){
	cerr<<"Fatal error: uncertainty for measurement is equal to 0"<<endl;
	exit (EXIT_FAILURE);
      }
      new ( (*fHits)[fNHits]) planeHit(fPlaneId, fPlanePosition, i, j, posX[i], posY[j], posXErr[i], posYErr[j]);
      fNHits++;
    }
  }
}

void tracker::init(UInt_t nHits, Double_t posX[], Double_t posY[], Double_t posXErr[], Double_t posYErr[]){
  fNHits=0;
  if(fHits->GetEntries()!=0) fHits->Clear();
  for(UInt_t i=0;i<nHits;i++){
    if(posXErr[i] == 0 || posYErr[i] ==0){
      cerr<<"Fatal error: uncertainty for measurement is equal to 0"<<endl;
      exit (EXIT_FAILURE);
    }
    new ( (*fHits)[fNHits]) planeHit(fPlaneId, fPlanePosition, i, i, posX[i], posY[i], posXErr[i], posYErr[i]);
    fNHits++;
  }
}

ClassImp(tracker);


