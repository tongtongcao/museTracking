#ifndef ROOT_MWPC
#define ROOT_MWPC

#include "TClonesArray.h"
#include "TVector3.h"

#include "utility.h"
#include "planeHit.h"
#include "plane.h"



class tracker: public plane{
 public:
  tracker(planeId id);
  ~tracker();

  Int_t getNHits() const{return fNHits;}
  void setNHits(Int_t nHits) { fNHits = nHits;}
  void init(UInt_t nHitsX, UInt_t nHitsY, Double_t posX[], Double_t posY[], Double_t posXErr[], Double_t posYErr[]);
  void init(UInt_t nHits, Double_t posX[], Double_t posY[], Double_t posXErr[], Double_t posYErr[]);

  TClonesArray* getHits() const{return fHits;}

 private:
  TClonesArray* fHits;
  Int_t fNHits;
  ClassDef(tracker,1)      

};

#endif
