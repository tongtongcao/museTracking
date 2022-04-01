#include <iostream>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "utility.h"
#include "matrix.h"
#include "planeHit.h"
#include "plane.h"
#include "tracker.h"
#include "trackState.h"
#include "trackSite.h"
#include "trackSystem.h"
#include "track.h"
#include "trackFinder.h"

using namespace std;


int main(){

  TFile *file=new TFile("../GEMglobal247.root");
  TTree *tree=(TTree*)file->Get("MUSEglobalMesurements");

  Int_t nHitsGem0, nHitsGem1, nHitsGem2, nHitsGem3;
  Double_t posXGem0[250], posXGem1[250], posXGem2[250], posXGem3[250];
  Double_t posYGem0[250], posYGem1[250], posYGem2[250], posYGem3[250];

  tree->SetBranchAddress("nHitsGem0", &nHitsGem0);
  tree->SetBranchAddress("nHitsGem1", &nHitsGem1);
  tree->SetBranchAddress("nHitsGem2", &nHitsGem2);
  tree->SetBranchAddress("nHitsGem3", &nHitsGem3);

  tree->SetBranchAddress("posXGem0", posXGem0);
  tree->SetBranchAddress("posXGem1", posXGem1);
  tree->SetBranchAddress("posXGem2", posXGem2);
  tree->SetBranchAddress("posXGem3", posXGem3);

  tree->SetBranchAddress("posYGem0", posYGem0);
  tree->SetBranchAddress("posYGem1", posYGem1);
  tree->SetBranchAddress("posYGem2", posYGem2);
  tree->SetBranchAddress("posYGem3", posYGem3);

  TFile *fileOut=new TFile("trackInformation.root","recreate");
  TTree *treeOut=new TTree("tracks","");

  Int_t nTracks=MAXNTRACKS;
  treeOut->Branch("nTracks",&nTracks,"nTracks/I");

  Int_t nHits[nTracks];
  Double_t chi2[nTracks];
  Double_t nDF[nTracks];
  Double_t chi2PerNDF[nTracks];
  treeOut->Branch("nHits",nHits,"nHits[nTracks]/I");
  treeOut->Branch("chi2",chi2,"chi2[nTracks]/D");
  treeOut->Branch("nDF",nDF,"nDF[nTracks]/D");
  treeOut->Branch("chi2PerNDF",chi2PerNDF,"chi2PerNDF[nTracks]/D");

  Double_t trackersMomentumGEM0[nTracks];
  Double_t trackersTimeGEM0[nTracks];
  Double_t trackersMXGEM0[nTracks];
  Double_t trackersMYGEM0[nTracks];
  Double_t trackersSXGEM0[nTracks];
  Double_t trackersSYGEM0[nTracks];
  Double_t trackersSTxGEM0[nTracks];
  Double_t trackersSTyGEM0[nTracks];

  treeOut->Branch("trackersMomentumGEM0",trackersMomentumGEM0,"trackersMomentumGEM0[nTracks]/D");
  treeOut->Branch("trackersTimeGEM0",trackersTimeGEM0,"trackersTimeGEM0[nTracks]/D");
  treeOut->Branch("trackersMXGEM0",trackersMXGEM0,"trackersMXGEM0[nTracks]/D");
  treeOut->Branch("trackersMYGEM0",trackersMYGEM0,"trackersMYGEM0[nTracks]/D");
  treeOut->Branch("trackersSXGEM0",trackersSXGEM0,"trackersSXGEM0[nTracks]/D");
  treeOut->Branch("trackersSYGEM0",trackersSYGEM0,"trackersSYGEM0[nTracks]/D");
  treeOut->Branch("trackersSTxGEM0",trackersSTxGEM0,"trackersSTxGEM0[nTracks]/D");
  treeOut->Branch("trackersSTyGEM0",trackersSTyGEM0,"trackersSTyGEM0[nTracks]/D");

  Double_t trackersMomentumGEM1[nTracks];
  Double_t trackersTimeGEM1[nTracks];
  Double_t trackersMXGEM1[nTracks];
  Double_t trackersMYGEM1[nTracks];
  Double_t trackersSXGEM1[nTracks];
  Double_t trackersSYGEM1[nTracks];
  Double_t trackersSTxGEM1[nTracks];
  Double_t trackersSTyGEM1[nTracks];

  treeOut->Branch("trackersMomentumGEM1",trackersMomentumGEM1,"trackersMomentumGEM1[nTracks]/D");
  treeOut->Branch("trackersTimeGEM1",trackersTimeGEM1,"trackersTimeGEM1[nTracks]/D");
  treeOut->Branch("trackersMXGEM1",trackersMXGEM1,"trackersMXGEM1[nTracks]/D");
  treeOut->Branch("trackersMYGEM1",trackersMYGEM1,"trackersMYGEM1[nTracks]/D");
  treeOut->Branch("trackersSXGEM1",trackersSXGEM1,"trackersSXGEM1[nTracks]/D");
  treeOut->Branch("trackersSYGEM1",trackersSYGEM1,"trackersSYGEM1[nTracks]/D");
  treeOut->Branch("trackersSTxGEM1",trackersSTxGEM1,"trackersSTxGEM1[nTracks]/D");
  treeOut->Branch("trackersSTyGEM1",trackersSTyGEM1,"trackersSTyGEM1[nTracks]/D");

  Double_t trackersMomentumGEM2[nTracks];
  Double_t trackersTimeGEM2[nTracks];
  Double_t trackersMXGEM2[nTracks];
  Double_t trackersMYGEM2[nTracks];
  Double_t trackersSXGEM2[nTracks];
  Double_t trackersSYGEM2[nTracks];
  Double_t trackersSTxGEM2[nTracks];
  Double_t trackersSTyGEM2[nTracks];

  treeOut->Branch("trackersMomentumGEM2",trackersMomentumGEM2,"trackersMomentumGEM2[nTracks]/D");
  treeOut->Branch("trackersTimeGEM2",trackersTimeGEM2,"trackersTimeGEM2[nTracks]/D");
  treeOut->Branch("trackersMXGEM2",trackersMXGEM2,"trackersMXGEM2[nTracks]/D");
  treeOut->Branch("trackersMYGEM2",trackersMYGEM2,"trackersMYGEM2[nTracks]/D");
  treeOut->Branch("trackersSXGEM2",trackersSXGEM2,"trackersSXGEM2[nTracks]/D");
  treeOut->Branch("trackersSYGEM2",trackersSYGEM2,"trackersSYGEM2[nTracks]/D");
  treeOut->Branch("trackersSTxGEM2",trackersSTxGEM2,"trackersSTxGEM2[nTracks]/D");
  treeOut->Branch("trackersSTyGEM2",trackersSTyGEM2,"trackersSTyGEM2[nTracks]/D");

  Double_t trackersMomentumGEM3[nTracks];
  Double_t trackersTimeGEM3[nTracks];
  Double_t trackersMXGEM3[nTracks];
  Double_t trackersMYGEM3[nTracks];
  Double_t trackersSXGEM3[nTracks];
  Double_t trackersSYGEM3[nTracks];
  Double_t trackersSTxGEM3[nTracks];
  Double_t trackersSTyGEM3[nTracks];

  treeOut->Branch("trackersMomentumGEM3",trackersMomentumGEM3,"trackersMomentumGEM3[nTracks]/D");
  treeOut->Branch("trackersTimeGEM3",trackersTimeGEM3,"trackersTimeGEM3[nTracks]/D");
  treeOut->Branch("trackersMXGEM3",trackersMXGEM3,"trackersMXGEM3[nTracks]/D");
  treeOut->Branch("trackersMYGEM3",trackersMYGEM3,"trackersMYGEM3[nTracks]/D");
  treeOut->Branch("trackersSXGEM3",trackersSXGEM3,"trackersSXGEM3[nTracks]/D");
  treeOut->Branch("trackersSYGEM3",trackersSYGEM3,"trackersSYGEM3[nTracks]/D");
  treeOut->Branch("trackersSTxGEM3",trackersSTxGEM3,"trackersSTxGEM3[nTracks]/D");
  treeOut->Branch("trackersSTyGEM3",trackersSTyGEM3,"trackersSTyGEM3[nTracks]/D");

  Double_t trackersMomentum[nTracks][kNTrackers];
  Double_t trackersTime[nTracks][kNTrackers];
  Double_t trackersMX[nTracks][kNTrackers];
  Double_t trackersMY[nTracks][kNTrackers];
  Double_t trackersSX[nTracks][kNTrackers];
  Double_t trackersSY[nTracks][kNTrackers];
  Double_t trackersSTx[nTracks][kNTrackers];
  Double_t trackersSTy[nTracks][kNTrackers];

  particleId partId;

  tracker *uTrackers[kNTrackers];

  Int_t nHitsMeas[kNTrackers];
  Double_t posTrackerX[kNTrackers][250], posTrackerY[kNTrackers][250], posTrackerXErr[kNTrackers][250], posTrackerYErr[kNTrackers][250];

  Double_t resolutionGEM=0.1; // mm

  for(Int_t iEvent=0; iEvent<tree->GetEntries(); iEvent++){
    if(iEvent%100 == 0) cout<<iEvent<<" events are processed"<<endl;

    tree->GetEntry(iEvent);

    nHitsMeas[0]=nHitsGem0;
    for(Int_t iTrackGEM0=0; iTrackGEM0<nHitsGem0; iTrackGEM0++){
      posTrackerX[0][iTrackGEM0]=posXGem0[iTrackGEM0];
      posTrackerY[0][iTrackGEM0]=posYGem0[iTrackGEM0];
      posTrackerXErr[0][iTrackGEM0]=resolutionGEM;
      posTrackerYErr[0][iTrackGEM0]=resolutionGEM;
    }

    nHitsMeas[1]=nHitsGem1;
    for(Int_t iTrackGEM1=0; iTrackGEM1<nHitsGem1; iTrackGEM1++){
      posTrackerX[1][iTrackGEM1]=posXGem1[iTrackGEM1];
      posTrackerY[1][iTrackGEM1]=posYGem1[iTrackGEM1];
      posTrackerXErr[1][iTrackGEM1]=resolutionGEM;
      posTrackerYErr[1][iTrackGEM1]=resolutionGEM;
    }

    nHitsMeas[2]=nHitsGem2;
    for(Int_t iTrackGEM2=0; iTrackGEM2<nHitsGem2; iTrackGEM2++){
      posTrackerX[2][iTrackGEM2]=posXGem2[iTrackGEM2];
      posTrackerY[2][iTrackGEM2]=posYGem2[iTrackGEM2];
      posTrackerXErr[2][iTrackGEM2]=resolutionGEM;
      posTrackerYErr[2][iTrackGEM2]=resolutionGEM;
    }

    nHitsMeas[3]=nHitsGem3;
    for(Int_t iTrackGEM3=0; iTrackGEM3<nHitsGem3; iTrackGEM3++){
      posTrackerX[3][iTrackGEM3]=posXGem3[iTrackGEM3];
      posTrackerY[3][iTrackGEM3]=posYGem3[iTrackGEM3];
      posTrackerXErr[3][iTrackGEM3]=resolutionGEM;
      posTrackerYErr[3][iTrackGEM3]=resolutionGEM;
    }

    for(Int_t iTracker=0; iTracker<kNTrackers; iTracker++){
      uTrackers[iTracker]=new tracker((planeId)kTrackerId[iTracker]);
      uTrackers[iTracker]->init(nHitsMeas[iTracker], posTrackerX[iTracker], posTrackerY[iTracker], posTrackerXErr[iTracker], posTrackerYErr[iTracker]);
    }

    partId=kPositronId;
    trackFinder *trackFind=new trackFinder(partId, kParticleMass[partId], kParticleCharge[partId], kInitMomentum, uTrackers);
    TClonesArray* theTracks=new TClonesArray("track", MAXNTRACKS, kTRUE);
    theTracks->Clear();

    if(!trackFind->processHits(theTracks)){
      nTracks=0;

      treeOut->Fill();

      theTracks->SetOwner(kTRUE);
      delete theTracks;
      theTracks=NULL;
      delete trackFind;
      trackFind=NULL;

      for(Int_t iTracker=0; iTracker<kNTrackers; iTracker++){
	delete uTrackers[iTracker];
	uTrackers[iTracker]=NULL;
      }
      continue;
   
    }
    else nTracks=theTracks->GetEntries();

    for(int iTrack=0; iTrack<nTracks; iTrack++){
      nHits[iTrack]=((track*)theTracks->At(iTrack))->getNHits();
      chi2[iTrack]=((track*)theTracks->At(iTrack))->getChi2();
      nDF[iTrack]=((track*)theTracks->At(iTrack))->getNDF();
      chi2PerNDF[iTrack]=((track*)theTracks->At(iTrack))->getChi2PerNDF();
      ((track*)theTracks->At(iTrack))->getTrackersPathMomentum(trackersMomentum[iTrack]);
      ((track*)theTracks->At(iTrack))->getTrackersPathTime(trackersTime[iTrack]);
      ((track*)theTracks->At(iTrack))->getTrackersMX(trackersMX[iTrack]);
      ((track*)theTracks->At(iTrack))->getTrackersMY(trackersMY[iTrack]);
      ((track*)theTracks->At(iTrack))->getTrackersSX(trackersSX[iTrack]);
      ((track*)theTracks->At(iTrack))->getTrackersSY(trackersSY[iTrack]);
      ((track*)theTracks->At(iTrack))->getTrackersSTx(trackersSTx[iTrack]);
      ((track*)theTracks->At(iTrack))->getTrackersSTy(trackersSTy[iTrack]);

      trackersMomentumGEM0[iTrack]=trackersMomentum[iTrack][0];
      trackersTimeGEM0[iTrack]=trackersTime[iTrack][0];
      trackersMXGEM0[iTrack]=trackersMX[iTrack][0];
      trackersMYGEM0[iTrack]=trackersMY[iTrack][0];
      trackersSXGEM0[iTrack]=trackersSX[iTrack][0];
      trackersSYGEM0[iTrack]=trackersSY[iTrack][0];
      trackersSTxGEM0[iTrack]=trackersSTx[iTrack][0];
      trackersSTyGEM0[iTrack]=trackersSTy[iTrack][0];

      trackersMomentumGEM1[iTrack]=trackersMomentum[iTrack][1];
      trackersTimeGEM1[iTrack]=trackersTime[iTrack][1];
      trackersMXGEM1[iTrack]=trackersMX[iTrack][1];
      trackersMYGEM1[iTrack]=trackersMY[iTrack][1];
      trackersSXGEM1[iTrack]=trackersSX[iTrack][1];
      trackersSYGEM1[iTrack]=trackersSY[iTrack][1];
      trackersSTxGEM1[iTrack]=trackersSTx[iTrack][1];
      trackersSTyGEM1[iTrack]=trackersSTy[iTrack][1];

      trackersMomentumGEM2[iTrack]=trackersMomentum[iTrack][2];
      trackersTimeGEM2[iTrack]=trackersTime[iTrack][2];
      trackersMXGEM2[iTrack]=trackersMX[iTrack][2];
      trackersMYGEM2[iTrack]=trackersMY[iTrack][2];
      trackersSXGEM2[iTrack]=trackersSX[iTrack][2];
      trackersSYGEM2[iTrack]=trackersSY[iTrack][2];
      trackersSTxGEM2[iTrack]=trackersSTx[iTrack][2];
      trackersSTyGEM2[iTrack]=trackersSTy[iTrack][2];

      trackersMomentumGEM3[iTrack]=trackersMomentum[iTrack][3];
      trackersTimeGEM3[iTrack]=trackersTime[iTrack][3];
      trackersMXGEM3[iTrack]=trackersMX[iTrack][3];
      trackersMYGEM3[iTrack]=trackersMY[iTrack][3];
      trackersSXGEM3[iTrack]=trackersSX[iTrack][3];
      trackersSYGEM3[iTrack]=trackersSY[iTrack][3];
      trackersSTxGEM3[iTrack]=trackersSTx[iTrack][3];
      trackersSTyGEM3[iTrack]=trackersSTy[iTrack][3];
    }


    treeOut->Fill();
 
    theTracks->SetOwner(kTRUE);
    delete theTracks;
    theTracks=NULL;
    delete trackFind;
    trackFind=NULL;

    for(Int_t iTracker=0; iTracker<kNTrackers; iTracker++){
      delete uTrackers[iTracker];
      uTrackers[iTracker]=NULL;
    }

  }

  treeOut->Write();
  fileOut->Write();

  return 0;
}

