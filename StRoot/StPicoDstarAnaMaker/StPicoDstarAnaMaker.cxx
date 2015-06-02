#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoDstarAnaMaker.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoHFMaker/StHFPair.h"
#include "StPhysicalHelixD.hh"

#include "SystemOfUnits.h"
ClassImp(StPicoDstarAnaMaker)

StPicoDstarAnaMaker::StPicoDstarAnaMaker(char const * name,char const * inputFilesList, 
    char const * outName,StPicoDstMaker* picoDstMaker): 
  StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mOutFileName(outName), mInputFileList(inputFilesList),
  mOutputFile(NULL), mChain(NULL), mEventCounter(0), mHFCuts(NULL)
{}

Int_t StPicoDstarAnaMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFileList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoDstarAnaMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoDstarAnaMaker - Could not open list of files. ABORT!" << endm;
      return kStErr;
   }

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
   mOutputFile->cd();

   if (!mHFCuts)
    mHFCuts = new StHFCuts;   
   mHFCuts->init();

   // -------------- USER VARIABLES -------------------------

   const char *mEventCutNames[7]   = {"all", "good run", "trigger", "#it{v}_{z}", "#it{v}_{z}-#it{v}^{VPD}_{z}", "accepted", ""};

   ntp_DMeson = new TNtuple("ntp","DMeson Tree","flag:kaonPID:pionPID:softPionPID:"
			    "dcaKaon:dcaPion:dcaSoftPion:dcaDaugters_D0:kaonPt:pionPt:softPionPt:"
			    "theta_D0:decayL_D0:pt_D0:mass_D0:eta_D0:phi_D0:"
			    "pt_Dstar:mass_Dstar:eta_Dstar:phi_Dstar");
   D0_sig = new TH1F("D0_sig","D0 candidate mass",200,0.5,2.5);
   eventStat = new TH1F("eventStat","event Cut stat",7,0,7);
   for(int ii=0; ii<7; ii++)
     eventStat->GetXaxis()->SetBinLabel(ii+1, mEventCutNames[ii]);
   // -------------------------------------------------------
   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoDstarAnaMaker::~StPicoDstarAnaMaker()
{
   /*  */
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarAnaMaker::Finish()
{
   LOG_INFO << " StPicoDstarAnaMaker - writing data and closing output file " <<endm;
   mOutputFile->cd();
   // save user variables here
   ntp_DMeson->Write();
   D0_sig->Write();
   eventStat->Write();
   mOutputFile->Close();
   return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarAnaMaker::Make()
{
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoDstarAnaMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();

   if (!picoDst)
   {
      LOG_WARN << "StPicoDstarAnaMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if(mPicoD0Event->runId() != picoDst->event()->runId() ||
       mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
     LOG_ERROR <<" StPicoDstarAnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
     LOG_ERROR <<" StPicoDstarAnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
     exit(1);
   }

   // -------------- USER ANALYSIS -------------------------
   
   // check if good event (including bad run)
   int aEventStat[7];
   if(!mHFCuts->isGoodEvent(const_cast<const StPicoDst*>(picoDst), aEventStat)) 
     return kStOk;
   fillEventStats(aEventStat);
   TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();
   StPicoEvent * mEvent = picoDst->event();
   StThreeVectorF primVtx = mEvent->primaryVertex();
   int mBField = mEvent->bField();
   int const nTracks = picoDst -> numberOfTracks();
   
   vector < int > softPion;
   for (unsigned short iTrack = 0; iTrack < nTracks; iTrack++){
     StPicoTrack const * track = picoDst->track(iTrack);

     if(!isSoftPion( track, primVtx) ) continue;
     softPion.push_back(iTrack);
   } 
   
   for (int idx = 0; idx < aKaonPion->GetEntries(); idx++)
   {
      // this is an example of how to get the kaonPion pairs and their corresponsing tracks
      StKaonPion const* Dzero = (StKaonPion*)aKaonPion->At(idx);
      
      if(!isGoodPair(Dzero)) continue;
      //check theta star
      if(Dzero->cosThetaStar()>0.8 ) continue;
      StPicoTrack const * kaon= picoDst->track(Dzero->kaonIdx());
      if( !mHFCuts -> isTPCKaon(kaon) ) continue;
  
      StPicoTrack const * pion1= picoDst->track(Dzero->pionIdx());
      if( !mHFCuts -> isTPCPion(pion1) ) continue;
  
      for (unsigned short iTrack = 0; iTrack < softPion.size(); iTrack++){
	StPicoTrack const * pion2 = picoDst->track(softPion.at(iTrack));

	StPhysicalHelixD helix = pion2->dcaGeometry().helix();
	helix.moveOrigin(helix.pathLength(primVtx));
	StThreeVectorF const pion2Mom(pion2->pMom().x(), pion2->pMom().y(), pion2->pMom().z());
	StLorentzVectorF const pionFourMom(pion2Mom,pion2Mom.massHypothesis(.13957));
	StLorentzVectorF const Dstar = pionFourMom + (Dzero->lorentzVector());
	float const pT = sqrt(pow(Dstar.px(), 2.0) + pow(Dstar.py(), 2.0) );
	if( (Dstar.m()-Dzero->m() ) >0.175 || (Dstar.m()-Dzero->m() ) < 0.115) continue;
	float ntVar[30];
	int flag=0;
	int ii = 0;
	if(pion2->charge()>0) // Dplus
	  flag = 1;
	ntVar[ii++] = flag;
	//PID
	ntVar[ii++] = fabs(kaon->nSigmaKaon());
	ntVar[ii++] = fabs(pion1->nSigmaPion());
	ntVar[ii++] = fabs(pion2->nSigmaPion());
	//Track info
 	ntVar[ii++] = Dzero->kaonDca();
	ntVar[ii++] = Dzero->pionDca();
	ntVar[ii++] = (helix.origin()-primVtx).mag();
	ntVar[ii++] = Dzero->dcaDaughters();
	ntVar[ii++] = kaon->gPt();
	ntVar[ii++] = pion1->gPt();
	ntVar[ii++] = pion2->pMom().perp();
	//D0 information
	ntVar[ii++] = Dzero->pointingAngle();
	ntVar[ii++] = Dzero->decayLength();
	ntVar[ii++] = Dzero->pt();
	ntVar[ii++] = Dzero->m();
	ntVar[ii++] = Dzero->eta();
	ntVar[ii++] = Dzero->phi();
	//Dstar info
	ntVar[ii++] = pT;
	ntVar[ii++] = Dstar.m();
	ntVar[ii++] = Dstar.pseudoRapidity();
	ntVar[ii++] = Dstar.phi();
	//Fill
	ntp_DMeson->Fill(ntVar);
      }

   }

   return kStOK;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodPair(StKaonPion const* const kp) const
{
  if(!kp) return false;
  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());
  
  //  To be replaced by mHFCuts->isGoodSecondaryVertexPair(kp))
  bool d0Sig = kp->m() > 0.5 && kp->m() < 2.1 &&
    std::cos(kp->pointingAngle()) > mHFCuts->cutSecondaryPairCosThetaMin() &&
    kp->decayLength()  > mHFCuts->cutSecondaryPairDecayLengthMin() && 
    kp->decayLength()  < mHFCuts->cutSecondaryPairDecayLengthMax() &&
    kp->dcaDaughters() < mHFCuts->cutSecondaryPairDcaDaughtersMax();

  if( mHFCuts->isGoodTrack(kaon) && mHFCuts->isGoodTrack(pion) &&
      mHFCuts->isTPCKaon(kaon) && mHFCuts->isTPCPion(pion) && 
      d0Sig )
    D0_sig->Fill(kp->m());

  bool pairCuts = kp->m() > mHFCuts->cutSecondaryPairMassMin() && 
    kp->m() < mHFCuts->cutSecondaryPairMassMax() &&
    std::cos(kp->pointingAngle()) > mHFCuts->cutSecondaryPairCosThetaMin() &&
    kp->decayLength()  > mHFCuts->cutSecondaryPairDecayLengthMin() && 
    kp->decayLength()  < mHFCuts->cutSecondaryPairDecayLengthMax() &&
    kp->dcaDaughters() < mHFCuts->cutSecondaryPairDcaDaughtersMax();

  return (mHFCuts->isGoodTrack(kaon) && mHFCuts->isGoodTrack(pion) &&
	  mHFCuts->isTPCKaon(kaon) && mHFCuts->isTPCPion(pion) && 
	  pairCuts);
}
// _________________________________________________________
bool StPicoDstarAnaMaker::isPion(StPicoTrack const * const trk) const {
  if(trk->pMom().mag()<0.6 || fabs(trk->nSigmaProton())>3.0) return false;
  // -- good pion
  return true;
}
bool StPicoDstarAnaMaker::isSoftPion(StPicoTrack const * const trk, StThreeVectorF &vtx) const {
  // -- Good tof pion if availabe
  //if( !mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk),StHFCuts::kPion) ) return false;
  // -- good pion
  if( trk->pMom().perp()<0.2 || fabs(trk->nSigmaPion())>3.0 ) return false;
  //From vertex
  StPhysicalHelixD helix = trk->dcaGeometry().helix();
  helix.moveOrigin(helix.pathLength(vtx));
  if( (helix.origin()-vtx).mag() > 0.008 ) return false;
  return true;
}
//________________________________________________________________________
void StPicoDstarAnaMaker::fillEventStats(int *aEventStat) {
  // -- Fill event statistics   
  for (unsigned int idx = 0; idx < mEventStatsBin; ++idx) {
    if (aEventStat[idx])
      break;
    eventStat->Fill(idx);
  }
}
