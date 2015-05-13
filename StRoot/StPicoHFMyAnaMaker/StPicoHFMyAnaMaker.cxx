#include <vector>

#include "TClonesArray.h"

#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"

#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"

#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoHFMaker/StHFPair.h"
#include "StPicoHFMaker/StHFTriplet.h"
#include "StPicoHFMaker/StHFQuadruplet.h"

#include "StPicoHFMyAnaMaker.h"
#include "SystemOfUnits.h"
ClassImp(StPicoHFMyAnaMaker)

// _________________________________________________________
StPicoHFMyAnaMaker::StPicoHFMyAnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName,  
					   char const* inputHFListHFtree = "") :
  StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
  mDecayChannel(kChannel1) {
  // constructor
}

// _________________________________________________________
StPicoHFMyAnaMaker::~StPicoHFMyAnaMaker() {
  // destructor
}

// _________________________________________________________
int StPicoHFMyAnaMaker::InitHF() {
    ntp_DMeson = new TNtuple("ntp","DMeson Tree","flag:"
			     "dca1:dca2:dca3:dca4:dcaMax:"
			     "pt1:pt2:pt3:theta_hs:decayL_hs:"
			     "pt_hs:mass_hs:eta_hs:phi_hs:kaonTof");
  
  return kStOK;
}

// _________________________________________________________
void StPicoHFMyAnaMaker::ClearHF(Option_t *opt="") {
  return;
}

// _________________________________________________________
int StPicoHFMyAnaMaker::FinishHF() {
  ntp_DMeson->Write();
  return kStOK;
}

// _________________________________________________________
int StPicoHFMyAnaMaker::MakeHF() {
  // -- process event
  //    ADD YOUR PROCESSING CODE HERE
  //    ... it is usefull to use the methods below
  //     - createCandidates()
  //     - analyzeCandidates()

  if (isMakerMode() == StPicoHFMaker::kWrite) {
    createCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kRead) {
    // -- the reading back of the perviously written trees happens in the background
    analyzeCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kAnalyze) {
    createCandidates();
    analyzeCandidates();
  }

  return kStOK;
}

// _________________________________________________________
int StPicoHFMyAnaMaker::createCandidates() {
  // create candidate pairs/ triplet and fill them in arrays (in StPicoHFEvent)
  
  // -- ADD USER CODE TO CREATE PARTICLE CANDIDATES --------
  //    - vectors mIdxPicoKaons, mIdxPicoPions mIdxPicoProtons
  //      have been filled in the background using the cuts in HFCuts

  // -- Decay channel D0 4 decay --- EXAMPLE
  if (mDecayChannel == StPicoHFMyAnaMaker::kChannel1) {

    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
      StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
      if( !isPion(pion1) ) continue;
      if (! isVertexGoodTrack(pion1, mPrimVtx, mBField) ) continue;
	  
      for (unsigned short idxPion2 = 0; idxPion2 < mIdxPicoPions.size(); ++idxPion2) {
	StPicoTrack const *pion2 = mPicoDst->track(mIdxPicoPions[idxPion2]);
	if( !isPion(pion2) ) continue;
	if (! isVertexGoodTrack(pion2, mPrimVtx, mBField) ) continue;

	for (unsigned short idxPion3 = 0; idxPion3 < mIdxPicoPions.size(); ++idxPion3) {
	  StPicoTrack const *pion3 = mPicoDst->track(mIdxPicoPions[idxPion3]);
	  int const pi_charge = pion1->charge()*pion2->charge()*pion3->charge();
	  if( TMath::Abs(pi_charge)!=1) continue;
	  if( !isPion(pion3) ) continue;
	  if (! isVertexGoodTrack(pion3, mPrimVtx, mBField) ) continue;

	  for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
	    StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
	    if (kaon->charge() == pi_charge) continue;
	    if( !isKaon(kaon) ) continue;
	    if( !mHFCuts->isHybridTOFHadron(kaon, mHFCuts->getTofBetaBase(kaon), StHFCuts::kKaon) ) continue;
	    if (! isVertexGoodTrack(kaon, mPrimVtx, mBField) ) continue;
	    
	    StHFQuadruplet quad(pion1, pion2, pion3, kaon,
				mHFCuts->getHypotheticalMass(StHFCuts::kPion), mHFCuts->getHypotheticalMass(StHFCuts::kPion),
				mHFCuts->getHypotheticalMass(StHFCuts::kPion),mHFCuts->getHypotheticalMass(StHFCuts::kKaon),
				mIdxPicoPions[idxPion1],mIdxPicoPions[idxPion2],mIdxPicoPions[idxPion3],mIdxPicoKaons[idxKaon],
				mPrimVtx,mBField);
	    if (!mHFCuts->isGoodSecondaryVertexQuadruplet(&quad)) 
	      continue;
	    mPicoHFEvent->addHFSecondaryVertexQuadruplet(&quad);
	  }
	}
      }
    }
  }
  return kStOK;
}

// _________________________________________________________
int StPicoHFMyAnaMaker::analyzeCandidates() {  
  // -- Decay channel1
  if (mDecayChannel == StPicoHFMyAnaMaker::kChannel1) {
    
    // -- fill nTuple / hists for secondary pairs
    TClonesArray const * aCandidates= mPicoHFEvent->aHFSecondaryVertices();
    if ( mPicoHFEvent->nHFSecondaryVertices() <= 0 ) return kStOk;
    for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {

      StHFQuadruplet const* quadruplet = static_cast<StHFQuadruplet*>(aCandidates->At(idx));
      StPicoTrack const* pion1 = mPicoDst->track(quadruplet->particle1Idx());
      StPicoTrack const* pion2 = mPicoDst->track(quadruplet->particle2Idx());
      StPicoTrack const* pion3 = mPicoDst->track(quadruplet->particle3Idx());
      StPicoTrack const* kaon = mPicoDst->track(quadruplet->particle4Idx());
      int ii=0;
      float ntVar[30];
      float const dca1 = quadruplet->particle1Dca();
      float const dca2 = quadruplet->particle2Dca();
      float const dca3 = quadruplet->particle3Dca();
      float const dca4 = quadruplet->particle4Dca();
      
      float const dcaD12 = quadruplet->dcaDaughters12();
      float const dcaD13 = quadruplet->dcaDaughters13();
      float const dcaD14 = quadruplet->dcaDaughters14();
      float const dcaD23 = quadruplet->dcaDaughters23();
      float const dcaD24 = quadruplet->dcaDaughters24();
      float const dcaD34 = quadruplet->dcaDaughters34();

      float dcaMax = dcaD12 > dcaD13 ? dcaD12 : dcaD13;
      //Jsut saving max
      dcaMax = dcaMax > dcaD14 ? dcaMax : dcaD14;
      dcaMax = dcaMax > dcaD23 ? dcaMax : dcaD23;
      dcaMax = dcaMax > dcaD24 ? dcaMax : dcaD24;
      dcaMax = dcaMax > dcaD34 ? dcaMax : dcaD34;
      
      float kaonTOF = 0;
      if( mHFCuts->getTofBetaBase(kaon) > 0){
	float const tofBeta = mHFCuts->getTofBetaBase(kaon);
	float const ptot    = kaon->dcaGeometry().momentum().mag();
	float const betaInv = sqrt(ptot*ptot + pow(mHFCuts->getHypotheticalMass(StHFCuts::kKaon),2.0)) / ptot;
	kaonTOF = 1.0/tofBeta - betaInv;
      }

      float const pt1=sqrt(pow(pion1->gMom().x(),2.0)+pow(pion1->gMom().y(),2.0));
      float const pt2=sqrt(pow(pion2->gMom().x(),2.0)+pow(pion2->gMom().y(),2.0));
      float const pt3=sqrt(pow(pion3->gMom().x(),2.0)+pow(pion3->gMom().y(),2.0));
      float const pt4=sqrt(pow(kaon->gMom().x(),2.0)+pow(kaon->gMom().y(),2.0));
      float const pt=sqrt(pow(quadruplet->px(),2.0)+pow(quadruplet->py(),2.0));
      float flag = -999;
      if( kaon -> charge() < 0 ) flag = 0; // -- D0 ?
      if( kaon -> charge() > 0 ) flag = 1; // -- D0bar?
      
      // ---
      // Saving to NTUPLE
      ntVar[ii++] = flag;
      ntVar[ii++] = quadruplet->particle1Dca();
      ntVar[ii++] = quadruplet->particle2Dca();
      ntVar[ii++] = quadruplet->particle3Dca();
      ntVar[ii++] = quadruplet->particle4Dca();
      ntVar[ii++] = dcaMax;
      ntVar[ii++] = pt1;
      ntVar[ii++] = pt2;
      ntVar[ii++] = pt3;
      ntVar[ii++] = pt4;
      ntVar[ii++] = quadruplet->pointingAngle();
      ntVar[ii++] = quadruplet->decayLength();
      ntVar[ii++] = pt;
      ntVar[ii++] = quadruplet->m();
      ntVar[ii++] = quadruplet->eta();
      ntVar[ii++] = quadruplet->phi();
      ntVar[ii++] = kaonTOF;
      ntp_DMeson->Fill(ntVar);
    
    } // for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {
    // for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {
  } // else  if (mDecayChannel == StPicoHFMyAnaMaker::kChannel1) {

 return kStOK;
}

// _________________________________________________________
bool StPicoHFMyAnaMaker::isPion(StPicoTrack const * const trk) const {
  // -- good pion
  return true;
}

// _________________________________________________________
bool StPicoHFMyAnaMaker::isKaon(StPicoTrack const * const trk) const {
  // -- good kaon
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCKaon(trk));
} 

// _________________________________________________________
bool StPicoHFMyAnaMaker::isProton(StPicoTrack const * const trk) const {
  // -- good proton
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCProton(trk) && mHFCuts->isTOFProton(trk));
}
// _________________________________________________________
bool StPicoHFMyAnaMaker::isCloseTracks(StPicoTrack const * const trk1, StPicoTrack const * const trk2, StThreeVectorF const & vtx, float bField) const {
  
  StPhysicalHelixD p1Helix = trk1->dcaGeometry().helix();
  StPhysicalHelixD p2Helix = trk2->dcaGeometry().helix();
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));
  if( ( p1Helix.origin()-vtx ).mag()>0.2 || ( p2Helix.origin()-vtx ).mag()>0.2 ) return false;
  //Requires loading constants
  StThreeVectorF const p1Mom = p1Helix.momentum(bField * kilogauss);
  StThreeVectorF const p2Mom = p2Helix.momentum(bField * kilogauss);
  StPhysicalHelixD const p1StraightLine(p1Mom, p1Helix.origin(), 0, trk1->charge());
  StPhysicalHelixD const p2StraightLine(p2Mom, p2Helix.origin(), 0, trk2->charge());
  //DCA
  pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
  StThreeVectorF const p1AtDcaToP2 = p1StraightLine.at(ss.first);
  StThreeVectorF const p2AtDcaToP1 = p2StraightLine.at(ss.second);
  int const dca = (p1AtDcaToP2-p1AtDcaToP2).mag();
  if(dca > 0.01) return false;
// -- good pair
  return true;
}
// _________________________________________________________
bool StPicoHFMyAnaMaker::isVertexGoodTrack(StPicoTrack const * const trk1, StThreeVectorF const & vtx, float bField) const {
  StPhysicalHelixD p1Helix = trk1->dcaGeometry().helix();
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  if( ( p1Helix.origin()-vtx ).mag()>0.02) return false;
  return true;
}
