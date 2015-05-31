#ifndef StPicoDstarAnaMaker_h
#define StPicoDstarAnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "TNtuple.h"
#include "StThreeVectorF.hh"
class TString;
class TFile;
class TNtuple;
class StPicoTrack;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StHFCuts;
class StHFPair;

class StPicoDstarAnaMaker : public StMaker
{
  public:
    StPicoDstarAnaMaker(char const * name, char const * inputFilesList, 
        char const * outName,StPicoDstMaker* picoDstMaker);
    virtual ~StPicoDstarAnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setHFCuts(StHFCuts* cuts);    

  private:
    StPicoDstarAnaMaker() {}
    void readNextEvent();

    bool isGoodPair(StKaonPion const*) const;
    bool isPion(StPicoTrack const * const) const;
    bool isSoftPion(StPicoTrack const * const, StThreeVectorF &) const;
    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;

    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;

    StHFCuts* mHFCuts;

    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    TNtuple *ntp_DMeson;
    ClassDef(StPicoDstarAnaMaker, 1)
};

inline int StPicoDstarAnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoDstarAnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

inline void StPicoDstarAnaMaker::setHFCuts(StHFCuts* cuts)   
{ 
  mHFCuts = cuts; 
}

#endif
