void runPicoDStarAnaMaker(TString d0list, TString outFileName, TString badRunListFileName = "picoList_bad_MB.list")
{
   //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
   string SL_version = "SL15c";
   string env_SL = getenv("STAR");
   if (env_SL.find(SL_version) == string::npos)
   {
      cout << "Environment Star Library does not match the requested library in runPicoD0EventMaker.C. Exiting..." << endl;
      exit(1);
   }

   gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
   loadSharedLibraries();

   gSystem->Load("StPicoDstMaker");
   gSystem->Load("StPicoPrescales");
   gSystem->Load("StBTofUtil");
   gSystem->Load("StPicoD0EventMaker");
   gSystem->Load("StPicoHFMaker");
   gSystem->Load("StPicoDstarAnaMaker");



   
   chain = new StChain();

   // create list of picoDst files
   TString command = "sed 's/hft\\\/d0tree/picodsts/g' " + d0list + " >correspondingPico.list";
   gSystem->Exec(command.Data());
   command = "sed -i 's/picoD0/picoDst/g' correspondingPico.list";
   gSystem->Exec(command.Data());
   StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0, "correspondingPico.list", "picoDstMaker");
   StPicoDstarAnaMaker*  picoDstarAnaMaker = new StPicoDstarAnaMaker("picoDstarAnaMaker", d0list, outFileName.Data(), picoDstMaker);

   StHFCuts* dStarCuts = new StHFCuts("dStarCuts");

   /* -------------- USER variables -------------------------

   // -- File name of bad run list
   dStarCuts->setBadRunListFileName(badRunListFileName);

   // add your cuts here.

   // tracking
   dStarCuts->setCutNHitsFitMax(20);

   // pions
   dStarCuts->setCutTPCNSigmaPion(3.0);

   // kaons
   dStarCuts->setCutTPCNSigmaKaon(2.0);

   // kaonPion pair cuts
   float dcaDaughtersMax = 0.008;  // maximum
   float decayLengthMin  = 0.0030; // minimum
   float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
   float cosThetaMin     = 0.90;   // minimum
   float minMass         = 1.6;
   float maxMass         = 2.1;
   dStarCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass);
   */
   // ========================================================================================
   // HF cuts for Dstar only
   // ========================================================================================
   StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");
   picoDstarAnaMaker->setHFCuts(hfCuts);
   hfCuts->setBadRunListFileName(badRunListFileName);
   hfCuts->setCutNHitsFitMax(15); 
   hfCuts->setCutRequireHFT(true);
   hfCuts->setCutSecondaryPair(0.01, 0.003, 999999., 
			       0.99, 1.7, 2.0);

   hfCuts->setCutPtRange(0.6,20.0,StHFCuts::kPion);
   hfCuts->setCutPtRange(0.6,20.0,StHFCuts::kKaon);
   //TPC setters
   hfCuts->setCutTPCNSigmaPion(3.0);
   hfCuts->setCutTPCNSigmaKaon(3.0);
   //TOF setters, need to set pt range as well
  hfCuts->setCutTOFDeltaOneOverBeta(0.05, StHFCuts::kKaon);
  hfCuts->setCutPtotRangeHybridTOF(0.,20.0,StHFCuts::kKaon);
  hfCuts->setCutTOFDeltaOneOverBeta(0.05, StHFCuts::kPion);
  hfCuts->setCutPtotRangeHybridTOF(0.,20.0,StHFCuts::kPion);
  // ========================================================================================
   chain->Init();
   int nEntries = picoDstarAnaMaker->getEntries();
   cout<<"Loading "<<nEntries<<" events"<<endl;
   for (int iEvent = 0; iEvent < nEntries; ++iEvent)
     {
       if(iEvent%10000==0)
	 cout << "Working on eventNumber " << iEvent << endl;
      chain->Clear();
      int iret = chain->Make();
      if (iret)
      {
         cout << "Bad return code!" << iret << endl;
         break;
      }
   }

   chain->Finish();
   delete chain;

   // delete list of picos
   command = "rm -f correspondingPico.list";
   gSystem->Exec(command.Data());

}
