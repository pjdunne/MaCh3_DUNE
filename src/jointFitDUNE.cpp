#include <ctime>
#include <iostream>
#include <vector>
#include <sys/time.h>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStopwatch.h>

#include "samplePDFDUNE/samplePDFDUNEBase.h"
#include "mcmc/mcmc.h"

//#include "covariance/BANFF_to_xsec_mapping.h"

  //CAFAna binned Oscillation 
std::vector<double> get_default_CAFana_bins(){
   // From CAFana - probability binning -
   const int kNumTrueEnergyBins = 100; 
   // N+1 bin low edges
   std::vector<double> edges(kNumTrueEnergyBins+1);
  
   const double Emin = 0.5; // 500 MeV: there's really no events below there
  
   // How many edges to generate. Allow room for 0-Emin bi            const double N = kNumTrueEnergyBins-1;
   const double N = kNumTrueEnergyBins-1;
   const double A = N*Emin;
  
   edges[0] = 0;
  
   for(int i = 1; i <= N; ++i){
     edges[kNumTrueEnergyBins-i] = A/i;
   }
  
   edges[kNumTrueEnergyBins] = 120; // Replace the infinity that would be here
   return edges;
                                 
}


int main(int argc, char **argv)
{
  std::cout << "IM HERE" << std::endl;
  manager *fitMan = new manager(argv[1]);
  std::cout << "Using POT: " << fitMan->GetPOT()  << "(nu mode), " << fitMan->GetNubarPOT() << "(anti-nu mode)" << std::endl;

  // Fake data set
  bool fakedata = fitMan->GetFakeDataFitFlag();
  // Real data fit
  bool datafit = fitMan->GetRealDataFitFlag();
  // Toy fit
  bool toyfit = fitMan->GetToyFitFlag();
  // Asimov fit
  bool asimovfit = fitMan->GetAsimovFitFlag();

  // there's a check inside the manager class that does this; left here for demonstrative purposes
  if (fitMan->GetGoodConfig() == false) {
    std::cerr << "Didn't find a good config in input configuration" << std::endl;
    throw;
  }

  // Do a stats-only fit?
  // Note: only works with asimov fit (at the moment!)
  bool statsonly = fitMan->GetStatOnly();
  if (statsonly) {std::cout << "Doing a stats-only fit, ignoring ND280 data and using BANFF tuning" << std::endl;}
  

  // Bool about whether or not to use the near detector
  // If you're doing a stats-only fit or if you're using the BANFF matrix, you don't want to use
  // the near detector. (In the case of a BANFF matrix fit, using the ND will do absolutely nothing.
  // In the case of a stats-only fit it will add a constant value to the LLH, which won't affect
  // the fit results. In both cases it takes about an hour to load the data/MC and set up the ND
  // class, so that's just a massive waste of time and let's not do it).
  
  ///Let's ask the manager what are the file with covariance matrix
  TString fluxCovMatrixFile  = fitMan -> GetFluxCovMatrix() ;
  TString fluxCovMatrixName  = fitMan -> GetFluxCovName() ;
  TString xsecCovMatrixFile  = fitMan -> GetXsecCovMatrix() ;
  TString xsecCovMatrixName  = fitMan -> GetXsecCovName() ;

  //###########################################################################################################
  // Covariance Objects

  //covarianceNDDet_2019Poly* det;
  covarianceXsec *xsec;
  covarianceOsc *osc;

  xsec= new covarianceXsec(xsecCovMatrixName, xsecCovMatrixFile);


  // Setting flat priors based on XSECPARAMFLAT list in configuration file 
  std::vector<int> XsecFlatParams  = fitMan->GetXsecFlat();
  if (XsecFlatParams.size() == 1 && XsecFlatParams.at(0) == -1) {
    for (int j = 0; j < xsec->getSize(); j++) {
      xsec->setEvalLikelihood(j, false);
    }
  } else {
    for (unsigned j = 0; j < XsecFlatParams.size(); j++) {
      xsec->setEvalLikelihood(XsecFlatParams.at(j), false);
    }
  }


  std::vector<double> oscpars;
  // set up oscillation parameters (init params are in constructor)
  if (fitMan->GetUseBeta()) // if using beta in nue/nuebar oscillation probability
    {
      std::cerr << "[ERROR]: Not setup for beta in 2020." << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
      //osc = new covarianceOsc("osc_cov","inputs/oscillation_covariance_6par_nondouble_beta.root");
    }
  else // Default: normal oscillation covariance matrix (not including beta)
    {
      osc = new covarianceOsc("osc_cov","inputs/oscillation_covariance_6par_nondouble_PDG2019.root");
    }

  // oscpars from manager in order:
  // sin2th_12, sin2th_23, sin2th_13, delm2_12, delm2_23, delta_cp
  oscpars=fitMan->GetOscParameters();

  if (!(oscpars.size()==6 || oscpars.size()==7))
    {
      std::cout<<"Input osc pars not of right size, there should be six entries (or seven if setting beta)"<<std::endl;
      std::cout<<"oscpars.size() = " << oscpars.size() << std::endl;
      exit(1);
    }

  // Add beta (if set in config file)
  if (fitMan->GetUseBeta() && oscpars.size()==6)
    {
      oscpars.push_back(1);
      std::cout << "Using beta in nuebar appearance probability. Input value: beta = 1" << std::endl;
    }

  std::cout<<"Using these oscillation parameters: ";
  for(unsigned ipar=0;ipar<oscpars.size();ipar++)
    std::cout<<" "<<oscpars.at(ipar);
  std::cout << std::endl;

  osc->setFlipDeltaM23(true);

  // Ask config file whether to use reactor constraint
  bool useRC =  fitMan -> GetRC() ;
  std::cout << "use reactor prior is : " << useRC << std::endl ;
  osc->useReactorPrior(useRC); // this is hard coded inside, and is bad

  // Use prior for 12 parameters only
  //osc->setEvalLikelihood(0,false);
  osc->setEvalLikelihood(1,false);
  osc->setEvalLikelihood(2,false);
  //osc->setEvalLikelihood(3,false);
  osc->setEvalLikelihood(4,false);
  osc->setEvalLikelihood(5,false);

  osc->setParameters(oscpars);
  osc->acceptStep();


  xsec->setStepScale(fitMan->GetXsecStepScale());
  osc->setStepScale(fitMan->GetOscStepScale());

  //###########################################################################################################
  // SK Samples

  // Hardcoded binning_opt = 2 (erec-theta nue binning)
  std::cout << "Fitting electorn rings in Erec-theta. " << std::endl ;

  // make PDFs
  bool domaqeh=true;
  std::cout << "I AM CPU" << std::endl;
  std::vector<samplePDFDUNEBase*> pdfs;
  //!!add pdfs vector to make things easier
  samplePDFDUNEBase *numu_pdf = new samplePDFDUNEBase(fitMan->GetPOT(), "configs/SamplePDFDune_FHC_numuselec.cfg", xsec);
  samplePDFDUNEBase *numubar_pdf = new samplePDFDUNEBase(fitMan->GetNubarPOT(), "configs/SamplePDFDune_RHC_numuselec.cfg", xsec);
  samplePDFDUNEBase *nue_pdf = new samplePDFDUNEBase(fitMan->GetPOT(), "configs/SamplePDFDune_FHC_nueselec.cfg", xsec);
  samplePDFDUNEBase *nuebar_pdf = new samplePDFDUNEBase(fitMan->GetNubarPOT(), "configs/SamplePDFDune_RHC_nueselec.cfg", xsec);
  pdfs.push_back(numu_pdf);
  pdfs.push_back(numubar_pdf);
  pdfs.push_back(nue_pdf);
  pdfs.push_back(nuebar_pdf);
//#endif  

  for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
    pdfs[ipdf]->setUseBeta(fitMan->GetUseBeta());   // Set whether to use beta in oscillation probability calculation
    pdfs[ipdf]->setApplyBetaNue(fitMan->GetApplyBetaNue());   // Set option to apply beta to nue appearance probability instead of nuebar appearance probability
    pdfs[ipdf]->setApplyBetaDiag(fitMan->GetApplyBetaDiag());        // Set option to apply (1/beta) to nue appearance probability and beta to nuebar appearance probability
    pdfs[ipdf]->useNonDoubledAngles(true);
  }


  
  // tell PDF where covariances are
  //for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
    //pdfs[ipdf]->setXsecCov(xsec);
  //}

  //###########################################################################################################
  // Set covariance objects equal to output of previous chain
  
  std::vector<double> oscparstarts;
  std::map<TString,std::vector<double> > parstarts;
  int lastStep = 0;

 /*

  if(fitMan->GetStartFromPosterior()){//Start from values at the end of an already run chain
    //Read in paramter names and values from file
    std::cout << "MCMC getting starting position from " << fitMan->GetPosteriorFiles() << std::endl;
    TFile *infile = new TFile(fitMan->GetPosteriorFiles().c_str(), "READ");
    TTree *posts = (TTree*)infile->Get("posteriors");
    TObjArray* brlis = (TObjArray*)posts->GetListOfBranches();
    int nbr = brlis->GetEntries();
    TString branch_names[nbr];
    double branch_vals[nbr];
    int step_val;
    for (int i = 0; i < nbr; ++i) {
      TBranch *br = (TBranch*)brlis->At(i);
      TString bname = br->GetName();
      branch_names[i] = bname;
      if(branch_names[i] == "step") {
        posts->SetBranchAddress("step",&step_val);
        continue;
      }
      std::cout << " * Loading " << bname << std::endl;
      posts->SetBranchAddress(branch_names[i], &branch_vals[i]);
    }
    posts->GetEntry(posts->GetEntries()-1);
    std::map<TString,double> init_pars;
    for (int i = 0; i < nbr; ++i) {
      init_pars.insert( std::pair<TString, double>(branch_names[i], branch_vals[i]));
    }
    infile->Close();
    delete infile;
    
    //Make vectors of parameter value for each covariance
    std::vector<TString> covtypes;
    else{
      covtypes.push_back("xsec");
      covtypes.push_back("b");
      covtypes.push_back("ndd");
    }
    covtypes.push_back("skd_joint");
 
    for(unsigned icov=0;icov<covtypes.size();icov++){
      std::vector<double> covparstarts;
      std::map<TString, double>::const_iterator it;
      int iPar=0;
      while(it != init_pars.end()){
  	it = init_pars.find(covtypes[icov]+"_"+TString::Format("%d",iPar));
  	if (it != init_pars.end()) {
  	  covparstarts.push_back(it->second);
  	}
  	iPar++;
      }
      if(covparstarts.size()!=0) parstarts.insert(std::pair<TString,std::vector<double> >(covtypes[icov],covparstarts));
      else std::cout<<"Did not find any parameters in posterior tree for: "<<covtypes[icov]<<std::endl<<"assuming previous chain didn't use them"<<std::endl;
    }

    std::map<TString, double>::const_iterator itt;

    // set the oscillation parameters
    itt = init_pars.find("sin2th_12");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("sin2th_23");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("sin2th_13");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("delm2_12");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("delm2_23");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("delta_cp");
    oscparstarts.push_back(itt->second);

    lastStep = step_val;

  }
  */
  // set to nominal
  xsec->setParameters();

  //###########################################################################################################
  // Apply reweight and find event rates

  std::vector<double> CAFana_default_edges = get_default_CAFana_bins();
  for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
    pdfs[ipdf]->useBinnedOscReweighting(true, CAFana_default_edges.size()-1, &CAFana_default_edges[0]);
    pdfs[ipdf]->reweight(osc->getPropPars(), osc->getPropPars());     // Get nominal predictions
  }

  TH1D *numu_nominal = (TH1D*)numu_pdf->get1DHist()->Clone("numu_nominal");
  TH1D *numubar_nominal = (TH1D*)numubar_pdf->get1DHist()->Clone("numubar_nominal");
  TH1D *nue_nominal = (TH1D*)nue_pdf->get1DHist()->Clone("nue_nominal");
  TH1D *nuebar_nominal = (TH1D*)nuebar_pdf->get1DHist()->Clone("nuebar_nominal");

  std::cout << "Nominal event rates: " << std::endl;
  std::cout << "numu: " << numu_nominal->Integral() << std::endl;
  std::cout << "numu-bar: " << numubar_nominal->Integral() << std::endl;
  std::cout << "nue: " << nue_nominal->Integral() << std::endl;
  std::cout << "nue-bar: " << nuebar_nominal->Integral() << std::endl;

  // Check what type of fit you are doing
  std::vector<double> *numudat = 0;
  std::vector<double> *numubardat = 0;
  std::vector<std::vector<double> > *nuedat = 0;
  std::vector<std::vector<double> > *nuebardat = 0;

    
  //###########################################################################################################
  // Determine the type of fit we want to do

  TH1D *numu_data;
  TH1D *numubar_data;
  TH1D *nue_data;
  TH1D *nuebar_data;

  if (toyfit && !datafit && !fakedata && !asimovfit) // Toy fit
  {
    std::cout << "Loading toy data from " << fitMan->GetToyFilename().c_str() << std::endl;
    std::cout<< "Warning not implemented for CC1pi yet, just doing 4 sample!!"<<std::endl;

    int n_seeds = 4;
    std::cout << "- Using seed " << fitMan->GetSeed() << std::endl;
    TRandom3 rnd_machine(fitMan->GetSeed());

    // waste a few throws
    for (int h = 0; h < 100; ++h)
      double waste = rnd_machine.Rndm();

    std::vector<int> seed_array;

    for (int i = 0; i < n_seeds; ++i)
    {
      double rr = 1000000 * rnd_machine.Rndm();
      seed_array.push_back(rr);
      std::cout << "- Number " << i << " seed " << rr << std::endl;
    }

    // set the covariance objects to use nominal values thrown from the covariance matrix
    xsec->throwNominal(false, seed_array[1]);
    // if (useND280) {
    //   det->throwNominal(false, seed_array[2]);
    // }
    std::cout << "- Toy number " << fitMan->GetNtoy() << " from file " << fitMan->GetToyFilename().c_str() << std::endl;
    TFile *dat = new TFile(fitMan->GetToyFilename().c_str(), "READ");
    TTree *dattree = (TTree*)dat->Get("toyMC_2015");

    // Change this bit for new names
    dattree->SetBranchAddress("sk_numu", &numudat);
    dattree->SetBranchAddress("sk_nue", &nuedat);
    dattree->SetBranchAddress("sk_numubar", &numubardat);
    dattree->SetBranchAddress("sk_nuebar", &nuebardat);

    dattree->GetEntry(fitMan->GetNtoy());

    numu_pdf->addData(*numudat);
    nue_pdf->addData(*nuedat);
    numubar_pdf->addData(*numubardat);
    nuebar_pdf->addData(*nuebardat);

    
    delete dattree;
    delete dat; //dat->Close();
  }
  else if (fakedata && !datafit && !asimovfit) // Fake data
  {
    //Find out from Leila what to do about these files, maybe take straight from MC
    std::cout << "Loading fake data set from " << fitMan->GetDataFilename() << std::endl;
    TFile *f = new TFile(fitMan->GetDataFilename().c_str(), "OPEN");
    numu_data = (TH1D*)f->Get("numu");
    numu_data->SetDirectory(0);
    numubar_data = (TH1D*)f->Get("numubar");
    numubar_data->SetDirectory(0); // this keeps the TH1D alive after the file is closed
    nue_data = (TH1D*)f->Get("nue");
    nue_data->SetDirectory(0);
    nuebar_data = (TH1D*)f->Get("nuebar");
    nuebar_data->SetDirectory(0); // this keeps the TH1D alive after the file is closed

    f->Close();

    // SK
    //check these behave nicely in 2D
    numu_pdf->addData(numu_data);
    nue_pdf->addData(nue_data);
    numubar_pdf->addData(numubar_data);
    nuebar_pdf->addData(nuebar_data);

  }
  else if (asimovfit && !datafit) // Asimov fit
  {

    // -------------------------------------------------------------------- //
    // APPLY BANFF TUNING TO ASIMOV INPUT (Don't apply as a prior for fit, just use for Asimov "data")
    // Note: this is hardcoded for BANFF 2017b

      std::cout << "-----------------------" << std::endl; 

      // set nominal
      xsec->setParameters();


    numu_pdf->reweight(osc->getPropPars(), osc->getPropPars());
    TH1D *numu_asimov = (TH1D*)numu_pdf->get1DHist()->Clone("numu_asimov");
    nue_pdf->reweight(osc->getPropPars(), osc->getPropPars());
    TH1D *nue_asimov = (TH1D*)nue_pdf->get1DHist()->Clone("nue_asimov");
    numubar_pdf->reweight(osc->getPropPars(), osc->getPropPars());
    TH1D *numubar_asimov = (TH1D*)numubar_pdf->get1DHist()->Clone("numubar_asimov");
    nuebar_pdf->reweight(osc->getPropPars(), osc->getPropPars());
    TH1D *nuebar_asimov = (TH1D*)nuebar_pdf->get1DHist()->Clone("nuebar_asimov");

    // Print event rates to check
    std::cout << "-------- SK event rates for Asimov fit (Asimov fake data) ------------" << std::endl;
    std::cout << "FHC 1Rmu:   " << numu_asimov->Integral() << std::endl;
    std::cout << "FHC 1Re:    " << nue_asimov->Integral() << std::endl;
    std::cout << "RHC 1Rmu:   " << numubar_asimov->Integral() << std::endl;
    std::cout << "RHC 1Re:    " << nuebar_asimov->Integral() << std::endl;

    numu_pdf->addData(numu_asimov); 
    nue_pdf->addData(nue_asimov); 
    numubar_pdf->addData(numubar_asimov);  
    nuebar_pdf->addData(nuebar_asimov); 

    
  }
  else if (datafit)
  {
    std::cout << "Performing a real data fit from " << fitMan->GetDataFilename() << std::endl;

    TFile *f = new TFile(fitMan->GetDataFilename().c_str(), "OPEN");
    numu_data = (TH1D*)f->Get("numu");
    numu_data->SetDirectory(0);
    numubar_data = (TH1D*)f->Get("numubar");
    numubar_data->SetDirectory(0); 
    nue_data = (TH1D*)f->Get("nue_erecTheta");
    nue_data->SetDirectory(0);
    nuebar_data = (TH1D*)f->Get("nuebar_erecTheta");
    nuebar_data->SetDirectory(0); 
    delete f;

    // SK
    //check this behaves nicely with 2D
    numu_pdf->addData(numu_data);
    nue_pdf->addData(nue_data);
    numubar_pdf->addData(numubar_data);
    nuebar_pdf->addData(nuebar_data);

    // ND (don't have lines with 'sample->addData(xxx)', and it will set the data as the data automatically)
  }
  else
    std::cerr << "No valid option given for fit type" << std::endl;


  //###########################################################################################################

  // Back to actual nominal for fit
  // If starting from end values of a previous chain set everything to the values from there
  // and do acceptStep() to update fParCurr with these values
  if(fitMan->GetStartFromPosterior()) {
    if(parstarts.find("xsec")!=parstarts.end()) {
	xsec->setParameters(parstarts["xsec"]);
	xsec->acceptStep();
      }
      else xsec->setParameters();
  }

  else {
    xsec->setParameters();
  }

  //###########################################################################################################
  // MCMC

  //mcmc *markovChain = new mcmc(fitMan->getOutputFilename(), true);
  mcmc *markovChain = new mcmc(fitMan);

  //numu_fds->Write();

  // set up
  markovChain->setChainLength(fitMan->GetNSteps());
  markovChain->addOscHandler(osc, osc);
  if(lastStep > 0) markovChain->setInitialStepNumber(lastStep+1);

  // add samples
  // if (useND280) {
  //   markovChain->addSamplePDF(sample);
  // }
  markovChain->addSamplePDF(numu_pdf);
  markovChain->addSamplePDF(numubar_pdf);
  markovChain->addSamplePDF(nue_pdf);
  markovChain->addSamplePDF(nuebar_pdf);

  //!!add cc1pi

  // add systematic objects
  if (!statsonly) {
    markovChain->addSystObj(xsec);
    // if (useND280) {
    //   markovChain->addSystObj(det);
    // }
  }
  
  // run!
  markovChain->runMCMC();

  return 0;
}
