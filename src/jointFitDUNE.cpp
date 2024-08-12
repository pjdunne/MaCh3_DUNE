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

#include "samplePDFDUNE/samplePDFDUNEAtmBase.h"
#include "mcmc/mcmc.h"

int main(int argc, char **argv)
{

  #ifdef MULTITHREAD
    std::cout << "MUTLI THREADING IS ON" << std::endl;
  #else
    std::cout << "MUTLI THREADING IS OFF" << std::endl;
  #endif

  #ifdef CPU_ONLY
    std::cout << "GPU OSC CALC IS ON" << std::endl;
  #else
    std::cout << "ONLY CPU OSC CALC" << std::endl;
  #endif

  // Read config
  manager *fitMan = new manager(argv[1]);

  std::string  XsecMatrixFile = fitMan->raw()["General"]["Systematics"]["XsecCovFile"].as<std::string>(); 
  std::string  XsecMatrixName = fitMan->raw()["General"]["Systematics"]["XsecCovName"].as<std::string>();
  std::string  OscMatrixFile = fitMan->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
  std::string  OscMatrixName = fitMan->raw()["General"]["Systematics"]["OscCovName"].as<std::string>();
 
  std::string Outfile = fitMan->raw()["General"]["Output"]["FileName"].as<std::string>();

  // Bool about whether or not to use the near detector
  // If you're doing a stats-only fit or if you're using the BANFF matrix, you don't want to use
  // the near detector. (In the case of a BANFF matrix fit, using the ND will do absolutely nothing.
  // In the case of a stats-only fit it will add a constant value to the LLH, which won't affect
  // the fit results. In both cases it takes about an hour to load the data/MC and set up the ND
  // class, so that's just a massive waste of time and let's not do it).
  
  //###########################################################################################################
  // Covariance Objects

  //covarianceNDDet_2019Poly* det;
  covarianceXsec *xsec;
  covarianceOsc *osc;

  xsec = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str());


  // Setting flat priors based on XSECPARAMFLAT list in configuration file 

  bool statsonly = false;

  //std::vector<double> oscpars{0.307,0.528,0.0218,7.53e-5, 2.509e-3,-1.601}; // Asimov A
  std::vector<double> oscpars{0.310,0.582,0.0224,7.39e-5,2.525e-3,-2.498}; //NuFit 4.0 NH
  //
  // set up oscillation parameters (init params are in constructor)
  osc = new covarianceOsc(OscMatrixName.c_str(), OscMatrixFile.c_str());

  // oscpars from manager in order:
  // sin2th_12, sin2th_23, sin2th_13, delm2_12, delm2_23, delta_cp

  if (!(oscpars.size()==6 || oscpars.size()==7))
    {
      std::cout<<"Input osc pars not of right size, there should be six entries (or seven if setting beta)"<<std::endl;
      std::cout<<"oscpars.size() = " << oscpars.size() << std::endl;
      exit(1);
    }

  // Add beta (if set in config file)

  std::cout<<"Using these oscillation parameters: ";
  for(unsigned ipar=0;ipar<oscpars.size();ipar++)
    std::cout<<" "<<oscpars.at(ipar);
  std::cout << std::endl;

  osc->setFlipDeltaM23(true);

  // Ask config file whether to use reactor constraint

  // Use prior for 12 parameters only
  //osc->setEvalLikelihood(0,false);
  osc->setEvalLikelihood(1,false);
  osc->setEvalLikelihood(2,false);
  //osc->setEvalLikelihood(3,false);
  osc->setEvalLikelihood(4,false);
  osc->setEvalLikelihood(5,false);

  osc->setParameters(oscpars);
  osc->acceptStep();


  xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());
  osc->setStepScale(0.045);

  //###########################################################################################################
  // SK Samples

  // Hardcoded binning_opt = 2 (erec-theta nue binning)
  std::cout << "Fitting electorn rings in Erec-theta. " << std::endl ;

  // make PDFs
  bool domaqeh=true;
  std::cout << "I AM CPU" << std::endl;
  std::vector<samplePDFDUNEAtmBase*> pdfs;
  //!!add pdfs vector to make things easier
  samplePDFDUNEAtmBase *numu_pdf = new samplePDFDUNEAtmBase(1.3628319e+23, "configs/SamplePDFDune_FHC_numuselec.yaml", xsec);
  samplePDFDUNEAtmBase *numubar_pdf = new samplePDFDUNEAtmBase(1.3628319e+23, "configs/SamplePDFDune_RHC_numuselec.yaml", xsec);
  samplePDFDUNEAtmBase *nue_pdf = new samplePDFDUNEAtmBase(1.3628319e+23, "configs/SamplePDFDune_FHC_nueselec.yaml", xsec);
  samplePDFDUNEAtmBase *nuebar_pdf = new samplePDFDUNEAtmBase(1.3628319e+23, "configs/SamplePDFDune_RHC_nueselec.yaml", xsec);
  pdfs.push_back(numu_pdf);
  pdfs.push_back(numubar_pdf);
  pdfs.push_back(nue_pdf);
  pdfs.push_back(nuebar_pdf);

  for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
    pdfs[ipdf]->UseNonDoubledAngles(true);
    pdfs[ipdf] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
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

 

  if(fitMan->raw()["General"]["StartFromPos"].as<bool>()) {//Start from values at the end of an already run chain
    //Read in paramter names and values from file
    std::cout << "MCMC getting starting position from " << fitMan->raw()["General"]["PosFileName"].as<std::string>() << std::endl;
    TFile *infile = new TFile(fitMan->raw()["General"]["PosFileName"].as<std::string>().c_str(), "READ");
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
    covtypes.push_back("xsec");
 
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
  
  // set to nominal
  xsec->setParameters();

  //###########################################################################################################
  // Apply reweight and find event rates

  std::vector<double> CAFana_default_edges = get_default_CAFana_bins();
  for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
    //pdfs[ipdf]->useBinnedOscReweighting(true, CAFana_default_edges.size()-1, &CAFana_default_edges[0]);
    pdfs[ipdf]->reweight(osc->getPropPars());     // Get nominal predictions
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


  
  // APPLY BANFF TUNING TO ASIMOV INPUT (Don't apply as a prior for fit, just use for Asimov "data")
  // Note: this is hardcoded for BANFF 2017b

  std::cout << "-----------------------" << std::endl; 

  // set nominal
  xsec->setParameters();


  numu_pdf->reweight(osc->getPropPars());
  TH1D *numu_asimov = (TH1D*)numu_pdf->get1DHist()->Clone("numu_asimov");
  nue_pdf->reweight(osc->getPropPars());
  TH1D *nue_asimov = (TH1D*)nue_pdf->get1DHist()->Clone("nue_asimov");
  numubar_pdf->reweight(osc->getPropPars());
  TH1D *numubar_asimov = (TH1D*)numubar_pdf->get1DHist()->Clone("numubar_asimov");
  nuebar_pdf->reweight(osc->getPropPars());
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

  


  //###########################################################################################################

  // Back to actual nominal for fit
  // If starting from end values of a previous chain set everything to the values from there
  // and do acceptStep() to update fParCurr with these values
  //
  if(fitMan->raw()["General"]["StartFromPos"].as<bool>()) {
    if(parstarts.find("xsec")!=parstarts.end()) {
      xsec->setParameters(parstarts["xsec"]);
      xsec->acceptStep();
    }
    else {xsec->setParameters();}
  }
  else {
    xsec->setParameters();
  }

  //###########################################################################################################
  // MCMC

  mcmc *markovChain = new mcmc(fitMan);

  //numu_fds->Write();

  // set up
  int NSteps = fitMan->raw()["General"]["NSTEPS"].as<int>(); 
  markovChain->setChainLength(NSteps);
  markovChain->addOscHandler(osc, osc);
  if(lastStep > 0) markovChain->setInitialStepNumber(lastStep+1);

  // add samples
  markovChain->addSamplePDF(numu_pdf);
  markovChain->addSamplePDF(numubar_pdf);
  markovChain->addSamplePDF(nue_pdf);
  markovChain->addSamplePDF(nuebar_pdf);

  //start chain from random position
  xsec->throwParameters();
  osc->throwParameters();

  // add systematic objects
  if (!statsonly) {
    markovChain->addSystObj(xsec);
  }
  else {std::cout << "STATS ONLY!" << std::endl;}
  
  // run!
  markovChain->runMCMC();

  return 0;
}
