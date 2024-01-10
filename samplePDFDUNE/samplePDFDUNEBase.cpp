#include <TROOT.h>
#include "samplePDFDUNEBase.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "manager/manager.h"


//#define DEBUG

// Constructors for erec-binned errors

//!!rewrite execs to give arguments in new order
samplePDFDUNEBase::samplePDFDUNEBase(double pot, std::string mc_version, covarianceXsec* xsec_cov)
  : samplePDFBase(pot)
{ 
  std::cout << "- Using DUNE sample config in this file " << mc_version << std::endl;
  //ETA - safety feature so you can't pass a NULL xsec_cov
  if(xsec_cov == NULL){std::cerr << "[ERROR:] You've passed me a NULL xsec covariance matrix... I need this to setup splines!" << std::endl; throw;}
  init(pot, mc_version, xsec_cov);          

}

samplePDFDUNEBase::~samplePDFDUNEBase()
{
}

void samplePDFDUNEBase::init(double pot, std::string samplecfgfile, covarianceXsec *xsec_cov)
{

  Beta=1;
  useBeta=false;
  applyBetaNue=false;
  applyBetaDiag=false;

  //doubled_angle =true ;
  UseNonDoubledAngles(true);
  if (doubled_angle) std::cout << "- Using non doubled angles for oscillation parameters" << std::endl;

  osc_binned = false;
  if (osc_binned) std::cout << "- Using binned oscillation weights" << std::endl;

  modes = new TH1D("modes","",120,-60,60);

  // Create the libconfig object
  //libconfig::Config samplecfg;
  // Could turn on automatic type conversion in the future, but could be a bad idea..
  // cfg.setAutoConvert(true);

  std::string mtupleprefix;
  std::string mtuplesuffix;
  std::string splineprefix;
  std::string splinesuffix;

  char* sample_char = (char*)samplecfgfile.c_str();
  //ETA - trying to read stuff from yaml file
  manager* SampleManager = new manager(sample_char);

  //Bools
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  //Inputs
  mtupleprefix = SampleManager->raw()["InputFiles"]["mtupleprefix"].as<std::string>();
  mtuplesuffix = SampleManager->raw()["InputFiles"]["mtuplesuffix"].as<std::string>();
  splineprefix = SampleManager->raw()["InputFiles"]["splineprefix"].as<std::string>();
  splinesuffix = SampleManager->raw()["InputFiles"]["splinesuffix"].as<std::string>();

  //Binning
  BinningOpt = SampleManager->raw()["Binning"]["BinningOpt"].as<int>();  
  std::vector<double> sample_erec_bins = SampleManager->raw()["Binning"]["XVarBins"].as<std::vector<double>>();
  std::vector<double> sample_theta_bins = SampleManager->raw()["Binning"]["YVarBins"].as<std::vector<double>>();

  samplename = SampleManager->raw()["SampleName"].as<std::string>();

  std::vector<std::string> mtuple_files;
  std::vector<std::string> spline_files;
  std::vector<int> sample_vecno;
  std::vector<int> sample_oscnutype;
  std::vector<int> sample_nutype;
  std::vector<bool> sample_signal;
  
  //Loop over all the sub-samples
  for (auto const &osc_channel : SampleManager->raw()["SubSamples"]) {
	std::cout << "Found sub sample" << std::endl;
	mtuple_files.push_back(osc_channel["mtuplefile"].as<std::string>());
	spline_files.push_back(osc_channel["splinefile"].as<std::string>());
	sample_vecno.push_back(osc_channel["samplevecno"].as<int>());
	sample_nutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["nutype"].as<int>())));
	sample_oscnutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["oscnutype"].as<int>())));
	sample_signal.push_back(osc_channel["signal"].as<bool>());
  }


  //Now loop over the kinematic cuts
  for ( auto const &SelectionCuts : SampleManager->raw()["SelectionCuts"]) {
	std::cout << "Looping over selection cuts " << std::endl;
	SelectionStr.push_back(SelectionCuts["KinematicStr"].as<std::string>());

	SelectionBounds.push_back(SelectionCuts["Bounds"].as<std::vector<double>>());
	std::cout << "Found cut on string " << SelectionCuts["KinematicStr"].as<std::string>() << std::endl;
	std::cout << "With bounds " << SelectionCuts["Bounds"].as<std::vector<double>>()[0] << " to " << SelectionCuts["Bounds"].as<std::vector<double>>()[1] << std::endl;
  }
  NSelections = SelectionStr.size();

  //Make some arrays so we can initialise _hPDF1D and _hPDF2D with these
  double erec_bin_edges[sample_erec_bins.size()];
  double theta_bin_edges[sample_theta_bins.size()];
  for(unsigned erec_i = 0 ; erec_i < sample_erec_bins.size() ; erec_i++){erec_bin_edges[erec_i] = sample_erec_bins[erec_i];}
  for(unsigned theta_i = 0 ; theta_i < sample_theta_bins.size() ; theta_i++){theta_bin_edges[theta_i] = sample_theta_bins[theta_i];}

  // create dunemc storage
  int nSamples = SampleManager->raw()["NSubSamples"].as<int>();
  for (int i=0;i<nSamples;i++) {
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
  }

  //Now down with yaml file for sample
  delete SampleManager;
  std::cout << "Oscnutype size: " << sample_oscnutype.size() << ", dunemcSamples size: " << dunemcSamples.size() << endl;  
  if(sample_oscnutype.size() != dunemcSamples.size()){std::cerr << "[ERROR:] samplePDFDUNEBase::samplePDFDUNEBase() - something went wrong either getting information from sample config" << std::endl; throw;}

  for(unsigned iSample=0 ; iSample < dunemcSamples.size() ; iSample++){
    setupDUNEMC((mtupleprefix+mtuple_files[iSample]+mtuplesuffix).c_str(), &dunemcSamples[sample_vecno[iSample]], pot, sample_nutype[iSample], sample_oscnutype[iSample], sample_signal[iSample]);
  }

  for (int i=0;i<nSamples;i++) {
    struct fdmc_base obj = fdmc_base();
    MCSamples.push_back(obj);
  }

  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    setupFDMC(&dunemcSamples[sample_vecno[iSample]], &MCSamples[sample_vecno[iSample]]);
  }

  std::cout << "################" << std::endl;
  std::cout << "Setup FD MC   " << std::endl;
  std::cout << "################" << std::endl;

  // ETA - If xsec_cov hasn't been passed to the samplePDFDUNEBase constructor then it's NULL
  // and the old funcitonality is kepy
  // this calls this function in the core code
  // this needs to come after setupFDMC as otherwise MCSamples.splinefile will be NULL
  std::cout << "xsec_cov is: " << xsec_cov << std::endl;
  SetXsecCov(xsec_cov); 

  std::vector<std::string> spline_filepaths;

  std::cout << "Now setting up Splines" << std::endl;
  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
	setupSplines(&MCSamples[sample_vecno[iSample]] , (splineprefix+spline_files[iSample]+splinesuffix).c_str(), MCSamples[iSample].nutype);
  spline_filepaths.push_back(splineprefix+spline_files[iSample]+splinesuffix);
  }

  splineFile = new splinesDUNE(xsec_cov);

  //////////////////////////////////
  // Now add samples to spline monolith
  //////////////////////////////////

  //ETA - do we need to do this here?
  //Can't we have one splineFile object for all samples as Dan does in atmospherics fit?
  //Then just add spline files to monolith?
  std::cout<<"Adding samples to spline monolith"<<std::endl;
  std::cout << "samplename is " << samplename << std::endl;
  std::cout << "BinningOpt is " << BinningOpt << std::endl;
  std::cout << "SampleDetID is " << SampleDetID << std::endl;
  std::cout << "spline_filepaths is of size " << spline_filepaths.size() << std::endl;

  splineFile->AddSample(samplename, BinningOpt, SampleDetID, spline_filepaths);

  // Print statements for debugging
  splineFile->PrintArrayDimension();
  splineFile->CountNumberOfLoadedSplines(false, 1);
  splineFile->TransferToMonolith();
  std::cout << "--------------------------------" <<std::endl;

  std::cout << "################" << std::endl;
  std::cout << "Setup FD splines   " << std::endl;
  std::cout << "################" << std::endl;

  SetupWeightPointers();

  fillSplineBins();

  #ifdef USE_PROB3
  std::cout << "- Setup Prob3++" << std::endl;
  #else
  std::cout << "- Setup CUDAProb3" << std::endl;
  #endif
  
  _sampleFile->Close();
  char *histname = (char*)"blah";
  char *histtitle = (char*)"blahblah";

  std::cout << "-------------------------------------------------------------------" <<std::endl;

  //The binning here is arbitrary, now we get info from cfg so the 
  //set1DBinning and set2Dbinning calls below will make the binning
  //to be what we actually want
  _hPDF1D = new TH1D("hErec_nue", "Reconstructed Energy", 200, 0 , 50.0);
  dathist = new TH1D("dat_nue","",200,0, 50.0); 
  _hPDF2D = new TH2D(histname,histtitle,15,0,50.0*1000,15,0,150);
  dathist2d = new TH2D("dat2d_nue","",15,0,1500,15,0,150);

  //ETA Don't forget the -1 on the size here, as it's number of bins not bin edges
  set1DBinning(sample_erec_bins.size()-1, erec_bin_edges);
  set2DBinning(sample_erec_bins.size()-1, erec_bin_edges, sample_theta_bins.size()-1, theta_bin_edges); 

  return;
}


void samplePDFDUNEBase::SetupWeightPointers() {
  
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
	for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
	  //DB Setting total weight pointers
	  MCSamples[i].ntotal_weight_pointers[j] = 6;
	  MCSamples[i].total_weight_pointers[j] = new double*[MCSamples[i].ntotal_weight_pointers[j]];
	  MCSamples[i].total_weight_pointers[j][0] = &(dunemcSamples[i].pot_s);
	  MCSamples[i].total_weight_pointers[j][1] = &(dunemcSamples[i].norm_s);
	  MCSamples[i].total_weight_pointers[j][2] = &(MCSamples[i].osc_w[j]);
	  MCSamples[i].total_weight_pointers[j][3] = &(dunemcSamples[i].rw_berpaacvwgt[j]);
	  MCSamples[i].total_weight_pointers[j][4] = &(MCSamples[i].flux_w[j]);
	  MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
	}
  }

  return;
}


void samplePDFDUNEBase::setupDUNEMC(const char *sampleFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats)
{
  
  // set up splines
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;
  
  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree*)_sampleFile->Get("cafTree");
  // _data = (TTree*)_sampleFile->Get("caf");

  _isND = false;

  if(_data){
    std::cout << "Found mtuple tree is " << sampleFile << std::endl;
    std::cout << "N of entries: " << _data->GetEntries() << std::endl;
  }
  
  if (_isND) {
    _data->SetBranchStatus("*", 0);
    _data->SetBranchStatus("Ev", 1);
    _data->SetBranchAddress("Ev", &_ev);
    _data->SetBranchStatus("Ev_reco", 1);
    _data->SetBranchAddress("Ev_reco", &_erec);
    _data->SetBranchStatus("Ev_reco", 1);
    _data->SetBranchAddress("Ev_reco", &_erec_nue);
    _data->SetBranchStatus("mode",1);
    _data->SetBranchAddress("mode",&_mode);
    _data->SetBranchAddress("NuMomX", &_NuMomX);
    _data->SetBranchAddress("NuMomY", &_NuMomY);
    _data->SetBranchAddress("NuMomZ", &_NuMomZ);
    _data->SetBranchAddress("LepMomX", &_LepMomX);
    _data->SetBranchAddress("LepMomY", &_LepMomY);
    _data->SetBranchAddress("LepMomZ", &_LepMomZ);
    _data->SetBranchStatus("cvnnumu",1);
    _data->SetBranchAddress("cvnnumu", &_cvnnumu);
    _data->SetBranchStatus("cvnnue",1);
    _data->SetBranchAddress("cvnnue", &_cvnnue);
    _data->SetBranchStatus("isCC", 1);
    _data->SetBranchAddress("isCC", &_isCC);
    _data->SetBranchStatus("nuPDGunosc", 1);
    _data->SetBranchAddress("nuPDGunosc", &_nuPDGunosc);
    _data->SetBranchStatus("nuPDG", 1);
    _data->SetBranchAddress("nuPDG", &_nuPDGunosc);
    _data->SetBranchStatus("run", 1);
    _data->SetBranchAddress("run", &_run);
    _data->SetBranchStatus("isFHC", 1);
    _data->SetBranchAddress("isFHC", &_isFHC);
    _data->SetBranchStatus("BeRPA_A_cvwgt", 1);
    _data->SetBranchAddress("BeRPA_A_cvwgt", &_BeRPA_cvwgt);
    _data->SetBranchStatus("vtx_x", 1);
    _data->SetBranchAddress("vtx_x", &_vtx_x);
    _data->SetBranchStatus("vtx_y", 1);
    _data->SetBranchAddress("vtx_y", &_vtx_y);
    _data->SetBranchStatus("vtx_z", 1);
    _data->SetBranchAddress("vtx_z", &_vtx_z);  
    _data->SetBranchStatus("LepPDG",1); 
    _data->SetBranchAddress("LepPDG",&_ipnu); 
    _data->SetBranchStatus("LepNuAngle",1); 
    _data->SetBranchAddress("LepNuAngle",&_LepTheta); 
    _data->SetBranchStatus("Q2",1); 
    _data->SetBranchAddress("Q2",&_Q2); 
  }
  else {	
    _data->SetBranchStatus("*", 0);
    _data->SetBranchStatus("Ev", 1);
    _data->SetBranchAddress("Ev", &_ev);
    _data->SetBranchStatus("Ev_reco_numu", 1);
    _data->SetBranchAddress("Ev_reco_numu", &_erec);
    _data->SetBranchStatus("Ev_reco_nue", 1);
    _data->SetBranchAddress("Ev_reco_nue", &_erec_nue);
    _data->SetBranchStatus("mode",1);
    _data->SetBranchAddress("mode",&_mode);
    _data->SetBranchAddress("NuMomX", &_NuMomX);
    _data->SetBranchStatus("NuMomY", 1);
    _data->SetBranchAddress("NuMomY", &_NuMomY);
    _data->SetBranchAddress("NuMomZ", &_NuMomZ);
    _data->SetBranchStatus("LepMomX", 1);
    _data->SetBranchAddress("LepMomX", &_LepMomX);
    _data->SetBranchStatus("LepMomY", 1);
    _data->SetBranchAddress("LepMomY", &_LepMomY);
    _data->SetBranchStatus("LepMomZ", 1);
    _data->SetBranchAddress("LepMomZ", &_LepMomZ);
    _data->SetBranchStatus("cvnnumu",1);
    _data->SetBranchAddress("cvnnumu", &_cvnnumu);
    _data->SetBranchStatus("cvnnue",1);
    _data->SetBranchAddress("cvnnue", &_cvnnue);
    _data->SetBranchStatus("isCC", 1);
    _data->SetBranchAddress("isCC", &_isCC);
    _data->SetBranchStatus("nuPDGunosc", 1);
    _data->SetBranchAddress("nuPDGunosc", &_nuPDGunosc);
    _data->SetBranchStatus("nuPDG", 1);
    _data->SetBranchAddress("nuPDG", &_nuPDGunosc);
    _data->SetBranchStatus("run", 1);
    _data->SetBranchAddress("run", &_run);
    _data->SetBranchStatus("isFHC", 1);
    _data->SetBranchAddress("isFHC", &_isFHC);
    _data->SetBranchStatus("BeRPA_A_cvwgt", 1);
    _data->SetBranchAddress("BeRPA_A_cvwgt", &_BeRPA_cvwgt);
    _data->SetBranchStatus("vtx_x", 1);
    _data->SetBranchAddress("vtx_x", &_vtx_x);
    _data->SetBranchStatus("vtx_y", 1);
    _data->SetBranchAddress("vtx_y", &_vtx_y);
    _data->SetBranchStatus("vtx_z", 1);
    _data->SetBranchAddress("vtx_z", &_vtx_z);  
    _data->SetBranchStatus("LepPDG",1); 
    _data->SetBranchAddress("LepPDG",&_ipnu); 
    _data->SetBranchStatus("LepNuAngle",1); 
    _data->SetBranchAddress("LepNuAngle",&_LepTheta); 
    _data->SetBranchStatus("Q2",1); 
    _data->SetBranchAddress("Q2",&_Q2);
    _data->SetBranchStatus("weight",1); 
    _data->SetBranchAddress("weight",&_weight);
 }

  //FIX Commenting out for now
  TH1D* norm = (TH1D*)_sampleFile->Get("norm");
  if(!norm){
    std::cout<< "Add a norm KEY to the root file using MakeNormHists.cxx"<<std::endl;
    std::cout << "Ignoring for now" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  } 

  // now fill the actual variables
  // duneobj->norm_s = norm->GetBinContent(1);
  // duneobj->pot_s = pot/norm->GetBinContent(2);
  duneobj->norm_s = 1;
  duneobj->pot_s = 1;
  //duneobj->norm_s = 1.0;
  //duneobj->pot_s = 1.0;
  std::cout<< "pot_s = " << duneobj->pot_s << std::endl;
  std::cout<< "norm_s = " << duneobj->norm_s << std::endl;
  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;
 
  
  std::cout << "signal: " << duneobj->signal << std::endl;
  std::cout << "nevents: " << duneobj->nEvents << std::endl;

  // allocate memory for dunemc variables
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->rw_Q2 = new double[duneobj->nEvents];
  duneobj->beam_w = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->xsec_w = new double[duneobj->nEvents];
  duneobj->rw_cvnnumu = new double[duneobj->nEvents];
  duneobj->rw_cvnnue = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_run = new int[duneobj->nEvents];
  duneobj->rw_isFHC = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_truecz = new double[duneobj->nEvents];


  duneobj->energyscale_w = new double[duneobj->nEvents];
  duneobj->mode = new int[duneobj->nEvents];
  duneobj->rw_lower_erec_1d = new double[duneobj->nEvents]; //lower erec bound for bin
  duneobj->rw_upper_erec_1d = new double[duneobj->nEvents]; //upper erec bound for bin
  duneobj->rw_lower_erec_2d = new double[duneobj->nEvents]; //lower erec bound for bin
  duneobj->rw_upper_erec_2d = new double[duneobj->nEvents]; //upper erec bound for bin
  duneobj->rw_ipnu = new int*[duneobj->nEvents]; 

  //These spline bins get filled in fillSplineBins
  duneobj->enu_s_bin = new unsigned int[duneobj->nEvents];
  duneobj->erec_s_bin = new unsigned int[duneobj->nEvents];
  //duneobj->theta_s_bin = new unsigned int[duneobj->nEvents];
  duneobj->flux_bin = new int[duneobj->nEvents];
  duneobj->xsec_norms_bins = new std::list< int >[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];

  _data->GetEntry(0);
  
  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) // Loop through tree
    {
      _data->GetEntry(i);
      //      cout << "Event[" << i << "]: erec " << _erec << endl;
      if (iselike) {
        duneobj->rw_erec[i] = _erec_nue;}
        //std::cout << "NUE EREC Energy" << std::endl;
      else {duneobj->rw_erec[i] = _erec;} // in GeV
      duneobj->rw_etru[i] = _ev; // in GeV
      if ( i == 1 ) {std::cout << "Etrue = " << duneobj->rw_etru[i] << std::endl;}
      duneobj->rw_cvnnumu[i] = _cvnnumu; 
      duneobj->rw_cvnnue[i] = _cvnnue;
      // duneobj->rw_theta[i] = _LepNuAngle;
      duneobj->rw_theta[i] = _LepMomY/sqrt(pow(_LepMomZ, 2) + pow(_LepMomY, 2) + pow(_LepMomX, 2));
      duneobj->rw_Q2[i] = _Q2;
      duneobj->rw_isCC[i] = _isCC;
      duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
      duneobj->rw_nuPDG[i] = _nuPDG;
      duneobj->rw_run[i] = _run;
      duneobj->rw_isFHC[i] = _isFHC;
      duneobj->rw_berpaacvwgt[i] = _BeRPA_cvwgt;
      duneobj->rw_vtx_x[i] = _vtx_x;
      duneobj->rw_vtx_y[i] = _vtx_y;
      duneobj->rw_vtx_z[i] = _vtx_z;
      duneobj->rw_truecz[i] = _NuMomY/_ev;

	  //Assume everything is on Argon for now....
	  duneobj->Target[i] = 40;
 
      duneobj->beam_w[i] = 1.0;
      duneobj->xsec_w[i] = 1.0;

      // fill modes
      modes->Fill(_mode);
      //!!possible cc1pi exception might need to be 11
      int mode= TMath::Abs(_mode);       

      duneobj->mode[i]=SIMBMode_ToMaCh3Mode(mode, _isCC);
 
      duneobj->energyscale_w[i] = 1.0;
      
      duneobj->flux_w[i] = _weight;
    }
  
  //std::vector<double> etruVector(duneobj->rw_etru, duneobj->rw_etru + duneobj->nEvents);
  //duneobj->kNu = new cudaprob3::BeamCpuPropagator<double>(duneobj->nEvents, 1); 
  //duneobj->kNu = new cudaprob3::BeamCudaPropagatorSingle(0, duneobj->nEvents); 
  //duneobj->kNu->setEnergyList(etruVector);
  _sampleFile->Close();
  std::cout << "Sample set up OK" << std::endl;
  
}

double samplePDFDUNEBase::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent){

 KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter)); 
 double KinematicValue = -999;

 switch(KinPar){
   case kTrueNeutrinoEnergy:
	 KinematicValue = dunemcSamples[iSample].rw_etru[iEvent]; 
	 break;
   case kTrueXPos:
	 KinematicValue = dunemcSamples[iSample].rw_vtx_x[iEvent];
	   break;
   case kTrueYPos:
	 KinematicValue = dunemcSamples[iSample].rw_vtx_y[iEvent];
	   break;
   case kTrueZPos:
	 KinematicValue = dunemcSamples[iSample].rw_vtx_z[iEvent];
	   break;
   default:
	 std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type..." << std::endl;
	 throw;
 }

  return KinematicValue;
}

double samplePDFDUNEBase::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent){
  KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
  double KinematicValue = -999;

  switch(KinPar){
   case kTrueNeutrinoEnergy:
	 KinematicValue = dunemcSamples[iSample].rw_etru[iEvent]; 
	 break;
   case kTrueXPos:
	 KinematicValue = dunemcSamples[iSample].rw_vtx_x[iEvent];
	   break;
   case kTrueYPos:
	 KinematicValue = dunemcSamples[iSample].rw_vtx_y[iEvent];
	   break;
   case kTrueZPos:
	 KinematicValue = dunemcSamples[iSample].rw_vtx_z[iEvent];
	   break;
   default:
	 std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type..." << std::endl;
	 throw;
 }
  return KinematicValue;
}

void samplePDFDUNEBase::setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj) 
{

  fdobj->nEvents = duneobj->nEvents;
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->x_var = new double*[fdobj->nEvents];
  fdobj->y_var = new double*[fdobj->nEvents];
  // fdobj->enu_s_bin = new unsigned int[fdobj->nEvents];
  // fdobj->xvar_s_bin = new unsigned int[fdobj->nEvents];
  // fdobj->yvar_s_bin = new unsigned int[fdobj->nEvents];
  fdobj->rw_etru = new double*[fdobj->nEvents];
  fdobj->XBin = new int[fdobj->nEvents];
  fdobj->YBin = new int[fdobj->nEvents];    
  fdobj->NomXBin = new int[fdobj->nEvents];
  fdobj->NomYBin = new int[fdobj->nEvents];
  fdobj->XBin = new int [fdobj->nEvents];
  fdobj->YBin = new int [fdobj->nEvents];;	 
  fdobj->rw_lower_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_lower_lower_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_upper_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_upper_upper_xbinedge = new double [fdobj->nEvents];
  fdobj->mode = new int*[fdobj->nEvents];
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents]; 
  fdobj->xsec_spline_pointers = new const double**[fdobj->nEvents];
  fdobj->nxsec_norm_pointers = new int[fdobj->nEvents];
  fdobj->xsec_norm_pointers = new const double**[fdobj->nEvents];
  fdobj->xsec_norms_bins = new std::list< int >[fdobj->nEvents];
  fdobj->xsec_w = new double[fdobj->nEvents];
  fdobj->flux_w = new double[fdobj->nEvents];
  fdobj->osc_w = new double[fdobj->nEvents];
  fdobj->isNC = new bool[fdobj->nEvents];
  fdobj->rw_etru = new double*[fdobj->nEvents];
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents];
  fdobj->xsec_spline_pointers = new const double**[fdobj->nEvents];
  fdobj->ntotal_weight_pointers = new int[fdobj->nEvents];
  fdobj->total_weight_pointers = new double**[fdobj->nEvents];
  fdobj->Target = new int*[fdobj->nEvents];

  //DB Example Atmospheric Implementation
  fdobj->Unity = 1.;
  fdobj->osc_w_pointer = new const double*[fdobj->nEvents];
  fdobj->rw_truecz = new double[fdobj->nEvents];
  TRandom3* CZ = new TRandom3(1);

  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
	fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
	fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
	fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
	//ETA - set these to a dummy value to begin with, these get set in set1DBinning or set2DBinning
	fdobj->NomXBin[iEvent] = -1;
	fdobj->NomYBin[iEvent] = -1;
	fdobj->XBin[iEvent] = -1;
	fdobj->YBin[iEvent] = -1;	 
	fdobj->rw_lower_xbinedge[iEvent] = -1;
	fdobj->rw_lower_lower_xbinedge[iEvent] = -1;
	fdobj->rw_upper_xbinedge[iEvent] = -1;
	fdobj->rw_upper_upper_xbinedge[iEvent] = -1;
	fdobj->xsec_w[iEvent] = 1.0;
	fdobj->osc_w[iEvent] = 1.0;
	fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
	fdobj->flux_w[iEvent] = duneobj->flux_w[iEvent];
	fdobj->SampleDetID = SampleDetID;

	//DB Example Atmospheric Implementation
	fdobj->osc_w_pointer[iEvent] = &(fdobj->Unity);
	// fdobj->rw_truecz[iEvent] = CZ->Uniform(-1.,1.);
	fdobj->rw_truecz[iEvent] = duneobj->rw_truecz[iEvent];

	//ETA - this is where the variables that you want to bin your samples in are defined
	//If you want to bin in different variables this is where you put it for now
	switch(BinningOpt){
	  case 0:
	  case 1:
		//Just point to xvar to the address of the variable you want to bin in
		//This way we don't have to update both fdmc and skmc when we apply shifts
		//to variables we're binning in
		fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
		// fdobj->y_var[iEvent] = &(duneobj->dummy_y);//ETA - don't think we even need this as if we have a 1D sample we never need this, just not sure I like an unitialised variable in fdmc struct? 
    fdobj->y_var[iEvent] = &(duneobj->rw_truecz[iEvent]);
    break;
	  case 2:
		//Just point to xvar to the address of the variable you want to bin in
		//This way we don't have to update both fdmc and skmc when we apply shifts
		//to variables we're binning in
		fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
		fdobj->y_var[iEvent] = &(duneobj->rw_theta[iEvent]);
		break;
	  default:
		std::cout << "[ERROR:] " << __FILE__ << ":" << __LINE__ << " unrecognised binning option" << BinningOpt << std::endl;
		throw;
		break;
	}

  }

  return;
}

void samplePDFDUNEBase::setupSplines(fdmc_base *fdobj, const char *splineFile, int nutype) {

  int nevents = fdobj->nEvents;
  std::cout << "##################" << std::endl;
  std::cout << "Initialising splines from file: " << (splineFile) << std::endl;
  std::cout << "##################" << std::endl;

  switch (BinningOpt){
	case 0:
	case 1:
	  fdobj->splineFile = new splinesDUNE(XsecCov);
	  if (!(nutype==1 || nutype==-1 || nutype==2 || nutype==-2)){
		std::cerr << "problem setting up splines in erec" << std::endl;
	  }
	  break;
	case 2:
	  std::cout << "Creating splineDUNEBase" << std::endl;
	  fdobj->splineFile = new splinesDUNE(XsecCov);
    std::cout << "splineFile is -> " << fdobj->splineFile << std::endl;
	  if (!(nutype==1 || nutype==-1 || nutype==2 || nutype==-2)) {
		std::cerr << "problem setting up splines in erec" << std::endl;
	  } 
	  break;
    default:
	  break;
  }

  // ETA - Moved SetupSplineInfoArrays to be here
  // fdobj->splineFile->SetupSplineInfoArray(xsecCov);
  // fdobj->splineFile->SetSplineInfoArrays();

  return;
}

//This is currently here just for show. We'll implement functional parameters soon!
double samplePDFDUNEBase::CalcXsecWeightFunc(int iSample, int iEvent) 
{
  return 1.0;
}

std::vector<double> samplePDFDUNEBase::ReturnKinematicParameterBinning(std::string KinematicParameterStr) 
{
  std::cout << "ReturnKinematicVarBinning" << std::endl;
  std::vector<double> binningVector;
  return binningVector;
}

double samplePDFDUNEBase::getDiscVar(int iSample, int iEvent, int varindx) 
{
  std::cout << "getDiscVar" << std::endl;
  return 0.0;
}

double samplePDFDUNEBase::getCovLikelihood() 
{
  std::cout << "getCovLikelihood" << std::endl;
  return 0.0;
}

void samplePDFDUNEBase::printPosteriors()
{
  std::cout << "printPosteriors" << std::endl;
}


int samplePDFDUNEBase::getNMCSamples()
{
  // std::cout << "getNMCSamples" << std::endl;
  // return 0;
  return dunemcSamples.size();
}

int samplePDFDUNEBase::getNEventsInSample(int sample)
{
  // std::cout << "getNEventsInSample" << std::endl;
  // return 0;
  return dunemcSamples[sample].nEvents;
}
