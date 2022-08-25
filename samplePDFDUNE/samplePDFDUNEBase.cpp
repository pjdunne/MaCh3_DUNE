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
  //flux = NULL;
  //skdet_joint = NULL;
  // ETA - If xsec_cov hasn't been passed to the samplePDFDUNEBase constructor then it's NULL
  // and the old funcitonality is kepy
  setXsecCov(xsec_cov); 
  // Initialise xsec model variable, gets set properly when xsec cov is set
  // What cross-section model are we using? (only applies to 2015 onwards)
  // 2015a = 0, 2015b = 1, 2015c = 2, 2016a = 3, see Structs.h
  // Undefined -> old 2012 model no longer supported

  //doubled_angle =true ;
  useNonDoubledAngles(true);
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

  //ETA - try initialising a manager with the samplecfg
  //MANAGER NOT READY YET
  /*
  manager* sample_manager = new manager(sample_char, false);
  std::cout << "isRHC is " << sample_manager->getIsRHC() << std::endl;
  IsRHC = sample_manager->getIsRHC();
  iselike = sample_manager->getiselike();
  up_bnd = sample_manager->getUpBnd();
  mtupleprefix = (std::string)sample_manager->getMtuplePrefix();
  mtuplesuffix = (std::string)sample_manager->getMtupleSuffix();
  splineprefix = (std::string)sample_manager->getSplinePrefix();
  splinesuffix = (std::string)sample_manager->getSplineSuffix();
  std::vector<std::string> mtuple_files = sample_manager->getMtupleFiles();
  std::vector<std::string> spline_files = sample_manager->getSplineFiles();
  std::vector<int> sample_vecno = sample_manager->getSamplevecno();
  std::vector<int> sample_oscnutype = sample_manager->getOscNuType();
  std::vector<int> sample_nutype = sample_manager->getNuType();
  std::vector<bool> sample_signal = sample_manager->getSignal(); */


  manager* SampleManager = new manager(sample_char);

  //Bools
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  //SampleDetID = SampleManager->raw()[""][""];
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  //Cuts
  up_bnd = SampleManager->raw()["SampleCuts"]["up_bnd"].as<float>();

  //Inputs
  mtupleprefix = SampleManager->raw()["InputFiles"]["mtupleprefix"].as<std::string>();
  mtuplesuffix = SampleManager->raw()["InputFiles"]["mtuplesuffix"].as<std::string>();
  splineprefix = SampleManager->raw()["InputFiles"]["splineprefix"].as<std::string>();
  splinesuffix = SampleManager->raw()["InputFiles"]["splinesuffix"].as<std::string>();



  //IsRHC = false;
  SampleDetID = 25;
  //iselike = false;
  //up_bnd = 30.;
  //mtupleprefix = "inputs/DUNE_CAF_files/FD_FHC_ger_";
  //mtuplesuffix = "_numuselec.root";
  //splineprefix = "inputs/DUNE_spline_files/FD_FHC_ger_";
  //splinesuffix = "_numuselec_splines.root";;
  std::vector<std::string> mtuple_files; mtuple_files.push_back("numu_x_numu");
  std::vector<std::string> spline_files; spline_files.push_back("numu_x_numu");
  std::vector<int> sample_vecno; sample_vecno.push_back(0);
  std::vector<int> sample_oscnutype; sample_oscnutype.push_back(2);
  std::vector<int> sample_nutype; sample_nutype.push_back(2);
  std::vector<bool> sample_signal; sample_signal.push_back(false);

  //ETA - only have support for erec and erec-theta binning options.
  //if you want to add in the option to do p-theta binning and enu-Q2
  //then just add in the binning option to the manager class and
  //add a switch statement below
  std::vector<double> sample_erec_bins {0.,  0.5,  1.,  1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4.,  5.,   6.,  10.};
  std::vector<double> sample_theta_bins {0.,  0.5,  1.,  1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4.,  5.,   6.,  10.};
  /* if(iselike){
    sample_erec_bins = sample_manager->getNueErecBins();
    sample_theta_bins = sample_manager->getNueThetaBins();
  }
  else{
	sample_erec_bins = sample_manager->getNumuErecBins();
    sample_theta_bins = sample_manager->getNumuThetaBins();
  } */

  //Make some arrays so we can initialise _hPDF1D and _hPDF2D with these
  double erec_bin_edges[sample_erec_bins.size()];
  double theta_bin_edges[sample_theta_bins.size()];
  for(unsigned erec_i = 0 ; erec_i < sample_erec_bins.size() ; erec_i++){erec_bin_edges[erec_i] = sample_erec_bins[erec_i];}
  for(unsigned theta_i = 0 ; theta_i < sample_theta_bins.size() ; theta_i++){theta_bin_edges[theta_i] = sample_theta_bins[theta_i];}


  // create dunemc storage
  int nSamples = 1;
  for (int i=0;i<nSamples;i++) {
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
  }


  //delete sample_manager;

  cout << "Oscnutype size: " << sample_oscnutype.size() << ", dunemcSamples size: " << dunemcSamples.size() << endl;  
  if(sample_oscnutype.size() != dunemcSamples.size()){std::cerr << "[ERROR:] samplePDFDUNEBase::samplePDFDUNEBase() - something went wrong either getting information from sample config" << std::endl; throw;}

  for(unsigned iSample=0 ; iSample < dunemcSamples.size() ; iSample++){
    setupDUNEMC((mtupleprefix+mtuple_files[iSample]+mtuplesuffix).c_str(), &dunemcSamples[sample_vecno[iSample]], pot, sample_nutype[iSample], sample_oscnutype[iSample], sample_signal[iSample]);
  }

  for (int i=0;i<nSamples;i++) {
    struct fdmc_base obj = fdmc_base();
    MCSamples.push_back(obj);
  }

  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    setupFDMC(&dunemcSamples[sample_vecno[iSample]], &MCSamples[sample_vecno[iSample]], (splineprefix+spline_files[iSample]+splinesuffix).c_str());
  }


  fillSplineBins();
  
  _sampleFile->Close();
  char *histname = (char*)"blah";
  char *histtitle = (char*)"blahblah";

  std::cout << "-------------------------------------------------------------------" <<std::endl;

  //The binning here is arbitrary, now we get info from cfg so the 
  //set1DBinning and set2Dbinning calls below will make the binning
  //to be what we actually want
  _hPDF1D = new TH1D("hErec_nue", "Reconstructed Energy", 200, 0 , up_bnd);
  dathist = new TH1D("dat_nue","",200,0,up_bnd); 
  _hPDF2D = new TH2D(histname,histtitle,15,0,up_bnd*1000,15,0,150);
  dathist2d = new TH2D("dat2d_nue","",15,0,1500,15,0,150);

  //ETA Don't forget the -1 on the size here, as it's number of bins not bin edges
  set1DBinning(sample_erec_bins.size()-1, erec_bin_edges);
  set2DBinning(sample_erec_bins.size()-1, erec_bin_edges, sample_theta_bins.size()-1, theta_bin_edges); 

}

void samplePDFDUNEBase::setupDUNEMC(const char *sampleFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats)
{
  
  // set up splines
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;

  
  std::cout << "Setup splines! " << std::endl;
  
  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree*)_sampleFile->Get("caf");

  if(_data){
    std::cout << "Found mtuple tree is " << sampleFile << std::endl;
    std::cout << "N of entries: " << _data->GetEntries() << std::endl;
  }
  


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
  _data->SetBranchStatus("isFD", 1);
  _data->SetBranchAddress("isFD", &_isFD);
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
  _data->SetBranchAddress("LepPDG",&_LepTheta); 
  _data->SetBranchStatus("Q2",1); 
  _data->SetBranchAddress("Q2",&_Q2); 


  TH1D* norm = (TH1D*)_sampleFile->Get("norm");
  if(!norm){
    std::cout<< "Add a norm KEY to the root file using MakeNormHists.cxx"<<std::endl;
    std::cout << "Ignoring for now" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  } 


  // now fill the actual variables
  duneobj->norm_s = norm->GetBinContent(1);
  duneobj->pot_s = pot/norm->GetBinContent(2);
  //duneobj->norm_s = 0.08369969415;
  //duneobj->pot_s = 1.0;
  std::cout<< "pot_s = " << duneobj->pot_s << std::endl;
  std::cout<< "norm_s = " << duneobj->norm_s << std::endl;
  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;
  
  std::cout << "nevents: " << duneobj->nEvents << std::endl;


  // allocate memory for dunemc variables
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->osc_w = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->rw_Q2 = new double[duneobj->nEvents];
  duneobj->beam_w = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->skdet_w = new double[duneobj->nEvents];
  duneobj->xsec_w = new double[duneobj->nEvents];
  duneobj->rw_cvnnumu = new double[duneobj->nEvents];
  duneobj->rw_cvnnue = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_run = new int[duneobj->nEvents];
  duneobj->rw_isFHC = new int[duneobj->nEvents];
  duneobj->rw_isFD = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];



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
  duneobj->Eb_bin = new int[duneobj->nEvents];

  _data->GetEntry(0);
  
  double dir_beam[3]     = {MaCh3Utils::SKNuDir[0], MaCh3Utils::SKNuDir[1], MaCh3Utils::SKNuDir[2]};// T2K beam direction at SK
  

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
      duneobj->rw_theta[i] = _LepNuAngle;
      duneobj->rw_Q2[i] = _Q2;
      duneobj->rw_isCC[i] = _isCC;
      duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
      duneobj->rw_nuPDG[i] = _nuPDG;
      duneobj->rw_run[i] = _run;
      duneobj->rw_isFHC[i] = _isFHC;
      duneobj->rw_isFD[i] = _isFD;
      duneobj->rw_berpaacvwgt[i] = _BeRPA_cvwgt;
      duneobj->rw_vtx_x[i] = _vtx_x;
      duneobj->rw_vtx_y[i] = _vtx_y;
      duneobj->rw_vtx_z[i] = _vtx_z;

      
      duneobj->osc_w[i] = 1.0; // not oscillated
      duneobj->beam_w[i] = 1.0;
      duneobj->xsec_w[i] = 1.0;
      duneobj->skdet_w[i] = 1.0;


      // fill modes
      modes->Fill(_mode);
      //!!possible cc1pi exception might need to be 11
      int mode= TMath::Abs(_mode);

      //Adjust SIMB modes into CC and NC
       

      //duneobj->mode[i]=SIMBMode_ToMaCh3Mode(mode, _isCC);
      duneobj->mode[i]= mode;
 
      duneobj->energyscale_w[i] = 1.0;
      
      duneobj->flux_w[i] = 1.0;

      //ETA - Now Find the correct global bin number which then gets used in getDiscVar	  
      //doing this here saves repeating this many times when running a chain
      duneobj->Eb_bin[i] = -999;
 
    }
  




  //std::vector<double> etruVector(duneobj->rw_etru, duneobj->rw_etru + duneobj->nEvents);
  //duneobj->kNu = new cudaprob3::BeamCpuPropagator<double>(duneobj->nEvents, 1); 
  //duneobj->kNu = new cudaprob3::BeamCudaPropagatorSingle(0, duneobj->nEvents); 
  //duneobj->kNu->setEnergyList(etruVector);
  _sampleFile->Close();
  std::cout << "Sample set up OK" << std::endl;
  
}

void samplePDFDUNEBase::setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj, const char *splineFile) 
{
  fdobj->nEvents = duneobj->nEvents;
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;

  fdobj->x_var = new double*[fdobj->nEvents];
  fdobj->y_var = new double*[fdobj->nEvents];
  fdobj->enu_s_bin = new unsigned int[fdobj->nEvents];
  fdobj->xvar_s_bin = new unsigned int[fdobj->nEvents];
  fdobj->yvar_s_bin = new unsigned int[fdobj->nEvents];
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

  switch(BinningOpt){
	case 0:
	case 1:
	  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
		//Just point to xvar to the address of the variable you want to bin in
		//This way we don't have to update both fdmc and skmc when we apply shifts
		//to variables we're binning in
		fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
		fdobj->y_var[iEvent] = &(duneobj->dummy_y);//ETA - don't think we even need this as if we have a 1D sample we never need this, just not sure I like an unitialised variable in fdmc struct? 
		fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
		fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
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
	  }
	  break;
	case 2:
	  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
		//Just point to xvar to the address of the variable you want to bin in
		//This way we don't have to update both fdmc and skmc when we apply shifts
		//to variables we're binning in
		fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
		fdobj->y_var[iEvent] = &(duneobj->rw_theta[iEvent]);
		fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
		fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
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
		//std::cout << "Set MCSamples flux_w to " << skmcSamples[iSample].flux_w[iEvent] << std::endl;
	  }
	  break;
	default:
	  std::cout << "[ERROR:] ___LINE___ unrecognised binning option" << std::endl;
	  break;
  }
  

  std::cout << "Setting up Splines" << std::endl;
  setupSplines(fdobj, splineFile, fdobj->nutype, fdobj->signal);
  
  std::cout << "~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "Done is setupFDMC" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~" << std::endl;

  return;
  }
  


double samplePDFDUNEBase::calcFuncSystWeight(int iSample, int iEvent) 
{
  std::cout << "calcFuncSystWeight" << std::endl;
  return 0.0;
}


double samplePDFDUNEBase::ReturnKinematicParameter(KinematicTypes Var, int i, int j) 
{
  std::cout << "ReturnKinematicVar" << std::endl;
  return 0.0;
}

std::vector<double> samplePDFDUNEBase::ReturnKinematicParameterBinning(KinematicTypes Var) 
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
  std::cout << "getNMCSamples" << std::endl;
  return 0;
}

int samplePDFDUNEBase::getNEventsInSample(int sample)
{
  std::cout << "getNEventsInSample" << std::endl;
  return 0;
}
