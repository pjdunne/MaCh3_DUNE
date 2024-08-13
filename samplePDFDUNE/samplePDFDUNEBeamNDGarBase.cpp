#include <TROOT.h>

#include "samplePDFDUNEBeamNDGarBase.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "manager/manager.h"
#include "TLeaf.h"

samplePDFDUNEBeamNDGarBase::samplePDFDUNEBeamNDGarBase(double pot, std::string mc_version, covarianceXsec* xsec_cov) : samplePDFBase(pot) { 
  if(xsec_cov == NULL){std::cerr << "[ERROR:] You've passed me a NULL xsec covariance matrix... I need this to setup splines!" << std::endl; throw;}
  init(pot, mc_version, xsec_cov);          
}

samplePDFDUNEBeamNDGarBase::~samplePDFDUNEBeamNDGarBase() {
  delete sr;
  delete _sampleFile;
}

void samplePDFDUNEBeamNDGarBase::init(double pot, std::string samplecfgfile, covarianceXsec *xsec_cov) {
  char* sample_char = (char*)samplecfgfile.c_str();
  //ETA - trying to read stuff from yaml file
  manager* SampleManager = new manager(sample_char);

  //Bools
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();

  iscalo_reco = SampleManager->raw()["SampleBools"]["iscalo_reco"].as<bool>(); //NK determine what reco used
  muonscore_threshold = SampleManager->raw()["SampleCuts"]["muonscore_threshold"].as<float>(); //NK determine what muon score threshold to use

  //Inputs
  std::string mtupleprefix = SampleManager->raw()["InputFiles"]["mtupleprefix"].as<std::string>();
  std::string mtuplesuffix = SampleManager->raw()["InputFiles"]["mtuplesuffix"].as<std::string>();
  std::string splineprefix = SampleManager->raw()["InputFiles"]["splineprefix"].as<std::string>();
  std::string splinesuffix = SampleManager->raw()["InputFiles"]["splinesuffix"].as<std::string>();

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
  
  std::vector<double> tempselection(3);
  for(int iSelec=0; iSelec<NSelections; iSelec++){
    tempselection[0] = ReturnKinematicParameterFromString(SelectionStr[iSelec]);
    tempselection[1] = SelectionBounds[iSelec][0];
    tempselection[2] = SelectionBounds[iSelec][1];
    StoredSelection.push_back(tempselection);
  }
  
  //Make some arrays so we can initialise _hPDF1D and _hPDF2D with these
  double erec_bin_edges[sample_erec_bins.size()];
  double theta_bin_edges[sample_theta_bins.size()];
  for(unsigned erec_i = 0 ; erec_i < sample_erec_bins.size() ; erec_i++){erec_bin_edges[erec_i] = sample_erec_bins[erec_i];}
  for(unsigned theta_i = 0 ; theta_i < sample_theta_bins.size() ; theta_i++){theta_bin_edges[theta_i] = sample_theta_bins[theta_i];}
  
  // create dunendgarmc storage
  int nSamples = SampleManager->raw()["NSubSamples"].as<int>();
  for (int i=0;i<nSamples;i++) {
    struct dunemc_base obj = dunemc_base();
    dunendgarmcSamples.push_back(obj);
  }
  //Now down with yaml file for sample
  delete SampleManager;
  for(unsigned iSample=0 ; iSample < dunendgarmcSamples.size() ; iSample++){
    setupDUNEMC((mtupleprefix+mtuple_files[iSample]+mtuplesuffix).c_str(), &dunendgarmcSamples[sample_vecno[iSample]], pot, sample_nutype[iSample], sample_oscnutype[iSample], sample_signal[iSample]);
  }

  for (int i=0;i<nSamples;i++) {
    struct fdmc_base obj = fdmc_base();
    MCSamples.push_back(obj);
  }
  
  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    setupFDMC(&dunendgarmcSamples[sample_vecno[iSample]], &MCSamples[sample_vecno[iSample]], (splineprefix+spline_files[iSample]+splinesuffix).c_str());
  }
  
  std::cout << "################" << std::endl;
  std::cout << "Setup FD MC   " << std::endl;
  std::cout << "################" << std::endl;

  // ETA - If xsec_cov hasn't been passed to the samplePDFDUNEBaseND constructor then it's NULL
  // and the old funcitonality is kept
  // this calls this function in the core code
  // this needs to come after setupFDMC as otherwise MCSamples.splinefile will be NULL
  SetXsecCov(xsec_cov); 
  
  std::vector<std::string> spline_filepaths;
  
  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
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
  
  SetupNormParameters();
  SetupWeightPointers();
  
  fillSplineBins();
  
  _sampleFile->Close();
  
  std::cout << "-------------------------------------------------------------------" <<std::endl;
  
  //The binning here is arbitrary, now we get info from cfg so the 
  //set1DBinning and set2Dbinning calls below will make the binning
  //to be what we actually want
  _hPDF1D = new TH1D("hErec_nue", "Reconstructed Energy", 200, 0 , 50.0);
  dathist = new TH1D("dat_nue","",200,0, 50.0); 
  _hPDF2D = new TH2D("blah","blah",15,0,50.0*1000,15,0,150);
  dathist2d = new TH2D("dat2d_nue","",15,0,1500,15,0,150);

  //ETA Don't forget the -1 on the size here, as it's number of bins not bin edges
  set1DBinning(sample_erec_bins.size()-1, erec_bin_edges);
  set2DBinning(sample_erec_bins.size()-1, erec_bin_edges, sample_theta_bins.size()-1, theta_bin_edges); 
}

void samplePDFDUNEBeamNDGarBase::SetupWeightPointers() {
  for (int i = 0; i < (int)dunendgarmcSamples.size(); ++i) {
    for (int j = 0; j < dunendgarmcSamples[i].nEvents; ++j) {
      //DB Setting total weight pointers
      MCSamples[i].ntotal_weight_pointers[j] = 6;
      MCSamples[i].total_weight_pointers[j] = new double*[MCSamples[i].ntotal_weight_pointers[j]];
      MCSamples[i].total_weight_pointers[j][0] = &(dunendgarmcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunendgarmcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = &(MCSamples[i].osc_w[j]);
      MCSamples[i].total_weight_pointers[j][3] = &(dunendgarmcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(MCSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
    }
  }
}

void samplePDFDUNEBeamNDGarBase::setupDUNEMC(const char *sampleFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats) {
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;
  
  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree*)_sampleFile->Get("cafTree");
  if(_data){
    std::cout << "Found mtuple tree is " << sampleFile << std::endl;
    std::cout << "N of entries: " << _data->GetEntries() << std::endl;
  }
  
  _data->SetBranchStatus("*", 1);
  _data->SetBranchAddress("rec", &sr);

  duneobj->norm_s = 1.0;
  duneobj->pot_s = (pot)/1e21;

  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;
 
  std::cout << "signal: " << duneobj->signal << std::endl;
  std::cout << "nevents: " << duneobj->nEvents << std::endl;

  // allocate memory for dunendgarmc variables
  duneobj->rw_yrec = new double[duneobj->nEvents];
  duneobj->rw_elep_reco = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 

  duneobj->mode = new int[duneobj->nEvents];

  duneobj->nproton = new int[duneobj->nEvents];
  duneobj->nneutron = new int[duneobj->nEvents];
  duneobj->npip = new int[duneobj->nEvents];
  duneobj->npim = new int[duneobj->nEvents];
  duneobj->npi0 = new int[duneobj->nEvents];

  duneobj->nrecomuon = new int[duneobj->nEvents];
  duneobj->ntruemuon = new int[duneobj->nEvents];
  duneobj->nmuonsratio = new double[duneobj->nEvents];
  duneobj->ntruemuonprim = new int[duneobj->nEvents];

  duneobj->nrecoparticles = new int[duneobj->nEvents];
  duneobj->in_fdv = new bool[duneobj->nEvents];
  duneobj->rw_elep_true = new double[duneobj->nEvents];

  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_rad = new double[duneobj->nEvents];

  duneobj->rw_lep_pT = new double[duneobj->nEvents];
  duneobj->rw_lep_pZ = new double[duneobj->nEvents];

  duneobj->rw_reco_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_reco_rad = new double[duneobj->nEvents];

  duneobj->Target = new int[duneobj->nEvents];

  int num_no_ixns =0;
  int num_no_recparticles = 0;
  int num_in_fdv = 0;
  int num_in_fdv_noreco = 0;
  int num_notin_fdv =0;
  int num_nanenergy =0;
  int num_nanparticles =0;

  //FILL DUNE STRUCT
  for (int i = 0; i < (duneobj->nEvents); ++i) { // Loop through tree
    _data->GetEntry(i);
    double radius = pow((pow((sr->mc.nu[0].vtx.y+150),2) + pow((sr->mc.nu[0].vtx.z-1486),2)),0.5);
    if(std::abs(sr->mc.nu[0].vtx.x)<=209.0 &&  radius<=227.02){
      //std::cout<<"this event is within the fiducial volume"<<std::endl;
      num_in_fdv++;
      duneobj->in_fdv[i] = 1;
    }
    else{
      //std::cout<<"this event is NOT within the fiducial volume"<<std::endl;
      num_notin_fdv++;
      duneobj->in_fdv[i] = 0;
    }
    
    if(sr->common.ixn.ngsft == 0){
      //duneobj->rw_erec[i] = (double)(0);
      float erec_total =0;
      duneobj->rw_elep_reco[i] = double(0);
      duneobj->rw_yrec[i] = (double)(0);
      num_no_ixns++;
      duneobj->nrecoparticles[i] = (int)(0);
    }
    else{
      duneobj->nrecoparticles[i] = (int)(0);
      float erec_total =0;
      float elep_reco =0;
      float muonscore = muonscore_threshold;
      int nixns = (int)(sr->common.ixn.ngsft);
      for(int i_ixn =0; i_ixn<nixns; i_ixn++){
	int nrecoparticles = (int)(sr->common.ixn.gsft[i_ixn].part.ngsft);
	duneobj->nrecoparticles[i] += (int)(sr->common.ixn.gsft[i_ixn].part.ngsft);
	int nanparticles = 0;
	if(nrecoparticles ==0){
	  double radius = pow((pow((sr->mc.nu[0].vtx.y+150.),2) + pow((sr->mc.nu[0].vtx.z-1486.),2)),0.5);
	  if(std::abs(sr->mc.nu[0].vtx.x)<=209.0 || radius<=227.02){
	    //std::cout<<"within fd"<<std::endl;
	    num_in_fdv_noreco++;
	  }   
	  num_no_recparticles++;}
	for(int i_part =0; i_part<nrecoparticles; i_part++){
	  float erec_part = (float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].E);
	  //std::cout<<"erec_part: "<<erec_part<<std::endl;
	  if(std::isnan(erec_part)){nanparticles++;}
	  erec_total+=erec_part;
	  if((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].score.gsft_pid.muon_score>muonscore)){
	    if(erec_part>elep_reco){
	      elep_reco = erec_part;
	      duneobj->rw_reco_vtx_x[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.x));
	      duneobj->rw_reco_vtx_y[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.y));
	      duneobj->rw_reco_vtx_z[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.z));
	      duneobj->rw_lep_pT[i] = (double)(pow(pow((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.x), 2) + pow((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.y), 2), 0.5));
	      duneobj->rw_lep_pZ[i] = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.z);
	    }
	    duneobj->nrecomuon[i]++; 
	  }
	}
	num_nanparticles = num_nanparticles + (nanparticles/nrecoparticles);
      } //ADD PRIMARY LEPTON ENERGY ELEP_RECO
      if(std::isnan(erec_total)){std::cout<<"nan energy"<<std::endl; num_nanenergy++; erec_total = (float)(sr->common.ixn.gsft[0].Enu.lep_calo);}
      if(iscalo_reco){duneobj->rw_erec[i]=(double)(sr->common.ixn.gsft[0].Enu.lep_calo);  /*std::cout<<"calo erec: "<<(double)(sr->common.ixn.gsft[0].Enu.lep_calo)<<std::endl;*/}
      else{duneobj->rw_erec[i]=(double)(erec_total);}
      duneobj->rw_elep_reco[i] = (double)(elep_reco);
    }
    
    if(duneobj->rw_erec[i] != 0){duneobj->rw_yrec[i] = (double)(((duneobj->rw_erec[i])-(duneobj->rw_elep_reco[i]))/(duneobj->rw_erec[i]));}
    else{duneobj->rw_yrec[i] = (double)(0);}
    duneobj->rw_etru[i] = (double)(sr->mc.nu[0].E); // in GeV
    duneobj->rw_isCC[i] = (int)(sr->mc.nu[0].iscc);
    duneobj->rw_nuPDGunosc[i] = sr->mc.nu[0].pdgorig;
    duneobj->rw_nuPDG[i] = sr->mc.nu[0].pdg;
    duneobj->rw_berpaacvwgt[i] = _BeRPA_cvwgt;
    
    int ntrueparticles = (int)(sr->mc.nu[0].nprim);
    for(int i_truepart =0; i_truepart<ntrueparticles; i_truepart++){
      if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 13){
	duneobj->ntruemuon[i]++;
	duneobj->ntruemuonprim[i]++;
      }
    }
    int ntruesecparticles = (int)(sr->mc.nu[0].nsec);
    for(int i_truepart =0; i_truepart<ntruesecparticles; i_truepart++){
      if(std::abs(sr->mc.nu[0].sec[i_truepart].pdg) == 13){
	duneobj->ntruemuon[i]++;
      }
    }
    
    duneobj->nproton[i] = sr->mc.nu[0].nproton;
    duneobj->nneutron[i] = sr->mc.nu[0].nneutron;
    duneobj->npip[i] = sr->mc.nu[0].npip;
    duneobj->npim[i] = sr->mc.nu[0].npim;
    duneobj->npi0[i] = sr->mc.nu[0].npi0;
    
    duneobj->nmuonsratio[i] = (double)(duneobj->nrecomuon[i])/(double)(duneobj->ntruemuonprim[i]);
    duneobj->rw_vtx_x[i] = (double)(sr->mc.nu[0].vtx.x);
    duneobj->rw_vtx_y[i] = (double)(sr->mc.nu[0].vtx.y);
    duneobj->rw_vtx_z[i] = (double)(sr->mc.nu[0].vtx.z);
    
    duneobj->rw_rad[i] = (double)(pow((pow((duneobj->rw_vtx_y[i]+150),2) + pow((duneobj->rw_vtx_z[i]-1486),2)),0.5)); 
    duneobj->rw_reco_rad[i] = (double)(pow(pow((duneobj->rw_reco_vtx_y[i]+150),2) + pow((duneobj->rw_reco_vtx_z[i]-1486), 2), 0.5));
    duneobj->rw_elep_true[i] = (double)(sr->mc.nu[0].prim[0].p.E);
    
    //Assume everything is on Argon for now....
    duneobj->Target[i] = 40;
    
    _mode = sr->mc.nu[0].mode;
    _isCC = (int)(sr->mc.nu[0].iscc);
    
    int mode= TMath::Abs(_mode);       
    duneobj->mode[i]=GENIEMode_ToMaCh3Mode(mode, _isCC);
    
    duneobj->flux_w[i] = 1.0;
  }
  
  std::cout<<"num evts in fdv: "<<num_in_fdv<<std::endl;
  std::cout<<"num evts not in fdv: "<<num_notin_fdv<<std::endl;
  std::cout<<"num no reco particles in fdv"<< num_in_fdv_noreco<<std::endl;
  std::cout<<"num no ixns: "<<num_no_ixns<<std::endl;
  std::cout<<"num no rec particles: "<<num_no_recparticles<<std::endl;
  std::cout<<"num nan energy: "<<num_nanenergy<<std::endl;
  std::cout<<"num nan particles: "<<num_nanparticles<<std::endl;
  std::cout << "Sample set up OK" << std::endl;
}

double samplePDFDUNEBeamNDGarBase::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
 KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
 return ReturnKinematicParameter(KinPar, iSample, iEvent);
}

double samplePDFDUNEBeamNDGarBase::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
  return ReturnKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEBeamNDGarBase::ReturnKinematicParameter(KinematicTypes KinematicParameter, int iSample, int iEvent) {
 double KinematicValue = -999;
 
 switch(KinematicParameter) {
 case kTrueNeutrinoEnergy:
   KinematicValue = dunendgarmcSamples[iSample].rw_etru[iEvent]; 
   break;
 case kPionMultiplicity:
   KinematicValue = dunendgarmcSamples[iSample].npip[iEvent]+dunendgarmcSamples[iSample].npim[iEvent+dunendgarmcSamples[iSample].npi0[iEvent]];
   break;
 case kRecoNeutrinoEnergy:
   KinematicValue = dunendgarmcSamples[iSample].rw_erec[iEvent];
   break;
 case kTrueXPos:
   KinematicValue = dunendgarmcSamples[iSample].rw_vtx_x[iEvent];
   break;
 case kTrueYPos:
   KinematicValue = dunendgarmcSamples[iSample].rw_vtx_y[iEvent];
   break;
 case kTrueZPos:
   KinematicValue = dunendgarmcSamples[iSample].rw_vtx_z[iEvent];
   break;
 case kTrueRad:
   KinematicValue = dunendgarmcSamples[iSample].rw_rad[iEvent];
   break;
 case kNRecoParticles:
   KinematicValue = dunendgarmcSamples[iSample].nrecoparticles[iEvent];
   break; 
 case kInFDV:
   KinematicValue = dunendgarmcSamples[iSample].in_fdv[iEvent];
   break;
 case kTrueMinusRecoEnergyRatio:
   KinematicValue = (dunendgarmcSamples[iSample].rw_etru[iEvent]-dunendgarmcSamples[iSample].rw_erec[iEvent])/dunendgarmcSamples[iSample].rw_etru[iEvent];
   break;
 case kTrueMinusRecoEnergy:
   KinematicValue = (dunendgarmcSamples[iSample].rw_etru[iEvent]-dunendgarmcSamples[iSample].rw_erec[iEvent]);
   break;
 case kNRecoMuons:
   KinematicValue = dunendgarmcSamples[iSample].nrecomuon[iEvent];
   break;
 case kNTrueMuons:
   KinematicValue = dunendgarmcSamples[iSample].ntruemuon[iEvent];
   break;
 case kNMuonsRecoOverTruth:
   KinematicValue = dunendgarmcSamples[iSample].nmuonsratio[iEvent];
   break;
 case kRecoLepEnergy:
   KinematicValue = dunendgarmcSamples[iSample].rw_elep_reco[iEvent];
   break;
 case kTrueLepEnergy:
   KinematicValue = dunendgarmcSamples[iSample].rw_elep_true[iEvent];
   break;
 case kRecoXPos:
   KinematicValue = dunendgarmcSamples[iSample].rw_reco_vtx_x[iEvent];
   break;
 case kRecoYPos:
   KinematicValue = dunendgarmcSamples[iSample].rw_reco_vtx_y[iEvent];
   break;
 case kRecoZPos:
   KinematicValue = dunendgarmcSamples[iSample].rw_reco_vtx_z[iEvent];
   break;
 case kRecoRad:
   KinematicValue = dunendgarmcSamples[iSample].rw_reco_rad[iEvent];
   break;
 case kLepPT:
   KinematicValue = dunendgarmcSamples[iSample].rw_lep_pT[iEvent];
   break;
 case kLepPZ:
   KinematicValue = dunendgarmcSamples[iSample].rw_lep_pZ[iEvent];
   break;
 case kM3Mode:
   KinematicValue = dunendgarmcSamples[iSample].mode[iEvent];
   break;
 default:
   std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type..." << std::endl;
   throw;
 }
 
 return KinematicValue;
}

void samplePDFDUNEBeamNDGarBase::setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj, const char *splineFile) {
  fdobj->nEvents = duneobj->nEvents;
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->x_var = new double*[fdobj->nEvents];
  fdobj->y_var = new double*[fdobj->nEvents];
  fdobj->rw_etru = new double*[fdobj->nEvents];
  fdobj->XBin = new int[fdobj->nEvents];
  fdobj->YBin = new int[fdobj->nEvents];    
  fdobj->NomXBin = new int[fdobj->nEvents];
  fdobj->NomYBin = new int[fdobj->nEvents];
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
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents];
  fdobj->xsec_spline_pointers = new const double**[fdobj->nEvents];
  fdobj->ntotal_weight_pointers = new int[fdobj->nEvents];
  fdobj->total_weight_pointers = new double**[fdobj->nEvents];
  fdobj->Target = new int*[fdobj->nEvents];
  
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
    
    //ETA - this is where the variables that you want to bin your samples in are defined
    //If you want to bin in different variables this is where you put it for now
    switch(BinningOpt){
    case 0:
    case 1:
      //Just point to xvar to the address of the variable you want to bin in
      //This way we don't have to update both fdmc and skmc when we apply shifts
      //to variables we're binning in
      fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
      fdobj->y_var[iEvent] = &(duneobj->dummy_y);//ETA - don't think we even need this as if we have a 1D sample we never need this, just not sure I like an unitialised variable in fdmc struct? 
      break;
    case 2:
      //Just point to xvar to the address of the variable you want to bin in
      //This way we don't have to update both fdmc and skmc when we apply shifts
      //to variables we're binning in
      fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
      fdobj->y_var[iEvent] = &(duneobj->rw_yrec[iEvent]);
      break;
    default:
      std::cout << "[ERROR:] " << __FILE__ << ":" << __LINE__ << " unrecognised binning option" << BinningOpt << std::endl;
      throw;
      break;
    }
  }
  
}

std::vector<double> samplePDFDUNEBeamNDGarBase::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameterStr));
  return ReturnKinematicParameterBinning(KinPar);
}

std::vector<double> samplePDFDUNEBeamNDGarBase::ReturnKinematicParameterBinning(KinematicTypes KinPar) {
  std::vector<double> binningVector;
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    for(double ibins =0; ibins<10*10; ibins++){
      double binval = ibins/10;
      binningVector.push_back(binval);
    }
    break;
  case kRecoNeutrinoEnergy:
    for(double ibins =0; ibins<10*10; ibins++){
      double binval = ibins/10;
      binningVector.push_back(binval);
    } 
    break;
  case kRecoXPos:
  case kTrueXPos:
    for(double ibins =0; ibins<259*2; ibins++){
      binningVector.push_back(ibins-259);
    }
    break;
  case kRecoYPos:
  case kTrueYPos:
    for(double ibins =0; ibins<277*2; ibins++){
      binningVector.push_back(ibins-277-150);
    }
    break;
  case kRecoZPos:
  case kTrueZPos:
    for(double ibins =0; ibins<277*2; ibins++){
      binningVector.push_back(ibins-277+1486);
    }
    break;
  case kPionMultiplicity:
    for(double ibins =0; ibins<10; ibins++){
      binningVector.push_back(ibins);
    }
    break;
  case kNRecoParticles:
    for(double ibins =0; ibins<50; ibins++){
      binningVector.push_back(ibins);
    }
    break; 
  case kInFDV:
    for(double ibins =0; ibins<3; ibins++){
      binningVector.push_back(ibins);
    }
    break;
  case kNMuonsRecoOverTruth:
  case kTrueMinusRecoEnergyRatio:
    for(double ibins =0; ibins<20*10; ibins++){
      binningVector.push_back(-10+(double)(ibins)/10);
    }
    break;
  case kTrueMinusRecoEnergy:
    for(double ibins =0; ibins<20*10; ibins++){
      binningVector.push_back(-10+(double)(ibins)/10);
    }
    break;
  case kNTrueMuons:
  case kNRecoMuons:
    for(double ibins =0; ibins<10; ibins++){
      binningVector.push_back(ibins);
    }
    break;
  case kRecoLepEnergy:
    for(double ibins =0; ibins<10*10; ibins++){
      binningVector.push_back((double)(ibins)/10);
    } 
    break;
  case kTrueLepEnergy:
    for(double ibins =0; ibins<10*10; ibins++){
      binningVector.push_back((double)(ibins)/10);
    } 
    break;
  case kTrueRad:
  case kRecoRad:
    for(double ibins =0; ibins<300; ibins++){
      binningVector.push_back(ibins);
    }
    break;
  case kLepPT:
  case kLepPZ:
    for(double ibins =0; ibins<10*10; ibins++){
      binningVector.push_back((double)(ibins)/10);
    }
    break;
  default:
    for(double ibins =0; ibins<10*100; ibins++){
      binningVector.push_back(ibins/100);
    }
    break;
  }
  
  return binningVector;
}
