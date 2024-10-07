#include <TROOT.h>
#include "samplePDFDUNEAtmBase.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "manager/manager.h"

samplePDFDUNEAtmBase::samplePDFDUNEAtmBase(double pot_, std::string mc_version_, covarianceXsec* xsec_cov_) : samplePDFFDBase(pot_, mc_version_, xsec_cov_) {
  // create dunemc storage
  for (int i=0;i<nSamples;i++) {
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
  }

  Initialise();
}

samplePDFDUNEAtmBase::~samplePDFDUNEAtmBase() {
}

void samplePDFDUNEAtmBase::Init() {
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();
  
  std::cout << "-------------------------------------------------------------------" <<std::endl;
}

void samplePDFDUNEAtmBase::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splinesDUNE* DUNESplines = new splinesDUNE(XsecCov);
    splineFile = (splineFDBase*)DUNESplines;
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splineFile = nullptr;
  }
  
  return;
}

void samplePDFDUNEAtmBase::SetupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
    for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
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
}

int samplePDFDUNEAtmBase::setupExperimentMC(int iSample) {
  const char *sampleFile = (mtupleprefix+mtuple_files[iSample]+mtuplesuffix).c_str();
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  int nutype = sample_nutype[iSample];
  int oscnutype = sample_oscnutype[iSample];
  bool signal = sample_signal[iSample];

  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;
  
  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree*)_sampleFile->Get("cafTree");

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

  duneobj->norm_s = 1;
  duneobj->pot_s = 1;

  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;

  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->rw_Q2 = new double[duneobj->nEvents];
  duneobj->beam_w = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_cvnnumu = new double[duneobj->nEvents];
  duneobj->rw_cvnnue = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_run = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_truecz = new double[duneobj->nEvents];

  duneobj->mode = new int[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];

  _data->GetEntry(0);
  for (int i = 0; i < duneobj->nEvents; ++i) { // Loop through tree
    _data->GetEntry(i);
    
    if (iselike) {
      duneobj->rw_erec[i] = _erec_nue;
    } else {
      duneobj->rw_erec[i] = _erec;
    }
    
    duneobj->rw_etru[i] = _ev; // in GeV
    
    duneobj->rw_cvnnumu[i] = _cvnnumu; 
    duneobj->rw_cvnnue[i] = _cvnnue;
    duneobj->rw_theta[i] = _LepMomY/sqrt(pow(_LepMomZ, 2) + pow(_LepMomY, 2) + pow(_LepMomX, 2));
    duneobj->rw_Q2[i] = _Q2;
    duneobj->rw_isCC[i] = _isCC;
    duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj->rw_nuPDG[i] = _nuPDG;
    duneobj->rw_run[i] = _run;
    duneobj->rw_berpaacvwgt[i] = _BeRPA_cvwgt;
    duneobj->rw_vtx_x[i] = _vtx_x;
    duneobj->rw_vtx_y[i] = _vtx_y;
    duneobj->rw_vtx_z[i] = _vtx_z;
    duneobj->rw_truecz[i] = _NuMomY/_ev;
    
    //Assume everything is on Argon for now....
    duneobj->Target[i] = 40;
    
    duneobj->beam_w[i] = 1.0;
    int mode= TMath::Abs(_mode);       
    
    duneobj->mode[i]=SIMBMode_ToMaCh3Mode(mode, _isCC);
    duneobj->flux_w[i] = _weight;
  }
  
  _sampleFile->Close();
  return duneobj->nEvents;
}

double* samplePDFDUNEAtmBase::ReturnKinematicParameterByReference(KinematicTypes KinPar, int iSample, int iEvent) {
  double* KinematicValue;
  
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &dunemcSamples[iSample].rw_etru[iEvent]; 
    break;
  case kTrueXPos:
    KinematicValue = &dunemcSamples[iSample].rw_vtx_x[iEvent];
    break;
  case kTrueYPos:
    KinematicValue = &dunemcSamples[iSample].rw_vtx_y[iEvent];
    break;
  case kTrueZPos:
    KinematicValue = &dunemcSamples[iSample].rw_vtx_z[iEvent];
    break;
  case kTrueCosZ:
    KinematicValue = &dunemcSamples[iSample].rw_truecz[iEvent];
    break;
  default:
    std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type '" << KinPar << "'" << std::endl;
    throw;
  }
  
  return KinematicValue;
}

double* samplePDFDUNEAtmBase::ReturnKinematicParameterByReference(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
  return ReturnKinematicParameterByReference(KinPar,iSample,iEvent);
}

double* samplePDFDUNEAtmBase::ReturnKinematicParameterByReference(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return ReturnKinematicParameterByReference(KinPar,iSample,iEvent);
}

double samplePDFDUNEAtmBase::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  return *ReturnKinematicParameterByReference(KinematicVariable, iSample, iEvent);
}

double samplePDFDUNEAtmBase::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  return *ReturnKinematicParameterByReference(KinematicParameter, iSample, iEvent);
}

void samplePDFDUNEAtmBase::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  fdmc_base *fdobj = &(MCSamples[iSample]);
  
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->SampleDetID = SampleDetID;
  
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent) {
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
    fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
    fdobj->flux_w[iEvent] = duneobj->flux_w[iEvent];
    fdobj->rw_truecz[iEvent] = duneobj->rw_truecz[iEvent];
    
    //ETA - this is where the variables that you want to bin your samples in are defined
    //If you want to bin in different variables this is where you put it for now
    switch(nDimensions){
    case 0:
    case 1:
      //Just point to xvar to the address of the variable you want to bin in
      //This way we don't have to update both fdmc and skmc when we apply shifts
      //to variables we're binning in
      fdobj->x_var[iEvent] = &(duneobj->rw_erec[iEvent]);
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
      std::cout << "[ERROR:] " << __FILE__ << ":" << __LINE__ << " unrecognised binning option" << nDimensions << std::endl;
      throw;
      break;
    }
  }
}

std::vector<double> samplePDFDUNEAtmBase::ReturnKinematicParameterBinning(std::string KinematicParameterStr)  {
  std::cout << "ReturnKinematicVarBinning" << std::endl;
  std::vector<double> binningVector;
  return binningVector;
}
