#include <TROOT.h>

#include "samplePDFDUNEBeamFDBase.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "manager/manager.h"

samplePDFDUNEBeamFDBase::samplePDFDUNEBeamFDBase(double pot_, std::string mc_version_, covarianceXsec* xsec_cov_) : samplePDFFDBase(pot_, mc_version_, xsec_cov_) {
  // create dunemc storage
  for (int i=0;i<nSamples;i++) {
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
  }
  
  Initialise();
}

samplePDFDUNEBeamFDBase::~samplePDFDUNEBeamFDBase() {
}

void samplePDFDUNEBeamFDBase::Init() {
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  tot_escale_fd_pos = -999;
  tot_escale_sqrt_fd_pos = -999;
  tot_escale_invsqrt_fd_pos = -999;
  had_escale_fd_pos = -999;
  had_escale_sqrt_fd_pos = -999;
  had_escale_invsqrt_fd_pos = -999;
  mu_escale_fd_pos = -999;
  mu_escale_sqrt_fd_pos = -999;
  mu_escale_invsqrt_fd_pos = -999;
  n_escale_fd_pos = -999;
  n_escale_sqrt_fd_pos = -999;
  n_escale_invsqrt_fd_pos = -999;
  em_escale_fd_pos = -999;
  em_escale_sqrt_fd_pos = -999;
  em_escale_invsqrt_fd_pos = -999;
  had_res_fd_pos = -999;
  mu_res_fd_pos = -999;
  n_res_fd_pos = -999;
  em_res_fd_pos = -999;
  cvn_numu_fd_pos = -999;
  cvn_nue_fd_pos = -999;

  nFDDetectorSystPointers = funcParsIndex.size();
  FDDetectorSystPointers = std::vector<const double*>(nFDDetectorSystPointers);

  int func_it = 0;
  for (std::vector<int>::iterator it = funcParsIndex.begin(); it != funcParsIndex.end(); ++it, ++func_it) {
    std::string name = funcParsNames.at(func_it);
    
    if (name == "TotalEScaleFD") {
      tot_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_fd_pos);
    }
    else if (name == "TotalEScaleSqrtFD") {
      tot_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_sqrt_fd_pos);
    }
    else if (name == "TotalEScaleInvSqrtFD") {
      tot_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_invsqrt_fd_pos);
    }
    else if (name == "HadEScaleFD") {
      had_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_fd_pos);
    }
    else if (name == "HadEScaleSqrtFD") {
      had_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_sqrt_fd_pos);
    }
    else if (name == "HadEScaleInvSqrtFD") {
      had_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_invsqrt_fd_pos);
    }
    else if (name == "MuEScaleFD") {
      mu_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_fd_pos);
    }
    else if (name == "MuEScaleSqrtFD") {
      mu_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_sqrt_fd_pos);
    }
    else if (name == "MuEScaleInvSqrtFD") {
      mu_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_invsqrt_fd_pos);
    }
    else if (name == "NEScaleFD") {
      n_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_fd_pos);
    }
    else if (name == "NEScaleSqrtFD") {
      n_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_sqrt_fd_pos);
    }
    else if (name == "NEScaleInvSqrtFD") {
      n_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_invsqrt_fd_pos);
    }
    else if (name == "EMEScaleFD") {
      em_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_fd_pos);
    }
    else if (name == "EMEScaleSqrtFD") {
      em_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_sqrt_fd_pos);
    }
    else if (name == "EMEScaleInvSqrtFD") {
      em_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_invsqrt_fd_pos);
    }
    else if (name == "HadResFD") {
      had_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_res_fd_pos);
    }
    else if (name == "MuResFD") {
      mu_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_res_fd_pos);
    }
    else if (name == "NResFD") {
      n_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_res_fd_pos);
    }
    else if (name == "EMResFD") {
      em_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_res_fd_pos);
    }
    else if (name == "CVNNumuFD") {
      cvn_numu_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(cvn_numu_fd_pos);
    }
    else if (name == "CVNNueFD") {
      cvn_nue_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(cvn_nue_fd_pos);
    }
    
    else { 
      std::cerr << "Found a functional parameter which wasn't specified in the xml | samplePDFDUNEBeamFDBase:" << name << std::endl;
      throw;
    }
  }

  splinesDUNE* DUNESplines = new splinesDUNE(XsecCov);
  splineFile = (splineFDBase*)DUNESplines;
  InitialiseSplineObject();

  std::cout << "-------------------------------------------------------------------" <<std::endl;
}

void samplePDFDUNEBeamFDBase::SetupWeightPointers() {
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


void samplePDFDUNEBeamFDBase::setupExperimentMC(int iSample) {
  const char *sampleFile = (mtupleprefix+mtuple_files[iSample]+mtuplesuffix).c_str();
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  int nutype = sample_nutype[iSample];
  int oscnutype = sample_oscnutype[iSample];
  bool signal = sample_signal[iSample];
  
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;
  
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
  _data->SetBranchStatus("RecoHadEnNumu", 1);
  _data->SetBranchAddress("RecoHadEnNumu", &_erec_had);
  _data->SetBranchStatus("RecoHadEnNue", 1);
  _data->SetBranchAddress("RecoHadEnNue", &_erec_had_nue);
  _data->SetBranchStatus("RecoLepEnNumu", 1);
  _data->SetBranchAddress("RecoLepEnNumu", &_erec_lep);
  _data->SetBranchStatus("RecoLepEnNue", 1);
  _data->SetBranchAddress("RecoLepEnNue", &_erec_lep_nue);

  _data->SetBranchStatus("eRecoP", 1);
  _data->SetBranchAddress("eRecoP", &_eRecoP);
  _data->SetBranchStatus("eRecoPip", 1);
  _data->SetBranchAddress("eRecoPip", &_eRecoPip);
  _data->SetBranchStatus("eRecoPim", 1);
  _data->SetBranchAddress("eRecoPim", &_eRecoPim);
  _data->SetBranchStatus("eRecoPi0", 1);
  _data->SetBranchAddress("eRecoPi0", &_eRecoPi0);
  _data->SetBranchStatus("eRecoN", 1);
  _data->SetBranchAddress("eRecoN", &_eRecoN);

  _data->SetBranchStatus("LepE", 1);
  _data->SetBranchAddress("LepE", &_LepE);
  _data->SetBranchStatus("eP", 1);
  _data->SetBranchAddress("eP", &_eP);
  _data->SetBranchStatus("ePip", 1);
  _data->SetBranchAddress("ePip", &_ePip);
  _data->SetBranchStatus("ePim", 1);
  _data->SetBranchAddress("ePim", &_ePim);
  _data->SetBranchStatus("ePi0", 1);
  _data->SetBranchAddress("ePi0", &_ePi0);
  _data->SetBranchStatus("eN", 1);
  _data->SetBranchAddress("eN", &_eN);

  _data->SetBranchStatus("mode",1);
  _data->SetBranchAddress("mode",&_mode);
  _data->SetBranchStatus("cvnnumu",1);
  _data->SetBranchAddress("cvnnumu", &_cvnnumu);
  _data->SetBranchStatus("cvnnue",1);
  _data->SetBranchAddress("cvnnue", &_cvnnue);
  _data->SetBranchStatus("isCC", 1);
  _data->SetBranchAddress("isCC", &_isCC);
  _data->SetBranchStatus("nuPDGunosc", 1);
  _data->SetBranchAddress("nuPDGunosc", &_nuPDGunosc);
  _data->SetBranchStatus("nuPDG", 1);
  _data->SetBranchAddress("nuPDG", &_nuPDG);
  _data->SetBranchStatus("BeRPA_A_cvwgt", 1);
  _data->SetBranchAddress("BeRPA_A_cvwgt", &_BeRPA_cvwgt);
  _data->SetBranchStatus("vtx_x", 1);
  _data->SetBranchAddress("vtx_x", &_vtx_x);
  _data->SetBranchStatus("vtx_y", 1);
  _data->SetBranchAddress("vtx_y", &_vtx_y);
  _data->SetBranchStatus("vtx_z", 1);
  _data->SetBranchAddress("vtx_z", &_vtx_z);  

  TH1D* norm = (TH1D*)_sampleFile->Get("norm");
  if(!norm){
    std::cout<< "Add a norm KEY to the root file using MakeNormHists.cxx"<<std::endl;
    std::cout << "Ignoring for now" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // now fill the actual variables
  duneobj->norm_s = norm->GetBinContent(1);
  duneobj->pot_s = (pot)/norm->GetBinContent(2);

  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;

  // allocate memory for dunemc variables
  duneobj->rw_cvnnumu = new double[duneobj->nEvents];
  duneobj->rw_cvnnue = new double[duneobj->nEvents];
  duneobj->rw_cvnnumu_shifted = new double[duneobj->nEvents];
  duneobj->rw_cvnnue_shifted = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_erec_had = new double[duneobj->nEvents];
  duneobj->rw_erec_lep = new double[duneobj->nEvents];

  duneobj->rw_eRecoP = new double[duneobj->nEvents];
  duneobj->rw_eRecoPip = new double[duneobj->nEvents];
  duneobj->rw_eRecoPim = new double[duneobj->nEvents];
  duneobj->rw_eRecoPi0 = new double[duneobj->nEvents];
  duneobj->rw_eRecoN = new double[duneobj->nEvents];

  duneobj->rw_LepE = new double[duneobj->nEvents];
  duneobj->rw_eP = new double[duneobj->nEvents];
  duneobj->rw_ePip = new double[duneobj->nEvents];
  duneobj->rw_ePim = new double[duneobj->nEvents];
  duneobj->rw_ePi0 = new double[duneobj->nEvents];
  duneobj->rw_eN = new double[duneobj->nEvents];

  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];

  duneobj->mode = new int[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];

  _data->GetEntry(0);

  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) { // Loop through tree
    _data->GetEntry(i);
    duneobj->rw_cvnnumu[i] = (double)_cvnnumu;
    duneobj->rw_cvnnue[i] = (double)_cvnnue;
    duneobj->rw_cvnnumu_shifted[i] = (double)_cvnnumu; 
    duneobj->rw_cvnnue_shifted[i] = (double)_cvnnue;
    if (iselike) {
      duneobj->rw_erec[i] = (double)_erec_nue;
      duneobj->rw_erec_shifted[i] = (double)_erec_nue; 
      duneobj->rw_erec_had[i] = (double)_erec_had_nue;
      duneobj->rw_erec_lep[i] = (double)_erec_lep_nue;
    } else {
      duneobj->rw_erec[i] = (double)_erec; 
      duneobj->rw_erec_shifted[i] = (double)_erec; 
      duneobj->rw_erec_had[i] = (double)_erec_had; 
      duneobj->rw_erec_lep[i] = (double)_erec_lep; 
    }
    
    duneobj->rw_eRecoP[i] = (double)_eRecoP; 
    duneobj->rw_eRecoPip[i] = (double)_eRecoPip; 
    duneobj->rw_eRecoPim[i] = (double)_eRecoPim; 
    duneobj->rw_eRecoPi0[i] = (double)_eRecoPi0; 
    duneobj->rw_eRecoN[i] = (double)_eRecoN; 
    
    duneobj->rw_LepE[i] = (double)_LepE; 
    duneobj->rw_eP[i] = (double)_eP; 
    duneobj->rw_ePip[i] = (double)_ePip; 
    duneobj->rw_ePim[i] = (double)_ePim; 
    duneobj->rw_ePi0[i] = (double)_ePi0; 
    duneobj->rw_eN[i] = (double)_eN; 
    
    duneobj->rw_etru[i] = (double)_ev;
    duneobj->rw_theta[i] = (double)_LepNuAngle;
    duneobj->rw_isCC[i] = _isCC;
    duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj->rw_nuPDG[i] = _nuPDG;
    duneobj->rw_berpaacvwgt[i] = (double)_BeRPA_cvwgt;
    duneobj->rw_vtx_x[i] = (double)_vtx_x;
    duneobj->rw_vtx_y[i] = (double)_vtx_y;
    duneobj->rw_vtx_z[i] = (double)_vtx_z;
    
    //Assume everything is on Argon for now....
    duneobj->Target[i] = 40;
    
    int mode= TMath::Abs(_mode);       
    duneobj->mode[i]=SIMBMode_ToMaCh3Mode(mode, _isCC);
    
    duneobj->flux_w[i] = 1.0;
  }
  
  _sampleFile->Close();
  return duneobj->nEvents;
}

double samplePDFDUNEBeamFDBase::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
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
  case kCVNNumu:
    KinematicValue = dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
    break;
  case kCVNNue:
    KinematicValue = dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
    break;
  default:
    std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type..." << std::endl;
    throw;
  }
  
  return KinematicValue;
}

double samplePDFDUNEBeamFDBase::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
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
 case kCVNNumu:
   KinematicValue = dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
   break;
 case kCVNNue:
   KinematicValue = dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
   break;
 default:
   std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type..." << std::endl;
   throw;
 }
 
 return KinematicValue;
}

void samplePDFDUNEBeamFDBase::setupFDMC(int iSample) {
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
    
    //ETA - this is where the variables that you want to bin your samples in are defined
    //If you want to bin in different variables this is where you put it for now
    switch(BinningOpt){
    case 0:
    case 1:
      //Just point to xvar to the address of the variable you want to bin in
      //This way we don't have to update both fdmc and skmc when we apply shifts
      //to variables we're binning in
      fdobj->x_var[iEvent] = &(duneobj->rw_erec_shifted[iEvent]);
      fdobj->y_var[iEvent] = &(duneobj->dummy_y);//ETA - don't think we even need this as if we have a 1D sample we never need this, just not sure I like an unitialised variable in fdmc struct? 
      break;
    case 2:
      //Just point to xvar to the address of the variable you want to bin in
      //This way we don't have to update both fdmc and skmc when we apply shifts
      //to variables we're binning in
      fdobj->x_var[iEvent] = &(duneobj->rw_erec_shifted[iEvent]);
      fdobj->y_var[iEvent] = &(duneobj->rw_theta[iEvent]);
      break;
    default:
      std::cout << "[ERROR:] " << __FILE__ << ":" << __LINE__ << " unrecognised binning option" << BinningOpt << std::endl;
      throw;
      break;
    }
    
  }
}
 
void samplePDFDUNEBeamFDBase::applyShifts(int iSample, int iEvent) {
   
   // reset erec back to original value
  dunemcSamples[iSample].rw_erec_shifted[iEvent] = dunemcSamples[iSample].rw_erec[iEvent];

  // reset cvnnumu back to original value
  dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent] = dunemcSamples[iSample].rw_cvnnumu[iEvent];

  // reset cvnnue back to original value
  dunemcSamples[iSample].rw_cvnnue_shifted[iEvent] = dunemcSamples[iSample].rw_cvnnue[iEvent];

  //Calculate values needed
  double sqrtErecHad =  sqrt(dunemcSamples[iSample].rw_erec_had[iEvent]);
  double sqrtErecLep =  sqrt(dunemcSamples[iSample].rw_erec_lep[iEvent]);
  double sqrteRecoPi0 = sqrt(dunemcSamples[iSample].rw_eRecoPi0[iEvent]);
  double sqrteRecoN = sqrt(dunemcSamples[iSample].rw_eRecoN[iEvent]);
  double sumEhad = dunemcSamples[iSample].rw_eRecoP[iEvent] + dunemcSamples[iSample].rw_eRecoPip[iEvent] + dunemcSamples[iSample].rw_eRecoPim[iEvent];
  double sqrtSumEhad = sqrt(sumEhad);

  double invSqrtErecHad =  1/(sqrtErecHad+0.1);
  double invSqrtErecLep =  1/(sqrtErecLep+0.1);
  double invSqrteRecoPi0 =  1/(sqrteRecoPi0+0.1);
  double invSqrteRecoN =  1/(sqrteRecoN+0.1);
  double invSqrtSumEhad =  1/(sqrtSumEhad+0.1);

  bool CCnumu {dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==14) && dunemcSamples[iSample].nutype==2};
  bool CCnue {dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==12) && dunemcSamples[iSample].nutype==1};
  bool NotCCnumu {!(dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==14)) && dunemcSamples[iSample].nutype==2};


  TotalEScaleFD(FDDetectorSystPointers[0], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_had[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], NotCCnumu);

  TotalEScaleSqrtFD(FDDetectorSystPointers[1], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_had[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecHad, sqrtErecLep, NotCCnumu);

  TotalEScaleInvSqrtFD(FDDetectorSystPointers[2], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_had[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecHad, invSqrtErecLep, NotCCnumu);

  HadEScaleFD(FDDetectorSystPointers[3], &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad);

  HadEScaleSqrtFD(FDDetectorSystPointers[4], &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, sqrtSumEhad);

  HadEScaleInvSqrtFD(FDDetectorSystPointers[5], &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, invSqrtSumEhad);

  MuEScaleFD(FDDetectorSystPointers[6], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], CCnumu);

  MuEScaleSqrtFD(FDDetectorSystPointers[7], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, CCnumu);

  MuEScaleInvSqrtFD(FDDetectorSystPointers[8], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, CCnumu);

  NEScaleFD(FDDetectorSystPointers[9], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent]);

  NEScaleSqrtFD(FDDetectorSystPointers[10], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent], sqrteRecoN);

  NEScaleInvSqrtFD(FDDetectorSystPointers[11], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent], invSqrteRecoN);

  EMEScaleFD(FDDetectorSystPointers[12], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], CCnue);

  EMEScaleSqrtFD(FDDetectorSystPointers[13], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, sqrteRecoPi0, CCnue);

  EMEScaleInvSqrtFD(FDDetectorSystPointers[14], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, invSqrteRecoPi0, CCnue);

  HadResFD(FDDetectorSystPointers[15], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoP[iEvent], dunemcSamples[iSample].rw_eRecoPip[iEvent], dunemcSamples[iSample].rw_eRecoPim[iEvent], dunemcSamples[iSample].rw_eP[iEvent], dunemcSamples[iSample].rw_ePip[iEvent], dunemcSamples[iSample].rw_ePim[iEvent]);

  MuResFD(FDDetectorSystPointers[16], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], dunemcSamples[iSample].rw_LepE[iEvent], CCnumu);

  NResFD(FDDetectorSystPointers[17], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent], dunemcSamples[iSample].rw_eN[iEvent]);

  EMResFD(FDDetectorSystPointers[18], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_ePi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], dunemcSamples[iSample].rw_LepE[iEvent], CCnue);

  CVNNumuFD(FDDetectorSystPointers[19], &dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent]);

  CVNNueFD(FDDetectorSystPointers[20], &dunemcSamples[iSample].rw_cvnnue_shifted[iEvent]);
}

std::vector<double> samplePDFDUNEBeamFDBase::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  std::cout << "ReturnKinematicVarBinning" << std::endl;
  std::vector<double> binningVector;
  return binningVector;
}
