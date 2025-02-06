#include "samplePDFDUNEBeamND.h"

//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"

//ROOT includes
#include "TError.h"

samplePDFDUNEBeamND::samplePDFDUNEBeamND(std::string mc_version_, covarianceXsec* xsec_cov_, covarianceOsc* osc_cov_) : samplePDFFDBase(mc_version_, xsec_cov_, osc_cov_) {  
  Initialise();
}

samplePDFDUNEBeamND::~samplePDFDUNEBeamND() {
}

void samplePDFDUNEBeamND::Init() {
  dunendmcSamples.resize(nSamples,dunemc_base());
  
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  IsELike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  tot_escale_nd_pos = -999;
  tot_escale_sqrt_nd_pos = -999;
  tot_escale_invsqrt_nd_pos = -999;
  had_escale_nd_pos = -999;
  had_escale_sqrt_nd_pos = -999;
  had_escale_invsqrt_nd_pos = -999;
  mu_escale_nd_pos = -999;
  mu_escale_sqrt_nd_pos = -999;
  mu_escale_invsqrt_nd_pos = -999;
  n_escale_nd_pos = -999;
  n_escale_sqrt_nd_pos = -999;
  n_escale_invsqrt_nd_pos = -999;
  em_escale_nd_pos = -999;
  em_escale_sqrt_nd_pos = -999;
  em_escale_invsqrt_nd_pos = -999;
  had_res_nd_pos = -999;
  mu_res_nd_pos = -999;
  n_res_nd_pos = -999;
  em_res_nd_pos = -999;
  IsELike = SampleManager->raw()["SampleBools"]["IsELike"].as<bool>();
 
  std::cout << "-------------------------------------------------------------------" <<std::endl;
}

void samplePDFDUNEBeamND::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    SplineHandler = std::unique_ptr<splineFDBase>(new splinesDUNE(XsecCov));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    SplineHandler = nullptr;
  }

  return;
}

void samplePDFDUNEBeamND::SetupWeightPointers() {
  for (int i = 0; i < (int)dunendmcSamples.size(); ++i) {
    for (int j = 0; j < dunendmcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 6;
      MCSamples[i].total_weight_pointers[j].resize(MCSamples[i].ntotal_weight_pointers[j]);
      MCSamples[i].total_weight_pointers[j][0] = &(dunendmcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunendmcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][3] = &(dunendmcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(dunendmcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
    }
  }
}

int samplePDFDUNEBeamND::setupExperimentMC(int iSample) {
  int CurrErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  caf::StandardRecord* sr = new caf::StandardRecord();
  dunemc_base* duneobj = &dunendmcSamples[iSample];


  std::string FileName = mc_files[iSample];
  MACH3LOG_INFO("Reading File: {}",FileName);
  TFile* File = TFile::Open(FileName.c_str());
  if (!File || File->IsZombie()) {
    MACH3LOG_ERROR("Did not find File: {}",FileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  TTree* Tree = File->Get<TTree>("cafTree");
  if (!Tree){
    MACH3LOG_ERROR("Did not find Tree::cafTree in File: {}",FileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  std::cout << "nEvents:" <<duneobj->nEvents<< std::endl;

  Tree->SetBranchStatus("*", 1);
  Tree->SetBranchAddress("rec", &sr);

  if (!IsRHC) { 
    duneobj->norm_s = (1e21/1.5e21);
  } else {
    duneobj->norm_s = (1e21/1.905e21);
  }
  duneobj->pot_s = (pot)/1e21;

  duneobj->nEvents = Tree->GetEntries();
  duneobj->nupdg = new int[duneobj->nEvents];
  duneobj->nupdgUnosc = new int[duneobj->nEvents];

  duneobj->rw_yrec = new double[duneobj->nEvents];
  duneobj->rw_erec_lep = new double[duneobj->nEvents];
  duneobj->rw_erec_had = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_reco_q = new double[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 

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

  duneobj->nupdg = new int[duneobj->nEvents];
  duneobj->nupdgUnosc = new int[duneobj->nEvents];
  duneobj->mode = new double[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];

  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];

  duneobj->rw_vtx_x_end = new double[duneobj->nEvents];
  duneobj->rw_vtx_y_end = new double[duneobj->nEvents];
  duneobj->rw_vtx_z_end = new double[duneobj->nEvents];

  duneobj->rw_reco_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_z = new double[duneobj->nEvents];

  duneobj->rw_isFHC = new bool[duneobj->nEvents];
  
  //FILL DUNE STRUCT
  for (int iEvent=0;iEvent<duneobj->nEvents; iEvent++) { // Loop through tree
    Tree->GetEntry(iEvent);

    int ndlp = sr->common.ixn.ndlp;

    for(int i=0; i<ndlp; i++) {
      int ntracks = sr->nd.lar.dlp[i].ntracks;
      for(int j=0; j<ntracks; j++) {
          duneobj->rw_erec[iEvent] += (double)(sr->nd.lar.dlp[i].tracks[j].E/1000.); // MeV to GeV
      }
    }
    duneobj->rw_reco_vtx_x[iEvent] = (double)(sr->nd.lar.dlp[0].tracks[0].start.X());
    duneobj->rw_reco_vtx_y[iEvent] = (double)(sr->nd.lar.dlp[0].tracks[0].start.Y());
    duneobj->rw_reco_vtx_z[iEvent] = (double)(sr->nd.lar.dlp[0].tracks[0].start.Z());
    

    int nnu = sr->mc.nnu;
    
    duneobj->rw_etru[iEvent] += (double)(sr->mc.nu[0].E);
    duneobj->rw_vtx_x[iEvent] = (double)(sr->mc.nu[0].vtx.X());
    duneobj->rw_vtx_y[iEvent] = (double)(sr->mc.nu[0].vtx.Y());
    duneobj->rw_vtx_z[iEvent] = (double)(sr->mc.nu[0].vtx.Z());
    
    // duneobj->rw_theta[iEvent] = RecoNuMomentumVector.Y();
    // std::cout << "rw_erec:" <<rw_erec[iEvent]<< std::endl;
    duneobj->rw_erec_lep[iEvent] = (double)(sr->common.ixn.dlp[0].Enu.lep_calo);
    duneobj->rw_erec_had[iEvent] = (double)(sr->common.ixn.dlp[0].Enu.calo);
    duneobj->rw_yrec[iEvent] = (double)((_erec-_erec_lep)/_erec);
    
    duneobj->nupdg[iEvent] = (int)(sample_nupdg[iSample]);
    std::cout << "nuPDG:" <<duneobj->nupdg[iEvent]<< std::endl;
    duneobj->nupdgUnosc[iEvent] = (int)(sample_nupdgunosc[iSample]);
    std::cout << "nuPDGunosc:" <<duneobj->nupdgUnosc[iEvent]<< std::endl;
    // duneobj->rw_theta[iEvent] = (double)_LepNuAngle;
    duneobj->rw_isCC[iEvent] = sr->mc.nu[0].iscc;
    // duneobj->rw_isFHC[iEvent] = (double)(1.0);
    // duneobj->rw_reco_q[iEvent] = _reco_q;
    // duneobj->rw_nuPDGunosc[iEvent] = _nuPDGunosc;
    // duneobj->rw_nuPDG[iEvent] = _nuPDG;
    // duneobj->rw_berpaacvwgt[iEvent] = (double)_BeRPA_cvwgt;
    
    // duneobj->rw_eRecoP[iEvent] = (double)_eRecoP; 
    // duneobj->rw_eRecoPip[iEvent] = (double)_eRecoPip; 
    // duneobj->rw_eRecoPim[iEvent] = (double)_eRecoPim; 
    // duneobj->rw_eRecoPi0[iEvent] = (double)_eRecoPi0; 
    // duneobj->rw_eRecoN[iEvent] = (double)_eRecoN; 
    
    // duneobj->rw_LepE[iEvent] = (double)_LepE; 
    // duneobj->rw_eP[iEvent] = (double)_eP; 
    // duneobj->rw_ePip[iEvent] = (double)_ePip; 
    // duneobj->rw_ePim[iEvent] = (double)_ePim; 
    // duneobj->rw_ePi0[iEvent] = (double)_ePi0; 
    // duneobj->rw_eN[iEvent] = (double)_eN; 
    
    //Assume everything is on Argon for now....
    duneobj->Target[iEvent] = 40;
    
    int mode= TMath::Abs(_mode);       
    duneobj->mode[iEvent]=(double)GENIEMode_ToMaCh3Mode(mode, _isCC);
    


    duneobj->flux_w[iEvent] = (double)(sr->mc.nu[0].genweight);
  }

  delete Tree;
  delete File;

  gErrorIgnoreLevel = CurrErrorLevel;

  return duneobj->nEvents;
}


const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  double* KinematicValue;
  
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &dunendmcSamples[iSample].rw_etru[iEvent]; 
    break;
  case kRecoQ:
    KinematicValue = &dunendmcSamples[iSample].rw_reco_q[iEvent];
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &dunendmcSamples[iSample].rw_erec[iEvent];
    break;
  case kIsFHC:
    KinematicValue = static_cast<double*>(static_cast<void*>(&dunendmcSamples[iSample].rw_isFHC[iEvent]));
    break;
  case kOscChannel:
    KinematicValue = static_cast<double*>(static_cast<void*>(&dunendmcSamples[iSample].nupdgUnosc[iEvent]));
    break;
  case kMode:
    KinematicValue = &dunendmcSamples[iSample].mode[iEvent];
    break;
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    std::cout << KinPar << ReturnStringFromKinematicParameter(KinPar) << std::endl;
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}

const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEBeamND::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicVariable, iSample, iEvent);
}

double samplePDFDUNEBeamND::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicParameter, iSample, iEvent);
}

void samplePDFDUNEBeamND::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);
  
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
    std::cout << "This is Event "<<iEvent<< std::endl;
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    std::cout << "rw_etru:" <<duneobj->rw_etru[iEvent]<< std::endl;
    std::cout << "rw_erec:" <<duneobj->rw_erec[iEvent]<< std::endl;
    std::cout << "rw_erec_lep:" <<duneobj->rw_erec_lep[iEvent]<< std::endl;
    std::cout << "rw_erec_had:" <<duneobj->rw_erec_had[iEvent]<< std::endl;
    std::cout << "rw_reco_vtx_x:" <<duneobj->rw_reco_vtx_x[iEvent]<< std::endl;
    std::cout << "rw_vtx_x:" <<duneobj->rw_vtx_x[iEvent]<< std::endl;    
    std::cout << "rw_reco_vtx_y:" <<duneobj->rw_reco_vtx_y[iEvent]<< std::endl;
    std::cout << "rw_vtx_y:" <<duneobj->rw_vtx_y[iEvent]<< std::endl;
    std::cout << "rw_reco_vtx_z:" <<duneobj->rw_reco_vtx_z[iEvent]<< std::endl;
    std::cout << "rw_vtx_z:" <<duneobj->rw_vtx_z[iEvent]<< std::endl;
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    std::cout << "mode:" <<duneobj->mode[iEvent]<< std::endl;
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
    std::cout << "Target:" <<duneobj->Target[iEvent]<< std::endl;
    fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
    std::cout << "isCC:" <<duneobj->rw_isCC[iEvent]<< std::endl;
    std::cout << "isNC:" <<fdobj->isNC[iEvent]<< std::endl;
    fdobj->nupdgUnosc[iEvent] = &(duneobj->nupdgUnosc[iEvent]);
    std::cout << "nupdgUnosc:" <<duneobj->nupdgUnosc[iEvent]<< std::endl;
    fdobj->nupdg[iEvent] = &(duneobj->nupdg[iEvent]);
    std::cout << "nupdg:" <<duneobj->nupdg[iEvent]<< std::endl;
    std::cout << "-------------------------------------------------------------------" <<std::endl;
  }
}

void samplePDFDUNEBeamND::applyShifts(int iSample, int iEvent) {
  // reset erec back to original value
  
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] = dunendmcSamples[iSample].rw_erec[iEvent];

  //Calculate values needed
  double sqrtErecHad =  sqrt(dunendmcSamples[iSample].rw_erec_had[iEvent]);
  double sqrtErecLep =  sqrt(dunendmcSamples[iSample].rw_erec_lep[iEvent]);
  double sqrteRecoPi0 = sqrt(dunendmcSamples[iSample].rw_eRecoPi0[iEvent]);
  double sqrteRecoN = sqrt(dunendmcSamples[iSample].rw_eRecoN[iEvent]);
  double sumEhad = dunendmcSamples[iSample].rw_eRecoP[iEvent] + dunendmcSamples[iSample].rw_eRecoPip[iEvent] + dunendmcSamples[iSample].rw_eRecoPim[iEvent];
  double sqrtSumEhad = sqrt(sumEhad);

  double invSqrtErecHad =  1/(sqrtErecHad+0.1);
  double invSqrtErecLep =  1/(sqrtErecLep+0.1);
  double invSqrteRecoPi0 =  1/(sqrteRecoPi0+0.1);
  double invSqrteRecoN =  1/(sqrteRecoN+0.1);
  double invSqrtSumEhad =  1/(sqrtSumEhad+0.1);

  bool CCnumu {dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].nupdg[iEvent])==14 && dunendmcSamples[iSample].nupdgUnosc[iEvent]==2};
  bool CCnue {dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].nupdg[iEvent])==12 && dunendmcSamples[iSample].nupdgUnosc[iEvent]==1};
  bool NotCCnumu {!(dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].nupdg[iEvent])==14) && dunendmcSamples[iSample].nupdgUnosc[iEvent]==2};

}

std::vector<double> samplePDFDUNEBeamND::ReturnKinematicParameterBinning(std::string KinematicParameterStr) 
{
  std::vector<double> binningVector;
  return binningVector;
}

int samplePDFDUNEBeamND::ReturnKinematicParameterFromString(std::string KinematicParameterStr) {

  if(KinematicParameterStr == "TrueNeutrinoEnergy") return kTrueNeutrinoEnergy;
  if(KinematicParameterStr == "RecoQ") return kRecoQ;
  if(KinematicParameterStr == "RecoNeutrinoEnergy") return kRecoNeutrinoEnergy;
  if(KinematicParameterStr == "IsFHC") return kRecoNeutrinoEnergy;
  if(KinematicParameterStr == "OscChannel") return kOscChannel;
  if(KinematicParameterStr == "Mode") return kMode;
  std::cout << "HELP ME " <<KinematicParameterStr<<std::endl;
  return -1;
}

std::string samplePDFDUNEBeamND::ReturnStringFromKinematicParameter(int KinematicParameter) {
  switch(KinematicParameter){
    case kTrueNeutrinoEnergy: return "TrueNeutrinoEnergy";
    case kRecoQ: return "RecoQ";
    case kRecoNeutrinoEnergy: return "RecoNeutrinoEnergy";
    case kIsFHC: return "IsFHC";
    case kOscChannel: return "OscChannel";
    case kMode: return "Mode";
    default: return "";
  }
}
