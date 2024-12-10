#include "samplePDFDUNEAtm.h"

//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"

//ROOT includes
#include "TError.h"

samplePDFDUNEAtm::samplePDFDUNEAtm(std::string mc_version_, covarianceXsec* xsec_cov_, covarianceOsc* osc_cov_) : samplePDFFDBase(mc_version_, xsec_cov_, osc_cov_) {
  Initialise();
}

samplePDFDUNEAtm::~samplePDFDUNEAtm() {
}

void samplePDFDUNEAtm::Init() {
  dunemcSamples.resize(nSamples,dunemc_base());
  
  IsELike = SampleManager->raw()["SampleBools"]["IsELike"].as<bool>();
}

void samplePDFDUNEAtm::SetupSplines() {
  SplineHandler = nullptr;
}

void samplePDFDUNEAtm::SetupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
    for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 3;
      MCSamples[i].total_weight_pointers[j].push_back(&(dunemcSamples[i].flux_w[j]));
      MCSamples[i].total_weight_pointers[j].push_back(MCSamples[i].osc_w_pointer[j]);
      MCSamples[i].total_weight_pointers[j].push_back(&(MCSamples[i].xsec_w[j]));
    }
  }
}

int samplePDFDUNEAtm::setupExperimentMC(int iSample) {
  int CurrErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  caf::StandardRecord* sr = new caf::StandardRecord();
  dunemc_base* duneobj = &dunemcSamples[iSample];

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
  
  Tree->SetBranchStatus("*", 1);
  Tree->SetBranchAddress("rec", &sr);

  duneobj->nEvents = Tree->GetEntries();
  duneobj->norm_s = 1;
  duneobj->pot_s = 1;

  duneobj->signal = sample_signal[iSample];
  
  duneobj->mode = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];
  
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_truecz = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
 
  for (int iEvent=0;iEvent<duneobj->nEvents;iEvent++) {
    Tree->GetEntry(iEvent);

    if ((iEvent % (duneobj->nEvents/10))==0) {
      MACH3LOG_INFO("\tProcessing event: {}/{}",iEvent,duneobj->nEvents);
    }

    duneobj->mode[iEvent] = (double)SIMBMode_ToMaCh3Mode(TMath::Abs(sr->mc.nu[0].mode),(int)sr->mc.nu[0].iscc);
    duneobj->rw_isCC[iEvent] = sr->mc.nu[0].iscc;
    duneobj->Target[iEvent] = 40;
    
    duneobj->rw_etru[iEvent] = (double)(sr->mc.nu[0].E);

    TVector3 TrueNuMomentumVector = (TVector3(sr->mc.nu[0].momentum.X(),sr->mc.nu[0].momentum.Y(),sr->mc.nu[0].momentum.Z())).Unit();
    duneobj->rw_truecz[iEvent] = TrueNuMomentumVector.Y();

    duneobj->flux_w[iEvent] = sr->mc.nu[0].genweight;

    TVector3 RecoNuMomentumVector;
    if (IsELike) {
      duneobj->rw_erec[iEvent] = sr->common.ixn.pandora[0].Enu.e_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.heshw.X(),sr->common.ixn.pandora[0].dir.heshw.Y(),sr->common.ixn.pandora[0].dir.heshw.Z())).Unit();
    } else {
      duneobj->rw_erec[iEvent] = sr->common.ixn.pandora[0].Enu.lep_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.lngtrk.X(),sr->common.ixn.pandora[0].dir.lngtrk.Y(),sr->common.ixn.pandora[0].dir.lngtrk.Z())).Unit();      
    }
    duneobj->rw_theta[iEvent] = RecoNuMomentumVector.Y();
    
  }

  delete Tree;
  delete File;

  gErrorIgnoreLevel = CurrErrorLevel;
  
  return duneobj->nEvents;
}

void samplePDFDUNEAtm::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);  

  //Make sure that this is only set if you're an atmoshperic object
  fdobj->rw_truecz.resize(fdobj->nEvents);
  
  fdobj->signal = duneobj->signal;
  fdobj->SampleDetID = SampleDetID;
  fdobj->ChannelIndex = iSample;
  
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent) {
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->isNC[iEvent] = !duneobj->rw_isCC[iEvent];
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]);
    
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->rw_truecz[iEvent] = &(duneobj->rw_truecz[iEvent]);

    fdobj->nupdgUnosc[iEvent] = &(sample_nupdgunosc[iSample]);
    fdobj->nupdg[iEvent] = &(sample_nupdg[iSample]);
  }
}

const double* samplePDFDUNEAtm::GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  double* KinematicValue;

  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iSample].rw_etru[iEvent]);
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iSample].rw_erec[iEvent]);
    break;
  case kTrueCosZ:
    KinematicValue = &(dunemcSamples[iSample].rw_truecz[iEvent]);
    break;
  case kRecoCosZ:
    KinematicValue = &(dunemcSamples[iSample].rw_theta[iEvent]);
    break;
  case kOscChannel:
    KinematicValue = &(MCSamples[iSample].ChannelIndex);
    break;
  case kMode:
    KinematicValue = MCSamples[iSample].mode[iEvent];
    break;
  default:
    MACH3LOG_ERROR("Unknown KinPar: {}",KinPar);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}

const double* samplePDFDUNEAtm::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

const double* samplePDFDUNEAtm::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEAtm::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicVariable, iSample, iEvent);
}

double samplePDFDUNEAtm::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicParameter, iSample, iEvent);
}

std::vector<double> samplePDFDUNEAtm::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameterStr));
  return ReturnKinematicParameterBinning(KinPar);
}

std::vector<double> samplePDFDUNEAtm::ReturnKinematicParameterBinning(KinematicTypes KinPar)  {
  std::vector<double> binningVector;
  
  switch (KinPar) {

  case kRecoNeutrinoEnergy:
  case kTrueNeutrinoEnergy:
    for (int i=0;i<20;i++) {
      binningVector.emplace_back((double)i);
    }
    binningVector.emplace_back(100.);
    binningVector.emplace_back(1000.);
    break;
    
  default:
    MACH3LOG_ERROR("Unknown KinPar: {}",KinPar);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return binningVector;
}

int samplePDFDUNEAtm::ReturnKinematicParameterFromString(std::string KinematicParameterStr) {
  int ReturnVal;

  if (KinematicParameterStr == "TrueNeutrinoEnergy") {ReturnVal = kTrueNeutrinoEnergy;}
  else if (KinematicParameterStr == "RecoNeutrinoEnergy") {ReturnVal = kRecoNeutrinoEnergy;}
  else if (KinematicParameterStr == "TrueCosineZ") {ReturnVal = kTrueCosZ;}
  else if (KinematicParameterStr == "RecoCosineZ") {ReturnVal = kRecoCosZ;}
  else if (KinematicParameterStr == "OscChannel") {ReturnVal = kOscChannel;}
  else if (KinematicParameterStr == "Mode") {ReturnVal = kMode;}
  else {
    MACH3LOG_ERROR("KinematicParameterStr: {} not found",KinematicParameterStr);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return ReturnVal;
}

std::string samplePDFDUNEAtm::ReturnStringFromKinematicParameter(int KinematicParameter) {
  return "";
}
