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
 
  // std::cout << "-------------------------------------------------------------------" <<std::endl;
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
      MCSamples[i].ntotal_weight_pointers[j] = 5;
      MCSamples[i].total_weight_pointers[j].resize(MCSamples[i].ntotal_weight_pointers[j]);
      MCSamples[i].total_weight_pointers[j][0] = &(dunendmcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunendmcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][3] = &(dunendmcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(MCSamples[i].xsec_w[j]);
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
  
  /*
  -----------------------------------------------------------------------------------------------------
  -----------------------------------------------------------------------------------------------------
  Initialising Branch Variables
  -----------------------------------------------------------------------------------------------------
  -----------------------------------------------------------------------------------------------------
  */ 
  Tree->SetBranchStatus("*", 1);
  Tree->SetBranchAddress("rec", &sr);

  if (!IsRHC) { 
    duneobj->norm_s = (1e21/1.5e21);
  } else {
    duneobj->norm_s = (1e21/1.905e21);
  }
  duneobj->pot_s = 1;

  duneobj->nEvents = Tree->GetEntries();
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_iscontained = new int[duneobj->nEvents];
  duneobj->nupdg = new int[duneobj->nEvents];
  duneobj->nupdgUnosc = new int[duneobj->nEvents];
  duneobj->mode = new double[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];
  duneobj->rw_isFHC = new bool[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];

  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_yrec = new double[duneobj->nEvents]; // Not filled in MicroProdN3p1
  duneobj->rw_LepE = new double[duneobj->nEvents]; // Reco E for now!!!!!

  duneobj->rw_eP = new double[duneobj->nEvents];
  duneobj->rw_ePip = new double[duneobj->nEvents];
  duneobj->rw_ePim = new double[duneobj->nEvents];
  duneobj->rw_ePi0 = new double[duneobj->nEvents];
  duneobj->rw_eN = new double[duneobj->nEvents];
  duneobj->rw_eMuon = new double[duneobj->nEvents];

  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_vtx_end_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_end_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_end_z = new double[duneobj->nEvents];
  duneobj->rw_px = new double[duneobj->nEvents];
  duneobj->rw_py = new double[duneobj->nEvents];
  duneobj->rw_pz = new double[duneobj->nEvents];

  duneobj->rw_MuMom = new double[duneobj->nEvents];
  duneobj->rw_PipMom = new double[duneobj->nEvents];
  duneobj->rw_PimMom = new double[duneobj->nEvents];
  duneobj->rw_Pi0Mom = new double[duneobj->nEvents];
  duneobj->rw_PMom = new double[duneobj->nEvents];
  duneobj->rw_NMom = new double[duneobj->nEvents];

  duneobj->rw_MuTheta = new double[duneobj->nEvents];
  duneobj->rw_PipTheta = new double[duneobj->nEvents];
  duneobj->rw_PimTheta = new double[duneobj->nEvents];
  duneobj->rw_Pi0Theta = new double[duneobj->nEvents];
  duneobj->rw_PTheta = new double[duneobj->nEvents];
  duneobj->rw_NTheta = new double[duneobj->nEvents];

  duneobj->rw_erec_lep = new double[duneobj->nEvents]; // Not filled in MicroProdN3p1
  duneobj->rw_erec_had = new double[duneobj->nEvents]; // Not filled in MicroProdN3p1
  duneobj->rw_erec = new double[duneobj->nEvents]; 
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->rw_reco_q = new double[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 

  duneobj->rw_eRecoP = new double[duneobj->nEvents];
  duneobj->rw_eRecoPip = new double[duneobj->nEvents];
  duneobj->rw_eRecoPim = new double[duneobj->nEvents];
  duneobj->rw_eRecoPi0 = new double[duneobj->nEvents];
  duneobj->rw_eRecoN = new double[duneobj->nEvents];
  duneobj->rw_eRecoMuon = new double[duneobj->nEvents];

  duneobj->rw_MuMomReco = new double[duneobj->nEvents];

  duneobj->rw_reco_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_end_x = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_end_y = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_end_z = new double[duneobj->nEvents];
  duneobj->rw_reco_px = new double[duneobj->nEvents];
  duneobj->rw_reco_py = new double[duneobj->nEvents];
  duneobj->rw_reco_pz = new double[duneobj->nEvents];

  duneobj->rw_reco_pid = new double[duneobj->nEvents];

  duneobj->rw_E_diff = new double[duneobj->nEvents];
  duneobj->rw_E_diff_Muon = new double[duneobj->nEvents];
  duneobj->rw_E_diff_Pip = new double[duneobj->nEvents];
  duneobj->rw_E_diff_Pim = new double[duneobj->nEvents];
  duneobj->rw_E_diff_Pi0 = new double[duneobj->nEvents];
  duneobj->rw_E_diff_P = new double[duneobj->nEvents];
  duneobj->rw_E_diff_N = new double[duneobj->nEvents]; 

  duneobj->rw_reco_px = new double[duneobj->nEvents];
  duneobj->rw_reco_py = new double[duneobj->nEvents];
  duneobj->rw_reco_pz = new double[duneobj->nEvents];
    /*
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
    Initialising Base Variables
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
    */   
  for (int iEvent=0;iEvent<duneobj->nEvents; iEvent++) { // Loop through tree
    Tree->GetEntry(iEvent);




    duneobj->rw_isCC[iEvent] = sr->mc.nu[0].iscc;
    duneobj->rw_isFHC[iEvent] = (double)(1.0);
    duneobj->nupdg[iEvent] = (int)(sr->mc.nu[0].pdg);
    duneobj->nupdgUnosc[iEvent] = (int)(sr->mc.nu[0].pdgorig);    
    duneobj->Target[iEvent] = 40;
    int mode= TMath::Abs(_mode);       
    duneobj->mode[iEvent]=(double)GENIEMode_ToMaCh3Mode(mode, _isCC);
    duneobj->flux_w[iEvent] = (double)(sr->mc.nu[0].genweight);

    /*
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
    ND-LAr Truth Information
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
    */    
    duneobj->rw_etru[iEvent] = (double)(sr->mc.nu[0].E);
   
    int nprim = sr->mc.nu[0].nprim;
    for (int i = 0; i < nprim; i++) {
      int pdg = sr->mc.nu[0].prim[i].pdg;
      double energy = (double)(sr->mc.nu[0].prim[i].p.E);
      double px = (double)(sr->mc.nu[0].prim[i].p.px);
      double py = (double)(sr->mc.nu[0].prim[i].p.py);
      double pz = (double)(sr->mc.nu[0].prim[i].p.pz);
      double momentum = sqrt(px * px + py * py + pz * pz);
      double theta = acosl(pz / momentum);

      switch (pdg) {
      case 2212:
        duneobj->rw_eP[iEvent] = energy;
        duneobj->rw_PMom[iEvent] = momentum;
        duneobj->rw_PTheta[iEvent] = theta;
        break;
      case 211:
        duneobj->rw_ePip[iEvent] = energy;
        duneobj->rw_PipMom[iEvent] = momentum;
        duneobj->rw_PipTheta[iEvent] = theta;
        break;
      case -211:
        duneobj->rw_ePim[iEvent] = energy;
        duneobj->rw_PimMom[iEvent] = momentum;
        duneobj->rw_PimTheta[iEvent] = theta;
        break;
      case 111:
        duneobj->rw_ePi0[iEvent] = energy;
        duneobj->rw_Pi0Mom[iEvent] = momentum;
        duneobj->rw_Pi0Theta[iEvent] = theta;
        break;
      case 2112:
        duneobj->rw_eN[iEvent] = energy;
        duneobj->rw_NMom[iEvent] = momentum;
        duneobj->rw_NTheta[iEvent] = theta;
        break;
      case 13:
        duneobj->rw_eMuon[iEvent] = energy;
        duneobj->rw_MuMom[iEvent] = momentum;
        duneobj->rw_MuTheta[iEvent] = theta;
        break;
      }
    }

    /*
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
    ND-LAr Reco Information
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
    */
  //  int ndlar_ndlp = sr->nd.lar.ndlp;
  //  for (int i = 0; i < ndlar_ndlp; i++) {
  //    int ntracks = sr->nd.lar.dlp[i].ntracks;
  //    for (int j = 0; j < ntracks; j++) {
  //      erec_total+=sr->nd.lar.dlp[i].tracks[j].E;
  //    }
  //    int nshowers = sr->nd.lar.dlp[i].nshowers;
  //    for (int j = 0; j < nshowers; j++) {
  //      erec_total+=sr->nd.lar.dlp[i].showers[j].Evis;
  //    }
  //    if (duneobj->rw_iscontained[iEvent] && duneobj->rw_isCC[iEvent]){
  //      duneobj->rw_erec[iEvent] = erec_total;
  //      duneobj->rw_E_diff[iEvent] = duneobj->rw_etru[iEvent]-duneobj->rw_erec[iEvent];
  //    }
  //    else {
  //      duneobj->rw_E_diff[iEvent] = -9999;
  //      duneobj->rw_erec[iEvent] = -9999;
  //    }
  //  }

    double erec_total = 0.0;
    int nmuons =0;
    int common_ndlp = sr->common.ixn.ndlp;
    for(int i=0; i<common_ndlp; i++) {
      int part_ndlp = sr->common.ixn.dlp[i].part.ndlp;
      for(int j=0; j<part_ndlp; j++) {
        //Total E_rec
        // if(sr->common.ixn.dlp[i].part.dlp[j].E > 50.0){
        //   erec_total += (double)(0);
        //   std::cout << "PDG: " << sr->common.ixn.dlp[i].part.dlp[j].pdg << ", E: " << sr->common.ixn.dlp[i].part.dlp[j].E << std::endl;
        // } // Fixing a Kaon bug for MicroProdN3p1 
        duneobj->rw_reco_vtx_x[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].start.X());
        duneobj->rw_reco_vtx_y[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].start.Y());
        duneobj->rw_reco_vtx_z[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].start.Z());
        
        duneobj->rw_reco_vtx_end_x[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].end.X());
        duneobj->rw_reco_vtx_end_y[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].end.Y());
        duneobj->rw_reco_vtx_end_z[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].end.Z());

        duneobj->rw_reco_px[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].p.x);
        duneobj->rw_reco_py[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].p.y);
        duneobj->rw_reco_pz[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].p.z);
        
        duneobj->rw_reco_pid[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].pdg);
        duneobj->rw_iscontained[iEvent] = (double)(sr->common.ixn.dlp[i].part.dlp[j].contained);

        std::map<int, double*> particle_energy_map = {
          {2212, &duneobj->rw_eRecoP[iEvent]},
          {211, &duneobj->rw_eRecoPip[iEvent]},
          {-211, &duneobj->rw_eRecoPim[iEvent]},
          {111, &duneobj->rw_eRecoPi0[iEvent]},
          {2112, &duneobj->rw_eRecoN[iEvent]},
          {13, &duneobj->rw_eRecoMuon[iEvent]}
          // {-13, &duneobj->rw_LepE[iEvent]} //How to do muon +/-
        };

        auto it = particle_energy_map.find(sr->common.ixn.dlp[i].part.dlp[j].pdg);
        if (it != particle_energy_map.end()) {
          if (duneobj->rw_iscontained[iEvent] && duneobj->rw_isCC[iEvent]) {
            *(it->second) += (double)(sr->common.ixn.dlp[i].part.dlp[j].E);
          } else {
            *(it->second) += -9999;
          }
        }
        
      }
    }

    std::vector<std::pair<double, double>> energy_pairs = {
      {duneobj->rw_eMuon[iEvent], duneobj->rw_eRecoMuon[iEvent]},
      {duneobj->rw_ePip[iEvent], duneobj->rw_eRecoPip[iEvent]},
      {duneobj->rw_ePim[iEvent], duneobj->rw_eRecoPim[iEvent]},
      {duneobj->rw_ePi0[iEvent], duneobj->rw_eRecoPi0[iEvent]},
      {duneobj->rw_eP[iEvent], duneobj->rw_eRecoP[iEvent]},
      {duneobj->rw_eN[iEvent], duneobj->rw_eRecoN[iEvent]}
    };
    std::vector<double*> energy_diffs = {
      &duneobj->rw_E_diff_Muon[iEvent],
      &duneobj->rw_E_diff_Pip[iEvent],
      &duneobj->rw_E_diff_Pim[iEvent],
      &duneobj->rw_E_diff_Pi0[iEvent],
      &duneobj->rw_E_diff_P[iEvent],
      &duneobj->rw_E_diff_N[iEvent]
    };
    for (size_t i = 0; i < energy_pairs.size(); ++i) {
      if (energy_pairs[i].first == 0 || energy_pairs[i].second == 0) {
      *energy_diffs[i] = -9999;
      } else {
      *energy_diffs[i] = energy_pairs[i].first - energy_pairs[i].second;
      }
    }


  }
  
  delete Tree;
  delete File;

  gErrorIgnoreLevel = CurrErrorLevel;

  return duneobj->nEvents;
}


const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  double* KinematicValue;
  static double tempValue;

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
    tempValue = static_cast<double>(dunendmcSamples[iSample].rw_isFHC[iEvent]);
    KinematicValue = &tempValue;
    break;
  case kOscChannel:
    tempValue = static_cast<double>(dunendmcSamples[iSample].nupdgUnosc[iEvent]);
    KinematicValue = &tempValue;
    break;
  case kMode:
    KinematicValue = &dunendmcSamples[iSample].mode[iEvent];
    break;
  case kMuonMom:
    KinematicValue = &dunendmcSamples[iSample].rw_MuMom[iEvent];
    break;
  case kMuonEnergy:
    KinematicValue = &dunendmcSamples[iSample].rw_eMuon[iEvent];
    break;
  case kRecoMuonEnergy:
    KinematicValue = &dunendmcSamples[iSample].rw_eRecoMuon[iEvent];
    break;
  case kMuonTheta:
    KinematicValue = &dunendmcSamples[iSample].rw_MuTheta[iEvent];
    break;
  case kPipMom:
    KinematicValue = &dunendmcSamples[iSample].rw_PipMom[iEvent];
    break;
  case kPipTheta:
    KinematicValue = &dunendmcSamples[iSample].rw_PipTheta[iEvent];
    break;
  case kMuonEDiff:
    KinematicValue = &dunendmcSamples[iSample].rw_E_diff_Muon[iEvent];
    break;
  case kPipEDiff:
    KinematicValue = &dunendmcSamples[iSample].rw_E_diff_Pip[iEvent];
    break;
  case kEDiff:
    KinematicValue = &dunendmcSamples[iSample].rw_E_diff[iEvent];
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
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
    fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
    fdobj->nupdgUnosc[iEvent] = &(duneobj->nupdgUnosc[iEvent]);
    fdobj->nupdg[iEvent] = &(duneobj->nupdg[iEvent]);

    // std::cout << "This is Event "<<iEvent<< std::endl;
    // std::cout << "rw_etru:" <<duneobj->rw_etru[iEvent]<< std::endl;
    // std::cout << "rw_erec:" <<duneobj->rw_erec[iEvent]<< std::endl;
    // int common_ndlp = duneobj->common_ndlp[iEvent];
    // std::cout << "rw_reco_vtx_x:" <<duneobj->rw_reco_vtx_x[iEvent]<< std::endl;
    // std::cout << "rw_vtx_x:" <<duneobj->rw_vtx_x[iEvent]<< std::endl;    
    // std::cout << "rw_reco_vtx_y:" <<duneobj->rw_reco_vtx_y[iEvent]<< std::endl;
    // std::cout << "rw_vtx_y:" <<duneobj->rw_vtx_y[iEvent]<< std::endl;
    // std::cout << "rw_reco_vtx_z:" <<duneobj->rw_reco_vtx_z[iEvent]<< std::endl;
    // std::cout << "rw_vtx_z:" <<duneobj->rw_vtx_z[iEvent]<< std::endl;
    // std::cout << "rw_eRecoP:" <<duneobj->rw_eRecoP[iEvent]<< std::endl;
    // std::cout << "rw_eP:" <<duneobj->rw_eP[iEvent]<< std::endl;
    // std::cout << "rw_eRecoPip:" <<duneobj->rw_eRecoPip[iEvent]<< std::endl;
    // std::cout << "rw_ePip:" <<duneobj->rw_ePip[iEvent]<< std::endl;
    // std::cout << "rw_eRecoPim:" <<duneobj->rw_eRecoPim[iEvent]<< std::endl;
    // std::cout << "rw_ePim:" <<duneobj->rw_ePim[iEvent]<< std::endl;
    // std::cout << "rw_eRecoPi0:" <<duneobj->rw_eRecoPi0[iEvent]<< std::endl;
    // std::cout << "rw_ePi0:" <<duneobj->rw_ePi0[iEvent]<< std::endl;
    // std::cout << "rw_eRecoN:" <<duneobj->rw_eRecoN[iEvent]<< std::endl;
    // std::cout << "rw_eN:" <<duneobj->rw_eN[iEvent]<< std::endl;
    // std::cout << "mode:" <<duneobj->mode[iEvent]<< std::endl;
    // std::cout << "Target:" <<duneobj->Target[iEvent]<< std::endl;
    // std::cout << "isCC:" <<duneobj->rw_isCC[iEvent]<< std::endl;
    // std::cout << "isNC:" <<fdobj->isNC[iEvent]<< std::endl;
    // for(int i = 0; i < common_ndlp; i++) {
    //   int part_ndlp = duneobj->part_ndlp[iEvent];
    //   for(int j = 0; j < part_ndlp; j++) {
    //     std::cout << "rw_reco_px:" <<duneobj->rw_reco_px[iEvent][i][j]<< std::endl;
    //     std::cout << "rw_reco_py:" <<duneobj->rw_reco_py[iEvent][i][j]<< std::endl;
    //     std::cout << "rw_reco_pz:" <<duneobj->rw_reco_pz[iEvent][i][j]<< std::endl;
    //   }
    // }
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
  if(KinematicParameterStr == "IsFHC") return kIsFHC;
  if(KinematicParameterStr == "OscChannel") return kOscChannel;
  if(KinematicParameterStr == "Mode") return kMode;
  if(KinematicParameterStr == "MuonMom") return kMuonMom;
  if(KinematicParameterStr == "MuonEnergy") return kMuonEnergy;
  if(KinematicParameterStr == "MuonTheta") return kMuonTheta;
  if (KinematicParameterStr == "RecoMuonEnergy") return kRecoMuonEnergy;
  if(KinematicParameterStr == "PipMom") return kPipMom;
  if(KinematicParameterStr == "PipTheta") return kPipTheta;
  if(KinematicParameterStr == "MuonEDiff") return kMuonEDiff;
  if(KinematicParameterStr == "PipEDiff") return kPipEDiff;
  if(KinematicParameterStr == "EDiff") return kEDiff;

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
    case kMuonMom: return "MuonMom";
    case kMuonEnergy: return "MuonEnergy";
    case kRecoMuonEnergy: return "RecoMuonEnergy";
    case kMuonTheta: return "MuonTheta";
    case kPipMom: return "PipMom";
    case kPipTheta: return "PipTheta";
    case kMuonEDiff: return "MuonEDiff";
    case kPipEDiff: return "PipEDiff";
    case kEDiff: return "EDiff";
    default: return "";
  }
}
