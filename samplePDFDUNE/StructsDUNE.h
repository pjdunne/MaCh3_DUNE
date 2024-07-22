#ifndef _StructsDUNE_h_
#define _StructsDUNE_h_

// Run low or high memory versions of structs
// N.B. for 64 bit systems sizeof(float) == sizeof(double) so not a huge effect
#define __LOW__MEMORY_STRUCTS__

#ifdef __LOW_MEMORY_STRUCTS__
#define __float__ float
#define __int__ short int
#else
#define __float__ double
#define __int__ int
#endif

// Include some healthy defines for constructors
#define __BAD_DOUBLE__ -999.99
#define __BAD_INT__ -999

#define __LARGE_WEIGHT__ 100



#include <iostream>
#include <vector>
#include <TSpline.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <TH2Poly.h>
 // end MaCh3Utils namespace


enum KinematicTypes {

  kTrueNeutrinoEnergy = 0,
  kRecoNeutrinoEnergy = 1,
  kTrueQ2 = 2,
  kRecoQ2 = 3,
  kTrueQ0 = 4,
  kRecoQ0 = 5,
  kTrueXPos = 6,
  kTrueYPos = 7,
  kTrueZPos = 8,
  kPionMultiplicity = 9,
  kNRecoParticles = 10,
  kInFDV = 11,
  kTrueMinusRecoEnergy = 12,
  kTrueMinusRecoEnergyRatio = 13, 
  kNMuons = 14,
  kRecoLepEnergy = 15,
  kM3Mode = 16,
  kOscChannel = 17,
  kNKinematicParams
};

inline int ReturnKinematicParameterFromString(std::string KinematicParameterStr){

  if (KinematicParameterStr.find("TrueNeutrinoEnergy") != std::string::npos) {return kTrueNeutrinoEnergy;}
  if (KinematicParameterStr.find("RecoNeutrinoEnergy") != std::string::npos) {return kRecoNeutrinoEnergy;}
  if (KinematicParameterStr.find("TrueQ2") != std::string::npos) {return kTrueQ2;}
  if (KinematicParameterStr.find("RecoQ2") != std::string::npos) {return kRecoQ2;}
  if (KinematicParameterStr.find("TrueXPos") != std::string::npos) {return kTrueXPos;}
  if (KinematicParameterStr.find("TrueYPos") != std::string::npos) {return kTrueYPos;}
  if (KinematicParameterStr.find("TrueZPos") != std::string::npos) {return kTrueZPos;}
  if (KinematicParameterStr.find("PionMultiplicity") != std::string::npos) {return kPionMultiplicity;}
  if (KinematicParameterStr.find("NRecoParticles") != std::string::npos) {return kNRecoParticles;}
  if (KinematicParameterStr.find("InFDV") != std::string::npos) {return kInFDV;}
  if (KinematicParameterStr.find("TrueMinusRecoEnergyRatio") != std::string::npos) {return kTrueMinusRecoEnergyRatio;}
  else if (KinematicParameterStr.find("TrueMinusRecoEnergy") != std::string::npos) {return kTrueMinusRecoEnergy;}
  if (KinematicParameterStr.find("NMuons") != std::string::npos) {return kNMuons;}
  if (KinematicParameterStr.find("RecoLepEnergy") != std::string::npos) {return kRecoLepEnergy;}
  if (KinematicParameterStr.find("M3Mode") != std::string::npos) {return kM3Mode;}
  if (KinematicParameterStr.find("OscChannel") != std::string::npos) {return kOscChannel;}

  return kNKinematicParams; 
}

inline std::string ReturnKinematicParameterStringFromEnum(KinematicTypes KinematicParameter){

  if (KinematicParameter == kTrueNeutrinoEnergy) {return "TrueNeutrinoEnergy";}
  if (KinematicParameter == kRecoNeutrinoEnergy) {return "RecoNeutrinoEnergy";}
  if (KinematicParameter == kTrueQ2) {return "TrueQ2";}
  if (KinematicParameter == kRecoQ2) {return "RecoQ2";}
  if (KinematicParameter == kTrueXPos) {return "TrueXPos";}
  if (KinematicParameter == kTrueYPos) {return "TrueYPos";}
  if (KinematicParameter == kTrueZPos) {return "TrueZPos";}
  if (KinematicParameter == kPionMultiplicity) {return "PionMultiplicity";}
  if (KinematicParameter == kNRecoParticles) {return "NRecoParticles";}
  if (KinematicParameter == kInFDV) {return "InFDV";}
  if (KinematicParameter == kTrueMinusRecoEnergy) {return "TrueMinusRecoEnergy";}
  if (KinematicParameter == kTrueMinusRecoEnergyRatio) {return "TrueMinusRecoEnergyRatio";}
  if (KinematicParameter == kNMuons) {return "NMuons";}
  if (KinematicParameter == kRecoLepEnergy) {return "NRecoLepEnergy";}
  if (KinematicParameter == kM3Mode) {return "M3Mode";}
  if (KinematicParameter == kOscChannel) {return "OscChannel";}

  return "NULL"; 
}


// ********************************
// FD Detector Systematic Functions
// ********************************

// -------------------------------------------------------------------------
// Global FD Energy Scale Systematics - Essentially Calibration Uncertainty
// Don't shift muon energy since that isn't calculated calorimetrically
// -------------------------------------------------------------------------


// Total Energy Scale
inline void TotalEScaleFD(const double * par, double * erec, double erecHad, double erecLep, int isCC, int nuPDG, int nutype) {

  (*erec) += (*par) * erecHad;

  //if not true CC numu event AND reco nue event
  if ( !( isCC==1 && abs(nuPDG) == 14 ) && nutype == 1 )
  {
    (*erec) += (*par) * erecLep;
  }

}

// Total Energy Scale Sqrt
inline void TotalEScaleSqrtFD(const double * par, double * erec, double erecHad, double erecLep, int isCC, int nuPDG, int nutype) {

  (*erec) += (*par) * pow(erecHad, 0.5) * erecHad ;

  //if not true CC numu AND reco nue event
  if ( !( isCC==1 && abs(nuPDG) == 14 ) && nutype == 1 )
  {
    (*erec) += (*par) * pow(erecLep, 0.5) * erecLep;
  }

}

// Total Energy Scale Inverse Sqrt
inline void TotalEScaleInvSqrtFD(const double * par, double * erec, double erecHad, double erecLep, int isCC, int nuPDG, int nutype) {

  (*erec) += (*par) * pow(erecHad+0.1, -0.5) * erecHad ;

  //if not true CC numu AND reco nue event
  if ( !( isCC==1 && abs(nuPDG) == 14 ) && nutype == 1 )
  {
    (*erec) += (*par) * pow(erecLep+0.1, -0.5) * erecLep;
  }

}

// ---------------------------------------------------------------
// Particle-Specific Uncertainties - Essentially Particle Response
// ---------------------------------------------------------------


// ---------------------------------------------------------------
// CHARGED HADRONS
// ---------------------------------------------------------------


// Charged Hadron Energy Scale 
inline void HadEScaleFD(const double * par, double * erec, double eRecoP, double eRecoPip, double eRecoPim ) {

  // Protons + Positive Pions + Negative Pions
  double sumE = eRecoP + eRecoPip + eRecoPim;
  (*erec) += (*par) * sumE;
  
}

// Charged Hadron Energy Scale Sqrt
inline void HadEScaleSqrtFD(const double * par, double * erec, double eRecoP, double eRecoPip, double eRecoPim ) {

  // Protons + Positive Pions + Negative Pions
  double sumE = eRecoP + eRecoPip + eRecoPim;
  (*erec) += (*par) * pow(sumE, 0.5) * sumE;
  
}

// Charged Hadron Energy Scale Inv Sqrt
inline void HadEScaleInvSqrtFD(const double * par, double * erec, double eRecoP, double eRecoPip, double eRecoPim ) {

  // Protons + Positive Pions + Negative Pions
  double sumE = eRecoP + eRecoPip + eRecoPim;
  (*erec) += (*par) * pow(sumE+0.1, -0.5) * sumE;
  
}


// ---------------------------------------------------------------
// Muons
// ---------------------------------------------------------------


// Muon Energy Scale
inline void MuEScaleFD(const double * par, double * erec, double erecLep, int isCC, int nuPDG, int nutype) {

  //if true CC numu AND reco numu event
  if ( isCC==1 && abs(nuPDG) == 14 && nutype == 2 )
  {
    (*erec) += (*par) * erecLep;
  }
  
}

// Muon Energy Scale Sqrt
inline void MuEScaleSqrtFD(const double * par, double * erec, double erecLep, int isCC, int nuPDG, int nutype) {

  //if true CC numu AND reco numu event
  if ( isCC==1 && abs(nuPDG) == 14 && nutype == 2 )
  {
    (*erec) += (*par) * pow(erecLep, 0.5) * erecLep;
  }
  
}

// Muon Energy Scale Inverse Sqrt
inline void MuEScaleInvSqrtFD(const double * par, double * erec, double erecLep, int isCC, int nuPDG, int nutype) {

  //if true CC numu AND reco numu event
  if ( isCC==1 && abs(nuPDG) == 14 && nutype == 2 )
  {
    (*erec) += (*par) * pow(erecLep+0.1, -0.5) * erecLep;
  }
  
}

// ---------------------------------------------------------------
// Neutrons
// ---------------------------------------------------------------


// Neutron Energy Scale
inline void NEScaleFD(const double * par, double * erec, double eRecoN) {

  (*erec) += (*par) * eRecoN;
  
}

// Neutron Energy Scale Sqrt
inline void NEScaleSqrtFD(const double * par, double * erec, double eRecoN) {

  (*erec) += (*par) * pow(eRecoN, 0.5) * eRecoN;
  
}

// Neutron Energy Scale Inverse Sqrt
inline void NEScaleInvSqrtFD(const double * par, double * erec, double eRecoN) {

  (*erec) += (*par) * pow(eRecoN+0.1, -0.5) * eRecoN;
  
}

// ---------------------------------------------------------------
// Electromagnetic Shower
// ---------------------------------------------------------------


// Electromagnetic Shower Energy Scale
inline void EMEScaleFD(const double * par, double * erec, double eRecoPi0, double erecLep, int isCC, int nuPDG, int nutype) {

  (*erec) += (*par) * eRecoPi0;

  //if true CC nue AND reco nue event
  if ( isCC==1 && abs(nuPDG) == 12 && nutype == 1 )
  {
    (*erec) += (*par) * erecLep;
  }
 
}

// Electromagnetic Shower Energy Scale Sqrt
inline void EMEScaleSqrtFD(const double * par, double * erec, double eRecoPi0, double erecLep, int isCC, int nuPDG, int nutype) {

  (*erec) += (*par) * pow(eRecoPi0, 0.5) * eRecoPi0;

  //if true CC nue AND reco nue event
  if ( isCC==1 && abs(nuPDG) == 12 && nutype == 1 )
  {
    (*erec) += (*par) * pow(erecLep, 0.5) * erecLep;
  }
 
}

// Electromagnetic Shower Energy Scale Inverse Sqrt
inline void EMEScaleInvSqrtFD(const double * par, double * erec, double eRecoPi0, double erecLep, int isCC, int nuPDG, int nutype) {

  (*erec) += (*par) * pow(eRecoPi0+0.1, -0.5) * eRecoPi0;

  //if true CC nue AND reco nue event
  if ( isCC==1 && abs(nuPDG) == 12 && nutype == 1 )
  {
    (*erec) += (*par) * pow(erecLep+0.1, -0.5) * erecLep;
  }
 
}

// ---------------------------------------------------------------
// Resolution Uncertainties
// ---------------------------------------------------------------
 
// ---------------------------------------------------------------
// CHARGED HADRONS
// ---------------------------------------------------------------
inline void HadResFD(const double * par, double * erec, double eRecoP, double eRecoPip, double eRecoPim, double eP, double ePip, double ePim) {

  // Reco Sum: Protons + Positive Pions + Negative Pions
  double recoSum = eRecoP + eRecoPip + eRecoPim;

  // True Sum: Protons + Positive Pions + Negative Pions
  double trueSum = eP + ePip + ePim;

  (*erec) += (*par) * (trueSum - recoSum);

}

// ---------------------------------------------------------------
// Muons
// ---------------------------------------------------------------
inline void MuResFD(const double * par, double * erec, double erecLep, double LepE, int isCC, int nuPDG, int nutype) {

  //if true CC numu AND reco numu event
  if ( isCC==1 && abs(nuPDG) == 14 && nutype == 2 )
  {
    (*erec) += (*par) * (LepE - erecLep);
  }

}


// ---------------------------------------------------------------
// Neutron
// ---------------------------------------------------------------
inline void NResFD(const double * par, double * erec, double eRecoN, double eN) {

  (*erec) += (*par) * (eN - eRecoN);
  
}

// ---------------------------------------------------------------
// Electromagnetic Shower
// ---------------------------------------------------------------

inline void EMResFD(const double * par, double * erec, double eRecoPi0, double ePi0, double erecLep, double LepE, int isCC, int nuPDG, int nutype) {

  (*erec) += (*par) * (ePi0 - eRecoPi0);

  //if true CC nue AND reco nue event
  if ( isCC==1 && abs(nuPDG) == 12 && nutype == 1 )
  {
    (*erec) += (*par) * (LepE - erecLep);
  }
 
}


// ***************************
// Struct to describe the GENIE interaction modes that will be used during the fit
// Note: these modes can be found in https://github.com/GENIE-MC/Generator/blob/R-2_12_10/src/Interaction/ScatteringType.h
enum MaCh3_Mode {
  // ***************************

  // CCQE
  kMaCh3_CCQE              = 0,
  // 1 Kaon
  kMaCh3_CC_Single_Kaon     = 1,
  // CC DIS
  kMaCh3_CC_DIS            = 2,
  // CC RES
  kMaCh3_CC_RES             = 3,
  // CC COH
  kMaCh3_CC_COH             = 4,
  // CC Diffractive
  kMaCh3_CC_Diffractive     = 5,
  // CC Electron EL
  kMaCh3_CC_Nue_EL          = 6,
  // CC Inverse Muon Decay
  kMaCh3_CC_IMD             = 7,
  // CC Atmospheric Neutrino Gamma
  kMaCh3_CC_AMnuGamma       = 8,
  // CC MEC (aka 2p2h)
  kMaCh3_CC_MEC             = 9,
  // CC Coherent Elastic
  kMaCh3_CC_COHEL           = 10,
  // CC Inverse Beta Decay
  kMaCh3_CC_IBD             = 11,
  // CC GlashowRES 
  kMaCh3_CC_GlashowRES      = 12,
  // CC IMD Annihilation
  kMaCh3_CC_IMDAnnihalation = 13,
  
  // NCQE
  kMaCh3_NCQE              = 14,
  // NC DIS
  kMaCh3_NC_DIS            = 15,
  // NC RES
  kMaCh3_NC_RES             = 16,
  // NC COH
  kMaCh3_NC_COH             = 17,
  // CC Diffractive
  kMaCh3_NC_Diffractive     = 18,
  // NC Electron EL
  kMaCh3_NC_Nue_EL          = 19,
  // NC Inverse Muon Decay
  kMaCh3_NC_IMD             = 20,
  // NC Atmospheric Neutrino Gamma
  kMaCh3_NC_AMnuGamma       = 21,
  // NC MEC (aka 2p2h)
  kMaCh3_NC_MEC             = 22,
  // NC Coherent Elastic
  kMaCh3_NC_COHEL           = 23,
  // NC Inverse Beta Decay
  kMaCh3_NC_IBD             = 24,
  // NC GlashowRES 
  kMaCh3_NC_GlashowRES      = 25,
  // NC IMD Annihilation
  kMaCh3_NC_IMDAnnihalation = 26,

  // Keep a counter of the number of MaCh3 modes we have
  kMaCh3_nModes
};


// ***************************
// Struct to describe the SIMB interaction modes which are found inside the FD CAF files
// Note: this modes can be found in https://nusoft.fnal.gov/larsoft/doxsvn/html/namespacesimb.html
enum SIMB_Mode {
  // ***************************

  // Unknown
  kUnknownInteraction        =   -1,
  // QE
  kQE                        =    0,
  // Resonant
  kRes                       =    1,
  // DIS
  kDIS                       =    2,
  // Coherent
  kCoh                       =    3,
  // Coherent Elastic
  kCohElastic                =    4,
  // Electron Scattering
  kElectronScattering        =    5,
  // Inverse Muon Decay Annihliation
  kIMDAnnihilation           =    6,
  // Inverse Beta Decay
  kInverseBetaDecay          =    7, 
  // Glasgow Resonance
  kGlashowResonance          =    8,
  // Atmospheric Muon Nu Gamma
  kAMNuGamma                 =    9,
  // MEC aka 2p2h
  kMEC                       =   10,
  // Diffractive
  kDiffractive               =   11,
  //
  kEM                        =   12,
  //
  kWeakMix                   =   13,
  // Just keep a counter of the number of modes
  kSIMB_nModes
};

// ***************************
// Struct to describe the GENIE interaction modes which are found inside the ND CAF files
// Note: this modes can be found in https://github.com/GENIE-MC/Generator/blob/R-2_12_10/src/Interaction/ScatteringType.N
enum GENIE_Mode {
  // ***************************

  // Unknown
  gUnknownInteraction        =   -1,
  // QE
  gQE                        =    1,
  // Single Kaon
  gSingleKaon                =    2,
  // DIS
  gDIS                       =    3,
  // RES
  gRes                       =    4,
  // Coherent
  gCoh                       =    5,
  // Diffractive
  gDiffractive                =    6,
  // Nu-e Elastic
  gElectronScattering        =    7, 
  // Inverse Muon Decay
  gIMD                       =    8,
  // Atmospheric Muon Nu Gamma
  gAMNuGamma                 =    9,
  // MEC aka 2p2h
  gMEC                       =   10,
  // Coherent Elastic
  gCohElastic                =   11,
  // Inverse Beta Decay
  gInverseBetaDecay          =   12, 
  // Glasgow Resonance
  gGlashowResonance          =   13,
  // Inverse Muon Decay Annihliation
  gIMDAnnihilation           =   14,
  // Just keep a counter of the number of modes
  gGENIE_nModes
};

// **********************
// Convert an input SIMB mode to the MaCh3 (GENIE) modes
inline int SIMBMode_ToMaCh3Mode(int SIMB_mode, int isCC) {
  // **********************

  int ReturnMode = kMaCh3_nModes;
  
  if (isCC == 1) {
    switch (SIMB_mode) {
        //Unknown
      case kUnknownInteraction:
        ReturnMode = kMaCh3_nModes;
      	    break;
        //QE
      case kQE:
        ReturnMode = kMaCh3_CCQE; // CC QE in MaCh3 DUNE
        break;
        // DIS
      case kDIS:
        ReturnMode = kMaCh3_CC_DIS; //CC DIS in MaCh3 DUNE
        break;
        // RES
      case kRes:
        ReturnMode = kMaCh3_CC_RES; // CC RES in MaCh3 DUNE
        break;
        // Coherent
      case kCoh:
        ReturnMode = kMaCh3_CC_COH; // CC Coh in MaCh3 DUNE
        break;
        // Diffractive
      case kDiffractive:
        ReturnMode = kMaCh3_CC_Diffractive; // CC multi-pi in MaCh3
        break;
        // Nue Elastic
      case kElectronScattering:
        ReturnMode = kMaCh3_CC_Nue_EL; // CC Nue scattering in MaCh3 DUNE
        break;
        // Atmospheric Mu Gamma
      case kAMNuGamma:
        ReturnMode = kMaCh3_CC_AMnuGamma; // CC Am Nu Mu in MaCh3 DUNE
        break;
        // MEC
      case kMEC:
        ReturnMode = kMaCh3_CC_MEC; // CC MEC in MaCh3 DUNE
        break;
        // Coherent Elastic
      case kCohElastic:
        ReturnMode = kMaCh3_CC_COHEL; // CC Coherent Elastic in MaCh3 DUNE
        break;
        // Inverse Beta Decay
      case kInverseBetaDecay:
        ReturnMode = kMaCh3_CC_IBD; // CC Inverse Beta Decay in MaCh3 DUNE
        break;
        // Glashow Resonance
      case kGlashowResonance:
        ReturnMode = kMaCh3_CC_GlashowRES; // CC Glashow Reaonance in DUNE
        break;
        // IMD Annihilation
      case kIMDAnnihilation:
        ReturnMode = kMaCh3_CC_IMDAnnihalation; // CC Inverse Muon Decay Annihalation in DUNE
        break;
      default:
        ReturnMode = kMaCh3_nModes; // Something else in MaCh3 (sand?)
        break;
    }

  }


  if (isCC == 0) {
    switch (SIMB_mode) {
        //Unknown
      case kUnknownInteraction:
        ReturnMode = kMaCh3_nModes;
      	    break;
        //QE
      case kQE:
        ReturnMode = kMaCh3_NCQE; // NC QE in MaCh3 DUNE
        break;
        // DIS
      case kDIS:
        ReturnMode = kMaCh3_NC_DIS; // NC DIS in MaCh3 DUNE
        break;
        // RES
      case kRes:
        ReturnMode = kMaCh3_NC_RES; // NC RES in MaCh3 DUNE
        break;
        // Coherent
      case kCoh:
        ReturnMode = kMaCh3_NC_COH; // NC Coh in MaCh3 DUNE
        break;
        // Diffractive
      case kDiffractive:
        ReturnMode = kMaCh3_NC_Diffractive; // CC multi-pi in MaCh3
        break;
        // Nue Elastic
      case kElectronScattering:
        ReturnMode = kMaCh3_NC_Nue_EL; // NC Nue scattering in MaCh3 DUNE
        break;
        // Atmospheric Mu Gamma
      case kAMNuGamma:
        ReturnMode = kMaCh3_NC_AMnuGamma; // NC Am Nu Mu in MaCh3 DUNE
        break;
        // MEC
      case kMEC:
        ReturnMode = kMaCh3_NC_MEC; // NC MEC in MaCh3 DUNE
        break;
        // Coherent Elastic
      case kCohElastic:
        ReturnMode = kMaCh3_NC_COHEL; // NC Coherent Elastic in MaCh3 DUNE
        break;
        // Inverse Beta Decay
      case kInverseBetaDecay:
        ReturnMode = kMaCh3_NC_IBD; // Inverse Beta Decay in MaCh3 DUNE
        break;
        // Glashow Resonance
      case kGlashowResonance:
        ReturnMode = kMaCh3_NC_GlashowRES; 
        break;
        // IMD Annihilation
      case kIMDAnnihilation:
        ReturnMode = kMaCh3_NC_IMDAnnihalation; 
        break;
      default:
        ReturnMode = kMaCh3_nModes; // Something else in MaCh3 (sand?)
        break;
    }

  }

  return ReturnMode;
};

// **********************
// Convert an input SIMB mode to the MaCh3 (GENIE) modes
inline int GENIEMode_ToMaCh3Mode(int GENIE_mode, int isCC) {
  // **********************

  int ReturnMode = kMaCh3_nModes;
  
  if (isCC == 1) {
    switch (GENIE_mode) {
        //Unknown
      case gUnknownInteraction:
        ReturnMode = kMaCh3_nModes;
      	    break;
        //QE
      case gQE:
        ReturnMode = kMaCh3_CCQE; // CC QE in MaCh3 DUNE
        break;
        // DIS
      case gDIS:
        ReturnMode = kMaCh3_CC_DIS; //CC DIS in MaCh3 DUNE
        break;
        // RES
      case gRes:
        ReturnMode = kMaCh3_CC_RES; // CC RES in MaCh3 DUNE
        break;
        // Coherent
      case gCoh:
        ReturnMode = kMaCh3_CC_COH; // CC Coh in MaCh3 DUNE
        break;
        // Diffractive
      case gDiffractive:
        ReturnMode = kMaCh3_CC_Diffractive; // CC multi-pi in MaCh3
        break;
        // Nue Elastic
      case gElectronScattering:
        ReturnMode = kMaCh3_CC_Nue_EL; // CC Nue scattering in MaCh3 DUNE
        break;
        // Atmospheric Mu Gamma
      case gAMNuGamma:
        ReturnMode = kMaCh3_CC_AMnuGamma; // CC Am Nu Mu in MaCh3 DUNE
        break;
        // MEC
      case gMEC:
        ReturnMode = kMaCh3_CC_MEC; // CC MEC in MaCh3 DUNE
        break;
        // Coherent Elastic
      case gCohElastic:
        ReturnMode = kMaCh3_CC_COHEL; // CC Coherent Elastic in MaCh3 DUNE
        break;
        // Inverse Beta Decay
      case gInverseBetaDecay:
        ReturnMode = kMaCh3_CC_IBD; // CC Inverse Beta Decay in MaCh3 DUNE
        break;
        // Glashow Resonance
      case gGlashowResonance:
        ReturnMode = kMaCh3_CC_GlashowRES; // NC 1 gamma
        break;
        // IMD Annihilation
      case gIMDAnnihilation:
        ReturnMode = kMaCh3_CC_IMDAnnihalation; // NC other in MaCh3
        break;
      case gIMD:
        ReturnMode = kMaCh3_CC_IMDAnnihalation; // Stick Inverse Muon Decay into above
        break;
      case gSingleKaon:
        ReturnMode = kMaCh3_CC_IMDAnnihalation; // Stick Single Kaon into above
        break;
      default:
        std::cout<< "kMaCh3_nModes"<<kMaCh3_nModes<< " from " << GENIE_mode <<  std::endl;
        ReturnMode = kMaCh3_nModes; // Something else in MaCh3 (sand?)
        break;
    }

  }


  if (isCC == 0) {
    switch (GENIE_mode) {
        //Unknown
      case kUnknownInteraction:
        ReturnMode = kMaCh3_nModes;
      	    break;
        //QE
      case gQE:
        ReturnMode = kMaCh3_NCQE; // NC QE in MaCh3 DUNE
        break;
        // DIS
      case gDIS:
        ReturnMode = kMaCh3_NC_DIS; // NC DIS in MaCh3 DUNE
        break;
        // RES
      case gRes:
        ReturnMode = kMaCh3_NC_RES; // NC RES in MaCh3 DUNE
        break;
        // Coherent
      case gCoh:
        ReturnMode = kMaCh3_NC_COH; // NC Coh in MaCh3 DUNE
        break;
        // Diffractive
      case gDiffractive:
        ReturnMode = kMaCh3_NC_Diffractive; // CC multi-pi in MaCh3
        break;
        // Nue Elastic
      case gElectronScattering:
        ReturnMode = kMaCh3_NC_Nue_EL; // NC Nue scattering in MaCh3 DUNE
        break;
        // Atmospheric Mu Gamma
      case gAMNuGamma:
        ReturnMode = kMaCh3_NC_AMnuGamma; // NC Am Nu Mu in MaCh3 DUNE
        break;
        // MEC
      case gMEC:
        ReturnMode = kMaCh3_NC_MEC; // NC MEC in MaCh3 DUNE
        break;
        // Coherent Elastic
      case gCohElastic:
        ReturnMode = kMaCh3_NC_COHEL; // NC Coherent Elastic in MaCh3 DUNE
        break;
        // Inverse Beta Decay
      case gInverseBetaDecay:
        ReturnMode = kMaCh3_NC_IBD; // Inverse Beta Decay in MaCh3 DUNE
        break;
        // Glashow Resonance
      case gGlashowResonance:
        ReturnMode = kMaCh3_NC_GlashowRES;
        break;
        // IMD Annihilation
      case gIMDAnnihilation:
        ReturnMode = kMaCh3_NC_IMDAnnihalation;
        break;
      case gIMD:
        ReturnMode = kMaCh3_NC_IMDAnnihalation; 
        break;
      case gSingleKaon:
        ReturnMode = kMaCh3_NC_IMDAnnihalation; 
        break;
      default:
        ReturnMode = kMaCh3_nModes; // Something else in MaCh3 (sand?)
        break;
    }

  }

  return ReturnMode;
};

// *****************
// Enum for detector, used in flux->getBin
enum Detector_enum {
  // *****************
  kND280 = 0,
  kSK = 1
};

// *******************************
// Convert String to MaCh3Mode
inline int DUNEString_ToMaCh3Mode(std::string Mode){
  // *******************************
  if (Mode.find("ccqe") != std::string::npos) {return kMaCh3_CCQE;}
  if (Mode.find("ccdis") != std::string::npos) {return kMaCh3_CC_DIS;}
  if (Mode.find("ccres") != std::string::npos) {return kMaCh3_CC_RES;}
  if (Mode.find("cccohel") != std::string::npos) {return kMaCh3_CC_COHEL;}
  else if (Mode.find("cccoh") != std::string::npos) {return kMaCh3_CC_COH;}
  if (Mode.find("ccdiff") != std::string::npos) {return kMaCh3_CC_Diffractive;}
  if (Mode.find("ccnueel") != std::string::npos) {return kMaCh3_CC_Nue_EL;}
  if (Mode.find("ccmec") != std::string::npos) {return kMaCh3_CC_MEC;}
  if (Mode.find("ccIMD") != std::string::npos) {return kMaCh3_CC_IMD;}
  if (Mode.find("ccglasres") != std::string::npos) {return kMaCh3_CC_GlashowRES;}
  if (Mode.find("ccimdannihilation") != std::string::npos) {return kMaCh3_CC_IMDAnnihalation;}
  if (Mode.find("ccamnugamma") != std::string::npos) {return kMaCh3_CC_AMnuGamma;}
  if (Mode.find("ccibd") != std::string::npos) {return kMaCh3_CC_IBD;}
  if (Mode.find("ncqe") != std::string::npos) {return kMaCh3_NCQE;}
  if (Mode.find("ncdis") != std::string::npos) {return kMaCh3_NC_DIS;}
  if (Mode.find("ncres") != std::string::npos) {return kMaCh3_NC_RES;}
  if (Mode.find("nccohel") != std::string::npos) {return kMaCh3_NC_COHEL;}
  else if (Mode.find("nccoh") != std::string::npos) {return kMaCh3_NC_COH;}
  if (Mode.find("ncdiff") != std::string::npos) {return kMaCh3_NC_Diffractive;}
  if (Mode.find("ncnueel") != std::string::npos) {return kMaCh3_NC_Nue_EL;}
  if (Mode.find("ncmec") != std::string::npos) {return kMaCh3_NC_MEC;}
  if (Mode.find("ncIMD") != std::string::npos) {return kMaCh3_NC_IMD;}
  if (Mode.find("ncglasres") != std::string::npos) {return kMaCh3_NC_GlashowRES;}
  if (Mode.find("ncimdannihilation") != std::string::npos) {return kMaCh3_NC_IMDAnnihalation;}
  if (Mode.find("ncamnugamma") != std::string::npos) {return kMaCh3_NC_AMnuGamma;}
  if (Mode.find("ncibd") != std::string::npos) {return kMaCh3_NC_IBD;}
  return kMaCh3_nModes;
}
// *******************************
// Convert MaCh3 mode to SK name
inline std::string MaCh3mode_ToDUNEString(MaCh3_Mode i) {
  // *******************************
  std::string name;

  switch(i) {
    case kMaCh3_CCQE:
      name = "ccqe";
      break;
	case kMaCh3_CC_Single_Kaon:
	  name = "unknown";
	  break;
	case kMaCh3_CC_DIS:
	  name = "ccdis";
	  break;
    case kMaCh3_CC_RES:
      name = "ccres";
      break;
	case kMaCh3_CC_COH:
	  name = "cccoh";
	  break;
	case kMaCh3_CC_COHEL:
	  name = "cccohel";
	  break;
	case kMaCh3_CC_Diffractive:
	  name = "ccdiff";
	  break;
	case kMaCh3_CC_Nue_EL:
	  name = "ccnueel";
	  break;
    case kMaCh3_CC_MEC:
      name = "ccmec";
      break;
	case kMaCh3_CC_IMD:
	  name = "ccIMD";
	  break;
	case kMaCh3_CC_GlashowRES:
	  name = "ccglasres";
	  break;
	case kMaCh3_CC_IMDAnnihalation:
	  name = "ccimdannihilation";
	  break;
	case kMaCh3_CC_AMnuGamma:
	  name = "ccamnugamma";
	  break;
	case kMaCh3_CC_IBD:
	  name = "ccibd";
	  break;
    case kMaCh3_NCQE:
      name = "ncqe";
      break;
	case kMaCh3_NC_DIS:
	  name = "ncdis";
	  break;
    case kMaCh3_NC_RES:
      name = "ncres";
      break;
	case kMaCh3_NC_COH:
	  name = "nccoh";
	  break;
	case kMaCh3_NC_COHEL:
	  name = "nccohel";
	  break;
	case kMaCh3_NC_Diffractive:
	  name = "ncdiff";
	  break;
	case kMaCh3_NC_Nue_EL:
	  name = "ncnueel";
	  break;
    case kMaCh3_NC_MEC:
      name = "ncmec";
      break;
	case kMaCh3_NC_IMD:
	  name = "ncIMD";
	  break;
	case kMaCh3_NC_GlashowRES:
	  name = "ncglasres";
	  break;
	case kMaCh3_NC_IMDAnnihalation:
	  name = "ncimdannihilation";
	  break;
	case kMaCh3_NC_AMnuGamma:
	  name = "ncamnugamma";
	  break;
	case kMaCh3_NC_IBD:
	  name = "ncibd";
	  break;
	case kMaCh3_nModes:
	  name = "unknown";
	  break;
    default:
      std::cerr << "Did not find the MaCh3 mode you specified " << i << std::endl;
      name = "UNKNOWN_BAD";
  }

  return name;
}

enum MaCh3_Spline_Modes {

  // ***************************

  // CCQE
  kMaCh3_Spline_CCQE              = 0,
  // 1 Kaon
  kMaCh3_Spline_CC_Single_Kaon     = 1,
  // CC DIS
  kMaCh3_Spline_CC_DIS            = 2,
  // CC RES
  kMaCh3_Spline_CC_RES             = 3,
  // CC COH
  kMaCh3_Spline_CC_COH             = 4,
  // CC Diffractive
  kMaCh3_Spline_CC_Diffractive     = 5,
  // CC Electron EL
  kMaCh3_Spline_CC_Nue_EL          = 6,
  // CC Inverse Muon Decay
  kMaCh3_Spline_CC_IMD             = 7,
  // CC Atmospheric Neutrino Gamma
  kMaCh3_Spline_CC_AMnuGamma       = 8,
  // CC MEC (aka 2p2h)
  kMaCh3_Spline_CC_MEC             = 9,
  // CC Coherent Elastic
  kMaCh3_Spline_CC_COHEL           = 10,
  // CC Inverse Beta Decay
  kMaCh3_Spline_CC_IBD             = 11,
  // CC GlashowRES 
  kMaCh3_Spline_CC_GlashowRES      = 12,
  // CC IMD Annihilation
  kMaCh3_Spline_CC_IMDAnnihalation = 13,
  
  // NCQE
  kMaCh3_Spline_NCQE              = 14,
  // NC DIS
  kMaCh3_Spline_NC_DIS            = 15,
  // NC RES
  kMaCh3_Spline_NC_RES             = 16,
  // NC COH
  kMaCh3_Spline_NC_COH             = 17,
  // CC Diffractive
  kMaCh3_Spline_NC_Diffractive     = 18,
  // NC Electron EL
  kMaCh3_Spline_NC_Nue_EL          = 19,
  // NC Inverse Muon Decay
  kMaCh3_Spline_NC_IMD             = 20,
  // NC Atmospheric Neutrino Gamma
  kMaCh3_Spline_NC_AMnuGamma       = 21,
  // NC MEC (aka 2p2h)
  kMaCh3_Spline_NC_MEC             = 22,
  // NC Coherent Elastic
  kMaCh3_Spline_NC_COHEL           = 23,
  // NC Inverse Beta Decay
  kMaCh3_Spline_NC_IBD             = 24,
  // NC GlashowRES 
  kMaCh3_Spline_NC_GlashowRES      = 25,
  // NC IMD Annihilation
  kMaCh3_Spline_NC_IMDAnnihalation = 26,

  // Keep a counter of the number of MaCh3 modes we have
  kMaCh3_Spline_nModes

};

inline int MaCh3Mode_to_SplineMode(int iMode){

  //No grouping of modes in MaCh3
  return iMode;
}

// Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
double OverflowIntegral(TH2Poly*);

// Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly*);

// Poly Projectors
TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins);
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins);

#endif
