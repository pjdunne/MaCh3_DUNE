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


 //CAFAna binned Oscillation 
std::vector<double> get_default_CAFana_bins(){
  // From CAFana - probability binning -
  const int kNumTrueEnergyBins = 100;

  // N+1 bin low edges
  std::vector<double> edges(kNumTrueEnergyBins+1);

  const double Emin = 0.5; // 500 MeV: there's really no events below there

  // How many edges to generate. Allow room for 0-Emin bi            const double N = kNumTrueEnergyBins-1;
  const double N = kNumTrueEnergyBins-1;
  const double A = N*Emin;

  edges[0] = 0;

  for(int i = 1; i <= N; ++i){
	edges[kNumTrueEnergyBins-i] = A/i;
  }

  edges[kNumTrueEnergyBins] = 120; // Replace the infinity that would be here
  return edges;

}

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
  kTrueCosZ = 9,
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
  if (KinematicParameterStr.find("TrueCosZ") != std::string::npos) {return kTrueCosZ;}

  return kNKinematicParams; 
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
  kUnknownInteraction      =   -1,
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
        ReturnMode = kMaCh3_CC_GlashowRES; // NC 1 gamma
        break;
        // IMD Annihilation
      case kIMDAnnihilation:
        ReturnMode = kMaCh3_CC_IMDAnnihalation; // NC other in MaCh3
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
        ReturnMode = kMaCh3_NC_GlashowRES; // NC 1 gamma
        break;
        // IMD Annihilation
      case kIMDAnnihilation:
        ReturnMode = kMaCh3_NC_IMDAnnihalation; // NC other in MaCh3
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


int MaCh3Mode_to_SplineMode(int iMode){

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
