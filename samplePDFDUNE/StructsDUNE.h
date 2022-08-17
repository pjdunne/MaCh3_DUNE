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


// ***************************
// A handy namespace for ND280 psyche extraction
namespace MaCh3Utils {
  // ***************************

  // ***************************
  // Return mass for given PDG
  double GetMassFromPDG(int PDG);
  // ***************************



  // Neutrino direction
  extern const double ND280NuDir[3];

  extern const double SKNuDir[3];

} // end MaCh3Utils namespace



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
// Struct to describe the NEUT interactionkWeakMix                   =   13 modes
// Note: this modes can be found in neut/src/nemodsel.F
enum NEUT_Mode {
  // ***************************

  // CCQE
  kNEUT_CCQE      = 1,
  // (Nieves) 2p2h aka MEC
  kNEUT_2p2h      = 2,
  // CC 1pi+ 1p
  kNEUT_CC1pip1p  = 11,
  // CC 1pi0
  kNEUT_CC1pi0    = 12,
  // CC 1pi+ 1n
  kNEUT_CC1pip1n  = 13,
  // CC 1pi diffractive
  kNEUT_CC1pidiff = 15,
  // CC coherent (1pi+ for neutrino)
  kNEUT_CCcoh     = 16,
  // CC 1 gamma
  kNEUT_CC1gam    = 17,
  // CC multi pi (1.3 < W < 2.0)
  kNEUT_CCMpi     = 21,
  // CC 1 eta
  kNEUT_CC1eta    = 22,
  // CC 1 Kaon
  kNEUT_CC1K      = 23,
  // CC DIS (PYTHIA/JETSET)
  kNEUT_CCDIS     = 26,

  // NC 1pi0 1n
  kNEUT_NC1pi01n  = 31,
  // NC 1pi0 1p
  kNEUT_NC1pi01p  = 32,
  // NC 1pi- 1p
  kNEUT_NC1pim1p  = 33,
  // NC 1pi+ 1n
  kNEUT_NC1pip1n  = 34,
  // NC 1pi diffractive
  kNEUT_NC1pidiff = 35,
  // NC coherent (1pi0)
  kNEUT_NCcoh     = 36,
  // NC 1 gamma 1n
  kNEUT_NC1gam1n  = 38,
  // NC 1 gamma 1p
  kNEUT_NC1gam1p  = 39,

  // NC multi pi
  kNEUT_NCMpi     = 41,
  // NC 1 eta 1n from Delta
  kNEUT_NC1eta1n  = 42,
  // NC 1 eta 1p from Delta
  kNEUT_NC1eta1p  = 43,
  // NC 1 K0 1 Lambda from Delta
  kNEUT_NC1K01Lam = 44,
  // NC 1 K+ 1 Lambda from Delta
  kNEUT_NC1Kp1Lam = 45,
  // NC DIS
  kNEUT_NCDIS     = 46,
  // NC elastic 1p
  kNEUT_NCEL1p    = 51,
  // NC elastic 1n
  kNEUT_NCEL1n    = 52,

  // Keep a counter of the number of NEUT modes we have
  kNEUT_nModes
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

// Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
double OverflowIntegral(TH2Poly*);

// Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly*);

// Poly Projectors
TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins);
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins);

#endif
