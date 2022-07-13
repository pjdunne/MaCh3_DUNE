#ifndef _Structs_h_
#define _Structs_h_

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

// *******************
// Template to make vector out of an array of any length
template< typename T, size_t N >
std::vector<T> MakeVector( const T (&data)[N] ) {
// *******************
  return std::vector<T>(data, data+N);
}

// *******************
// Normalisations for cross-section parameters
class XsecNorms2 {
  // *******************
  public:
    // Bins for normalisation parameter
    TAxis *ebins;
    // Name of parameters
    std::string name;
    // Mode which parameter applies to
    int mode;
    // PDG which parameter applies to
    std::vector<int> pdgs;
    // Targets which parameter applies to
    std::vector<int> targets;
    // Parameter number of this normalisation in current NIWG parameterisation
    int startbin;
};

// *******************
// Normalisations for cross-section parameters
class XsecNorms3 {
  // *******************
 public:
  // Bins for normalisation parameter
  TAxis *ebins;
  // Name of parameters
  std::string name;
  // Mode which parameter applies to
  std::vector<int> modes;
  // Horn currents which parameter applies to
  std::vector<int> horncurrents;
  // PDG which parameter applies to
  std::vector<int> pdgs;
  // Preosc PDG which parameter applies to
  std::vector<int> preoscpdgs;
  // Targets which parameter applies to
  std::vector<int> targets;
  //Does this parameter have kinematic bounds
  bool hasKinBounds;
  //Etrue bounds
  double etru_bnd_low;
  double etru_bnd_high;
  //Etrue bounds
  double q2_true_bnd_low;
  double q2_true_bnd_high;
  // Parameter number of this normalisation in current NIWG parameterisation
  int index;
};


// *******************
// Add a struct to hold info about the splinified xsec parameters
// Used in ND280 code to accelerate CPU TSpline3 evaluation
// Adaptable for SK
struct FastSplineInfo {
// *******************
  // Number of points in spline
  int nPts;

  // Array of the knots positions
  double *xPts;

  // Array of what segment of spline we're currently interested in
  // Gets updated once per MCMC iteration
  int CurrSegment;

  //ETA trying to change things so we don't read in all flat splines
  int flat;
};

// *******************
// The kinematic types at ND280 (projections of the event)
enum ND280KinematicTypes {
  // *******************
  kLeptonMomentum = 0,
  kLeptonCosTheta = 1,
  kNeutrinoEnergy = 2,
  kQ2             = 3,
  kErecQE         = 4,
  kQ2QE           = 5,
  kLeptonTheta    = 6,
  kq0             = 7,
  kq3             = 8,
  kPionMomentum   = 9,
  kPionCosTheta   = 10
};

// *******************
// Get the name of a kinematic type
// Useful for setting x and y axis of TH2D samples and so on
inline std::string ND280Kinematic_ToString(ND280KinematicTypes type) {
  // *******************
  std::string ReturnString = "";
  switch(type) {
    case kLeptonMomentum:
      ReturnString = "p_{#mu} (MeV)";
      break;
    case kLeptonCosTheta:
      ReturnString = "cos #theta_{#mu}";
      break;
    case kNeutrinoEnergy:
      ReturnString = "E_{#nu}^{true} (GeV)";
      break;
    case kQ2:
      ReturnString = "Q^{2}_{true} (GeV^{2})";
      break;
    case kErecQE:
      ReturnString = "E_{#nu}^{QE} (GeV)";
      break;
    case kQ2QE:
      ReturnString = "Q^{2}_{QE} (GeV^{2})";
      break;
    case kLeptonTheta:
      ReturnString = "#theta_{#mu} (degrees)";
      break;
    case kq0:
      ReturnString = "|q_{0}| (Mev)";
      break;
    case kq3:
      ReturnString = "|q_{3}| (Mev)";
      break;
    default:
      std::cerr << "Error, did not find ToString conversion for " << type << std::endl;
      break;
  }
  return ReturnString;
};

// ***************************
// The base for the cross-section class
class xsecBase {
  // ***************************
  public:
    // Virtual functions
    xsecBase();
    ~xsecBase() {};
    void Print();
    // Mode
    __int__ mode;
    // Which species
    __int__ species;
    // Which normalisation parameter applies to this event
    __int__ normbin;
    // Which target the event is on
    __int__ target;
    // Q2 of event
    __float__ Q2;
    // Enu of event
    __float__ Enu;
    // One time weight to apply to event, e.g. NC1gamma
    __float__ weight;
};


// ********************************************
// NIWG 2015 struct
class xsec2015 : public xsecBase {
  // ********************************************

  public:
    // And the Struct's constructor
    xsec2015();
    ~xsec2015();
    void Print();

    // Still good for Winter 2016
    TSpline3* splMAQE;
    TSpline3* splCA5;
    TSpline3* splMARES;
    TSpline3* splBGSCLLOWPPI;
    TSpline3* splBGSCL;
    TSpline3* splBYDIS;
    TSpline3* splBYMPI;
    TSpline3* splAGKYMULT;
    TSpline3* splFEFABS;
    TSpline3* splFEFCX;
    TSpline3* splFEFQE;
    TSpline3* splFEFINEL;
    TSpline3* splFEFCXH;
    TSpline3* splFEFQEH;

    // New for Summer 2017
    TSpline3* splPDDC;
    TSpline3* splPDDO;

    // New for 2020
    TSpline3* splMECENULOW;
    TSpline3* splMECENUHIGH;
    TSpline3* splMECENUBARLOW;
    TSpline3* splMECENUBARHIGH;

    // Deprecated for Summer 2017
    TSpline3* splEBC;
    TSpline3* splEBO;
    TSpline3* splSCCV;
    TSpline3* splSCCA;

    // New for 2020
    TSpline3* spl2P2HEDEP_LOWENU;
    TSpline3* spl2P2HEDEP_HIGHENU;
    TSpline3* spl2P2HEDEP_LOWENUBAR;
    TSpline3* spl2P2HEDEP_HIGHENUBAR;
    TSpline3* splISO_BKG_LOWPPI;
	TSpline3* splDIS_BY;
	TSpline3* splMPI_BY;
	TSpline3* splMPI_AGKY_XSEC;

    // Relativistic RPA weight of event
    __float__ relRPA;
};


// ********************************************
// Generic xsec2018 class
// Can use TF1 or TSpline3 or TSpline5 here, tjoho
template <class T>
class xsec2018 : public xsecBase {
  // ********************************************
  public:
    // The light constructor
    xsec2018(__int__ NumberOfSplines) { 
      nParams = NumberOfSplines;
      Func.reserve(nParams);
      for (int i = 0; i < nParams; ++i) {
        Func[i] = NULL;
      }
    }

    // The empty constructor
    xsec2018() {
      nParams = 0;
      Func = NULL;
    };

    // The light destructor
    ~xsec2018() {
      for (int i = 0; i < nParams; ++i) {
        if (Func[i]) delete Func[i];
      }
    }

    // Get number of splines
    inline __int__ GetNumberOfParams() { return nParams; }

    // The Printer
    inline void Print() {
      xsecBase::Print();
      std::cout << "    Splines: " << std::endl;
      for (int i = 0; i < nParams; ++i) {
        if (!Func[i]) continue;
        std::cout << "    " << std::left << std::setw(25) << Func[i]->GetName() << std::endl;
      }
    }

    // Set the number of splines for this event
    inline void SetSplineNumber(__int__ NumberOfSplines) {
      nParams = NumberOfSplines;
      Func = new T[nParams];
    }

    // Get the function for the nth spline
    inline T GetFunc(__int__ nSpline) { return Func[nSpline]; }
    // Set the function for the nth spline
    inline void SetFunc(__int__ nSpline, T Function) { Func[nSpline] = Function; }
    // Eval the current variation
    inline double Eval(__int__ nSpline, __float__ variation) { 
      // Some will be NULL, check this
      if (Func[nSpline]) {
        return Func[nSpline]->Eval(variation);
      } else {
        return 1.0;
      }
    }
  private:
    // Number of parameters
    __int__ nParams;
    // The function
    T* Func;
};

// ************************
// A reduced TF1 class only
// Only saves parameters for each TF1 and how many parameters each parameter set has
class TF1_red {
// ************************
  public:
    // Empty constructor
    TF1_red() {
      length = 0;
      Par = NULL;
    }

    // Empty destructor
    ~TF1_red() {
      if (Par != NULL) {
        delete[] Par;
        Par = NULL;
      }
    }

    // The useful constructor with deep copy
    TF1_red(__int__ nSize, __float__* Array, __int__ Parameter) {
      length = nSize;
      for (int i = 0; i < length; ++i) {
        Par[i] = Array[i];
      }
      ParamNo = Parameter;
    }

    // The TF1 constructor with deep copy
    TF1_red(TF1* &Function, int Param = -1) {
      Par = NULL;
      SetFunc(Function, Param);
    }

    // Get the number
    inline std::string GetName() {
      std::stringstream ss;
      ss << ParamNo;
      return ss.str();
    }

    // Set the function
    inline void SetFunc(TF1* &Func, int Param = -1) {
      length = Func->GetNpar();
      if (Par != NULL) delete[] Par;
      Par = new __float__[length];
      for (int i = 0; i < length; ++i) {
        Par[i] = Func->GetParameter(i);
      }
      ParamNo = Param;
      delete Func;
      Func = NULL;
    }

    // Evaluate a variation
    inline double Eval(__float__ var) {
      // If we have 5 parameters we're using a fifth order polynomial
      if (length == 5) {
        return 1+Par[0]*var+Par[1]*var*var+Par[2]*var*var*var+Par[3]*var*var*var*var+Par[4]*var*var*var*var*var;
      // If we have 2 parameters we're using two linear equations
      } else if (length == 2) {
        return (var<=0)*(1+Par[0]*var)+(var>0)*(1+Par[1]*var);
      } else {
        std::cerr << "*** Error in reduced TF1 class!" << std::endl;
        std::cerr << "    Class only knows about 5th order polynomial and two superposed linear function" << std::endl;
        std::cerr << "    You have tried something else than this, which remains unimplemented" << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
    }

    // Set a parameter to a value 
    inline void SetParameter(__int__ Parameter, __float__ Value) {
      Par[Parameter] = Value;
    }

    // Get a parameter value
    double GetParameter(__int__ Parameter) {
      if (Parameter > length) {
        std::cerr << "Error: you requested parameter number " << Parameter << " but length is " << length << " parameters" << std::endl;
        throw;
        return -999.999;
      }
      return Par[Parameter];
    }

    // Set the size
    inline void SetSize(__int__ nSpline) { 
      length = nSpline;
      Par = new __float__[length];
    }
    // Get the size
    inline int GetSize() { return length; }
    inline void Print() {
      std::cout << "Printing TF1_red: " << std::endl;
      std::cout << "  ParamNo = " << ParamNo << std::endl;
      std::cout << "  Length  = " << length << std::endl;
      std::cout << "  a       = " << Par[0] << std::endl;
      std::cout << "  b       = " << Par[1] << std::endl;
      if (length == 5) {
        std::cout << "  c       = " << Par[2] << std::endl;
        std::cout << "  d       = " << Par[3] << std::endl;
        std::cout << "  e       = " << Par[4] << std::endl;
      }
    }

  private:
    // The parameters
    __float__* Par;
    __int__ length;
    // Save the parameter number this spline applies to
    __int__ ParamNo;
};

// ************************
// Reduced TSpline3 class
class TSpline3_red {
// ************************
  public:
    // Empty constructor
    TSpline3_red() {
      nPoints = 0;
      Par = NULL;
      XPos = NULL;
      YResp = NULL;
    }

    // The constructor that takes a TSpline3 pointer and copies in to memory
    TSpline3_red(TSpline3* &spline, int Param = -1) {
      Par = NULL;
      XPos = NULL;
      YResp = NULL;
      SetFunc(spline, Param);
    }
    
    // Set a function
    inline void SetFunc(TSpline3* &spline, int Param = -1) {
      nPoints = spline->GetNp();
      ParamNo = Param;
      if (Par != NULL) {
        for (int i = 0; i < nPoints; ++i) {
          delete[] Par[i];
          Par[i] = NULL;
        }
        delete[] Par;
        Par = NULL;
      }
      if (XPos != NULL) delete[] XPos;
      if (YResp != NULL) delete[] YResp;
      // Save the parameters for each knot
      Par = new __float__*[nPoints];
      // Save the positions of the knots
      XPos = new __float__[nPoints];
      // Save the y response at each knot
      YResp = new __float__[nPoints];
      for (int i = 0; i < nPoints; ++i) {
        // 3 is the size of the TSpline3 coefficients
        Par[i] = new __float__[3];
        double x, y, b, c, d = -999.99;
        spline->GetCoeff(i, x, y, b, c, d);
        XPos[i]   = x;
        YResp[i]  = y;
        Par[i][0] = b;
        Par[i][1] = c;
        Par[i][2] = d;
      }
      delete spline;
      spline = NULL;
    }

    // Empty destructor
    ~TSpline3_red() {
      for (int i = 0; i < nPoints; ++i) {
        if (Par[i] != NULL) {
          delete[] Par[i];
        }
      }
      delete[] Par;
      delete[] XPos;
      delete[] YResp;
      Par = NULL;
      XPos = YResp = NULL;
    }

    // Find the segment relevant to this variation in x
    // See root/hist/hist/src/TSpline3::FindX(double) or samplePDFND....::FindSplineSegment
    inline int FindX(double x) {
      // The segment we're interested in (klow in ROOT code)
      int segment = 0;
      int kHigh = nPoints-1;
      // If the variation is below the lowest saved spline point
      if (x <= XPos[0]){
        segment = 0;
        // If the variation is above the highest saved spline point
      } else if (x >= XPos[nPoints-1]) {
        // Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
        segment = kHigh;
        // If the variation is between the maximum and minimum, perform a binary search
      } else {
        // The top point we've got
        int kHalf = 0;
        // While there is still a difference in the points (we haven't yet found the segment)
        // This is a binary search, incrementing segment and decrementing kHalf until we've found the segment
        while (kHigh - segment > 1) {
          // Increment the half-step 
          kHalf = (segment + kHigh)/2;
          // If our variation is above the kHalf, set the segment to kHalf
          if (x > XPos[kHalf]) {
            segment = kHalf;
            // Else move kHigh down
          } else {
            kHigh = kHalf;
          }
        } // End the while: we've now done our binary search
      } // End the else: we've now found our point
      if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;
      return segment;
    }

    // Evaluate the weight from a variation
    inline double Eval(double var) {
      // Get the segment for this variation
      int segment = FindX(var);
      // The get the coefficients for this variation
      double x, y, b, c, d = -999.99;
      GetCoeff(segment, x, y, b, c, d);
      double dx = var - x;
      // Evaluate the third order polynomial
      double weight = y+dx*(b+dx*(c+d*dx));
      return weight;
    }

    // Get the number of points
    inline int GetNp() { return nPoints; }
    // Get the ith knot's x and y position
    inline void GetKnot(int i, double &xtmp, double &ytmp) { 
      xtmp = XPos[i];
      ytmp = YResp[i];
    }

    // Get the coefficient of a given segment
    inline void GetCoeff(int segment, double &x, double &y, double &b, double &c, double &d) {
      b = Par[segment][0];
      c = Par[segment][1];
      d = Par[segment][2];
      x = XPos[segment];
      y = YResp[segment];
    }

    // Get the number
    inline std::string GetName() {
      std::stringstream ss;
      ss << ParamNo;
      return ss.str();
    }

    // Make a TSpline3 from the reduced splines
    inline TSpline3* ConstructTSpline3() {
      TSpline3 *spline = new TSpline3(GetName().c_str(), XPos, YResp, nPoints);
      return spline;
    }

  private:
    // Number of points/knot in TSpline3
    __int__ nPoints;
    // Always uses a third order polynomial, so hard-code the number of coefficients in implementation
    __float__ **Par;
    // Positions of each x for each knot
    __float__ *XPos;
    // y-value for each knot
    __float__ *YResp;
    // Parameter number (which parameter is this spline for)
    __int__ ParamNo;
};


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

  // Keep a counter of the number of GENIE modes we have
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
// Convert an input NEUT mode (e.g. from psyche or SK input) to the MaCh3 modes
/*
inline int NEUTMode_ToMaCh3Mode(int NEUT_mode) {
  // **********************

  int ReturnMode = kMaCh3_nModes;

  switch (NEUT_mode) {
    // QE
    case kNEUT_CCQE:
      ReturnMode = kMaCh3_CCQE; // CCQE in MaCh3
      break;
      // MEC
    case kNEUT_2p2h:
      ReturnMode = kMaCh3_2p2h; // MEC in MaCh3
      break;
      // CC1pi+ and CC1pi0
    case kNEUT_CC1pip1p:
    case kNEUT_CC1pi0:
    case kNEUT_CC1pip1n:
      ReturnMode = kMaCh3_CC1pi; // CC1pi in MaCh3
      break;
      // CC coherent
    case kNEUT_CCcoh:
      ReturnMode = kMaCh3_CCcoh; // CC coh in MaCh3
      break;
      // CC multi-pion
    case kNEUT_CCMpi:
      ReturnMode = kMaCh3_CCMpi; // CC multi-pi in MaCh3
      break;
      // CC1gamma
    case kNEUT_CC1gam:
      // CC1eta
    case kNEUT_CC1eta:
      // CC1K
    case kNEUT_CC1K:
	  // CC1pi diffractive
	case kNEUT_CC1pidiff:
	  ReturnMode = kMaCh3_CCMisc; // ETA - splitting CCOther and CC DIS for 2019 OA
	  break;
      // CCDIS
    case kNEUT_CCDIS:
      ReturnMode = kMaCh3_CCDIS; // CC DIS in MaCh3; probably CCOther really?
      break;

      // NC1pi0 1n
    case kNEUT_NC1pi01n:
      // NC1pi0 1p
    case kNEUT_NC1pi01p:
      ReturnMode = kMaCh3_NC1pi0; // NC1pi0 in MaCh3
      break;
      // NC1pi- 1p
    case kNEUT_NC1pim1p:
      // NC1pi+ 1n
    case kNEUT_NC1pip1n:
      ReturnMode = kMaCh3_NC1pipm; // NC1pi+/- in MaCh3
      break;
      // NC coherent (1pi0)
    case kNEUT_NCcoh:
      ReturnMode = kMaCh3_NCcoh; // NC coh in MaCh3
      break;
      // All of the NCOther modes
    case kNEUT_NC1gam1n:
    case kNEUT_NC1gam1p:
      ReturnMode = kMaCh3_NC1gam; // NC 1 gamma
      break;
    case kNEUT_NCMpi:
    case kNEUT_NC1eta1n:
    case kNEUT_NC1eta1p:
    case kNEUT_NC1K01Lam:
    case kNEUT_NC1Kp1Lam:
    case kNEUT_NCDIS:
    case kNEUT_NCEL1p:
    case kNEUT_NCEL1n:
	  // NC1pi+- diffractive
	case kNEUT_NC1pidiff:
      ReturnMode = kMaCh3_NCoth; // NC other in MaCh3
      break;
    default:
      ReturnMode = kMaCh3_nModes; // Something else in MaCh3 (sand?)
      break;
  }

  return ReturnMode;
};
*/


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
// Enum to track the target material
enum TargetMat {
  // *****************
  kTarget_C = 12,
  kTarget_O = 16
};

// *****************
// Enum to track the incoming neutrino species
enum NuPDG {
  // *****************
  kNue = 12,
  kNumu = 14,
  kNumu_bar = -14,
  kNue_bar = -12
};

// *****************
// Enum for detector, used in flux->getBin
enum Detector_enum {
  // *****************
  kND280 = 0,
  kSK = 1
};

// *****************
// Enum to track xsec2015model (because NIWG like changing parameterisation!)
enum NIWGmodel_enum {
  // *****************
  kNIWG_2015a = 0,
  kNIWG_2015b = 1,
  kNIWG_2015c = 2,
  kNIWG_2017a = 3,
  kNIWG_2017b = 4,
  kNIWG_2018a = 5,
  kNIWG_2020a = 6,

  kNIWG_undefined = -1
};

// *****************
// Converted the NIWG model to a string
inline std::string NIWGmodel_ToString(NIWGmodel_enum i) {
  // *****************
  std::string name;

  switch(i) {
    case kNIWG_2015a:
      name = "NIWG_2015a";
      break;
    case kNIWG_2015b:
      name = "NIWG_2015b";
      break;
    case kNIWG_2015c:
      name = "NIWG_2015c";
      break;
    case kNIWG_2017a:
      name = "NIWG_2017a";
      break;
    case kNIWG_2017b:
      name = "NIWG_2017b";
      break;
    case kNIWG_2018a:
      name = "NIWG_2018a";
      break;
    case kNIWG_2020a:
      name = "NIWG_2020a";
      break;
    case kNIWG_undefined:
      name = "NIWG_undefined";
      break;
    default:
      name = "NIWG_undefined_DEFAULT";
      break;
  }

  return name;
}

// *****************
// Convert a MaCh3 mode to a string
/*
inline std::string MaCh3mode_ToString(MaCh3_Mode i) {
  // *****************

  std::string name = "";

  switch(i) {
    case kMaCh3_CCQE:
      name = "CCQE";
      break;
    case kMaCh3_2p2h:
      name = "2p2h";
      break;
    case kMaCh3_CC1pi:
      name = "CC1pi";
      break;
    case kMaCh3_CCcoh:
      name = "CCcoh";
      break;
    case kMaCh3_CCMpi:
      name = "CCMpi";
      break;
    case kMaCh3_CCDIS:
      name = "CCDIS";
      break;
    case kMaCh3_CCMisc: //ETA splitting CCDIS and CCOther for 2019 OA
      name = "CCMisc";
      break;
    case kMaCh3_NC1pi0:
      name = "NC1pi0";
      break;
    case kMaCh3_NC1pipm:
      name = "NC1pipm";
      break;
    case kMaCh3_NCcoh:
      name = "NCcoh";
      break;
    case kMaCh3_NCoth:
      name = "NCoth";
      break;
    case kMaCh3_NC1gam:
      name = "NC1gam";
      break;

    default:
      std::cerr << "Did not find the MaCh3 mode you specified " << i << std::endl;
      name = "UNKNOWN_BAD";
  }

  return name;
}
*/
// *****************
// Convert a MaCh3 mode to a string
/*
inline std::string MaCh3mode_ToFancyString(MaCh3_Mode i) {
  // *****************

  std::string name = "";

  switch(i) {
    case kMaCh3_CCQE:
      name = "CCQE";
      break;
    case kMaCh3_2p2h:
      name = "2p2h";
      break;
    case kMaCh3_CC1pi:
      name = "CC 1#pi^{#pm,0}";
      break;
    case kMaCh3_CCcoh:
      name = "CC coherent";
      break;
    case kMaCh3_CCMpi:
      name = "CC multi-#pi";
      break;
    case kMaCh3_CCDIS:
      name = "CC DIS";
      break;
    case kMaCh3_CCMisc: //ETA splitting CCDIS and CCOther for 2019 OA
      name = "CC Misc";
      break;
    case kMaCh3_NC1pi0:
      name = "NC 1#pi^{0}";
      break;
    case kMaCh3_NC1pipm:
      name = "NC 1#pi^{#pm}";
      break;
    case kMaCh3_NCcoh:
      name = "NC coherent";
      break;
    case kMaCh3_NCoth:
      name = "NC other";
      break;
    case kMaCh3_NC1gam:
      name = "NC 1#gamma";
      break;
    case kMaCh3_nModes:
      name = "Sand #mu";
      break;

    default:
      std::cerr << "Did not find the MaCh3 mode you specified " << i << std::endl;
      name = "UNKNOWN_BAD";
  }

  return name;
}
*/


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



// *******************************
// Convert MaCh3 mode to SK name
/*
inline std::string MaCh3mode_ToSKString(MaCh3_Mode i) {
  // *******************************
  std::string name;

  switch(i) {
    case kMaCh3_CCQE:
      name = "ccqe";
      break;
    case kMaCh3_2p2h:
      name = "mec";
      break;
    case kMaCh3_CC1pi:
      name = "cc1pi";
      break;
    case kMaCh3_CCcoh:
      name = "cccoh";
      break;
    case kMaCh3_CCMpi:
      name = "ccmpi";
      break;
    case kMaCh3_CCDIS:
      name = "ccdis";
      break;
    case kMaCh3_CCMisc: //ETA splitting CCOther and CCDIS for 2019 OA
      name = "ccmisc";
      break;
    case kMaCh3_NC1pi0:
      name = "ncpiz";
      break;
    case kMaCh3_NC1pipm:
      name = "ncpipm";
      break;
    case kMaCh3_NCcoh:
      name = "nccoh";
      break;
    case kMaCh3_NCoth:
      name = "ncoth";
      break;
    case kMaCh3_NC1gam:
      name = "nc1gamma";
      break;

    default:
      std::cerr << "Did not find the MaCh3 mode you specified " << i << std::endl;
      name = "UNKNOWN_BAD";
  }

  return name;
}
*/

// Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
double OverflowIntegral(TH2Poly*);

// Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly*);

// Poly Projectors
TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins);
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins);

// Helper to check if files exist or not
inline std::string file_exists(std::string filename) {
  std::ifstream infile(filename.c_str());
  if (!infile.good()) {
    std::cerr << "*** ERROR ***" << std::endl;
    std::cerr << "File " << filename << " does not exist" << std::endl;
    std::cerr << "Please try again" << std::endl;
    std::cerr << "*************" << std::endl;
    throw;
  }

  return filename;
}

// To calculate the binned low q2 normalisations weights, assuming Q2 in [GeV/c] unit
inline double const calcQ2norm(const double Q2, const double Q2_norm0_w, const double Q2_norm1_w, const double Q2_norm2_w, const double Q2_norm3_w, const double Q2_norm4_w, const double Q2_norm5_w,  const double Q2_norm6_w, const double Q2_norm7_w){

  if(Q2 >= 0.0 && Q2 <=0.05){
    return Q2_norm0_w;
  }else if(Q2> 0.05 && Q2 <=0.1){
    return Q2_norm1_w;
  }else if(Q2> 0.1 && Q2 <=0.15){
    return Q2_norm2_w;
  }else if(Q2> 0.15 && Q2 <=0.2){
    return Q2_norm3_w;
  }else if(Q2> 0.2 && Q2 <=0.25){
    return Q2_norm4_w;
  }else if(Q2> 0.25 && Q2 <=0.5){
    return Q2_norm5_w;
  }else if(Q2> 0.5 && Q2 <=1.0){
    return Q2_norm6_w;
  }else if(Q2> 1.0){
    return Q2_norm7_w;
  }

  return 1.0;
}

double returnSKQ2(float _pnu[50], float _dirnu[50][3], int _mode, int oscnutype);
double returnSKW2(float _pnu[50], float _dirnu[50][3], int _mode, int _ipnu[50], double Q2);
double returnSKPiPlusMom(float _pnu[50], float _ipnu[50], int numnu);

#endif
