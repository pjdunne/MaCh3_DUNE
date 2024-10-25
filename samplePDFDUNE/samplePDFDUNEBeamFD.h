#ifndef _samplePDFDUNEBeamFD_h_
#define _samplePDFDUNEBeamFD_h_

#include "StructsDUNE.h"

#include "splines/splinesDUNE.h"

#include "samplePDF/samplePDFFDBase.h"

class samplePDFDUNEBeamFD : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamFD(std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBeamFD();

  enum KinematicTypes {
    kTrueNeutrinoEnergy = 0,
    kRecoNeutrinoEnergy,
    kTrueXPos,
    kTrueYPos,
    kTrueZPos,
    kCVNNumu,
    kCVNNue,
    kM3Mode,
    kGlobalBinNumber,
    kOscChannel,
    kq0,
    kq3,
    kERecQE,
    kELepRec,
    kEHadRec,
    kNDefaultProjections
  };

  //More robust getters to make plots in different variables, mode, osc channel, systematic weighting and with bin range 
  TH1D* get1DVarHist(KinematicTypes Var1, int fModeToFill=-1, int fSampleToFill=-1, int WeightStyle=0, TAxis* Axis=0);
  TH1D* get1DVarHist(KinematicTypes Var1, std::vector< std::vector<double> > Selection, int WeightStyle=0, TAxis* Axis=0);
    

  struct UserProjection {
    std::string name;
    std::function<double(dunemc_base const&, int)> proj;
  };

  // add a user projection, capture the return value if you don't want to look the 
  // KinematiceParameter value up by string
  static int AddProjection(std::string const &, std::function<double(dunemc_base const&, int)>);

 protected:

  static std::vector<UserProjection> user_projections;

  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();
  
  /// @todo extract this completely to core
  ///@brief Setup our spline file, this calls InitialseSplineObject() under the hood
  void SetupSplines();

  double const &ReturnKinematicParameterByReference(int KinematicParameter,
                                                    int iSample, int iEvent);
  double ReturnKinematicParameter(int KinematicParameter, int iSample,
                                  int iEvent);

  std::vector<double> ReturnKinematicParameterBinning(int KinematicParameter);
  
  int ReturnKinematicParameterFromString(std::string KinematicParameterStr);
  std::string ReturnStringFromKinematicParameter(int KinematicParameterStr);

  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
  void applyShifts(int iSample, int iEvent);

  // dunemc
  std::vector<dunemc_base> dunemcSamples;

  double pot;
  bool iselike;

};

#endif
