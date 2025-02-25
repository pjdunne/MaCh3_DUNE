#ifndef _samplePDFDUNEBeamFD_h_
#define _samplePDFDUNEBeamFD_h_

#include "splines/splinesDUNE.h"
#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"

class samplePDFDUNEBeamFD : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamFD(std::string mc_version, covarianceXsec* xsec_cov, covarianceOsc* osc_cov);
  ~samplePDFDUNEBeamFD();

  enum KinematicTypes {kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueXPos,kTrueYPos,kTrueZPos,kCVNNumu,kCVNNue,kM3Mode,kOscChannel};
  
 protected:
  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();
  void SetupSplines();
  
  double ReturnKinematicParameter (double KinematicVariable, int iSample, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent); 

  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  inline std::string ReturnStringFromKinematicParameter(int KinematicParameterStr);
  
  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.; (void)iSample; (void)iEvent;}
  void applyShifts(int iSample, int iEvent);

  // dunemc
  std::vector<struct dunemc_base> dunemcSamples;

  double pot;
  bool iselike;

  //Positions of FD Detector systematics
  double tot_escale_fd_pos;
  double tot_escale_sqrt_fd_pos;
  double tot_escale_invsqrt_fd_pos;
  double had_escale_fd_pos;
  double had_escale_sqrt_fd_pos;
  double had_escale_invsqrt_fd_pos;
  double mu_escale_fd_pos;
  double mu_escale_sqrt_fd_pos;
  double mu_escale_invsqrt_fd_pos;
  double n_escale_fd_pos;
  double n_escale_sqrt_fd_pos;
  double n_escale_invsqrt_fd_pos;
  double em_escale_fd_pos;
  double em_escale_sqrt_fd_pos;
  double em_escale_invsqrt_fd_pos;
  double had_res_fd_pos;
  double mu_res_fd_pos;
  double n_res_fd_pos;
  double em_res_fd_pos;
  double cvn_numu_fd_pos;
  double cvn_nue_fd_pos;
  
  std::vector<const double*> FDDetectorSystPointers;
  int nFDDetectorSystPointers;

  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    {"TrueNeutrinoEnergy",kTrueNeutrinoEnergy},
    {"RecoNeutrinoEnergy",kRecoNeutrinoEnergy},
    {"TrueXPos",kTrueXPos},
    {"TrueYPos",kTrueYPos},
    {"TrueZPos",kTrueZPos},
    {"CVNNumu",kCVNNumu},
    {"CVNNue",kCVNNue},
    {"Mode",kM3Mode},
    {"OscillationChannel",kOscChannel}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kRecoNeutrinoEnergy,"RecoNeutrinoEnergy"},
    {kTrueXPos,"TrueXPos"},
    {kTrueYPos,"TrueYPos"},
    {kTrueZPos,"TrueZPos"},
    {kCVNNumu,"CVNNumu"},
    {kCVNNue,"CVNNue"},
    {kM3Mode,"Mode"},
    {kOscChannel,"OscillationChannel"}
  };
};

#endif
