#ifndef _samplePDFDUNEAtm_h_
#define _samplePDFDUNEAtm_h_

#include "splines/splinesDUNE.h"
#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"

class samplePDFDUNEAtm : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEAtm(std::string mc_version, covarianceXsec* xsec_cov, covarianceOsc* osc_cov);
  ~samplePDFDUNEAtm();

  enum KinematicTypes{kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueCosZ,kRecoCosZ,kOscChannel,kMode};
  
protected:
  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();
  void SetupSplines();
  
  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.; (void)iSample; (void)iEvent;}
  void applyShifts(int iSample, int iEvent) {(void)iSample; (void)iEvent;}
  
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameterStr);
  std::vector<double> ReturnKinematicParameterBinning(KinematicTypes KinPar);

  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    {"TrueNeutrinoEnergy",kTrueNeutrinoEnergy},
    {"RecoNeutrinoEnergy",kRecoNeutrinoEnergy},
    {"TrueCosineZ",kTrueCosZ},
    {"RecoCosineZ",kRecoCosZ},
    {"OscillationChannel",kOscChannel},
    {"Mode",kMode}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kRecoNeutrinoEnergy,"RecoNeutrinoEnergy"},
    {kTrueCosZ,"TrueCosineZ"},    
    {kRecoCosZ,"RecoCosineZ"},
    {kOscChannel,"OscillationChannel"},
    {kMode,"Mode"}
  };
  
  std::vector<struct dunemc_base> dunemcSamples;
  bool IsELike;
};

#endif
