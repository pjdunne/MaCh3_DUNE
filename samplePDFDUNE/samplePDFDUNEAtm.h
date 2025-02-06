#ifndef _samplePDFDUNEAtm_h_
#define _samplePDFDUNEAtm_h_

#include "splines/splinesDUNE.h"
#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"
#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"
/// @brief Base class for handling atmospheric samples
class samplePDFDUNEAtm : virtual public samplePDFFDBase
{
public:
  /// @brief Constructor
  /// @param mc_version Configuration file
  /// @param xsec_cov cross-section covariance matrix
  /// @param osc_cov oscillation covariance matrix
  samplePDFDUNEAtm(std::string mc_version, covarianceXsec* xsec_cov, covarianceOsc* osc_cov);
  /// @brief destructor
  ~samplePDFDUNEAtm();

  /// @brief Enum to identify kinematics
  enum KinematicTypes{kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueCosZ,kRecoCosZ,kOscChannel,kMode};
  
protected:
  /// @brief Initialises object
  void Init();

  /// @brief Function to setup MC from file
  /// @param iSample sample ID
  /// @return Total number of events
  int setupExperimentMC(int iSample);

  /// @brief Tells FD base which variables to point to/be set to
  /// @param iSample Sample ID
  void setupFDMC(int iSample);

  /// @brief Sets up pointers weights for each event (oscillation/xsec/etc.)
  void SetupWeightPointers();

  /// @brief Sets up splines 
  void SetupSplines();
  
  //DB functions which could be initialised to do something which is non-trivial
  
  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iSample sample ID
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
  
  /// @brief NOT IMPLEMENTED: Apply kinematic shifts
  /// @param iSample Sample Number
  /// @param iEvent Event number
  void applyShifts(int iSample, int iEvent) {}
  
  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinPar Kinematic parameter enum val
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  /// @brief Gets binning for a given parameter
  /// @param KinematicParameterStr Parameter name
  /// @return Vector containing parameter bins
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameterStr);

  /// @brief Gets binning for a given parameter
  /// @param KinPar Parameter ID
  /// @return Vector containing parameter bins
  std::vector<double> ReturnKinematicParameterBinning(KinematicTypes KinPar);

  /// @brief Get kinematic parameter ID from string name
  /// @param KinematicStr 
  /// @return Parameter ID
  inline int ReturnKinematicParameterFromString(std::string KinematicStr);

  /// @brief Get kinematic parameter name from ID
  /// @param KinematicVariable  Parameter ID
  /// @return Parameter Name
  inline std::string ReturnStringFromKinematicParameter(int KinematicVariable);

  /// Array filled with MC samples for each oscillation channel
  std::vector<struct dunemc_base> dunemcSamples;

  /// Is the sample e-like
  bool IsELike;
};

#endif
