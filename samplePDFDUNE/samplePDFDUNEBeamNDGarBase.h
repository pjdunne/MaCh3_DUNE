#ifndef _samplePDFDUNEBeamNDGarBase_h_
#define _samplePDFDUNEBeamNDGarBase_h_

#include <iostream>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TGraph2DErrors.h>
#include <vector>
#include <omp.h>
#include <list>

#include "splines/splinesDUNE.h"
#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"
#include "samplePDF/samplePDFFDBase.h"
#include "StandardRecord.h"

#include "StructsDUNE.h"

class samplePDFDUNEBeamNDGarBase : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamNDGarBase(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBeamNDGarBase();

  //DB This should be removed once core-develop has the virtual function removed
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
 protected:
  void Init();
  void setupDUNEMC(const char *sampleInputFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats=false);
  void setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj);

  void SetupWeightPointers();

  double ReturnKinematicParameter (KinematicTypes KinematicParameter, int i, int j);
  double ReturnKinematicParameter (std::string KinematicParameter, int iSample, int iEvent);
  double ReturnKinematicParameter (double KinematicVariable, int iSample, int iEvent);
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  std::vector<double> ReturnKinematicParameterBinning(KinematicTypes KinematicParameter);

  // dunendmc
  std::vector<struct dunemc_base> dunendgarmcSamples;

  TFile *_sampleFile;
  TTree *_data;
  TString _nutype;
  int _mode;

  // dunendgarmc Variables
  double _ev;
  double _erec;
  double _erec_nue;
  double _elep_reco;
  double _LepNuAngle;
  int _reco_numu;
  int _reco_nue;
  double _BeRPA_cvwgt = 1;
  int _isCC;
  int _nuPDGunosc;
  int _nuPDG;
  int _run;
  int _isND;
  int _isFHC;
  double _vtx_x;
  double _vtx_y;
  double _vtx_z;
  double _LepTheta;
  double _Q2;

  bool iscalo_reco; //NK Added so we can easily change what energy reconstruction we are using
  float muonscore_threshold; //NK Added so we can optimise muon threshold

  caf::StandardRecord* sr = new caf::StandardRecord();
};

#endif
