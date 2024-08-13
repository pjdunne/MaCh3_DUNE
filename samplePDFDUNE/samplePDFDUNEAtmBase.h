#ifndef _samplePDFDUNEAtmBase_h_
#define _samplePDFDUNEAtmBase_h_

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

#include "StructsDUNE.h"

class samplePDFDUNEAtmBase : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEAtmBase(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEAtmBase();

  //DB This should be removed once core-develop has the virtual function removed
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
 protected:
  void init(double pot, std::string mc_version, covarianceXsec *xsec_cov);
  void setupDUNEMC(const char *sampleInputFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats=false);
  void setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj);

  void SetupWeightPointers() override;
  
  double ReturnKinematicParameter (std::string KinematicParameter, int iSample, int iEvent) override;
  double ReturnKinematicParameter (double KinematicVariable, int iSample, int iEvent) override;
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);

  // dunemc
  std::vector<struct dunemc_base> dunemcSamples;

  TFile *_sampleFile;
  TTree *_data;
  TString _nutype;
  int _mode;

  int _ipnu;

  // DUNEMC Variables
  double _ev;
  double _erec;
  double _erec_nue;
  double _NuMomX;
  double _NuMomY;
  double _NuMomZ;
  double _LepMomX;
  double _LepMomY;
  double _LepMomZ;
  double _LepNuAngle;
  double _cvnnumu;
  double _cvnnue;
  double _BeRPA_cvwgt;
  double _maccres_cvwgt;
  double _nunpcc1_cvwgt;
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
  double _weight;
  bool iselike;
};

#endif
