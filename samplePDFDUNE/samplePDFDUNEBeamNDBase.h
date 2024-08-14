#ifndef _samplePDFDUNEBeamNDBase_h_
#define _samplePDFDUNEBeamNDBase_h_

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

class samplePDFDUNEBeamNDBase : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamNDBase(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBeamNDBase();

  //DB This should be removed once core-develop has the virtual function removed
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
 protected:
  void Init();
  void setupDUNEMC(const char *sampleInputFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats=false);
  void setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj);

  void SetupWeightPointers();

  double ReturnKinematicParameter (std::string KinematicParameter, int iSample, int iEvent);
  double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);

  //Apply shifts from functional parameters
  void applyShifts(int iSample, int iEvent);

  std::vector<struct dunemc_base> dunendmcSamples;

  TFile *_sampleFile;
  TTree *_data;
  TString _nutype;
  int _mode;

  // dunendmc Variables
  double _ev;
  double _erec;
  double _erec_lep;
  double _erec_had;
  int _reco_numu;
  int _reco_nue;

  double _eRecoP;
  double _eRecoPip;
  double _eRecoPim;
  double _eRecoPi0;
  double _eRecoN;

  double _LepNuAngle;
  double _LepE;
  double _eP;
  double _ePip;
  double _ePim;
  double _ePi0;
  double _eN;
  double _BeRPA_cvwgt;
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
  int _reco_q;

  // configuration 
  bool iselike;
  bool isND;

  //Positions of ND Detector systematics
  double tot_escale_nd_pos;
  double tot_escale_sqrt_nd_pos;
  double tot_escale_invsqrt_nd_pos;
  double had_escale_nd_pos;
  double had_escale_sqrt_nd_pos;
  double had_escale_invsqrt_nd_pos;
  double mu_escale_nd_pos;
  double mu_escale_sqrt_nd_pos;
  double mu_escale_invsqrt_nd_pos;
  double n_escale_nd_pos;
  double n_escale_sqrt_nd_pos;
  double n_escale_invsqrt_nd_pos;
  double em_escale_nd_pos;
  double em_escale_sqrt_nd_pos;
  double em_escale_invsqrt_nd_pos;
  double had_res_nd_pos;
  double mu_res_nd_pos;
  double n_res_nd_pos;
  double em_res_nd_pos;

  std::vector<const double*> NDDetectorSystPointers;
  int nNDDetectorSystPointers;
};

#endif
