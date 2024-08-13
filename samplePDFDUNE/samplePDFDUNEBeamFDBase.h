#ifndef _samplePDFDUNEBeamFDBase_h_
#define _samplePDFDUNEBeamFDBase_h_

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

class samplePDFDUNEBeamFDBase : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamFDBase(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBeamFDBase();

  //DB This should be removed once core-develop has the virtual function removed
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
 protected:
  void init(double pot, std::string mc_version, covarianceXsec *xsec_cov);
  void setupDUNEMC(const char *sampleInputFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats=false);
  void setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj, const char *splineFile);

  void setupWeightPointers();

  double ReturnKinematicParameter (std::string KinematicParameter, int iSample, int iEvent);
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);

  //Apply shifts from functional parameters
  void applyShifts(int iSample, int iEvent);

  // dunemc
  std::vector<struct dunemc_base> dunemcSamples;

  TFile *_sampleFile;
  TTree *_data;
  TString _nutype;
  int _mode;

  //Reco Variables
  double _erec;
  double _erec_nue;
  double _erec_had;
  double _erec_had_nue;
  double _erec_lep;
  double _erec_lep_nue;

  double _eRecoP;
  double _eRecoPip;
  double _eRecoPim;
  double _eRecoPi0;
  double _eRecoN;

  double _cvnnumu;
  double _cvnnue;
  double _vtx_x;
  double _vtx_y;
  double _vtx_z;

  //Truth Variables
  double _ev;
  double _LepE;
  double _eP;
  double _ePip;
  double _ePim;
  double _ePi0;
  double _eN;
  double _NuMomX;
  double _NuMomY;
  double _NuMomZ;
  double _LepMomX;
  double _LepMomY;
  double _LepMomZ;
  double _LepNuAngle;
  double _BeRPA_cvwgt;
  double _maccres_cvwgt;
  double _nunpcc1_cvwgt;
  int _isCC;
  int _nuPDGunosc;
  int _nuPDG;
  int _run;
  int _isFHC;
  double _LepTheta;
  double _Q2;

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
};

#endif
