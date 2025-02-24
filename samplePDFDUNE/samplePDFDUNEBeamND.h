#ifndef _samplePDFDUNEBeamND_h_
#define _samplePDFDUNEBeamND_h_

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

class samplePDFDUNEBeamND : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamND(std::string mc_version, covarianceXsec* xsec_cov, covarianceOsc* osc_cov);
  ~samplePDFDUNEBeamND();

  enum KinematicTypes {
    kTrueNeutrinoEnergy, 
    kRecoQ, 
    kRecoNeutrinoEnergy, 
    kIsFHC, 
    kOscChannel, 
    kMode, 
    kMuonMom, 
    kMuonEnergy, 
    kRecoMuonEnergy, 
    kMuonTheta, 
    kPipMom, 
    kPipTheta, 
    kMuonEDiff, 
    kPipEDiff, 
    kEDiff };
  
 protected:
  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();
  void SetupSplines();
  
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  int ReturnKinematicParameterFromString(std::string KinematicParameterStr);
  std::string ReturnStringFromKinematicParameter(int KinematicParameter);
  
  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
  void applyShifts(int iSample, int iEvent);

  std::vector<struct dunemc_base> dunendmcSamples;

  TFile *_sampleFile;
  TTree *_data;
  TString _nutype;
  int _mode;

  double pot;

  // dunendmc Variables
  double _ev;
  double _erec;
  double _erec_lep;
  double _erec_had;
  int _reco_numu;
  int _reco_nue;

  int _nuPDG;
  int _nuPDGunosc;

  double _MuMom;
  double _PipMom;
  double _PimMom;
  double _Pi0Mom;
  double _PMom;
  double _NMom;

  double _MuMomReco;

  double _MuTheta;
  double _PipTheta;
  double _PimTheta;
  double _Pi0Theta;
  double _PTheta;
  double _NTheta;

  double _MuThetaReco;

  double _eRecoMuon;
  double _eRecoPip;
  double _eRecoPim;
  double _eRecoPi0;
  double _eRecoP;
  double _eRecoN;

  double _LepNuAngle;
  double _LepE;
  double _eP;
  double _ePip;
  double _ePim;
  double _ePi0;
  double _eN;
  double _eMuon;

  double _BeRPA_cvwgt;
  int _isCC;
  int nupdgUnosc;
  int nupdg;
  int _run;
  int _isND;
  int _isFHC;
  double _vtx_x;
  double _vtx_y;
  double _vtx_z;
  double _vtx_end_x;
  double _vtx_end_y;
  double _vtx_end_z;
  double _reco_vtx_x;
  double _reco_vtx_y;
  double _reco_vtx_z;
  double _reco_vtx_end_x;
  double _reco_vtx_end_y;
  double _reco_vtx_end_z;
  double _reco_px;
  double _reco_py;
  double _reco_pz;
  double _px;
  double _py;
  double _pz;
  double _LepTheta;
  double _Q2;
  int _reco_q;
  double _reco_pid;
  int _iscontained;

  double _E_diff;
  double _E_diff_Muon;
  double _E_diff_Pip;
  double _E_diff_Pim;
  double _E_diff_Pi0;
  double _E_diff_P;
  double _E_diff_N;
  // double* _mu_track_length_truth_true;
  // double* _mu_track_length_truth_reco;
  // double* _mu_track_length_reco;

  // configuration 
  bool IsELike;
  bool isND;
  bool IsRHC;

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
