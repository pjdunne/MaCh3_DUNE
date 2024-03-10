#ifndef _samplePDFDUNEBaseND_h_
#define _samplePDFDUNEBaseND_h_

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


//constructors are same for all three so put in here
struct dunendmc_base {

  int nutype;            
  int oscnutype;    
  bool signal; // true if signue                                              
  int nEvents; // how many MC events are there                              
  int *Target; //Target the interaction was on
  
  // spline bins                                                              
  unsigned int *enu_s_bin;
  unsigned int *erec_s_bin;
  unsigned int *yrec_s_bin;

  // flux bin
  int *flux_bin;

  // xsec bins  
  std::list< int > *xsec_norms_bins;

  // oscillation: CUDAprob3
#if defined (CUDA)
  //cudaprob3::BeamCudaPropagatorSingle *kNu;

#else
  //cudaprob3::BeamCpuPropagator<double> *kNu; 
#endif

  // histo pdf bins
  float *rw_erec;
  float *rw_erec_shifted;
  float *rw_erec_had;
  float *rw_erec_lep;
  float *rw_yrec;

  float *rw_eRecoP;
  float *rw_eRecoPip;
  float *rw_eRecoPim;
  float *rw_eRecoPi0;
  float *rw_eRecoN;

  float *rw_LepE;
  float *rw_eP;
  float *rw_ePip;
  float *rw_ePim;
  float *rw_ePi0;
  float *rw_eN;

  float *rw_etru;
  float *rw_mom;
  float *rw_theta;
  float *rw_Q2;

  float *rw_pdf_bin_1d; //global bin in the PDF so we can do a faster fill1Dhist 
  float *rw_lower_erec_1d; // lower to check if Eb has moved the erec bin
  float *rw_upper_erec_1d; // upper to check if Eb has moved the erec bin
  float *rw_pdf_bin_2d; //global bin in the PDF so we can do a faster fill2Dhist 
  float *rw_lower_erec_2d; // lower to check if Eb has moved the erec bin
  float *rw_upper_erec_2d; // upper to check if Eb has moved the erec bin
  int *rw_reco_nue;
  int *rw_reco_numu;
  float *rw_berpaacvwgt;
  int    *rw_isCC;
  int    *rw_nuPDGunosc;
  int    *rw_nuPDG;
  int    *rw_run;
  bool    *rw_isFHC;
  float *rw_vtx_x;
  float *rw_vtx_y;
  float *rw_vtx_z;
  int *rw_reco_q;
  float dummy_y;

  float *reco_numu;

  float **rw_dirlep;
  int *mode;
  int *isbound;
  int **rw_ipnu;


  float pot_s; // s is for scale                                             
  float norm_s;//    

  float *beam_w;
  float *flux_w; // not the same as beam systematics weight!                 
  float *xsec_w;
  float *energyscale_w;
  //float *relRPA_w;


};

class samplePDFDUNEBaseND : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBaseND(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBaseND();

  void printPosteriors();

  void setupSplines(fdmc_base *fdobj, const char *splineFile, int nutype, int signal);

  double GetLikelihood();

  void setNDCovMatrix(TMatrixDSym *cov);
 protected:
  void init(double pot, std::string mc_version, covarianceXsec *xsec_cov);
  void setupDUNEMC(const char *sampleInputFile, dunendmc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats=false);
  void setupFDMC(dunendmc_base *duneobj, fdmc_base *fdobj, const char *splineFile);

  void setupWeightPointers();

  // oscillation: Prob3++ 
  TH1D *modes;
  
  bool osc_binned;
  // an axis to set binned oscillation weights
  TAxis *osc_binned_axis ;

  //Generic Function applying shifts
  double CalcXsecWeightFunc(int iSample, int iEvent);
  //double CalcFuncSystWeight(int iSample, int iEvent);
  //double ReturnKinematicParameter (KinematicTypes Var, int i, int j);
  double ReturnKinematicParameter (std::string KinematicParameter, int iSample, int iEvent);
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);

  //Likelihood
  double getCovLikelihood();
  double getDiscVar(int sample , int event , int varindx);

  int getNMCSamples();
  int getNEventsInSample(int sample);

  //Apply shifts from functional parameters
  void applyShifts(int iSample, int iEvent);

  // dunendmc
  std::vector<struct dunendmc_base> dunendmcSamples;


  TFile *_sampleFile;
  TTree *_data;
  TString _nutype;
  int _mode;
  double _pnu[50];
  double _wgtflx;
  double _wgtosc;
  int _ipnu;

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



  //covarianceFlux *flux;
  //covarianceSkDet_joint *skdet_joint;
  
  // configuration 
  bool iselike;
  bool isND;
  bool iscc1pi;

  // Note: the following 4 variables shouldn't be used any more! (From 14/1/2015 - KD). Just kept in for backwards compatibility in compiling, but they have no effect.
  bool do_flux_rw;
  bool do_xsec_rw;
  bool do_det_rw;
  
  //stuff for xsec norms  
  int nxsec_norm_modes;
  int *xsec_norm_modes;
  std::vector<std::vector<int> > xsec_norm_pdgs;
  std::vector<std::vector<int> > xsec_norm_targets;
  int *xsec_norm_startbin;
  std::vector<std::string> xsec_norm_names;
  //std::vector<XsecNorms3> xsec_norms;

  //for nuebar appearance
  double Beta;
  bool useBeta;
  bool applyBetaNue;
  bool applyBetaDiag;

  vector< double > rwFracs;
  vector< int > ids;

  // Using the smarter xsec covariance matrix reader?
  bool DoItSmart;

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

  // The ND detector covariance matrix
  TMatrixDSym *NDCovMatrix;
  
  //Temp ND detector covariance with added stats error
  TMatrixDSym *tempNDCovMatrix;

  // The inverse ND detector covariance matrix
  TMatrixDSym *NDInvCovMatrix;


  double **NDInvertCovMatrix;

  std::vector<const double*> NDDetectorSystPointers;
  int nNDDetectorSystPointers;

  // Parameters used in the DoItSmart xsec setup
  // The old version hard-codes these instead
  std::vector<int> normParsModes;
  std::vector< std::vector<int> > normParsTarget;
  std::vector< std::vector<int> > normParsPDG;
  std::vector<int> normParsIndex;
  //int nFuncParams;
  //std::vector<int> funcParsIndex;
  //std::vector<std::string> funcParsNames;
  // spline params (don't technically need, just included for print-screen sanity check purposes)
  int nSplineParams;
  std::vector<int> splineParsIndex;
  std::vector<std::string> splineParsNames;

 private:
  void calcXsecNormsBins(dunendmc_base* dunendmc);
  //void calcHistoBins();
  //TH1D *_hPDF1Dtest;
};

#endif

//init same for all should just take a config with file names

//think about splines

//setcov functions should be doable

//calcoscweights all same

//reweight all the same

//getcovlikelihood 
