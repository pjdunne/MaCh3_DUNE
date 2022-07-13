#ifndef _samplePDFDUNEBase_h_
#define _samplePDFDUNEBase_h_

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
//#include "CUDAProb3/beamcudapropagator.cuh"
//#include "CUDAProb3/atmoscudapropagator.cuh"
#include "splines/splineFDBase.h"
#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"
#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"



//constructors are same for all three so put in here
struct dunemc_base {
  int nutype;            
  int oscnutype;    
  bool signal; // true if signue                                              
  int nEvents; // how many MC events are there                              
  
  // spline bins                                                              
  unsigned int *enu_s_bin;
  unsigned int *erec_s_bin;

  // Global bin number for Eb variation templates
  int *Eb_bin;

  // flux bin
  int *flux_bin;

  //CAFAna binned Oscillation 
  std::vector<double> get_default_CAFana_bins(){
     // From CAFana - probability binning -
     const int kNumTrueEnergyBins = 100;
  
     // N+1 bin low edges
     std::vector<double> edges(kNumTrueEnergyBins+1);
  
     const double Emin = 0.5; // 500 MeV: there's really no events below there
  
     // How many edges to generate. Allow room for 0-Emin bi            const double N = kNumTrueEnergyBins-1;
     const double N = kNumTrueEnergyBins-1;
     const double A = N*Emin;
  
     edges[0] = 0;
  
     for(int i = 1; i <= N; ++i){
       edges[kNumTrueEnergyBins-i] = A/i;
     }
  
     edges[kNumTrueEnergyBins] = 120; // Replace the infinity that would be here
     return edges;
                                 
    }
  

  // xsec bins  
  std::list< int > *xsec_norms_bins;

  // oscillation: CUDAprob3
#if defined (CUDA)
  //cudaprob3::BeamCudaPropagatorSingle *kNu;

#else
  //cudaprob3::BeamCpuPropagator<double> *kNu; 
#endif
  // histo pdf bins
  double *rw_erec;
  double *rw_etru;
  double *rw_mom;
  double *rw_theta;
  double *rw_Q2;
  double *rw_pdf_bin_1d; //global bin in the PDF so we can do a faster fill1Dhist 
  double *rw_lower_erec_1d; // lower to check if Eb has moved the erec bin
  double *rw_upper_erec_1d; // upper to check if Eb has moved the erec bin
  double *rw_pdf_bin_2d; //global bin in the PDF so we can do a faster fill2Dhist 
  double *rw_lower_erec_2d; // lower to check if Eb has moved the erec bin
  double *rw_upper_erec_2d; // upper to check if Eb has moved the erec bin
  double *rw_cvnnumu;
  double *rw_cvnnue;
  double *rw_berpaacvwgt;
  int    *rw_isCC;
  int    *rw_nuPDGunosc;
  int    *rw_nuPDG;
  int    *rw_run;
  int    *rw_isFD;
  int    *rw_isFHC;
  double *rw_vtx_x;
  double *rw_vtx_y;
  double *rw_vtx_z;
  double dummy_y;

  float **rw_dirlep;
  int *mode;
  int *isbound;
  int **rw_ipnu;
  // ---------- Pi+ rweighting test only (commented out by default) ----------- //
  // See http://www.t2k.org/asg/oagroup/meeting/2017/2017-02-07-banff-oa-pre-meeting/SKpionreweighting/
  /* int *rw_numnu;
  double **rw_pnu;
  int **rw_ipnu; */
  // ---------------------- End: Pi+ rweighting test only --------------------- //


  double pot_s; // s is for scale                                             
  double norm_s;//    

  //double *osc_w_point;                                                      
  double *osc_w; // oscillation weight                                        
  double *beam_w;
  double *flux_w; // not the same as beam systematics weight!                 
  double *skdet_w;
  double *xsec_w;
  double *energyscale_w;
  //float *relRPA_w;

  splineBase *splineFile; 

};

class samplePDFDUNEBase : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBase(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBase();

  void printPosteriors();

  //DUNE FD FV cut
  inline bool IsInFDFV(double pos_x_cm, double pos_y_cm, double pos_z_cm) {
    return (abs(pos_x_cm) < 310 && abs(pos_y_cm) < 550 && pos_z_cm > 50 &&
      pos_z_cm < 1244);
      	      }
 protected:
  void init(double pot, std::string mc_version, covarianceXsec *xsec_cov);
  void setupDUNEMC(const char *sampleInputFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats=false);
  void setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj, const char *splineFile);

  // oscillation: Prob3++ 
  TH1D *modes;
  
  bool osc_binned;
  // an axis to set binned oscillation weights
  TAxis *osc_binned_axis ;
  double calcOscWeights(int nutype, int oscnutype, double en, double *oscpar);
  double calcOscWeights(int nutype, int oscnutype, double en, double *oscpar_nub, double *oscpar_nu);
  void calcOscWeightsKProb(int sample, int nutype, int oscnutype,  double *w, double *oscpar_nub, double *oscpar_nu);
  void calcOscWeightsKProb(int sample, int nutype, int oscnutype,  double *w, double *oscpar_nu);
  void calcOscWeightsGPU(double *en, double *w, int num, int nutype, int oscnutype, double *oscpar);
  void calcOscWeightsGPU(double *en, double *w, int num, int nutype, int oscnutype, double *oscpar_nub, double *oscpar_nu);

  //Generic Function applying shifts
  double calcFuncSystWeight(int iSample, int iEvent);
  double ReturnKinematicParameter (KinematicTypes Var, int i, int j);
  std::vector<double> ReturnKinematicParameterBinning(KinematicTypes Var);


  //Likelihood
  double getCovLikelihood();
  double getDiscVar(int sample , int event , int varindx);

  int getNMCSamples();
  int getNEventsInSample(int sample);


  // dunemc
  vector<struct dunemc_base> dunemcSamples;

  //fdmc
  vector<struct fdmc_base> fdmcSamples;

  TFile *_sampleFile;
  TTree *_data;
  TString _nutype;
  int _mode;
  float _pnu[50];
  double _wgtflx;
  double _wgtosc;
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
  int _isFD;
  int _isFHC;
  double _vtx_x;
  double _vtx_y;
  double _vtx_z;
  double _LepTheta;
  double _Q2;


  //covarianceFlux *flux;
  //covarianceSkDet_joint *skdet_joint;
  
  // configuration 
  bool iselike;
  bool iscc1pi;

  // Note: the following 4 variables shouldn't be used any more! (From 14/1/2015 - KD). Just kept in for backwards compatibility in compiling, but they have no effect.
  bool do_flux_rw;
  bool do_xsec_rw;
  bool do_det_rw;

  bool IsRHC; // If true, this instance of the class is RHC
  int SampleDetID;
  
  // option for binning scheme (0 = 1D erec, 1 = 2D p-theta, 2 = 2D erec-theta, 3 = 2D erec-Q2)  
  int BinningOpt;

  //stuff for xsec norms  
  int nxsec_norm_modes;
  int *xsec_norm_modes;
  std::vector<std::vector<int> > xsec_norm_pdgs;
  std::vector<std::vector<int> > xsec_norm_targets;
  int *xsec_norm_startbin;
  std::vector<std::string> xsec_norm_names;
  std::vector<XsecNorms3> xsec_norms;

  //for nuebar appearance
  double Beta;
  bool useBeta;
  bool applyBetaNue;
  bool applyBetaDiag;

  vector< double > rwFracs;
  vector< int > ids;

  // Using the smarter xsec covariance matrix reader?
  bool DoItSmart;


  // Parameters used in the DoItSmart xsec setup
  // The old version hard-codes these instead
  std::vector<int> normParsModes;
  std::vector< std::vector<int> > normParsTarget;
  std::vector< std::vector<int> > normParsPDG;
  std::vector<int> normParsIndex;
  int nFuncParams;
  std::vector<int> funcParsIndex;
  std::vector<std::string> funcParsNames;
  // spline params (don't technically need, just included for print-screen sanity check purposes)
  int nSplineParams;
  std::vector<int> splineParsIndex;
  std::vector<std::string> splineParsNames;


 private:
  void calcXsecNormsBins(dunemc_base* dunemc);
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
