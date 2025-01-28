#ifndef _samplePDFDUNEBaseNDGAr_h_
#define _samplePDFDUNEBaseNDGAr_h_

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
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "StructsDUNE.h"


//constructors are same for all three so put in here
struct dunendgarmc_base {

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
  double *rw_erec;
  double *rw_etru;
  double *rw_etrurec;
  double *rw_etrurec_nopionthreshold;
  double *rw_elep_reco;
  double *rw_mom;
  double *rw_theta;
  double *rw_Q2;
  double *rw_yrec;
  double *rw_pdf_bin_1d; //global bin in the PDF so we can do a faster fill1Dhist 
  double *rw_lower_erec_1d; // lower to check if Eb has moved the erec bin
  double *rw_upper_erec_1d; // upper to check if Eb has moved the erec bin
  double *rw_pdf_bin_2d; //global bin in the PDF so we can do a faster fill2Dhist 
  double *rw_lower_erec_2d; // lower to check if Eb has moved the erec bin
  double *rw_upper_erec_2d; // upper to check if Eb has moved the erec bin
  int *rw_reco_nue;
  int *rw_reco_numu;
  double *rw_berpaacvwgt;
  int    *rw_isCC;
  int    *rw_nuPDGunosc;
  int    *rw_nuPDG;
  int    *rw_run;
  bool    *rw_isFHC;
  double *rw_vtx_x;
  double *rw_vtx_y;
  double *rw_vtx_z;

  double *rw_W;
  double *rw_Q0;
  double *rw_Q3;
  double *rw_reco_Q2;
  double *rw_reco_Q0;
  double *rw_reco_Q3;

  double dummy_y;

  double *reco_numu;

  float **rw_dirlep;
  int *mode;
  int *isbound;
  int **rw_ipnu;

  int *nproton; ///< number of (post-FSI) primary protons
  int *nneutron; ///< number of (post-FSI) primary neutrons
  int *npip; ///< number of (post-FSI) primary pi+
  int *npim; ///< number of (post-FSI) primary pi-
  int *npi0; ///< number of (post-FSI) primary pi0

  int *nrecopion; //number of reconstructed pion
  int *ntruemuon; //number of true muons
  int *ntruemuonprim; //number of true primary muons
  int *nrecomuon; //number of reconstructed muons
  double *nmuonsratio; //number of reco muons divided by number of true muons
  bool *isnumuCC; //is the true interaction numuCC

  double *rw_lep_pT;  //transverse most energetic lepton momentum
  double *rw_lep_pZ; //parallel most energetic lepton momentum
  double *rw_lep_pY;
  double *rw_lep_pX;
  double *rw_lep_pMag;
  double *rw_reco_lep_pT;  //transverse most energetic reconstructed lepton momentum
  double *rw_reco_lep_pZ; //parallel most energetic reconstructed lepton momentum
  double *rw_reco_lep_pY;
  double *rw_reco_lep_pX;
  double *rw_reco_lep_pMag;

  double *rw_pi_pT;  //transverse most energetic pion momentum
  double *rw_pi_pZ; //parallel most energetic pion momentum
  double *rw_pi_pY;
  double *rw_pi_pX;
  double *rw_pi_pMag;
  double *rw_reco_pi_pT;  //transverse reconstructed most energetic momentum
  double *rw_reco_pi_pZ; //parallel reconstructed most energetic pion momentum
  double *rw_reco_pi_pY;
  double *rw_reco_pi_pX;
  double *rw_reco_pi_pMag;
  
  double *rejectedpart_ptot;
  double *rejectedpart_pT;
  double *rejectedpart_theta_angle;
  double *rejectedpart_radcurvature;
  double *rejectedpart_sigmamom;
  double *rejectedpart_sigmatheta;
  double *highestpart_pT;
  double *highestpart_theta_angle;
  double *highestpart_lengthtrackx;
  double *highestpart_lengthtrackyz;
  double *rejectedpart_track_theta_angle;
  double *rejectedpart_ratioradcurvature;
  double *rejectedpart_beta;

  double *rw_reco_vtx_x;
  double *rw_reco_vtx_y;
  double *rw_reco_vtx_z;
  double *rw_reco_rad;
  double *rw_rad;

  double *rw_reco_lep_energy;
  double *rw_lep_energy;

  double *rw_reco_pi_energy;
  double *muon_pi_reco_angle;
  double *pi_z_reco_angle;
  double *rw_pi_energy;
  double *muon_pi_angle;
  double *pi_z_angle;
  double *muon_z_angle;
  double *rw_pi_min_energy; 

  double *rw_elep_true;

  int *nrecoparticles;
  bool *in_fdv;

  bool *is_accepted;
  double *momres_nonaccepted;
  int *pdg_nonaccepted;

  double *ecaldepositfraction;
  int *particleevent;
  int *particlepdg;
  double *particleenergy;
  double *particletheta;
  double *particlededx;
  double *particlemomentum;
  double pot_s; // s is for scale                                             
  double norm_s;//    

  double *beam_w;
  double *flux_w; // not the same as beam systematics weight!                 
  double *xsec_w;
  double *energyscale_w;
  //float *relRPA_w;

};

class samplePDFDUNEBaseNDGAr : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBaseNDGAr(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBaseNDGAr();

  void printPosteriors();

  void setupSplines(fdmc_base *fdobj, const char *splineFile, int nutype, int signal);

  caf::StandardRecord* sr = new caf::StandardRecord();

 protected:
  void init(double pot, std::string mc_version, covarianceXsec *xsec_cov);
  void setupDUNEMC(const char *sampleInputFile, dunendgarmc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats=false);
  void setupFDMC(dunendgarmc_base *duneobj, fdmc_base *fdobj, const char *splineFile);

  void setupWeightPointers();

  void makePixelGrid(float pixel_spacing_cm);
  double FindNHits(float pixel_spacing_cm, float centre_circle_y, float centre_circle_z, double rad_curvature);
  double CalcDeDx(double beta, double bg, double gamma);
  double CalcBeta(double p_mag, double& bg, double& gamma);
  void IsParticleAccepted(dunendgarmc_base *duneobj, int& i_truepart, int& i, int& isnotaccepted, double& highestpT, float pixel_spacing_cm, int& tot_particles);

  // oscillation: Prob3++ 
  TH1D *modes;
  
  bool osc_binned;
  // an axis to set binned oscillation weights
  TAxis *osc_binned_axis ;

  //Generic Function applying shifts
  double CalcXsecWeightFunc(int iSample, int iEvent);
  //double CalcFuncSystWeight(int iSample, int iEvent);
  double ReturnKinematicParameter (KinematicTypes KinematicParameter, int i, int j);
  double ReturnKinematicParameter (std::string KinematicParameter, int iSample, int iEvent);
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  std::vector<double> ReturnKinematicParameterBinning(KinematicTypes KinematicParameter);

  //Likelihood
  double getCovLikelihood();
  double getDiscVar(int sample , int event , int varindx);

  int getNMCSamples(){return dunendgarmcSamples.size();};
  int getNEventsInSample(int sample);


  // dunendmc
  std::vector<struct dunendgarmc_base> dunendgarmcSamples;

  //
  TH1D* get1DVarHist(std::string KinematicVar1, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis); 
  TH1D* get1DVarHist(std::string KinematicVar1,std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis); 
  TH2D* get2DVarHist(std::string KinematicVar1, std::string KinematicVar2, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis, TAxis* Axis2); 
  TH2D* get2DVarHist(std::string KinematicVar1,std::string KinematicVar2,std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis, TAxis* Axis2); 

  TFile *_sampleFile;
  TTree *_data;
  TString _nutype;

  TFile *_sampleFile_geant;
  TTree *_data_geant;

  int _mode;
  float _pnu[50];
  double _wgtflx;
  double _wgtosc;
  int _ipnu;


  //Standard Record object
//  caf::StandardRecord* sr = new caf::StandardRecord();
  
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

  //NK - adding geant vectors
  std::vector<float> *_MCPStartX=0;
  std::vector<float> *_MCPStartY=0;
  std::vector<float> *_MCPStartZ=0;
  std::vector<float> *_MCPEndX=0;
  std::vector<float> *_MCPEndY=0;
  std::vector<float> *_MCPEndZ=0;
  std::vector<float> *_MCPStartPX=0;
  std::vector<float> *_MCPStartPY=0;
  std::vector<float> *_MCPStartPZ=0;
  std::vector<float> *_MCPEndPX=0;
  std::vector<float> *_MCPEndPY=0;
  std::vector<float> *_MCPEndPZ=0;
  std::vector<int> *_PDG = 0;
  std::vector<int> *_MCPTrkID=0;
  std::vector<int> *_SimHitTrkID=0;
  std::vector<float> *_SimHitEnergy=0;

  double pdgmass;
  //particle masses in GeV
  double m_chargedpi = 0.13957039;
  double m_pi0 = 0.1349768;
  double m_e = 0.00051099895;
  double m_mu = 0.1056583755;
  double m_p = 0.93827208816;
  double m_n = 0.9395654205;
  double m_chargedk = 0.493677;

  double TPCFidLength;
  double TPCFidRadius;
  double TPCInstrumentedLength;
  double TPCInstrumentedRadius;
  double ECALInnerRadius;
  double ECALOuterRadius;
  double ECALEndCapStart;
  double ECALEndCapEnd;

  double TPC_centre_x =0.;
  double TPC_centre_y = -150.;
  double TPC_centre_z = 1486.;
  double K_const = 0.307075; //4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol)
  double sternheimer_A = 0.1956;
  double sternheimer_K = 3.0000;
  double sternheimer_X0 = 0.2000;
  double sternheimer_X1 = 3.0000;
  double sternheimer_Cbar = 5.2146;
  double excitationenergy = 188.0; //excitation energy for electrons in argon gas in eV
  double density = 0.0167; //in g/cm^3

  double X0 = 1193; //in cm From Federico's Kalman Filter Paper
//  double density = -0.00615*294.26 + 1.928; //in g/cm^3
  //covarianceFlux *flux;
  //covarianceSkDet_joint *skdet_joint;

  //pixel vars
  double pixelymin;
  double pixelymax;
  double pixelzmin;
  double pixelzmax;
  std::vector<double> yboundarypositions;
  std::vector<double> zboundarypositions;
  // configuration 
  bool iselike;
  bool isND;
  bool iscc1pi;

  bool incl_geant; //NK - Added so we can use GArAnaTrees
  bool iscalo_reco; //NK Added so we can easily change what energy reconstruction we are using
  bool ecal_containment; //NK Do we count containment if the particle stops in the ECAL?
  float muonscore_threshold; //NK Added so we can optimise muon threshold
  float protondEdxscore;
  float protontofscore;
  float recovertexradiusthreshold;
  float pionenergy_threshold; //NK Added so we can find pion energy threshold
  float B_field;
  float momentum_resolution_threshold;
  float pixel_spacing;
  float spatial_resolution;
  float adc_sampling_frequency;
  float drift_velocity;
//  float hits_per_mm;
  double average_gain; //in electrons per ADC
  float pi0_reco_efficiency;  //efficiency for pi0 reco in ECAL 
  float gamma_reco_efficiency;  //efficiency for gamma reco in ECAL
 
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
  void calcXsecNormsBins(dunendgarmc_base* dunendgarmc);
//  void makePixelGrid(float pixel_spacing_cm);
//  double FindNHits(float pixel_spacing_cm, float centre_circle_y, float centre_circle_z, double rad_curvature);
//  double CalcDeDx(double pdgmass, double beta, double bg, double gamma);
//  double CalcBeta(double pdgmass, double p_mag, double& bg, double& gamma);
//  void IsParticleAccepted(dunendgarmc_base *duneobj, int& i_truepart, int& i, int& isnotaccepted, double& highestpT, double pdgmass, float pixel_spacing_cm);
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
