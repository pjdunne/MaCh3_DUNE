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


//constructors are same for all three so put in here
struct duneatmmc_base {

  int nutype;            
  int oscnutype;    
  bool signal; // true if signue                                              
  int nEvents; // how many MC events are there                              
  int *Target; //Target the interaction was on
  
  int _ipnu;
  
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
  int    *rw_isFHC;
  double *rw_vtx_x;
  double *rw_vtx_y;
  double *rw_vtx_z;
  double dummy_y;

  float **rw_dirlep;
  int *mode;
  int *isbound;
  int **rw_ipnu;

  TH1D *modes;

  double pot_s; // s is for scale                                             
  double norm_s;//    

  double *beam_w;
  double *flux_w; // not the same as beam systematics weight!                 
  double *xsec_w;
  double *energyscale_w;
  //float *relRPA_w;

  double *rw_truecz;

};

class samplePDFDUNEAtmBase : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEAtmBase(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEAtmBase();

  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
 protected:
  void init(double pot, std::string mc_version, covarianceXsec *xsec_cov);
  void setupDUNEMC(const char *sampleInputFile, duneatmmc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats=false);
  void setupFDMC(duneatmmc_base *duneobj, fdmc_base *fdobj);

  void SetupWeightPointers() override;
  
  double ReturnKinematicParameter (std::string KinematicParameter, int iSample, int iEvent) override;
  double ReturnKinematicParameter (double KinematicVariable, int iSample, int iEvent) override;
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);

  inline int getNMCSamples() {return dunemcSamples.size();}
  inline int getNEventsInSample(int sample) {return dunemcSamples[sample].nEvents;}
  
  // dunemc
  std::vector<struct duneatmmc_base> dunemcSamples;

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


  // configuration 
  bool iselike;
  bool iscc1pi;
};

#endif
