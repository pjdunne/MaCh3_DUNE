//
// A parallel implementation of samplePDFDUNEBase
//

#ifndef _samplePDFDUNEBaseGPU_h_
#define _samplePDFDUNEBaseGPU_h_

#include <iostream>
//#include <TTree.h>
//#include <TH1D.h>
//#include <TH2D.h>
//#include <TMath.h>
//#include <TFile.h>
#include <vector>
#include "samplePDFDUNEBase.h"
#include "Structs.h"
#include "../Prob3++/BargerPropagator.h"

class samplePDFDUNEBaseGPU : public samplePDFDUNEBase
{
 public:
  samplePDFDUNEBaseGPU(double pot, std::string mc_version, covarianceXsec* xsec_cov = NULL);
  ~samplePDFDUNEBaseGPU();

  void reweight(double *oscpar);
  void reweight(double *oscpar_nub, double *oscpar_nu);                                                                

 protected:
  void calcOscWeightsGPU(double *en, double *w, int num, int nutype, int oscnutype, double *oscpar);
  void calcOscWeightsGPU(double *en, double *w, int num, int nutype, int oscnutype, double *oscpar_nub, double *oscpar_nu);

};


#endif
