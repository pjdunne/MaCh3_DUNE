#include "samplePDFDUNEBaseGPU.h"
#include <assert.h>

extern "C" void setMNS(double x12, double x13, double x23, double m21, double m23, double Delta, bool kSquared);
extern "C" void GetProb(int Alpha, int Beta, double Path, double Density, double *Energy, int n, double *oscw);

samplePDFDUNEBaseGPU::samplePDFDUNEBaseGPU(double pot, std::string mc_version, covarianceXsec* xsec_cov)
  : samplePDFDUNEBase(pot, mc_version, xsec_cov)
{
}

samplePDFDUNEBaseGPU::~samplePDFDUNEBaseGPU()
{
}

// ******************************************************************
// Calculate the oscillation weights on the GPU
// Pass the oscillation parameters (setMNS)
// Calculate the probability (GetProb)
//
// This is only for neutrino
//
// en     = energy
// w      = weights
// num    =
// nutype = neutrino flavour
// oscnutype = 
void samplePDFDUNEBaseGPU::calcOscWeightsGPU(double *en, double *w, int num, int nutype, int oscnutype, double *oscpar)
{
  if (nutype < 0) // if antinu
    {
      setMNS(oscpar[0], oscpar[2], oscpar[1], oscpar[3], oscpar[4], oscpar[5], doubled_angle);
      GetProb(nutype, oscnutype, oscpar[7], oscpar[8], en, num, w);
    }
  else // if nu
    {
      setMNS(oscpar[0], oscpar[2], oscpar[1], oscpar[3], oscpar[4], oscpar[5], doubled_angle);
      GetProb(nutype, oscnutype, oscpar[7], oscpar[8], en, num, w);
    }


  if (std::isnan(w[10]))
    {
      std::cerr << "WARNING: GPU oscillation weight returned NaN! " << w[10] << std::endl;
      std::cerr << "offending osc pars (nue): " << oscpar[0] << " " << oscpar[2] << " " << oscpar[1] << " " << oscpar[3] << " " << oscpar[4] << " " << oscpar[5] << std::endl;
      //std::exit(EXIT_FAILURE);
    }
}

void samplePDFDUNEBaseGPU::calcOscWeightsGPU(double *en, double *w, int num, int nutype, int oscnutype, double *oscpar_nub, double *oscpar_nu)
{
  if (useBeta)
    Beta=oscpar_nub[9];
  
  if (nutype < 0) // if antinu
    {
      //std::cout << "nutype = " << nutype << ", oscnutype = " << oscnutype << std::endl;
      setMNS(oscpar_nub[0], oscpar_nub[2], oscpar_nub[1], oscpar_nub[3], oscpar_nub[4], /*-1**/oscpar_nub[5], doubled_angle);
      GetProb(nutype, oscnutype, oscpar_nub[7], oscpar_nub[8], en, num, w);
      //sleep(.1);
    }
  else // is nu
    {
      setMNS(oscpar_nu[0], oscpar_nu[2], oscpar_nu[1], oscpar_nu[3], oscpar_nu[4], oscpar_nu[5], doubled_angle);
      GetProb(nutype, oscnutype, oscpar_nu[7], oscpar_nu[8], en, num, w);
      //sleep(.1);

      //std::cout << en[0] << std::endl;
      }
  


  if (std::isnan(w[10]))
    {
      std::cerr << "WARNING: GPU oscillation weight returned NaN! " << w[10] << std::endl;
      std::cerr << "offending osc pars (nuebar): " << oscpar_nub[0] << " " << oscpar_nub[2] << " " << oscpar_nub[1] << " " << oscpar_nub[3] << " " << oscpar_nub[4] << " " << oscpar_nub[5] << std::endl;
      //std::exit(EXIT_FAILURE);
      }
}


void samplePDFDUNEBaseGPU::reweight(double *oscpar)
{    
  // This is an awkward hack for now. If RHC is present, it should do 2-parameter reweighting, even with two sets of identical osc pars (to make sure it gets the nu and nubar osc pars right). This will also work for nuebar, because we want the nu and nubar oscpars the same apart from beta.
  if (skmcSamples.size()==6)
    {
      reweight(oscpar, oscpar);
      //std::cout << "doing 2-par reweighting with beta" << std::endl;
    }
  else
    {
      PrepReweight();
      
      for (int i=0; i<(int)skmcSamples.size(); ++i)
	{
	  calcOscWeightsGPU(skmcSamples[i].rw_etru, skmcSamples[i].osc_w, skmcSamples[i].nEvents, skmcSamples[i].nutype, skmcSamples[i].oscnutype, oscpar);
	  
	  for(int j = 0; j < skmcSamples[i].nEvents; ++j)   
	    {  
	      calcWeights(skmcSamples[i], i, j);      
	    }  
	}
    }
}

void samplePDFDUNEBaseGPU::reweight(double *oscpar_nub, double *oscpar_nu)
{
  PrepReweight();
  std::cout << "IN REWEIGHT" << std::endl;
  for(size_t i = 0; i < skmcSamples.size(); ++i)
    {
      //calcOscWeightsGPU(skmcSamples[i].rw_etru, skmcSamples[i].osc_w, skmcSamples[i].nEvents, skmcSamples[i].nutype, skmcSamples[i].oscnutype, oscpar_nub, oscpar_nu);

      for(int j = 0; j < skmcSamples[i].nEvents; ++j)
	{
	  calcWeights(skmcSamples[i], i, j);
	  //calcXsecWeightsGPU(skmcSamples[i], j);         
	}
      //auto start = std::chrono::high_resolution_clock::now();
      calcOscWeightsGPU(skmcSamples[i].rw_etru, skmcSamples[i].osc_w, skmcSamples[i].nEvents, skmcSamples[i].nutype, skmcSamples[i].oscnutype, oscpar_nub, oscpar_nu);
      //auto stop = std::chrono::high_resolution_clock::now();
      //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      //std::cout << "calcOscWeight took: " << duration.count() << " ms " << std::endl;
    }
}
