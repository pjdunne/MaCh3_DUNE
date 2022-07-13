#include <iostream>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "samplePDFDUNE/samplePDFDUNEBase.h"


std::string getNameNoExt(std::string name, std::string ext)  
{                                                            
  std::size_t pos ;                                          
  pos = name.find(ext);                                      
  name = name.substr(0,pos);                                 
  return name ;                                              
}                                                            
                                                             
void saveCanvas(TCanvas* canvas, std::string name, std::string legend)                                                                  
{                                                            
  name = getNameNoExt(name, ".root") ;                       
  name = name + legend + ".root" ;                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".root") ;                       
  name = name + ".png" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".png") ;                        
  name = name + ".pdf" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".pdf") ;                        
  name = name + ".eps" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
} 

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


int main(int argc, char * argv[]) {

  // ----------------------- OPTIONS ---------------------------------------- //

  if(argc == 1){
    std::cout << "Usage: bin/mini_MCMC config.cfg" << std::endl;
    return 1;
  }

  //manager *fitMan = new manager(argv[1]);

  /*

  // there's a check inside the manager class that does this; left here for demonstrative purposes
  if (fitMan->getGoodConfig() == false) {
    std::cerr << "Didn't find a good config in input configuration" << std::endl;
    throw;
  }

  std::string  fluxMatrixFile = fitMan -> getFluxCovMatrix();
  std::string  fluxMatrixName = fitMan -> getFluxCovName();
  std::string  xsecMatrixFile = fitMan -> getXsecCovMatrix();
  std::string  xsecMatrixName = fitMan -> getXsecCovName();
  std::string  skdetMatrixFile = fitMan -> getSKdetCovMatrix();
  std::string  skdetMatrixName = fitMan -> getSKdetCovName();

  // Asimov fit
  bool asimovfit = fitMan->getAsimovFitFlag();
  
  */


  // ----------------------- COVARIANCE AND SAMPLEPDF OBJECTS ---------------------------------------- //

  gStyle -> SetPalette(1);

  // make file to save plots
  TFile *Outfile = new TFile("Dummy_Hist.root" , "RECREATE");


  covarianceXsec *xsec = new covarianceXsec("xsec_cov", "inputs/dummy_xsec_covariance.root") ;


  std::cout << "---------- Printing off nominal parameter values ----------" << std::endl;
  std::cout << "Cross section parameters:" << std::endl;
  xsec->printNominal();

  std::cout << "---------- Finished printing nominal parameter values ----------" << std::endl;


  xsec->setParameters();

  //xsec->setStepScale(fitMan->getXsecStepScale());
  xsec->setStepScale(0.01);

  covarianceOsc *osc = new covarianceOsc("osc_cov","inputs/oscillation_covariance_6par_nondouble_PDG2019.root");

  // oscpars from manager in order:
  // sin2th_12, sin2th_23, sin2th_13, delm2_12, delm2_23, delta_cp
  //std::vector<double> oscpars =fitMan->getOscParameters();   
  std::vector<double> oscpars{0.307,0.528,0.0218,7.53e-5, 2.509e-3,-1.601}; // Asimov A

  std::cout<<"Using these oscillation parameters: ";
  for(unsigned ipar=0;ipar<oscpars.size();ipar++) std::cout<<" "<<oscpars.at(ipar);
  std::cout << std::endl;

  osc->setFlipDeltaM23(true);

  // Ask config file whether to use reactor constraint
  //std::cout << "use reactor prior is : " << useRC << std::endl ;

  // Use prior for 12 parameters only
  //osc->setEvalLikelihood(0,false);
  osc->setEvalLikelihood(1,false);
  osc->setEvalLikelihood(2,false);
  //osc->setEvalLikelihood(3,false);
  osc->setEvalLikelihood(4,false);
  osc->setEvalLikelihood(5,false);
  // This line gives a crash and stack trace...
  osc->setParameters(oscpars);
  osc->acceptStep();

  // Set some sample....
  samplePDFDUNEBase * numu_pdf = new samplePDFDUNEBase(1.3628e+23, "configs/SamplePDFDune_FHC_numuselec.cfg", xsec);

  // unosc
  std::vector<double> oscpars_un(oscpars);
  oscpars_un[0] = 0;
  oscpars_un[1] = 0;
  oscpars_un[2] = 0;

  // Oscillated
  osc -> setParameters(oscpars_un);
  std::cout << "oscpars[0] = " << (osc -> getPropPars())[0] << std::endl
	    << "oscpars[1] = " << (osc -> getPropPars())[1] << std::endl
	    << "oscpars[2] = " << (osc -> getPropPars())[2] << std::endl
	    << "oscpars[3] = " << (osc -> getPropPars())[3] << std::endl
	    << "oscpars[4] = " << (osc -> getPropPars())[4] << std::endl
	    << "oscpars[5] = " << (osc -> getPropPars())[5] << std::endl;

  //std::vector<double> CAFana_default_edges = get_default_CAFana_bins();
  //numu_pdf -> useBinnedOscReweighting(true, CAFana_default_edges.size()-1, &CAFana_default_edges[0]);

  //TH1D	*numu_nominal_hist	= (TH1D*)numu_pdf      -> get1DHist() -> Clone("numu_nom");


  // Unoscillated
  osc -> setParameters(oscpars_un);
  numu_pdf -> reweight(osc->getPropPars(), osc->getPropPars());

  TH1D *numu_unosc = (TH1D*)numu_pdf -> get1DHist() -> Clone("numu_plain");


  TCanvas *nomcanv = new TCanvas("nomcanv","",1200,600);
  //nomcanv->Divide(2,1);
  //nomcanv->cd(1);
  numu_unosc -> Draw("HIST");
  //nomcanv->cd(2);
  //numu_nominal_hist->Draw("HIST");

  std::string plotname = "Dummy Hist" ;
  saveCanvas(nomcanv, plotname,"_nominal_spectra") ;

  Outfile -> cd();
  numu_unosc			-> Write("numu_unosc");
  //numu_nominal_hist			-> Write("numu_nominal_hist");


  std::cout	<< "Integrals of nominal hists: " << std::endl
    		<< "Numu unosc:      " << numu_unosc			-> Integral() << std::endl; 
    		//<< "Numu osc:        " << numu_nominal_hist		-> Integral() << std::endl;


  return 0;
 }
