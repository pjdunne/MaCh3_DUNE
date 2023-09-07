#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "samplePDFDUNE/samplePDFDUNEBase.h"

#include "mcmc/mcmc.h"


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



int main(int argc, char * argv[]) {

  // ----------------------- OPTIONS ---------------------------------------- //

  if(argc == 1){
    std::cout << "Usage: bin/mini_MCMC config.cfg" << std::endl;
    return 1;
  }

  manager *fitMan = new manager(argv[1]);


  std::string  XsecMatrixFile = fitMan->raw()["General"]["Systematics"]["XsecCovFile"].as<std::string>(); 
  std::string  XsecMatrixName = fitMan->raw()["General"]["Systematics"]["XsecCovName"].as<std::string>();
  std::string  OscMatrixFile = fitMan->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
  std::string  OscMatrixName = fitMan->raw()["General"]["Systematics"]["OscCovName"].as<std::string>();
 
  // Asimov fit
  bool asimovfit = true;



  // ----------------------- COVARIANCE AND SAMPLEPDF OBJECTS ---------------------------------------- //

  gStyle -> SetPalette(1);

  // make file to save plots
  std::string odir = fitMan->raw()["General"]["Output"]["ODIR"].as<std::string>();
  std::string ofilename = fitMan->raw()["General"]["Output"]["OUTPUTNAME"].as<std::string>();
  TFile *Outfile = new TFile((odir + "/" + ofilename).c_str(),"RECREATE");


  covarianceXsec *xsec = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str()) ;


  //std::cout << "---------- Printing off nominal parameter values ----------" << std::endl;
  //std::cout << "Cross section parameters:" << std::endl;
  // xsec->printNominal();

  //std::cout << "---------- Finished printing nominal parameter values ----------" << std::endl;


  xsec->setParameters();

  xsec->setStepScale(0.01);

  // Oscillation covariance
  covarianceOsc *osc = new covarianceOsc("osc_cov","inputs/osc_covariance_DUNE_PDG2021_v1.root");
  osc->setStepScale(0.045);

  // oscpars from manager in order:
  // sin2th_12, sin2th_23, sin2th_13, delm2_12, delm2_23, delta_cp
  std::vector<double> oscpars{0.307,0.528,0.0218,7.53e-5, 2.509e-3,-1.601}; // Asimov A

  std::cout<<"Using these oscillation parameters: ";
  for(unsigned ipar=0;ipar<oscpars.size();ipar++) std::cout<<" "<<oscpars.at(ipar);
  std::cout << std::endl;

  // Ask config file whether to use reactor constraint
  //bool useRC =  fitMan -> getRC() ;
  //std::cout << "use reactor prior is : " << useRC << std::endl ;
  //  osc->useReactorPrior(useRC); // this is hard coded inside, and is bad

  // Use prior for 12 parameters only
  //osc->setEvalLikelihood(0,false);
  osc->setEvalLikelihood(1,false);
  osc->setEvalLikelihood(2,false);
  //osc->setEvalLikelihood(3,false);
  osc->setEvalLikelihood(4,false);
  osc->setEvalLikelihood(5,false);
  // This line gives a crash and stack trace...
  osc->setParameters(oscpars);
  //osc->acceptStep();
  //osc->setStepScale(fitMan->getOscStepScale());

  //  osc->setFlipDeltaM23(false);

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Creating Oscillator object" << std::endl;
  std::string OscillatorCfgName = fitMan->raw()["General"]["OscillatorConfigName"].as<std::string>();
  Oscillator* Oscill = new Oscillator(OscillatorCfgName);

  Oscill->FillOscillogram(osc->getPropPars(),25.0,0.5);
  
  std::vector<samplePDFDUNEBase*> pdfs;
  // samplePDFDUNEBase *numu_pdf = new samplePDFDUNEBase(1.3628319e+23, "configs/SamplePDFDune_FHC_numuselec.yaml", xsec);
  // samplePDFDUNEBase *numubar_pdf = new samplePDFDUNEBase(1.3628319e+23, "configs/SamplePDFDune_RHC_numuselec.yaml", xsec);
  // samplePDFDUNEBase *nue_pdf = new samplePDFDUNEBase(1.3628319e+23, "configs/SamplePDFDune_FHC_nueselec.yaml", xsec);
  // samplePDFDUNEBase *nuebar_pdf = new samplePDFDUNEBase(1.3628319e+23, "configs/SamplePDFDune_RHC_nueselec.yaml", xsec);
  samplePDFDUNEBase *numu_pdf = new samplePDFDUNEBase(1.3628319e+23, "configs/AtmSample_numuselec.yaml", xsec);
  samplePDFDUNEBase *nue_pdf = new samplePDFDUNEBase(1.3628319e+23, "configs/AtmSample_nueselec.yaml", xsec);
  numu_pdf->SetOscillator(Oscill);
  nue_pdf->SetOscillator(Oscill);
  
  pdfs.push_back(numu_pdf);
  // pdfs.push_back(numubar_pdf);
  pdfs.push_back(nue_pdf);
  // pdfs.push_back(nuebar_pdf);

  // From CAFana
  // const int kNumTrueEnergyBins = 100;
  
  // // N+1 bin low edges
  // std::vector<double> edges(kNumTrueEnergyBins+1);
  
  // const double Emin = .5; // 500 MeV: there's really no events below there
  
  // // How many edges to generate. Allow room for 0-Emin bin
  // const double N = kNumTrueEnergyBins-1;
  // const double A = N*Emin;
  
  // edges[0] = 0;
  
  // for(int i = 1; i <= N; ++i){
  //   edges[kNumTrueEnergyBins-i] = A/i;
  // }
  
  // edges[kNumTrueEnergyBins] = 120; // Replace the infinity that would be here
  
  // numu_pdf -> useBinnedOscReweighting(true, kNumTrueEnergyBins, &edges[0]);
  

  osc -> setParameters(oscpars);
  std::cout << "oscpars[0] = " << (osc -> getPropPars())[0] << std::endl
	    << "oscpars[1] = " << (osc -> getPropPars())[1] << std::endl
	    << "oscpars[2] = " << (osc -> getPropPars())[2] << std::endl
	    << "oscpars[3] = " << (osc -> getPropPars())[3] << std::endl
	    << "oscpars[4] = " << (osc -> getPropPars())[4] << std::endl
	    << "oscpars[5] = " << (osc -> getPropPars())[5] << std::endl;

  // Setup Oscillation Calculator
  // for (unsigned ipdf = 0 ; ipdf < pdfs.size() ; ++ipdf) {
  //   pdfs[ipdf] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
  // }

  // Set to nominal pars
  vector<double> xsecpar = xsec->getNominalArray();  
  xsec->setParameters(xsecpar);
  //  xsec->Print();

  numu_pdf->reweight(osc->getPropPars());
  TH1D *numu_asimov = (TH1D*)numu_pdf->get1DHist()->Clone("numu_asimov");
  nue_pdf->reweight(osc->getPropPars());
  TH1D *nue_asimov = (TH1D*)nue_pdf->get1DHist()->Clone("nue_asimov");
  // numubar_pdf->reweight(osc->getPropPars());
  // TH1D *numubar_asimov = (TH1D*)numubar_pdf->get1DHist()->Clone("numubar_asimov");
  // nuebar_pdf->reweight(osc->getPropPars());
  // TH1D *nuebar_asimov = (TH1D*)nuebar_pdf->get1DHist()->Clone("nuebar_asimov");

  // Print event rates to check
  std::cout << "-------- SK event rates for Asimov fit (Asimov fake data) ------------" << std::endl;
  std::cout << "FHC 1Rmu:   " << numu_asimov->Integral() << std::endl;
  std::cout << "FHC 1Re:    " << nue_asimov->Integral() << std::endl;
  // std::cout << "RHC 1Rmu:   " << numubar_asimov->Integral() << std::endl;
  // std::cout << "RHC 1Re:    " << nuebar_asimov->Integral() << std::endl;

  TH2D *numu_asimov_2d = (TH2D*)numu_pdf->get2DHist()->Clone("numu_asimov_2d");
  TH2D *nue_asimov_2d = (TH2D*)nue_pdf->get2DHist()->Clone("nue_asimov_2d");

  numu_pdf->addData(numu_asimov_2d); 
  nue_pdf->addData(nue_asimov_2d); 
  // numubar_pdf->addData(numubar_asimov);  
  // nuebar_pdf->addData(nuebar_asimov); 
  
  
  TCanvas *sys_canv = new TCanvas("sys_canv","",1200,600);
  
  // Save pars before modifying
  std::vector<double> tpar = xsecpar;

  int n_points = fitMan->raw()["General"]["LLH"]["nbins"].as<int>();
  int parId = fitMan->raw()["General"]["LLH"]["parId1"].as<int>();
  
  int i = parId;

  // for (int i=56; i < xsecpar.size(); i++) {

  xsecpar = xsec->getNominalArray();
  double lowerb = xsec->getNominal(i)-3.1*sqrt((*xsec->getCovMatrix())(i,i));
  double upperb = xsec->getNominal(i)+3.1*sqrt((*xsec->getCovMatrix())(i,i));

  std::string histname = "xsec_" + std::to_string(i) + "_llhscan";
  std::string histtitle = "xsec_" + std::to_string(i);
  
  TH1D *hScan = new TH1D(histname.c_str(), histtitle.c_str(), 40, lowerb, upperb);
  xsecpar[i] = xsec->getNominal(i)-3*sqrt((*xsec->getCovMatrix())(i,i));

  //double dsigma = (2*3)/(n_points-1)
  double dsigma = 0.07594936708;
  double totalllh = 0;
  double samplellh = 0;
  double penaltyllh = 0;
  std::cout << "######## Xsec Systematic " << i << std::endl;
  for (int j=0; j < n_points; j++) {
    xsec->setParameters(xsecpar);
    for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++) {
      pdfs[ipdf] -> reweight(osc -> getPropPars());
  std::cout << "Sample " << ipdf << " has sample llh " << pdfs[ipdf]->GetLikelihood() << " for variation " << j <<  std::endl; 
      samplellh += pdfs[ipdf]->GetLikelihood();
    }
    penaltyllh += (xsec->GetLikelihood());
    //std::cout << "xsec par = " << i << " || xsec par value = " << xsecpar[i] << " || sample llh = " << samplellh << " || penalty llh  = " <<  penaltyllh << std::endl;
  totalllh = samplellh ;//+ penaltyllh;
    hScan->Fill(xsecpar[i], 2*totalllh);
    samplellh = 0;
    penaltyllh = 0;
    totalllh = 0;
    xsecpar[i] += dsigma*sqrt((*xsec->getCovMatrix())(i,i));
  }
  Outfile->cd();
  hScan->Write();
  std::cout << "Finished xsec param " << i << std::endl;
  // }
  return 0;
}


