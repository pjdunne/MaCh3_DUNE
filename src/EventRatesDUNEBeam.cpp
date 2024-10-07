#include <iostream>
#include <chrono>
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

#include "samplePDFDUNE/MaCh3DUNEFactory.h"


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

void DrawModeBreakDown(samplePDFDUNEBeamFDBase*& Sample, std::unique_ptr<TFile> &File){

  File->cd();
  TH1D* ModeHistogram = nullptr;
  THStack* StackedModes = new THStack(Sample->GetName().c_str(), "stack");
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);

  std::ofstream outfile;
  std::string FileName=Sample->GetName()+"_EventRates.txt";
  outfile.open(FileName);

  outfile << "\\begin{table}[ht]" << std::endl;
  outfile << "\\begin{center}" << std::endl;
  outfile << "\\caption{Integral breakdown for sample: " << Sample->GetName() << "}" << std::endl;
  outfile << "\\label{" << Sample->GetName() << "-EventRate}" << std::endl;


  int space = 14;
  TString nColumns;
  nColumns += "|c|";
  outfile << "\\begin{tabular}{|l" << nColumns.Data() << "}" << std::endl;
  outfile << "\\hline" << std::endl;
  outfile << "&" << std::setw(space) << Sample->GetName() << " " << std::endl;

  double total = 0;
  double mode_rate = 0;
  for (unsigned int iMode=0;iMode<kMaCh3_nModes ;iMode++) {

	////
	std::string name = Sample->GetName();
	name += MaCh3mode_ToDUNEString(static_cast<MaCh3_Mode>(iMode));
	ModeHistogram = (TH1D*)Sample->get1DVarHist(kRecoNeutrinoEnergy, iMode);
	ModeHistogram->SetName(name.c_str());
	mode_rate = ModeHistogram->Integral();
	total += mode_rate; 

	outfile << MaCh3mode_ToDUNEString(static_cast<MaCh3_Mode>(iMode)).c_str() << std::setw(space) << " & " << mode_rate << " \\\\" << std::endl;

	//Add this to the stacked histogram
	if(ModeHistogram->Integral() > 0){
	  ModeHistogram->SetFillColor(kRed+iMode);
	  ModeHistogram->SetLineColor(kRed+iMode);
	  StackedModes->Add(ModeHistogram);	
	  leg->AddEntry(ModeHistogram, MaCh3mode_ToDUNEString(static_cast<MaCh3_Mode>(iMode)).c_str(), "f");
	  ModeHistogram->Write();
	}
  }
  outfile << "\\hline " << std::endl;
  outfile << "&" << std::setw(space) << "Total:" << total << "\\\\ \\hline" << std::endl;

  std::string PdfName = Sample->GetName();
  PdfName += "_ModeStack.pdf";
  StackedModes->Write();
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  StackedModes->Draw("HIST");
  leg->Draw("SAMES");
  c1->Print(PdfName.c_str());
  delete c1;
  

  return;
}

int main(int argc, char * argv[]) {

  // ----------------------- OPTIONS ---------------------------------------- //
  if(argc == 1){
    std::cout << "Usage: bin/mini_MCMC config.cfg" << std::endl;
    return 1;
  }

  manager *fitMan = new manager(argv[1]);
  auto OutputFileName = fitMan->raw()["General"]["OutputFile"].as<std::string>();

  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;

  //####################################################################################
  //Create samplePDFSKBase Objs
  std::cout << "Loading T2K samples.." << "\n" << std::endl;
  std::vector<samplePDFDUNEBeamFDBase*> DUNEPdfs;
  MakeMaCh3DuneBeamInstance(fitMan, DUNEPdfs, xsec, osc); 
  //Setup the cross-section parameters
  //This should get the prior values.
  std::vector<double> XsecParVals = xsec->getNominalArray();

  xsec->setParameters(XsecParVals);
  xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());

  //Some place to store the histograms
  std::vector<TH1D*> oscillated_hists;
  std::vector<TH1D*> unoscillated_hists;
  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(new TFile(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  for (unsigned sample_i = 0 ; sample_i < DUNEPdfs.size() ; ++sample_i) {
    	
	std::string name = DUNEPdfs[sample_i]->GetName();
	sample_names.push_back(name);
	TString NameTString = TString(name.c_str());

	DUNEPdfs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
	osc->setParameters();
	DUNEPdfs[sample_i] -> reweight();
	TH1D *SampleHistogram = (TH1D*)DUNEPdfs[sample_i] -> get1DHist() -> Clone(NameTString+"_unosc");
	unoscillated_hists.push_back(SampleHistogram);

	TCanvas *nomcanv = new TCanvas("nomcanv","",1200,600);
	SampleHistogram -> Draw("HIST");

	std::string plotname = "Dummy_Hist";
	saveCanvas(nomcanv, plotname,"_nominal_spectra");
	//DrawModeBreakDown(DUNEPdfs[sample_i], OutputFile);
  }
 
  //Now print out some event rates, we'll make a nice latex table at some point 
  for (unsigned iPDF = 0; iPDF < DUNEPdfs.size() ; ++iPDF) {
	std::cout << "Integrals of nominal hists: " << std::endl;
	std::cout << sample_names[iPDF].c_str() << ": " << unoscillated_hists[iPDF]-> Integral() << std::endl;
	std::cout << "~~~~~~~~~~~~~~~~" << std::endl;
  }

  OutputFile->Write();
  OutputFile->Close();
  return 0;
 }
