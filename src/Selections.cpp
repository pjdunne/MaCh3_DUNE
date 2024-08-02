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

#include "samplePDFDUNE/samplePDFDUNEBase.h"
#include "samplePDFDUNE/samplePDFDUNEBaseND.h"
#include "samplePDFDUNE/samplePDFDUNEBaseNDGAr.h"
#include "manager/manager.h"


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

  double FDPOT = fitMan->raw()["General"]["FDPOT"].as<double>(); 
  double NDPOT = fitMan->raw()["General"]["NDPOT"].as<double>(); 


  // Asimov fit
  bool EvalXsec = fitMan->raw()["SigmaVariations"]["EvalXsec"].as<bool>();
  //bool PlotByMode = fitMan->raw()["SigmaVariations"]["PlotByMode"].as<bool>();


  // ----------------------- COVARIANCE AND SAMPLEPDF OBJECTS ---------------------------------------- //

  gStyle -> SetPalette(1);

  // make file to save the variation of xsec params to
  std::string OutputFileName = fitMan->raw()["General"]["Output"]["FileName"].as<std::string>();
  TFile *OutputFile = new TFile(OutputFileName.c_str() , "RECREATE");

  //initialise xsec
  covarianceXsec *xsec = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str());

  std::cout << "---------- Printing off nominal parameter values ----------" << std::endl;
  std::cout << "Cross section parameters:" << std::endl;
  xsec->printNominal();

  bool XsecParsAtGen = fitMan->raw()["General"]["Systematics"]["XsecAtGen"].as<bool>();

  std::cout << "---------- Finished printing nominal parameter values ----------" << std::endl;

  covarianceOsc *osc = new covarianceOsc(OscMatrixName.c_str(), OscMatrixFile.c_str());

  std::vector<double> oscpars = fitMan->raw()["General"]["OscParams"].as<std::vector<double>>();

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

  bool addFD = fitMan->raw()["General"]["IncludeFD"].as<bool>();
  bool addND = fitMan->raw()["General"]["IncludeND"].as<bool>();
  bool addNDGAr = fitMan->raw()["General"]["IncludeNDGAr"].as<bool>();
  
  bool efficiency = fitMan->raw()["Selections"]["Efficiency"].as<bool>();
  bool efficiencyvsenergy = fitMan->raw()["Selections"]["EfficiencyVsEnergy"].as<bool>();
  bool twodimhist = fitMan->raw()["Selections"]["2DHist"].as<bool>();
//  int n_efficiencysamples = fitMan->raw()["Selections"]["NEfficiencySamples"].as<int>();

  if (!addFD && !addND && !addNDGAr) {std::cerr << "[ERROR:] You've chosen NOT to include FD or ND samples in the config file... you need to add something!" << std::endl; throw;}


  std::vector<samplePDFFDBase*> SamplePDFs;
  
  if(addFD) { 
    samplePDFDUNEBase *numu_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_FHC_numuselec.yaml", xsec);
    SamplePDFs.push_back(numu_pdf);
    samplePDFDUNEBase *nue_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_FHC_nueselec.yaml", xsec);
    SamplePDFs.push_back(nue_pdf);
    samplePDFDUNEBase *numubar_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_RHC_numuselec.yaml", xsec);
    SamplePDFs.push_back(numubar_pdf);
    samplePDFDUNEBase *nuebar_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_RHC_nueselec.yaml", xsec);
    SamplePDFs.push_back(nuebar_pdf);
  }
 
  if(addND) {
    samplePDFDUNEBaseND * FHC_numuCCND_pdf = new samplePDFDUNEBaseND(NDPOT, "configs/SamplePDFDuneND_FHC_CCnumuselec.yaml", xsec);
    SamplePDFs.push_back(FHC_numuCCND_pdf);
    samplePDFDUNEBaseND * RHC_numuCCND_pdf = new samplePDFDUNEBaseND(NDPOT, "configs/SamplePDFDuneND_RHC_CCnumuselec.yaml", xsec);
    SamplePDFs.push_back(RHC_numuCCND_pdf);
  }

  if(addNDGAr) {
    if(efficiency) {
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf);
    }
    else if(efficiencyvsenergy){
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_0 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_0.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_0);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_0 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_0.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_0);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_1 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_1.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_1);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_1 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_1.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_1);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_2 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_2.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_2);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_2 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_2.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_2);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_3 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_3.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_3);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_3 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_3.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_3);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_4 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_4.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_4);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_4 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_4.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_4);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_5 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_5.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_5);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_5 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_5.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_5);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_6 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_6.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_6);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_6 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_6.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_6);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_7 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_7.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_7);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_7 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_7.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_7);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_8 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_8.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_8);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_8 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_8.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_8);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_reco_pdf_9 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_RecoCuts_9.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_reco_pdf_9);
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_true_pdf_9 = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec_TruthCuts_9.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_true_pdf_9);
    }
    else{
      samplePDFDUNEBaseNDGAr * FHC_numuCCNDGAr_pdf = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml", xsec);
      SamplePDFs.push_back(FHC_numuCCNDGAr_pdf);
    }
//    samplePDFDUNEBaseNDGAr * RHC_numuCCNDGAr_pdf = new samplePDFDUNEBaseNDGAr(NDPOT, "configs/SamplePDFDuneND_RHC_CCnumuselec.yaml", xsec);
//    SamplePDFs.push_back(RHC_numuCCNDGAr_pdf);
  }


  // Oscillated
  osc -> setParameters(oscpars);
  std::cout << "oscpars[0] = " << (osc -> getPropPars())[0] << std::endl
	<< "oscpars[1] = " << (osc -> getPropPars())[1] << std::endl
	<< "oscpars[2] = " << (osc -> getPropPars())[2] << std::endl
	<< "oscpars[3] = " << (osc -> getPropPars())[3] << std::endl
	<< "oscpars[4] = " << (osc -> getPropPars())[4] << std::endl
	<< "oscpars[5] = " << (osc -> getPropPars())[5] << std::endl;

  // unosc
  std::vector<double> oscpars_un(oscpars);
  oscpars_un[0] = 0;
  oscpars_un[1] = 0;
  oscpars_un[2] = 0;
  //oscpars_un[3] = 0;
  //oscpars_un[4] = 0;

  //Setup the cross-section parameters
  //This should get the prior values.
  std::vector<double> XsecParVals = xsec->getNominalArray();

  if(XsecParsAtGen){
	TFile* XsecFile = new TFile(XsecMatrixFile.c_str(), "READ");
	TVectorD* XsecGeneratedParamArray = (TVectorD*)XsecFile->Get("xsec_param_nom");
	std::cout << "Setting xsec systs to their generated values " << std::endl;
	for (unsigned param_i = 0 ; param_i < XsecParVals.size() ; ++param_i) {
	  std::cout << "Generated value for param " << param_i << " is " << (*XsecGeneratedParamArray)(param_i) << std::endl;
	  XsecParVals[param_i] = (*XsecGeneratedParamArray)(param_i);
	  std::cout << "Set parameter " << param_i << " to value " << XsecParVals[param_i] << std::endl;
	}
        XsecFile->Close();
  }
  else{
	std::cout << "Keeping xsec parameters at their prior values" << std::endl;
  }

  xsec->setParameters(XsecParVals);
  xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());
  osc->setStepScale(fitMan->raw()["General"]["Systematics"]["OscStepScale"].as<double>());


  //Some place to store the histograms
  std::vector<TH1D*> oscillated_hists;
  std::vector<TH1D*> unoscillated_hists;
  std::vector<std::string> sample_names;

  for (unsigned sample_i = 0 ; sample_i < SamplePDFs.size() ; ++sample_i) {

	std::string name = SamplePDFs[sample_i]->GetSampleName();

	sample_names.push_back(name);
	TString NameTString = TString(name.c_str());
	// Unoscillated
	osc -> setParameters(oscpars_un);
	osc -> acceptStep();

	SamplePDFs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
	SamplePDFs[sample_i] -> reweight(osc->getPropPars());
	TH1D *sample_unosc = (TH1D*)SamplePDFs[sample_i] -> get1DHist() -> Clone(NameTString+"_unosc");
	unoscillated_hists.push_back(sample_unosc);

	OutputFile -> cd();
	sample_unosc			-> Write(NameTString+"_unosc");

	osc -> setParameters(oscpars);
	osc -> acceptStep();

	SamplePDFs[sample_i] -> reweight(osc->getPropPars());
	TH1D *numu_osc = (TH1D*)SamplePDFs[sample_i] -> get1DHist()->Clone(NameTString+"_osc");
	oscillated_hists.push_back(numu_osc);

  }

  //Now print out some event rates, we'll make a nice latex table at some point 
  std::cout << "Integrals of nominal hists: " << std::endl;
  for (unsigned sample_i = 0; sample_i < SamplePDFs.size() ; ++sample_i) {
	std::cout << sample_names[sample_i].c_str() << " unosc:      " << unoscillated_hists[sample_i]-> Integral() << std::endl;
	std::cout << sample_names[sample_i].c_str() << "   osc:      " << oscillated_hists[sample_i]-> Integral() << std::endl; 
	std::cout << " " << std::endl;
  }
  std::cout << "~~~~~~~~~~~~~~~~" << std::endl;

  /////////////////
  //
  // This is where the actual variation of systematics happens
  //
  //////////////////

  if (EvalXsec)
  {
	std::cout << "EVALUATING XSEC PARAMS" << std::endl;
	OutputFile->cd();
	for(unsigned ipdf = 0 ; ipdf < SamplePDFs.size() ; ipdf++){
	  unoscillated_hists[ipdf]->Write();
	  oscillated_hists[ipdf]->Write();
	}
  }

  // set back to oscillated spectra
  osc->setParameters(oscpars_un);

  for(unsigned ipdf=0;ipdf<SamplePDFs.size();ipdf++){
	SamplePDFs[ipdf]->reweight(osc->getPropPars());
  }


  const int nsubsamples = 12;

  //Get the nice names of the modes  
  int nmodes = kMaCh3_nModes;
  std::string mode_names[nmodes];
  for(unsigned i=0;i<kMaCh3_nModes;i++){
	mode_names[i] = MaCh3mode_ToDUNEString((MaCh3_Mode)i);
	//	std::cout << mode_names[i] << std::endl;
  }

  char houtname[100];

  std::vector<std::vector<std::vector<TH1D*>>> mode_nom;


  //Print out to two decimal places
  std::cout << setprecision(2);

  //If you want to evaluate xsec systs
  if (EvalXsec)
  { 

	OutputFile->cd();

	std::vector<char*> xsec_names;

	// ------ Do Xsec Variations ----- //
	vector<double> xsecpar = xsec->getNominalArray();
	//      for (int i=0; i<int(xsecpar.size()); i++)

	std::vector<TH1D*> varied_hists;

	//Now loop over all xsec systs
	for (int i=0; i<int(xsecpar.size()); i++)
	{

	  char xsec_name[100];
	  sprintf(xsec_name, "%i : %s", i, xsec->getParName(i));
	  xsec_names.push_back(xsec_name);
	  xsecpar = xsec->getNominalArray();

	  for(unsigned ipdf = 0 ; ipdf < SamplePDFs.size() ; ipdf++) {

		//Now loop over the sigma values you want to evaluate at (usual -3, -1, 0, +1, +3)
		for (int j=-3; j<=3; j+=2) {
		  char sign = (j>0) ? 'p' : 'n';

		  //Set xsec par to the value
		  xsecpar[i] = xsec->getNominal(i)+j*sqrt((*xsec->getCovMatrix())(i,i));
		  xsec->setParameters(xsecpar);
		  //Reweight prediction
		  SamplePDFs[ipdf]->reweight(osc->getPropPars());

		  sprintf( houtname, "%s_xsec_par_%i_sig_%c%i", sample_names[ipdf].c_str(), i, sign,abs(j));
		  TH1D * hist = (TH1D*)SamplePDFs[ipdf]->get1DHist()->Clone(houtname);
		  
		    //LW - Right now core can't plot by mode
			/*

		    //Loop over modes
		    if (PlotByMode) { 
		      for(int mode = 0 ; mode < kMaCh3_nModes ; mode++) {

		        //Loop over oscillaiton channels
		        for (int samp = 0 ; samp < 12 ; samp++) {

		          sprintf(houtname,"%s_xsec_par_%i_sm_%s_%s_sig_%c%i", names[ipdf],i,OscChannelNames[samp],mode_names[mode].c_str(),sign,abs(j));
		          TH1D * histmode = (TH1D*)pdfs[ipdf]->getModeHist1D(samp, mode, 2)->Clone(houtname);

		          if(histmode->Integral() == 0){continue;}
		          histmode->Write();
		        }//end of loop over osc channels 

		      }//end of loop over modes

		    }//end of PlotByMode
		    */

		  hist->Write();

		  delete hist;

		}// loop over sigma values

	  }//end of loop over pdfs

	}//end of loop over xsec systs

  }//end of EvalXsec
  std::vector<std::string> HistVariables = fitMan->raw()["Selections"]["KinematicParsToPlot"].as<std::vector<std::string>>();
  bool StackByMode = fitMan->raw()["Selections"]["StackByMode"].as<bool>();
  std::vector<std::string> PlotModes = fitMan->raw()["Selections"]["PlotModes"].as<std::vector<std::string>>();
  int Weighted = fitMan->raw()["Selections"]["Weighted"].as<int>();
  for(int i =0; i<PlotModes.size(); i++){std::cout<<"Plot Modes "<<i<<" : "<<PlotModes[i]<<" MaCh3 Mode: "<<(int)DUNEString_ToMaCh3Mode(PlotModes[i])<<std::endl;}
//  std::vector<std::string> HistVariables = {"RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "PionMultiplicity", "NRecoParticles", "InFDV"};

  if(efficiency){ 
    for(unsigned ipdf = 0 ; ipdf < SamplePDFs.size() ; ipdf++) {
      for(int ihist=0; ihist<HistVariables.size(); ihist++){
        for(int imode = 0; imode<PlotModes.size(); imode++){
          int mode;
          if(DUNEString_ToMaCh3Mode(PlotModes[0])==(int)(kMaCh3_nModes)){mode =-1;}
          else{mode = (int)DUNEString_ToMaCh3Mode(PlotModes[imode]); std::cout<<"mode: "<<mode<<std::endl;}
          TH1D * histogram = (TH1D*)SamplePDFs[ipdf]->get1DVarHist(HistVariables[ihist], mode, -1, Weighted, NULL);
          string name = histogram->GetName();
          string nameadd;
          if(ipdf ==0){nameadd = "_RECO";}
          else if(ipdf == 1){nameadd = "_TRUTH";}
          string newname = name + nameadd;
          histogram ->SetName(newname.c_str());
          histogram ->Write();
          delete histogram;
        }
      }
    }
  }
  else if(efficiencyvsenergy){ 
    for(unsigned ipdf = 0 ; ipdf < SamplePDFs.size() ; ipdf++) {
      for(int ihist=0; ihist<HistVariables.size(); ihist++){
        for(int imode = 0; imode<PlotModes.size(); imode++){
          int mode;
          if(DUNEString_ToMaCh3Mode(PlotModes[0])==(int)(kMaCh3_nModes)){mode =-1;}
          else{mode = (int)DUNEString_ToMaCh3Mode(PlotModes[imode]); std::cout<<"mode: "<<mode<<std::endl;}
          TH1D * histogram = (TH1D*)SamplePDFs[ipdf]->get1DVarHist(HistVariables[ihist], mode, -1, Weighted, NULL);
          string name = histogram->GetName();
          string nameadd;
          string nameadd2;
          if(ipdf % 2 ==0){nameadd = "_RECO_"; nameadd2 = std::to_string(ipdf/2).c_str();}
          else if(ipdf % 2 != 0){nameadd = "_TRUTH_"; nameadd2 = std::to_string((ipdf-1)/2).c_str();}
          string newname = name + nameadd + nameadd2;
          histogram ->SetName(newname.c_str());
          histogram ->Write();
          delete histogram;
        }
      }
    }
  }
  else if(twodimhist){ 
    for(unsigned ipdf = 0 ; ipdf < SamplePDFs.size() ; ipdf++) {
      for(int ihist=0; ihist<HistVariables.size(); ihist+2){
        for(int imode = 0; imode<PlotModes.size(); imode++){
          int mode;
          std::cout<<"ipdf: "<<ipdf<<" ihist: "<<ihist<<" imode: "<<imode<< " SamplePDFs.size()"<< SamplePDFs.size()<<std::endl;
          if(DUNEString_ToMaCh3Mode(PlotModes[0])==(int)(kMaCh3_nModes)){mode =-1; std::cout<<"mode: "<<mode<<std::endl;}
          else{mode = (int)DUNEString_ToMaCh3Mode(PlotModes[imode]); std::cout<<"mode: "<<mode<<std::endl;}
          TH2D * histogram = (TH2D*)SamplePDFs[ipdf]->get2DVarHist(HistVariables[ihist], HistVariables[ihist+1], mode, -1, Weighted, NULL, NULL);
          histogram ->Write();
          delete histogram;
        }
      }
    }
  }
  else{ 
    for(unsigned ipdf = 0 ; ipdf < SamplePDFs.size() ; ipdf++) {
      for(int ihist=0; ihist<HistVariables.size(); ihist++){
        for(int imode = 0; imode<PlotModes.size(); imode++){
          int mode;
          if(DUNEString_ToMaCh3Mode(PlotModes[0])==(int)(kMaCh3_nModes)){mode =-1; std::cout<<"mode: "<<mode<<std::endl;}
          else{mode = (int)DUNEString_ToMaCh3Mode(PlotModes[imode]); std::cout<<"mode: "<<mode<<std::endl;}
          TH1D * histogram = (TH1D*)SamplePDFs[ipdf]->get1DVarHist(HistVariables[ihist], mode, -1, Weighted, NULL);
          histogram ->Write();
          delete histogram;
        }
      }
    }
  } 
/*
  if(StackByMode){ 
    for(unsigned ipdf = 0 ; ipdf < SamplePDFs.size() ; ipdf++) {
      for(int ihist=0; ihist<HistVariables.size(); ihist++){
        std::string stackname = HistVariables[ihist];
        std::cout<<"stackname: "<<stackname<<std::endl;
        THStack* stack = new THStack(stackname.c_str(), HistVariables[ihist].c_str());
        TCanvas* canvas = new TCanvas(HistVariables[ihist].c_str(), HistVariables[ihist].c_str());
        for(int imode = 0; imode<PlotModes.size(); imode++){
          int mode;
          if(DUNEString_ToMaCh3Mode(PlotModes[0])==(int)(kMaCh3_nModes)){mode =-1;}
          else{mode = (int)DUNEString_ToMaCh3Mode(PlotModes[imode]); std::cout<<"mode: "<<mode<<std::endl;}
          TH1D * histogram = (TH1D*)SamplePDFs[ipdf]->get1DVarHist(HistVariables[ihist], mode, -1, 1, NULL);
          histogram->Write();
          stack->Add(histogram);
          delete histogram;
        }
        canvas->Draw();
        canvas->cd();
        stack->Draw();
        //stack->Write();
        delete stack;
        canvas->Write();
        delete canvas;
      }
    }
  }
  else{ 
    for(unsigned ipdf = 0 ; ipdf < SamplePDFs.size() ; ipdf++) {
      for(int ihist=0; ihist<HistVariables.size(); ihist++){
        for(int imode = 0; imode<PlotModes.size(); imode++){
          int mode;
          if(DUNEString_ToMaCh3Mode(PlotModes[0])==(int)(kMaCh3_nModes)){mode =-1;}
          else{mode = (int)DUNEString_ToMaCh3Mode(PlotModes[imode]); std::cout<<"mode: "<<mode<<std::endl;}
          TH1D * histogram = (TH1D*)SamplePDFs[ipdf]->get1DVarHist(HistVariables[ihist], mode, -1, Weighted, NULL);
          histogram ->Write();
          delete histogram;
        }
      }
    }
  }*/ 
  std::cout<<"before closing outputfile"<<std::endl;
  OutputFile->Close();
  std::cout<<"after closing outputfile"<<std::endl;
  return 0;
}
