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

#include "samplePDFDUNE/samplePDFDUNEAtmBase.h"
#include "manager/manager.h"

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


int main (int argc, char *argv[])
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(2);
  std::cout.setf(ios::showpos);

  bool eval_xsec = true;
  bool eval_flux = false;
  bool eval_skdet = false;

  manager *fitMan = new manager(argv[1]);

  // Get Plot by Mode boolean
  bool PlotByMode = false;


  // make covariance objects 
  std::string  XsecMatrixFile = fitMan->raw()["General"]["Systematics"]["XsecCovFile"].as<std::string>(); 
  std::string  XsecMatrixName = fitMan->raw()["General"]["Systematics"]["XsecCovName"].as<std::string>();
  std::string  OscMatrixFile = fitMan->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
  std::string  OscMatrixName = fitMan->raw()["General"]["Systematics"]["OscCovName"].as<std::string>();
 
  std::string OutfileName = fitMan->raw()["General"]["Output"]["OUTPUTNAME"].as<std::string>();

  // Print to output so can look up later                                                                             
  std::cout << "Covariance files used: " << std::endl
	<< XsecMatrixFile.c_str() << std::endl;


  //  covarianceFlux* flux              = new covarianceFlux("total_flux_cov", fluxCovMatrixFile) ;
  covarianceXsec *xsec          = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str()) ;
  //covarianceSkDet_joint* skdet      = new covarianceSkDet_joint("SKJointError_Erec_Total", skDetCovMatrixFile);
  //  covarianceSkDet_joint* skdet      = new covarianceSkDet_joint(skDetCovMatrixName, skDetCovMatrixFile);

  // Make PDFs
  std::vector<char*> names;
  std::vector<samplePDFDUNEAtmBase*>pdfs; 
  samplePDFDUNEAtmBase *fhc_nue_pdf = new samplePDFDUNEAtmBase(1.3628319e+23, "configs/SamplePDFDune_FHC_nueselec.yaml",xsec);
  samplePDFDUNEAtmBase *fhc_numu_pdf = new samplePDFDUNEAtmBase(1.3628319e+23, "configs/SamplePDFDune_FHC_numuselec.yaml",xsec);
  samplePDFDUNEAtmBase *rhc_nue_pdf = new samplePDFDUNEAtmBase(1.3628319e+23, "configs/SamplePDFDune_RHC_nueselec.yaml",xsec);
  samplePDFDUNEAtmBase *rhc_numu_pdf = new samplePDFDUNEAtmBase(1.3628319e+23, "configs/SamplePDFDune_RHC_numuselec.yaml",xsec);

  pdfs.push_back(fhc_nue_pdf);
  pdfs.push_back(fhc_numu_pdf);
  pdfs.push_back(rhc_nue_pdf);
  pdfs.push_back(rhc_numu_pdf);
  names.push_back("fhc_nue");
  names.push_back("fhc_numu");
  names.push_back("rhc_nue");
  names.push_back("rhc_numu");

  std::vector<double> flux_postNDfit;
  std::vector<double> xsec_postNDfit;
  std::vector<double> BANFF_untuned;

  //  flux->setParameters();
  xsec->setParameters(); 
  //  skdet->setParameters();

  std::cout << "---------- Printing off nominal parameter values ----------" << std::endl;
  std::cout << "Cross section parameters:" << std::endl;
  xsec->printNominal();
  //  std::cout << "Flux parameters:" << std::endl;
  //  flux->printNominal();
  //  std::cout << "SK detector parameters:" << std::endl;
  //  skdet->printNominal();
  //  std::cout << "---------- Finished printing nominal parameter values ----------" << std::endl;

  std::cout<<"Setting covariances in samplepdfs"<<std::endl;
  // tell PDFs where covariances are and start setting up oscillation parameters
  for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
	//	pdfs[ipdf]->setFluxCov(flux);
	//	pdfs[ipdf]->setSkDetCov(skdet);
	//pdfs[ipdf]->setXsecCov(xsec);
	//pdfs[ipdf]->useNonDoubledAngles(true);
  }



  ///////////
  // setup osc stuff
  ///////////

  // set up oscillation parameters (init params are in constructor)
  covarianceOsc *osc = new covarianceOsc("osc_cov","inputs/oscillation_covariance_6par_nondouble_PDG2019.root");
  // oscpars from manager in order:
  // sin2th_12, sin2th_23, sin2th_13, delm2_12, delm2_23, delta_cp
  std::vector<double> oscpars{0.307,0.528,0.0218,7.53e-5, 2.509e-3,-1.601}; // Asimov A

  std::cout<<"Using these oscillation parameters: ";
  for(unsigned ipar=0;ipar<oscpars.size();ipar++) std::cout<<" "<<oscpars.at(ipar);
  std::cout << std::endl;

  // Ask config file whether to use reactor constraint
  bool useRC =  false;
  std::cout << "use reactor prior is : " << useRC << std::endl ;
  osc->useReactorPrior(useRC); // this is hard coded inside, and is bad

  // Use prior for 12 parameters only
  //osc->setEvalLikelihood(0,false);
  osc->setEvalLikelihood(1,false);
  osc->setEvalLikelihood(2,false);
  //osc->setEvalLikelihood(3,false);
  osc->setEvalLikelihood(4,false);
  osc->setEvalLikelihood(5,false);
  osc->setFlipDeltaM23(true);
  // This line gives a crash and stack trace...
  osc->setParameters(oscpars);
  //  osc->acceptStep();
  //  osc->setStepScale(fitMan->getOscStepScale());

  osc -> setParameters(oscpars);
  std::cout << "oscpars[0] = " << (osc -> getPropPars())[0] << std::endl
	<< "oscpars[1] = " << (osc -> getPropPars())[1] << std::endl
	<< "oscpars[2] = " << (osc -> getPropPars())[2] << std::endl
	<< "oscpars[3] = " << (osc -> getPropPars())[3] << std::endl
	<< "oscpars[4] = " << (osc -> getPropPars())[4] << std::endl
	<< "oscpars[5] = " << (osc -> getPropPars())[5] << std::endl;

  osc->setFlipDeltaM23(true);

  // Unoscillated
  std::vector<double> oscpars_un(oscpars);
  oscpars_un[0] = 0.;
  oscpars_un[1] = 0.;
  oscpars_un[2] = 0.;
  

  osc->setParameters(oscpars_un);
  std::cout << "Oscillated event rates" << std::endl;
  std::vector<TH1D*> osc_hists;
  std::vector<double> CAFana_default_edges = get_default_CAFana_bins();
  for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
    pdfs[ipdf] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
        pdfs[ipdf]->useBinnedOscReweighting(true, CAFana_default_edges.size()-1, &CAFana_default_edges[0]);
	pdfs[ipdf]->reweight(osc->getPropPars(), osc->getPropPars());
	std::string names_osc = std::string(names[ipdf]) + std::string("_osc");
	osc_hists.push_back((TH1D*)pdfs[ipdf]->get1DHist()->Clone(names_osc.c_str()));
  }




  std::cout << "Now using setting osc pars to unoscillated" << std::endl;

  osc->setParameters(oscpars_un);

  for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
	pdfs[ipdf]->reweight(osc->getPropPars(), osc->getPropPars());
  }

  //Make output files
  std::string outfilename = OutfileName.c_str();
  std::size_t pos;
  pos = outfilename.find(".root");
  outfilename = outfilename.substr(0,pos);

  std::string xsec_filename = outfilename + "_xsec" + ".root";
  //  std::string flux_filename = outfilename + "_flux" + ".root";
  //  std::string skdet_filename = outfilename + "_skdet" + ".root";

  TFile *outputfile_xsec;
  //  TFile *outputfile_flux;
  // TFile *outputfile_skdet;

  if (eval_xsec) outputfile_xsec = new TFile(xsec_filename.c_str(),"RECREATE");
  //  if (eval_flux) outputfile_flux = new TFile(flux_filename.c_str(),"RECREATE");
  //  if (eval_skdet) outputfile_skdet = new TFile(skdet_filename.c_str(),"RECREATE");

  //std::cout << "how about some nominal spectra?" << std::endl;
  // Make and save nominal spectra
  std::vector<TH1D*> nominal_hists;
  for(unsigned ipdf = 0 ; ipdf < pdfs.size() ; ipdf++){
	nominal_hists.push_back((TH1D*)pdfs[ipdf]->get1DHist()->Clone(names[ipdf]));
  }


  if (eval_xsec)
  {
	std::cout << "EVALUATING XSEC PARAMS" << std::endl;
	outputfile_xsec->cd();
	for(unsigned ipdf = 0 ; ipdf < pdfs.size() ; ipdf++){
	  nominal_hists[ipdf]->Write();
	  osc_hists[ipdf]->Write();
	}
  }

  // set back to oscillated spectra
  osc->setParameters(oscpars_un);

  for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
	pdfs[ipdf]->reweight(osc->getPropPars(), osc->getPropPars());
  }


  const int nsubsamples = 12;


  int nmodes = kMaCh3_nModes;
  std::string mode_names[nmodes];
  for(unsigned i=0;i<kMaCh3_nModes;i++){
	mode_names[i] = MaCh3mode_ToDUNEString((MaCh3_Mode)i);
	//	std::cout << mode_names[i] << std::endl;
  }
  char *sample_names[nsubsamples];
  sample_names[0] = "numu";
  sample_names[1] = "nue";
  sample_names[2] = "numub";
  sample_names[3] = "nueb";
  sample_names[4] = "signue";
  sample_names[5] = "signueb";
  sample_names[6] = "signumu";
  sample_names[7] = "signumub";
  sample_names[8] = "nuenutau";
  sample_names[9] = "numunutau";
  sample_names[10] = "nuebnutaub";
  sample_names[11] = "numubnutaub";

  char houtname[100];

  std::vector<std::vector<std::vector<TH1D*>>> mode_nom;

  // std::vector<std::vector<TH1D*>> mode_vec;

  // std::vector<TH1D*> samp_vec;

  // std::cout << "mode name: " << mode_names[0] << std::endl;
  // sprintf(houtname,"%s_sm_%s_%s_nom",names[0], sample_names[1],mode_names[0].c_str());
  // std::cout << "Histo output name: " << houtname << std::endl;
  // TH1D * modehist = pdfs[0]->getModeHist1D(1, 0, 2);
  // modehist->SetName(houtname);
  // std::cout << "modehist integral: " << modehist->Integral() << ", " << modehist << std::endl;
  // samp_vec.push_back(modehist);
  // mode_vec.push_back(samp_vec);
  // mode_nom.push_back(mode_vec);

  // outputfile_xsec->cd();
  // mode_nom[0][0][0]->Write();
  // outputfile_xsec->Close();

  // Get all pdf/subsample/mode histograms
  /*
  for(unsigned ipdf = 0 ; ipdf < pdfs.size() ; ipdf++){
	//    std::cout << "Looping through pdfs " << std::endl;
	std::vector<std::vector<TH1D*>> mode_vec;

	for (int mode=0; mode<kMaCh3_nModes; mode++)
	{
	  //	std::cout << "Now looping through modes: " << mode_names[mode] << std::endl;
	  std::vector<TH1D*> samp_vec;

	  for (int samp=0; samp<nsubsamples; samp++)
	  { 
		//	    std::cout << "Looping through subsample: " << samp << std::endl;
		sprintf(houtname,"%s_sm_%s_%s_nom",names[ipdf], sample_names[samp],mode_names[mode].c_str());
		//	    std::cout << "histo: " << houtname << std::endl;
		TH1D * modehist = pdfs[ipdf]->getModeHist1D(samp, mode, 2);
		modehist->SetName(houtname);
		//	    std::cout << "modehist integral: " << modehist->Integral() << ", " << modehist << std::endl;
		samp_vec.push_back(modehist);

	  }
	  mode_vec.push_back(samp_vec);

	}
	mode_nom.push_back(mode_vec);

  } */

  //  Write them to output file
  /*
  for(unsigned ipdf = 0 ; ipdf < pdfs.size() ; ipdf++){
	for (int mode=0; mode<kMaCh3_nModes; mode++){
	  for (int samp=0; samp<nsubsamples; samp++){
		//  mode_nom[ipdf][mode][samp]->SetDirectory(0);  
		if(eval_xsec){
		  outputfile_xsec->cd();
		  //std::cout << "writing output: " << mode_nom[ipdf][mode][samp]->GetName() << std::endl;
		  mode_nom[ipdf][mode][samp]->Write();
		}
	  }
	}
  }
  */

  std::cout << setprecision(2);



  if (eval_xsec)
  { 

	outputfile_xsec->cd();

	// std::cout << "\\begin{table}[!h]" << std::endl;
	// //std::cout << "\\caption{Cross-section parameters}" << std::endl;
	// std::cout << "\\begin{tabular}{|c|c||c|c||c|c|c|c|} " << std::endl << "\\hline" << std::endl;
	// std::cout << "Parameter & 1$\\sigma$ value & Sample & N$_{asimov}$ & -3$\\sigma$ & -1$\\sigma$ & +1$\\sigma$ & +3$\\sigma$ \\\\ " << std::endl << "\\hline" << std::endl;

	std::vector<char*> xsec_names;

	// ------ Do Xsec Variations ----- //
	vector<double> xsecpar = xsec->getNominalArray();
	//      for (int i=0; i<int(xsecpar.size()); i++)

	std::vector<TH1D*> varied_hists;

	for (int i=0; i<int(xsecpar.size()); i++)
	{

	  char xsec_name[100];
	  sprintf(xsec_name, "%i : %s", i, xsec->getParName(i));
	  xsec_names.push_back(xsec_name);
	  xsecpar = xsec->getNominalArray();

	  std::cout << xsec_names[i] << " & " << sqrt((*xsec->getCovMatrix())(i,i)) << std::endl;

	  //std::vector<std::vector<std::vector<TH1D*>>> mode_var;
	  //vector to store varied histogram
	  //	  std::vector<TH1D*> var_hists;
	  //	  var_hists.reserve(pdfs.size());

	  for(int ipdf = 0 ; ipdf < pdfs.size() ; ipdf++){
		// 	std::vector<std::vector<TH1D*>> mode_vec;
		// 	mode_var.push_back(mode_vec);

		// 	double sample_val[4];
		// 	for(int v = 0 ; v < 4 ; v++){
		// 	  sample_val[4] = 0.;
		// 	}


		//Now loop over the sigma values you want to evaluate at (usual -3, -1, 0, +1, +3)
		for (int j=-3; j<=3; j+=2){
		  char sign = (j>0) ? 'p' : 'n';

		  //Set xsec par to the value
		  xsecpar[i] = xsec->getNominal(i)+j*sqrt((*xsec->getCovMatrix())(i,i));
		  xsec->setParameters(xsecpar);
		  //Reweight prediction
		  pdfs[ipdf]->reweight(osc->getPropPars(), osc->getPropPars());

		  sprintf(houtname,"%s_xsec_par_%i_sig_%c%i",names[ipdf],i,sign,abs(j));
		  TH1D * hist = (TH1D*)pdfs[ipdf]->get1DHist()->Clone(houtname);
		  //		var_hists[ipdf] = (TH1D*)pdfs[ipdf]->get1DHist()->Clone(houtname);
		  std::cout << "Setting xsec par to value: " << xsecpar[i] << " histo: " << houtname << std::endl;
/*
                  if (PlotByMode) { 
		    //Loop over modes
		    for(int mode = 0 ; mode < kMaCh3_nModes ; mode++){
			  // 		std::vector<TH1D*> samp_vec;
			  // 		mode_var[ipdf].push_back(samp_vec);

			  //Loop over oscillaiton channels
			  for (int samp = 0 ; samp < 12 ; samp++){
			    // 		  mode_var[ipdf][mode].push_back(pdfs[ipdf]->getModeHist1D(samp, mode, 2));
			    // 		  mode_var[ipdf][samp][mode]->Write(houtname);
			    sprintf(houtname,"%s_xsec_par_%i_sm_%s_%s_sig_%c%i", names[ipdf],i,sample_names[samp],mode_names[mode].c_str(),sign,abs(j));
			    TH1D * histmode = (TH1D*)pdfs[ipdf]->getModeHist1D(samp, mode, 2)->Clone(houtname);

			    if(histmode->Integral() == 0){continue;}
			    histmode->Write();
			  } 

		    }

		  }
*/
                  hist->Write();
		  //var_hists[ipdf]->Write();

		}// loop over sigma values

	  }//end of loop over pdfs

	}//end of loop over xsec systs

  }//end of eval_xsec

  /*


  
  double numu_val[4], nue_val[4], numub_val[4], nueb_val[4], nue1pi_val[4];
  for (int v=0; v<4; v++)
  {
  numu_val[v] = 0.0;
  nue_val[v] = 0.0;
  numub_val[v] = 0.0;
  nueb_val[v] = 0.0;
  nue1pi_val[v] = 0.0;
  }

  for (int j=-3; j<=3; j+=2)
  {


  xsecpar[i] = xsec->getNominal(i)+j*sqrt((*xsec->getCovMatrix())(i,i));

  //std::cout << xsec->getName() << " parameter " << i << " with value " << xsecpar[i] << std::endl;
  xsec->setParameters(xsecpar);
  for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
  pdfs[ipdf]->reweight(osc->getPropPars());
  }

  hsknumu_var=numu_pdf->get1DHist();
  hsknue_var=nue_pdf->get1DHist();
  hsknumubar_var=numubar_pdf->get1DHist();
  hsknuebar_var=nuebar_pdf->get1DHist();
  hsknue1pi_var=nue1pi_pdf->get1DHist();

  for (int mode=0; mode<kMaCh3_nModes; mode++)
  {
  for (int samp=0; samp<6; samp++)
  {
  mode_var_numu[samp][mode] = numu_pdf->getModeHist1D(samp,mode,2);
  mode_var_numubar[samp][mode] = numubar_pdf->getModeHist1D(samp,mode,2);
  mode_var_nue[samp][mode] = nue_pdf->getModeHist1D(samp,mode,2);
  mode_var_nuebar[samp][mode] = nuebar_pdf->getModeHist1D(samp,mode,2);
  mode_var_nue1pi[samp][mode] = nue1pi_pdf->getModeHist1D(samp,mode,2);
  }
  }

  hsknumu_var->SetDirectory(0);
  hsknue_var->SetDirectory(0);
  hsknumubar_var->SetDirectory(0);
  hsknuebar_var->SetDirectory(0);
  hsknue1pi_var->SetDirectory(0);

  char sign = (j>0) ? 'p' : 'n';

  sprintf(houtname,"sknumu_xsec_par_%i_sig_%c%i",i,sign,abs(j));
  hsknumu_var->Write(houtname);
  sprintf(houtname,"sknue_xsec_par_%i_sig_%c%i",i,sign,abs(j));
  hsknue_var->Write(houtname);
  sprintf(houtname,"anu_sknumubar_xsec_par_%i_sig_%c%i",i,sign,abs(j));
  hsknumubar_var->Write(houtname);
  sprintf(houtname,"anu_sknuebar_xsec_par_%i_sig_%c%i",i,sign,abs(j));
  hsknuebar_var->Write(houtname);
  sprintf(houtname,"sknue1pi_xsec_par_%i_sig_%c%i",i,sign,abs(j));
  hsknue1pi_var->Write(houtname);

  //std::cout << "kMaCh3_nModes = " << kMaCh3_nModes << std::endl;

  for (int mode=0; mode<kMaCh3_nModes; mode++)
  {
  for (int samp=0; samp<6; samp++)
  {
  //std::cout << "mode, samp = " << mode << ", " << samp << std::endl;
  sprintf(houtname_numu,"sknumu_xsec_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));
  sprintf(houtname_numubar,"anu_sknumubar_xsec_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));
  sprintf(houtname_nue,"sknue_xsec_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));
  sprintf(houtname_nuebar,"anu_sknuebar_xsec_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));
  sprintf(houtname_nue1pi,"sknue1pi_xsec_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));

  mode_var_numu[samp][mode]->SetDirectory(0);
  mode_var_numubar[samp][mode]->SetDirectory(0);
  mode_var_nue[samp][mode]->SetDirectory(0);
  mode_var_nuebar[samp][mode]->SetDirectory(0);
  mode_var_nue1pi[samp][mode]->SetDirectory(0);

  mode_var_numu[samp][mode]->Write(houtname_numu);
  mode_var_numubar[samp][mode]->Write(houtname_numubar);
  mode_var_nue[samp][mode]->Write(houtname_nue);
  mode_var_nuebar[samp][mode]->Write(houtname_nuebar);
  mode_var_nue1pi[samp][mode]->Write(houtname_nue1pi);
}
}

int v;
if (j==-3) v=0;
if (j==-1) v=1;
if (j==1) v=2;
if (j==3) v=3;
numu_val[v] = hsknumu_var->Integral();
nue_val[v] = hsknue_var->Integral();
numub_val[v] = hsknumubar_var->Integral();
nueb_val[v] = hsknuebar_var->Integral();
nue1pi_val[v] = hsknue1pi_var->Integral();
}


if (!(numu_val[1]==numu_val[2] && numu_val[1]==numu_val[3] && numu_val[2]==numu_val[3]))
  std::cout << " & & FHC 1R$_{\\mu}$ & " << n_numu << " & " << (numu_val[0]-n_numu)/n_numu*100 << "\\% & " << (numu_val[1]-n_numu)/n_numu*100.0 << "\\% & " << (numu_val[2]-n_numu)/n_numu*100.0 << "\\% & " << (numu_val[3]-n_numu)/n_numu*100.0 << "\\% \\\\ " << "\\hline" << std::endl;
if (!(nue_val[1]==nue_val[2] && nue_val[1]==nue_val[3] && nue_val[2]==nue_val[3]))
  std::cout << " & & FHC 1R$_e$ & " << n_nue << " & " << (nue_val[0]-n_nue)/n_nue*100.0 << "\\% & " << (nue_val[1]-n_nue)/n_nue*100.0 << "\\% & " << (nue_val[2]-n_nue)/n_nue*100.0 << "\\% & " << (nue_val[3]-n_nue)/n_nue*100.0 << "\\% \\\\ " << "\\hline" << std::endl;
if (!(numub_val[1]==numub_val[2] && numub_val[1]==numub_val[3] && numub_val[2]==numub_val[3]))
  std::cout << " & & RHC 1R$_{\\mu}$ & " << n_numubar << " & " << (numub_val[0]-n_numubar)/n_numubar*100.0 << "\\% & " << (numub_val[1]-n_numubar)/n_numubar*100.0 << "\\% & " << (numub_val[2]-n_numubar)/n_numubar*100.0 << "\\% & " << (numub_val[3]-n_numubar)/n_numubar*100.0 << "\\% \\\\ " << "\\hline" << std::endl;
if (!(nueb_val[1]==nueb_val[2] && nueb_val[1]==nueb_val[3] && nueb_val[2]==nueb_val[3]))  
  std::cout << " & & RHC 1R$_e$ & " << n_nuebar << " & "  << (nueb_val[0]-n_nuebar)/n_nuebar*100.0 << "\\% & " << (nueb_val[1]-n_nuebar)/n_nuebar*100.0 << "\\% & " << (nueb_val[2]-n_nuebar)/n_nuebar*100.0 << "\\% & " << (nueb_val[3]-n_nuebar)/n_nuebar*100.0 << "\\% \\\\ " << "\\hline" << std::endl;
if (!(nue1pi_val[1]==nue1pi_val[2] && nue1pi_val[1]==nue1pi_val[3] && nue1pi_val[2]==nue1pi_val[3]))
  std::cout << " & & FHC 1R$_e$ 1$\\pi$ & " << n_nue1pi << " & " << (nue1pi_val[0]-n_nue1pi)/n_nue1pi*100.0 << "\\% & " << (nue1pi_val[1]-n_nue1pi)/n_nue1pi*100.0 << "\\% & " << (nue1pi_val[2]-n_nue1pi)/n_nue1pi*100.0 << "\\% & " << (nue1pi_val[3]-n_nue1pi)/n_nue1pi*100.0 << "\\% \\\\ " << "\\hline" << std::endl;	  

  // Split tables so they fit on a page
if (i==9 || i==20 || i==32)
{
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{table}" << std::endl << std::endl << "\\begin{table}[!h]" << std::endl << "\\small" << std::endl << "\\centering" << std::endl;
  //std::cout << "\\caption{Cross-section parameters}" << std::endl;
  std::cout << "\\begin{tabular}{|c|c||c|c||c|c|c|c|} " << std::endl << "\\hline" << std::endl;
  std::cout << "Parameter & 1$\\sigma$ value & Sample & N$_{asimov}$ & -3$\\sigma$ & -1$\\sigma$ & +1$\\sigma$ & +3$\\sigma$ \\\\ " << std::endl << "\\hline" << std::endl;
}
}

std::cout << "\\end{tabular}" << std::endl;
std::cout << "\\end{table}" << std::endl << std::endl;

xsec->setParameters();
outputfile_xsec->Close();
}
*/

/*
   if (eval_flux) // ----- Do Flux Variations ----- //
   {     
   outputfile_flux->cd(); 
   std::cout << "\\begin{table}[!h]" << std::endl << "\\small" << std::endl << "\\centering" << std::endl;
   std::cout << "\\caption{Flux parameters}" << std::endl;
   std::cout << "\\begin{tabular}{|c|c||c|c||c|c|c|c|} " << std::endl << "\\hline" << std::endl;
   std::cout << "Parameter & 1$\\sigma$ value & Sample & N$_{asimov}$ & -3$\\sigma$ & -1$\\sigma$ & +1$\\sigma$ & +3$\\sigma$ \\\\ " << std::endl << "\\hline" << std::endl;

   vector<double> fluxpar = flux->getNominalArray();
   if(BANFFtuned){
   fluxpar = flux_postNDfit;
   }
   else{
   fluxpar= flux->getNominalArray();
   }
   for (int i=50; i<int(fluxpar.size()); i++)
   {
   if(BANFFtuned){
   fluxpar = flux_postNDfit;
   }
   else{
   fluxpar= flux->getNominalArray();
   }

   if(BANFFtuned){
   std::cout << "Flux " << i << " & " << sqrt((*NDpostfitcov)(i-50,i-50)) << " & & & & & & \\\\ " << "\\hline" << std::endl;
   }
   else{
   std::cout << "Flux " << i << " & " << sqrt((*flux->getCovMatrix())(i,i)) << " & & & & & & \\\\ " << "\\hline" << std::endl;
   }
   double numu_val[4], nue_val[4], numub_val[4], nueb_val[4], nue1pi_val[4];
   for (int v=0; v<4; v++)
   {
   numu_val[v] = 0.0;
   nue_val[v] = 0.0;
   numub_val[v] = 0.0;
   nueb_val[v] = 0.0;
   nue1pi_val[v] = 0.0;
   }

   for (int j=-3; j<=3; j+=2)
   {
   if(BANFFtuned){
   fluxpar[i] = flux_postNDfit[i-50]+j*sqrt((*NDpostfitcov)(i-50,i-50));
   }
   else{
   fluxpar[i] = flux->getNominal(i)+j*sqrt((*flux->getCovMatrix())(i,i));
   }
//std::cout << flux->getName() << " parameter " << i << " with value " << fluxpar[i] << std::endl;
flux->setParameters(fluxpar);
for(unsigned ipdf=0;ipdf<pdfs.size();ipdf++){
pdfs[ipdf]->reweight(osc->getPropPars());
}


hsknumu_var=numu_pdf->get1DHist();
hsknue_var=nue_pdf->get1DHist();
hsknumubar_var=numubar_pdf->get1DHist();
hsknuebar_var=nuebar_pdf->get1DHist();
hsknue1pi_var=nue1pi_pdf->get1DHist();

for (int mode=0; mode<kMaCh3_nModes; mode++)
{
for (int samp=0; samp<6; samp++)
{
mode_var_numu[samp][mode] = numu_pdf->getModeHist1D(samp,mode,2);
mode_var_numubar[samp][mode] = numubar_pdf->getModeHist1D(samp,mode,2);
mode_var_nue[samp][mode] = nue_pdf->getModeHist1D(samp,mode,2);
mode_var_nuebar[samp][mode] = nuebar_pdf->getModeHist1D(samp,mode,2);
mode_var_nue1pi[samp][mode] = nue1pi_pdf->getModeHist1D(samp,mode,2);
}
}


hsknumu_var->SetDirectory(0);
hsknue_var->SetDirectory(0);
hsknumubar_var->SetDirectory(0);
hsknuebar_var->SetDirectory(0);
hsknue1pi_var->SetDirectory(0);

char sign = (j>0) ? 'p' : 'n';

sprintf(houtname,"sknumu_flux_par_%i_sig_%c%i",i,sign,abs(j));
hsknumu_var->Write(houtname);
sprintf(houtname,"sknue_flux_par_%i_sig_%c%i",i,sign,abs(j));
hsknue_var->Write(houtname);
sprintf(houtname,"anu_sknumubar_flux_par_%i_sig_%c%i",i,sign,abs(j));
hsknumubar_var->Write(houtname);
sprintf(houtname,"anu_sknuebar_flux_par_%i_sig_%c%i",i,sign,abs(j));
hsknuebar_var->Write(houtname);
sprintf(houtname,"sknue1pi_flux_par_%i_sig_%c%i",i,sign,abs(j));
hsknue1pi_var->Write(houtname);

for (int mode=0; mode<kMaCh3_nModes; mode++)
{
  for (int samp=0; samp<6; samp++)
  {
	sprintf(houtname_numu,"sknumu_flux_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));
	sprintf(houtname_numubar,"anu_sknumubar_flux_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));
	sprintf(houtname_nue,"sknue_flux_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));
	sprintf(houtname_nuebar,"anu_sknuebar_flux_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));
	sprintf(houtname_nue1pi,"sknue1pi_flux_par_%i_sm_%s_%s_sig_%c%i",i,sample_names[samp],mode_names[mode],sign,abs(j));

	if (mode==kMaCh3_CCMpi) // combine CC multipi and CC other
	  continue;
	else
	{
	  if (mode==kMaCh3_CCDIS)
	  {
		mode_var_numu[samp][kMaCh3_CCDIS]->Add(mode_var_numu[samp][kMaCh3_CCMpi]);
		mode_var_numubar[samp][kMaCh3_CCDIS]->Add(mode_var_numubar[samp][kMaCh3_CCMpi]);
		mode_var_nue[samp][kMaCh3_CCDIS]->Add(mode_var_nue[samp][kMaCh3_CCMpi]);
		mode_var_nuebar[samp][kMaCh3_CCDIS]->Add(mode_var_nuebar[samp][kMaCh3_CCMpi]);
		mode_var_nue1pi[samp][kMaCh3_CCDIS]->Add(mode_var_nue1pi[samp][kMaCh3_CCMpi]);
	  }

	  mode_var_numu[samp][mode]->SetDirectory(0);
	  mode_var_numubar[samp][mode]->SetDirectory(0);
	  mode_var_nue[samp][mode]->SetDirectory(0);
	  mode_var_nuebar[samp][mode]->SetDirectory(0);
	  mode_var_nue1pi[samp][mode]->SetDirectory(0);

	  mode_var_numu[samp][mode]->Write(houtname_numu);
	  mode_var_numubar[samp][mode]->Write(houtname_numubar);
	  mode_var_nue[samp][mode]->Write(houtname_nue);
	  mode_var_nuebar[samp][mode]->Write(houtname_nuebar);
	  mode_var_nue1pi[samp][mode]->Write(houtname_nue1pi);
	}
  }
}


int v;
if (j==-3) v=0;
if (j==-1) v=1;
if (j==1) v=2;
if (j==3) v=3;
numu_val[v] = hsknumu_var->Integral();
nue_val[v] = hsknue_var->Integral();
numub_val[v] = hsknumubar_var->Integral();
nueb_val[v] = hsknuebar_var->Integral();
nue1pi_val[v] = hsknue1pi_var->Integral();
}

if (!(numu_val[1]==numu_val[2] && numu_val[1]==numu_val[3] && numu_val[2]==numu_val[3]))
  std::cout << " & & FHC 1R$_{\\mu}$ & " << n_numu << " & " << (numu_val[0]-n_numu)/n_numu*100 << "\\% & " << (numu_val[1]-n_numu)/n_numu*100.0 << "\\% & " << (numu_val[2]-n_numu)/n_numu*100.0<< "\\% & " << (numu_val[3]-n_numu)/n_numu*100.0 << "\\% \\\\ " << "\\hline" << std::endl;
if (!(nue_val[1]==nue_val[2] && nue_val[1]==nue_val[3] && nue_val[2]==nue_val[3]))
  std::cout << " & & FHC 1R$_e$ & " << n_nue << " & " << (nue_val[0]-n_nue)/n_nue*100.0 << "\\% & " << (nue_val[1]-n_nue)/n_nue*100.0 << "\\% & " << (nue_val[2]-n_nue)/n_nue*100.0 << "\\% & " << (nue_val[3]-n_nue)/n_nue*100.0 << "\\% \\\\ " << "\\hline" << std::endl;
if (!(numub_val[1]==numub_val[2] && numub_val[1]==numub_val[3] && numub_val[2]==numub_val[3]))
  std::cout << " & & RHC 1R$_{\\mu}$ & " << n_numubar << " & " << (numub_val[0]-n_numubar)/n_numubar*100.0 << "\\% & " << (numub_val[1]-n_numubar)/n_numubar*100.0 << "\\% & " << (numub_val[2]-n_numubar)/n_numubar*100.0 << "\\% & " << (numub_val[3]-n_numubar)/n_numubar*100.0 << "\\% \\\\ " << "\\hline" << std::endl;
if (!(nueb_val[1]==nueb_val[2] && nueb_val[1]==nueb_val[3] && nueb_val[2]==nueb_val[3]))
  std::cout << " & & RHC 1R$_e$ & " << n_nuebar << " & "  << (nueb_val[0]-n_nuebar)/n_nuebar*100.0 << "\\% & " << (nueb_val[1]-n_nuebar)/n_nuebar*100.0 << "\\% & " << (nueb_val[2]-n_nuebar)/n_nuebar*100.0 << "\\% & " << (nueb_val[3]-n_nuebar)/n_nuebar*100.0 << "\\% \\\\ " << "\\hline" << std::endl;
if (!(nue1pi_val[1]==nue1pi_val[2] && nue1pi_val[1]==nue1pi_val[3] && nue1pi_val[2]==nue1pi_val[3]))
  std::cout << " & & FHC 1R$_e$ 1$\\pi$ & " << n_nue1pi << " & " << (nue1pi_val[0]-n_nue1pi)/n_nue1pi*100.0 << "\\% & " << (nue1pi_val[1]-n_nue1pi)/n_nue1pi*100.0 << "\\% & " << (nue1pi_val[2]-n_nue1pi)/n_nue1pi*100.0 << "\\% & " << (nue1pi_val[3]-n_nue1pi)/n_nue1pi*100.0 << "\\% \\\\ " << "\\hline" << std::endl;

  // Split tables so they fit on a page
if (i==64 || i == 82)
{
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{table}" << std::endl << std::endl << "\\begin{table}[!h]" << std::endl << "\\small" << std::endl << "\\centering" << std::endl;
  std::cout << "\\caption{Flux parameters}" << std::endl;
  std::cout << "\\begin{tabular}{|c|c||c|c||c|c|c|c|} " << std::endl << "\\hline" << std::endl;
  std::cout << "Parameter & 1$\\sigma$ value & Sample & N$_{asimov}$ & -3$\\sigma$ & -1$\\sigma$ & +1$\\sigma$ & +3$\\sigma$ \\\\ " << std::endl << "\\hline" << std::endl;
}
}   

std::cout << "\\end{tabular}" << std::endl;
std::cout << "\\end{table}" << std::endl << std::endl;

flux->setParameters();
outputfile_flux->Close();
}
*/
outputfile_xsec->Close();
return 0.;
}
