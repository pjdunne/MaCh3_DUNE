// Compile with something like
//
// g++ `root-config --cflags` -g -o DrawComp DrawComp.cpp -I`root-config --incdir` `root-config --glibs --libs`
//
// Run on ND280 only fits to get xsec and beam posteriors prettified. Essentially a compiled makeTransferMatrixAll
//
// ./DrawComp will tell you how
//
// Should probably prettify this into something managable like classes, or at least some structs
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "TObjArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TColor.h"

// The input cross-section covariance matrix
std::string XsecCov = "";

// Check that xsecCov is set
void CheckXsecCov();

// Helper to put a name on the cross-section plots
std::string GetXsecName(int xsecParamNo);
// Helper to get the limits for each cross-section parameter
void GetXsecLimits(int param, double &central, double &prior, double &down_error, double &up_error);
// Helper to plot the prefit limits
TH1D* MakePrefit(int ndraw, int nFlux);

// Get the highest posterior density numbers from a 1D posterior
void GetHPD(TH1D * const post, double &central, double &error, double &error_pos, double &error_neg);
void GetArithmetic(TH1D * const hpost, double &mean, double &error);
void GetGaussian(TH1D *& hpost, TF1 *& gauss, double &central, double &error);

// The main comparison functions
void DrawComp(std::string inputfile, bool drawcorr = false);
void DrawComp(std::string inputFile1, std::string title1, std::string inputFile2, std::string title2);
void DrawComp(std::string inputFile1, std::string title1, std::string inputFile2, std::string title2, std::string inputFile3, std::string title3);

// Decalring global variables to avoid copy pasting it 3 times.
// Probably take these as user input instead!
int BurnInCut = 20; //We cut 1/BurnInCut steps, we only consider steps after stationary state is achieved
bool plotFlux = false; //Plot Xsec prefit postfit
bool plotXsec = true; // You almost certainly want this to be true
bool plotCovFull = false; //If true will plot full (xsec and flux) covariance 

int main(int argc, char *argv[]) {

  std::cout << "******************************************" << std::endl;
  std::cout << "Will plot Xsec prefit postfit?            " << plotXsec << std::endl;
  std::cout << "Will plot Flux prefit postfit?            " << plotFlux << std::endl;
  std::cout << "Will plot covariance with flux and xsec?  " << plotCovFull << std::endl;
  std::cout << "******************************************" << std::endl;
  
  if (argc != 2 && argc != 3 && argc !=5 && argc != 7) {
    std::cerr << "./DrawComp root_file_to_analyse.root drawcorr (root_file_to_compare_to1.root title root_file_to_compare_to2.root title)" << std::endl;
    exit(-1);
  }

  // Have implemented some good old hard-coding for comparing many fits; argc == 2 is 1 fit, argc == 3 is 2 fits, etc
  // Could additionally do fancy stuff; do we want beam+fsi+xsec params, do we want do draw a covariance matrix, etc
  if (argc == 2 || argc == 3) {

    std::cout << "Producing single fit output" << std::endl;
    std::string filename = argv[1];
    bool drawCorr = false;
    if (argc == 3) {
      drawCorr = true;
    }
    DrawComp(filename, drawCorr);

    // If we want to compare two fits (e.g. binning changes or introducing new params/priors)
  } else if (argc == 3+2*1) {

    std::cout << "Producing two fit comparison" << std::endl;
    std::string filename  = argv[1];
    std::string title1    = argv[2];
    std::string filename2 = argv[3];
    std::string title2    = argv[4];
    DrawComp(filename, title1, filename2, title2);

    // If we want to compare three fits (e.g. binning changes or introducing new params/priors)
  } else if (argc == 4+3*1) {

    std::cout << "Producing three fit comparison" << std::endl;
    std::string filename = argv[1];
    std::string title1    = argv[2];
    std::string filename2 = argv[3];
    std::string title2    = argv[4];
    std::string filename3 = argv[5];
    std::string title3    = argv[6];
    DrawComp(filename, title1, filename2, title2, filename3, title3);

  }
  return 0;
}


// All right, let's re-write this in a faster more efficient way!
void DrawComp(std::string inputFile, bool drawCorr) {

  // drawCorr decideds if we want to draw correlations or not
  std::cout << "File for study:       " << inputFile << std::endl;
  std::cout << "Draw correlations?    " << drawCorr << std::endl;

  // Open the chain
  TChain* chain = new TChain("posteriors","");
  chain->Add(inputFile.c_str());
  //MCMC Burn-in cut
  int cut = chain->GetMaximum("step")/BurnInCut;
  std::stringstream ss;
  ss << "step > " << cut;

  std::cout<<"MCMC has "<<chain->GetMaximum("step")<<" steps, burn-in is set to"<<cut<<std::endl;
  // What entries are we interested in?
  // Here specifically 20k steps in, and events that have a logL < 5000
  std::string stepcut = ss.str();

  // Get the list of branches
  TObjArray* brlis = (TObjArray*)chain->GetListOfBranches();
  // Get the number of branches
  int nbr = brlis->GetEntries();
  std::cout << "# of branches: " << nbr << std::endl;
  // Make an array of TStrings
  TString bnames[nbr];
  // Array of doubles containing the nominal values (first step for xsec)
  std::vector<double> nom;
  nom.resize(nbr);

  // Have a counter for how many of each systematic we have
  int ndraw = 0;
  int nflux = 0;
  int nxsec = 0;

  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  chain->SetBranchStatus("*", false);
  chain->SetBranchStatus("step", true);
  for (int i = 0; i < nbr; i++) {

    // Get the TBranch and its name
    TBranch* br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();

    // Get first entry of this branch (i.e. systematic); this is the first reconfigure so should be the fake-data which it was generated for (for xsec this is true...!)
    //if(bname.BeginsWith("LogL") || bname.BeginsWith("ndd_")) continue;
    if(bname.BeginsWith("LogL")) continue;

    // Chechk how many xsec and flux params we have
    if(bname.BeginsWith("xsec_")) {
	  nxsec++;
      //plotXsec = true;
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      ndraw++;
    }
      //Flux covariance is integrated in xsec covariance
      
    else if(bname.BeginsWith("b_")) {
      plotFlux = true;
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      ndraw++;
      nflux++;
    }
	/*
    } else if(bname.BeginsWith("ndd_")) {
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      ndraw++;
      } 
	  */
    }
    
  

  TFile *TempFile = new TFile(inputFile.c_str(), "open");
  TempFile->ls();
  //TTree* Settings = (TTree*)(TempFile->Get("Settings"));
  //Settings->ls();
  // Get the xsec covariance matrix
  //std::string *XsecInput = 0;
  std::string XsecInput = "/vols/t2k/users/ljw20/software/MaCh3_DUNE_datatype/MaCh3_DUNE/inputs/xsec_covariance_DUNE_systs_2023.root";
  //Settings->SetBranchAddress("XsecCov", &XsecInput);
  //Settings->GetEntry(0);
  
  //delete Settings;
  TempFile->Close();
  delete TempFile;

  //XsecCov = "../"+*XsecInput;
  XsecCov = XsecInput;

  // Get first entry in chain
  chain->GetEntry(0);
  nom.resize(nbr);

  /*
  //Count how many xsec and flux params we have
  //flux params name always starts with b_ 
  for(int i = 0; i < ndraw; ++i) {
	std::string tempString = GetXsecName(i);
	if (tempString.rfind("b_", 0) == 0) {
	  nflux++;
	}
	else nxsec++;
  }
*/
  
  std::cout << "# useful entries (flux, xsec, fsi): " << ndraw << std::endl;
  std::cout << "# useful entries (flux): " << nflux << std::endl;
  std::cout << "# useful entries (xsec): " << nxsec << std::endl;
  std::cout << "************************************************" << std::endl;

  gStyle->SetOptFit(111);

  // Open a TCanvas to write the posterior onto
  TCanvas* c0 = new TCanvas("c0", "c0", 0, 0, 1024, 1024);
  c0->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c0->SetTickx();
  c0->SetTicky();
  c0->SetBottomMargin(0.1);
  c0->SetTopMargin(0.1);
  c0->SetRightMargin(0.03);
  c0->SetLeftMargin(0.15);

  // Make sure we can read files located anywhere and strip the .root ending
  inputFile = inputFile.substr(0, inputFile.find(".root"));

  TString canvasname = inputFile;
  // Append if we're drawing correlations
  // Open bracket means we want to make a pdf file
  if (drawCorr) {
    canvasname += "_drawCorr.pdf[";
  } else {
    canvasname += "_drawPar.pdf[";
  }
  c0->Print(canvasname);

  // Once the pdf file is open no longer need to bracket
  canvasname.ReplaceAll("[","");

  // We fit with this Gaussian
  TF1 *gauss    = new TF1("gauss","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",-5,5);
  gauss->SetLineWidth(2);
  gauss->SetLineColor(kOrange-5);

  // Some TVectors with the means, errors and Gaussian parameters of the PDFs
  TVectorD* mean_vec  = new TVectorD(ndraw);
  TVectorD* err_vec   = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec = new TVectorD(ndraw);
  TVectorD* gaus_err_vec  = new TVectorD(ndraw); 
  TVectorD* HPD_mean_vec  = new TVectorD(ndraw);
  TVectorD* HPD_err_p_vec   = new TVectorD(ndraw); 
  TVectorD* HPD_err_vec   = new TVectorD(ndraw); 
  TVectorD* HPD_err_m_vec   = new TVectorD(ndraw); 

  TVectorD* xsec_nom  = new TVectorD(nxsec);

  int covBinning = ndraw;
  //If we only plot correlation/covariance between xsec (without flux)
  if(!plotCovFull && drawCorr)
  {
     covBinning = nxsec; 
  }
  // Only want to draw covariance of xsec parameters!
  TMatrixT<double>* covariance = new TMatrixT<double>(covBinning,covBinning);
  TMatrixT<double>* correlation = new TMatrixT<double>(covBinning,covBinning);
  for (int i = 0; i < covBinning; ++i) {
    for (int j = 0; j < covBinning; ++j) {
      (*covariance)(i,j) = 0.;
      (*correlation)(i,j) = 0.;
    }
  }

  // Output file to write to
  TString rootfilename = inputFile;
  if (drawCorr) {
    rootfilename += "_drawCorr.root";
  } else {
    rootfilename += "_drawPar.root";
  }

  // The output file
  TFile* file = new TFile(rootfilename, "RECREATE");
  file->cd();
  // Make a directory with the parameters
  TDirectory *params = file->mkdir("params");
  TDirectory *corr = gDirectory;
  if (drawCorr) {
    corr = file->mkdir("corr");
  }

  bool isXsec = false;

  // Remember which parameters we're varying
  bool *vary = new bool[ndraw]();
  for (int i=0; i < ndraw; ++i) {
    vary[i] = false;
  }

  // ndraw is number of draws we want to do
  for(int i = 0; i < ndraw; ++i) {
    isXsec = false;

    // If we're interested in drawing the correlations need to invoke another for loop
    if (drawCorr) {
      // Make an array of bools to book-keep what parameters are varied

      if(!plotCovFull && i>=nxsec) //Skip if we only want xsec
      {
         continue;
      }
      
      file->cd();
      int nbins = 70;

      double asimovLine = 0.0;

      std::string tempString = std::string(bnames[i]);
      if (i<ndraw-nflux) { //first go xsec then flux
      isXsec = true;

        tempString = GetXsecName(i);

        double central, prior, down, up;
        GetXsecLimits(i, central, prior, down, up);
        (*xsec_nom)(i) = central;
        asimovLine = prior;
      }
      else //now flux parameters
      {
        TString TempText= Form("b_%i", i-nxsec);
        tempString = TempText.Data();  
      }
      file->cd();

      // Get the maximum and minimum for the parameter
      double maximum = chain->GetMaximum(bnames[i]);
      double minimum = chain->GetMinimum(bnames[i]);

      // This holds the posterior density
      TH1D *hpost = new TH1D(bnames[i], bnames[i], nbins, minimum, maximum);
      hpost->SetMinimum(0);
      hpost->GetYaxis()->SetTitle("Steps");
      hpost->GetYaxis()->SetNoExponent(false);
      hpost->SetTitle(bnames[i]);

      // Project bnames[i] onto hpost, applying stepcut
      chain->Project(bnames[i], bnames[i], stepcut.c_str());

      // Apply one smoothing
      hpost->Smooth();

      // Get the characteristics of the hpost
      double mean, rms;
      GetArithmetic(hpost, mean, rms);
      double peakval, sigma_p, sigma_m, sigma_hpd;
      GetHPD(hpost, peakval, sigma_hpd, sigma_p, sigma_m);
      double gauss_mean, gauss_rms;
      GetGaussian(hpost, gauss, gauss_mean, gauss_rms);

      std::cout << mean << " +/- " << rms << " (" << peakval << "+/-" << sigma_hpd << " + " << sigma_p << " - " << sigma_m << ")" << " (" << gauss_mean << "+/-" << gauss_rms << ")" << std::endl;

      TLine *hpd = new TLine(peakval, hpost->GetMinimum(), peakval, hpost->GetMaximum());
      hpd->SetLineColor(kBlack);
      hpd->SetLineWidth(2);
      hpd->SetLineStyle(kSolid);

      // Make the legend
      TLegend *leg = new TLegend(0.12, 0.6, 0.5, 0.97);
      leg->SetTextSize(0.04);
      leg->AddEntry(hpost, Form("#splitline{PDF}{#mu = %.2f, #sigma = %.2f}", hpost->GetMean(), hpost->GetRMS()), "l");
      leg->AddEntry(gauss, Form("#splitline{Gauss}{#mu = %.2f, #sigma = %.2f}", gauss->GetParameter(1), gauss->GetParameter(2)), "l");
      leg->AddEntry(hpd, Form("#splitline{HPD}{#mu = %.2f, #sigma = %.2f (+%.2f-%.2f)}", peakval, sigma_hpd, sigma_p, sigma_m), "l");

      if (isXsec && ((*xsec_nom)(i)) != 0) {
        mean = mean / ((*xsec_nom)(i));
        rms = rms / ((*xsec_nom)(i));
        gauss_mean = gauss_mean / ((*xsec_nom)(i));
        gauss_rms = gauss_rms / ((*xsec_nom)(i));
        peakval = peakval / ((*xsec_nom)(i));
        sigma_hpd = sigma_hpd / ((*xsec_nom)(i));
        sigma_p = sigma_p / ((*xsec_nom)(i));
        sigma_m = sigma_m / ((*xsec_nom)(i));
      } else if (isXsec && ((*xsec_nom)(i)) == 0) {
        mean = mean + 1.0;
        gauss_mean = gauss_mean + 1.0;
        peakval = peakval + 1.0;
      }

      (*mean_vec)(i)      = mean;
      (*err_vec)(i)       = rms;
      (*gaus_mean_vec)(i) = gauss->GetParameter(1);
      (*gaus_err_vec)(i)  = gauss->GetParameter(2);
      (*HPD_mean_vec)(i)  = peakval;
      (*HPD_err_p_vec)(i) = sigma_p;
      (*HPD_err_vec)(i)   = sigma_hpd;
      (*HPD_err_m_vec)(i) = sigma_m;
      (*covariance)(i,i)  = rms*rms;
      (*correlation)(i,i) = 1.0;

      hpost->SetLineWidth(2);
      hpost->SetLineColor(kBlue-1);
      hpost->SetMaximum(hpost->GetMaximum()*1.5);
      hpost->SetTitle(tempString.c_str());
      hpost->GetXaxis()->SetTitle((std::string(hpost->GetTitle())+" rel. nom").c_str());

      // Now make the TLine for the asimov
      TLine *asimov = new TLine(asimovLine, hpost->GetMinimum(), asimovLine, hpost->GetMaximum());
      asimov->SetLineColor(kRed-3);
      asimov->SetLineWidth(2);
      asimov->SetLineStyle(kDashed);
      hpost->Draw();
      hpd->Draw("same");
      asimov->Draw("same");

      leg->AddEntry(asimov, Form("#splitline{Asimov}{x = %.2f}", asimovLine), "l");
      leg->SetLineColor(0);
      leg->SetLineStyle(0);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->Draw("same");

      // Don't plot if this is a fixed histogram (i.e. the peak is the whole integral)
      /*
      if (hpost->GetMaximum() == hpost->Integral()*1.5) {
        vary[i] = false;
        delete hpost;
        delete asimov;
        delete hpd;
        delete leg;
        continue;
      }
      */

      // Store that this parameter is indeed being varied
      vary[i] = true;

      // Write to file
      c0->SetName(hpost->GetName());
      c0->SetTitle(hpost->GetTitle());
      c0->Print(canvasname);

      // cd into params directory in root file
      file->cd();
      params->cd();
      c0->Write();

      // Loop over the other parameters to get the correlations
      for (int j = 0; j <= i; j++) {

      // Skip the diagonal elements which we've already done above
      if (j == i) continue;
	
        // If thie parameter isn't varied
        if (vary[j] == false) continue;

        std::string tempString = std::string(bnames[j]);
        // For xsec parameters get the parameter name
        if (j<ndraw-nflux) {
            tempString = GetXsecName(j);

            double central, prior, down, up;
            GetXsecLimits(j, central, prior, down, up);
            (*xsec_nom)(j) = central;
            asimovLine = prior;
        }
        else //now flux parameters
        {
            TString TempText= Form("b_%i", i-nxsec);
            tempString = TempText.Data();  
        }
        file->cd();

        int nbins = 70;

        TString drawcmd = bnames[j]+":"+bnames[i];
        std::cout << drawcmd << std::endl;
        double maximum2 = chain->GetMaximum(bnames[j]);
        double minimum2 = chain->GetMinimum(bnames[j]);

        // TH2F to hold the correlation 
        TH2F *hpost2 = new TH2F(drawcmd, drawcmd, hpost->GetNbinsX(), hpost->GetBinLowEdge(0), hpost->GetBinLowEdge(hpost->GetNbinsX()+1), nbins, minimum2, maximum2);
        hpost2->SetMinimum(0);
        hpost2->GetXaxis()->SetTitle(hpost->GetXaxis()->GetTitle());
        hpost2->GetYaxis()->SetTitle(tempString.c_str());
        hpost2->GetZaxis()->SetTitle("Steps");
        std::string plotTitle = hpost->GetXaxis()->GetTitle();
        plotTitle += " vs " + tempString;
        hpost2->SetTitle(plotTitle.c_str());

        // The draw command we want, i.e. draw param j vs param i
        chain->Project(drawcmd, drawcmd, stepcut.c_str());
        hpost2->Draw("colz");

        // Get the covariance for these two parameters
        (*covariance)(i,j) = hpost2->GetCovariance();
        (*covariance)(j,i) = (*covariance)(i,j);
        std::cout << std::setw(10) << "covariance:" << (*covariance)(i,j) << std::endl;
        (*correlation)(i,j) = hpost2->GetCorrelationFactor();
        (*correlation)(j,i) = (*correlation)(i,j);
        std::cout << std::setw(10) << "correlation:" << (*correlation)(i,j) << std::endl;

        c0->SetName(hpost2->GetName());
        c0->SetTitle(hpost2->GetTitle());
        //we only want to plot correlation between xsec parames
        if(j<nflux && i<nflux) c0->Print(canvasname);
        //if((j<25+nxsec && i<25+nxsec) || (j<25+nxsec && i>nflux) || (j<nflux && i<25+nxsec) || (j<nflux && i<nflux))
        

        // Write it to root file
        file->cd();
        corr->cd();
        hpost2->Write();

        delete hpost2;

      } // End for (j = 0; j <= i; ++j)

      delete hpost;
      delete asimov;
      delete hpd;
      delete leg;
    } 
    else {       // If we arent't interested in drawing the correlations
      file->cd();
      int nbins = 250;

      double asimovLine = 1.0;

      std::string tempString = std::string(bnames[i]);
      if (i<ndraw-nflux) { //first go xsec then flux
      isXsec = true;

        tempString = GetXsecName(i);

        double central, prior, down, up;
        GetXsecLimits(i, central, prior, down, up);
        (*xsec_nom)(i) = central;
        asimovLine = prior;
      }
      else //now flux parameters
      {
        TString TempText= Form("b_%i", i-nxsec);
        tempString = TempText.Data();  
      }
      file->cd();

      double maxi = chain->GetMaximum(bnames[i]);
      double mini = chain->GetMinimum(bnames[i]);
      // This holds the posterior density
      TH1D *hpost = new TH1D(bnames[i], bnames[i], nbins, mini, maxi);
      hpost->SetMinimum(0);
      hpost->GetYaxis()->SetTitle("Steps");
      hpost->GetYaxis()->SetNoExponent(false);

      // Project bnames[i] onto hpost, applying stepcut
      chain->Project(bnames[i], bnames[i], stepcut.c_str());

      hpost->Smooth();

      // Get the characteristics of the hpost
      double mean, rms;
      GetArithmetic(hpost, mean, rms);
      double peakval, sigma_p, sigma_m, sigma_hpd;
      GetHPD(hpost, peakval, sigma_hpd, sigma_p, sigma_m);
      double gauss_mean, gauss_rms;
      GetGaussian(hpost, gauss, gauss_mean, gauss_rms);

      std::cout << mean << " +/- " << rms << " (" << peakval << "+/-" << sigma_hpd << " + " << sigma_p << " - " << sigma_m << ")" << " (" << gauss_mean << "+/-" << gauss_rms << ")" << std::endl;

      TLine *hpd = new TLine(peakval, hpost->GetMinimum(), peakval, hpost->GetMaximum());
      hpd->SetLineColor(kBlack);
      hpd->SetLineWidth(2);
      hpd->SetLineStyle(kSolid);

      TLegend *leg = new TLegend(0.23, 0.57, 0.73, 0.87);
      leg->SetTextSize(0.03);
      leg->AddEntry(hpost, Form("#splitline{PDF}{#mu = %.2f, #sigma = %.2f}", hpost->GetMean(), hpost->GetRMS()), "l");
      leg->AddEntry(gauss, Form("#splitline{Gauss}{#mu = %.2f, #sigma = %.2f}", gauss->GetParameter(1), gauss->GetParameter(2)), "l");
      leg->AddEntry(hpd, Form("#splitline{HPD}{#mu = %.2f, #sigma = %.2f (+%.2f-%.2f)}", peakval, sigma_hpd,sigma_p, sigma_m), "l");

      if (isXsec && ((*xsec_nom)(i)) != 0) {
        mean = mean / ((*xsec_nom)(i));
        rms = rms / ((*xsec_nom)(i));
        gauss_mean = gauss_mean / ((*xsec_nom)(i));
        gauss_rms = gauss_rms / ((*xsec_nom)(i));
        peakval = peakval / ((*xsec_nom)(i));
        sigma_hpd = sigma_hpd / ((*xsec_nom)(i));
        sigma_p = sigma_p / ((*xsec_nom)(i));
        sigma_m = sigma_m / ((*xsec_nom)(i));
      } else if (isXsec && ((*xsec_nom)(i)) == 0) {
        mean = mean + 1.0;
        gauss_mean = gauss_mean + 1.0;
        peakval = peakval + 1.0;
      }

      (*mean_vec)(i)      = mean;
      (*err_vec)(i)       = rms;
      (*gaus_mean_vec)(i) = gauss_mean;
      (*gaus_err_vec)(i)  = gauss_rms;
      (*HPD_mean_vec)(i) = peakval;
      (*HPD_err_p_vec)(i)  = sigma_p;
      (*HPD_err_vec)(i)   = sigma_hpd;
      (*HPD_err_m_vec)(i)  = sigma_m;
      (*covariance)(i,i)  = rms*rms;
      (*correlation)(i,i)  = 1.0;

      hpost->SetLineWidth(2);
      hpost->SetMaximum(hpost->GetMaximum()*1.5);
      hpost->SetTitle(tempString.c_str());
      hpost->GetXaxis()->SetTitle((std::string(hpost->GetTitle())+" rel. nom").c_str());

      // Now make the TLine for the asimov
      TLine *asimov = new TLine(asimovLine, hpost->GetMinimum(), asimovLine, hpost->GetMaximum());
      asimov->SetLineColor(kRed-3);
      asimov->SetLineWidth(2);
      asimov->SetLineStyle(kDashed);
      hpost->Draw();
      hpd->Draw("same");
      asimov->Draw("same");

      // Make the legend
      leg->AddEntry(asimov, Form("#splitline{Input}{x = %.2f}", asimovLine), "l");
      leg->SetLineColor(0);
      leg->SetLineStyle(0);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->Draw("same");

      if (hpost->GetMaximum() == hpost->Integral()*1.5) {
        vary[i] = false;
        delete hpost;
        delete asimov;
        delete hpd;
        delete leg;
        continue;
      }
      vary[i] = true;

      // Write to file
      c0->SetName(hpost->GetName());
      c0->SetTitle(hpost->GetTitle());
      c0->Print(canvasname);

      // cd into params directory in root file
      file->cd();
      params->cd();
      c0->Write();

      delete hpost;
      delete asimov;
      delete hpd;
      delete leg;

    } // End the if(drawCorr) else

  } // End for ndraw


  // Make the prefit plot
  TH1D* prefit = MakePrefit(ndraw ,nflux); 
  
  // cd into the output file
  file->cd();

  // Make a TH1D of the central values and the errors
  TH1D *paramPlot = new TH1D("paramPlot", "paramPlot", ndraw, 0, ndraw);
  paramPlot->SetName("mach3params");
  paramPlot->SetTitle(stepcut.c_str());
  paramPlot->SetFillStyle(3001);
  paramPlot->SetFillColor(kBlue-1);
  paramPlot->SetMarkerColor(paramPlot->GetFillColor());
  paramPlot->SetMarkerStyle(20);
  paramPlot->SetLineColor(paramPlot->GetFillColor());
  paramPlot->SetMarkerSize(prefit->GetMarkerSize());

  // Same but with Gaussian output
  TH1D *paramPlot_gauss = (TH1D*)(paramPlot->Clone());
  paramPlot_gauss->SetMarkerColor(kOrange-5);
  paramPlot_gauss->SetMarkerStyle(23);
  paramPlot_gauss->SetLineWidth(2);
  paramPlot_gauss->SetMarkerSize((prefit->GetMarkerSize())*0.75);
  paramPlot_gauss->SetFillColor(paramPlot_gauss->GetMarkerColor());
  paramPlot_gauss->SetFillStyle(3244);
  paramPlot_gauss->SetLineColor(paramPlot_gauss->GetMarkerColor());

  // Same but with Gaussian output
  TH1D *paramPlot_HPD = (TH1D*)(paramPlot->Clone());
  paramPlot_HPD->SetMarkerColor(kBlack);
  paramPlot_HPD->SetMarkerStyle(25);
  paramPlot_HPD->SetLineWidth(2);
  paramPlot_HPD->SetMarkerSize((prefit->GetMarkerSize())*0.5);
  paramPlot_HPD->SetFillColor(0);
  paramPlot_HPD->SetFillStyle(0);
  paramPlot_HPD->SetLineColor(paramPlot_HPD->GetMarkerColor());

  // Set labels and data
  for (int i = 0; i < ndraw; ++i) {

    paramPlot->SetBinContent(i+1, (*mean_vec)(i));
    paramPlot->SetBinError(i+1, (*err_vec)(i));

    paramPlot_gauss->SetBinContent(i+1, (*gaus_mean_vec)(i));
    paramPlot_gauss->SetBinError(i+1, (*gaus_err_vec)(i));

    paramPlot_HPD->SetBinContent(i+1, (*HPD_mean_vec)(i));
    double error = (*HPD_err_vec)(i);
    paramPlot_HPD->SetBinError(i+1, error);

    // Don't care about labelling beam for now
    if (i < nflux) {
      paramPlot->GetXaxis()->SetBinLabel(i+1, GetXsecName(i).c_str());
      paramPlot_gauss->GetXaxis()->SetBinLabel(i+1, (GetXsecName(i )+"_G").c_str());
      paramPlot_HPD->GetXaxis()->SetBinLabel(i+1, (GetXsecName(i)+"_G").c_str());
    }
  }


  // Make a TLegend
  TLegend *CompLeg = new TLegend(0.33, 0.20, 0.80, 0.50);
  CompLeg->AddEntry(prefit, "Prefit", "fp");
  CompLeg->AddEntry(paramPlot, "Postfit PDF", "fp");
  CompLeg->AddEntry(paramPlot_gauss, "Postfit Gauss", "fp");
  CompLeg->AddEntry(paramPlot_HPD, "Postfit HPD", "lfep");
  CompLeg->SetFillColor(0);
  CompLeg->SetFillStyle(0);
  CompLeg->SetLineWidth(0);
  CompLeg->SetLineStyle(0);
  CompLeg->SetBorderSize(0);

  file->cd();
  c0->SetBottomMargin(0.2);
  
  // Plot the flux parameters (0 to 100) if enabled
  // Have already looked through the branches earlier
  if (plotFlux == true) {

    prefit->GetYaxis()->SetRangeUser(0.6, 1.3);
    prefit->GetXaxis()->SetTitle("");
    prefit->GetXaxis()->LabelsOption("v");
    paramPlot->GetYaxis()->SetRangeUser(0.6, 1.3);
    paramPlot->GetXaxis()->SetTitle("");
    paramPlot->SetTitle(stepcut.c_str());
    paramPlot->GetXaxis()->LabelsOption("v");
    paramPlot_gauss->GetYaxis()->SetRangeUser(0.6, 1.3);
    paramPlot_gauss->GetXaxis()->SetTitle("");
    paramPlot_gauss->SetTitle(stepcut.c_str());
    paramPlot_gauss->GetXaxis()->LabelsOption("v");
    paramPlot_HPD->GetYaxis()->SetRangeUser(0.6, 1.3);
    paramPlot_HPD->GetXaxis()->SetTitle("");
    paramPlot_HPD->SetTitle(stepcut.c_str());
    paramPlot_HPD->GetXaxis()->LabelsOption("v");

    int nplots = 4;
    for (int i = 0; i < nplots; ++i) {
      int bot = (double(i)/double(nplots))*nflux+nxsec;
      int top = (double(i+1)/double(nplots))*nflux+nxsec;

      prefit->GetXaxis()->SetRangeUser(bot, top);
      paramPlot->GetXaxis()->SetRangeUser(bot, top);
      paramPlot_gauss->GetXaxis()->SetRangeUser(bot, top);
      paramPlot_HPD->GetXaxis()->SetRangeUser(bot, top);

      prefit->Write(Form("param_flux_prefit_%i", i));
      paramPlot->Write(Form("param_flux_%i", i));
      paramPlot_gauss->Write(Form("param_flux_gaus_%i", i));
      paramPlot_HPD->Write(Form("param_flux_HPD_%i", i));

      prefit->Draw("e2");
      paramPlot->Draw("e2, same");
      paramPlot_gauss->Draw("e2, same");
      paramPlot_HPD->Draw("e1, same");

      CompLeg->Draw("same");
      c0->Write("param_flux_canv");
      c0->Print(canvasname);
      c0->Clear();
    }
  }

  // Plot the xsec parameters (100 to ~125)
  // Have already looked through the branches earlier
  file->cd();
  if (plotXsec == true) {
    prefit->GetYaxis()->SetRangeUser(-0.5 , 2.5);
    prefit->GetXaxis()->SetTitle("");
    prefit->GetXaxis()->SetRangeUser(0, nxsec);
    prefit->GetXaxis()->LabelsOption("v");

    paramPlot->GetXaxis()->SetRangeUser(0, nxsec);
    paramPlot_gauss->GetXaxis()->SetRangeUser(0, nxsec);
    paramPlot_HPD->GetXaxis()->SetRangeUser(0, nxsec);

    prefit->Write("param_xsec_prefit");
    paramPlot->Write("param_xsec");
    paramPlot_gauss->Write("param_xsec_gaus");
    paramPlot_HPD->Write("param_xsec_HPD");

    prefit->Draw("e2");
    paramPlot->Draw("e2, same");
    paramPlot_gauss->Draw("e2, same");
    paramPlot_HPD->Draw("same");

    CompLeg->SetX1NDC(0.33);
    CompLeg->SetX2NDC(0.80);
    CompLeg->SetY1NDC(0.20);
    CompLeg->SetY2NDC(0.50);

    CompLeg->Draw("same");
    c0->Write("param_xsec_canv");
    c0->Print(canvasname);
    c0->Clear();
  }

  delete CompLeg;

  c0->SetLeftMargin(0.1);
  c0->SetBottomMargin(0.1);

  TH2D* hCov, *hCovSq, *hCorr = NULL;
  if (drawCorr) {
    // The covariance matrix from the fit
    hCov = new TH2D("hCov", "hCov", covBinning, 0, covBinning, covBinning, 0, covBinning);
    hCov->GetZaxis()->SetTitle("Covariance");
    // The covariance matrix square root, with correct sign
    hCovSq = new TH2D("hCovSq", "hCovSq", covBinning, 0, covBinning, covBinning, 0, covBinning);
    // The correlation
    hCorr = new TH2D("hCorr", "hCorr", covBinning, 0, covBinning, covBinning, 0, covBinning);
    hCorr->GetZaxis()->SetTitle("Correlation");
    hCorr->SetMinimum(-1);
    hCorr->SetMaximum(1);
    hCov->GetXaxis()->SetLabelSize(0.015);
    hCov->GetYaxis()->SetLabelSize(0.015);
    hCovSq->GetXaxis()->SetLabelSize(0.015);
    hCovSq->GetYaxis()->SetLabelSize(0.015);
    hCorr->GetXaxis()->SetLabelSize(0.015);
    hCorr->GetYaxis()->SetLabelSize(0.015);

    // Loop over the covariance matrix entries
    for(int i=0; i < ndraw; i++) {

      if (drawCorr && vary[i] == false) {
        continue;
      }

      if (i < nxsec) {
        hCov->GetXaxis()->SetBinLabel(i+1, GetXsecName(i).c_str());
        hCovSq->GetXaxis()->SetBinLabel(i+1, GetXsecName(i).c_str());
        hCorr->GetXaxis()->SetBinLabel(i+1, GetXsecName(i).c_str());
      }

      for (int j=0; j<ndraw; j++) {

        if (drawCorr && vary[j] == false) {
          continue;
        }

        if (j < nxsec) {
          hCov->GetYaxis()->SetBinLabel(j+1, GetXsecName(j).c_str());
          hCovSq->GetYaxis()->SetBinLabel(j+1, GetXsecName(j).c_str());
          hCorr->GetYaxis()->SetBinLabel(j+1, GetXsecName(j).c_str());
        }

        // The value of the covariance
        double cov = (*covariance)(i,j);
        double corr = (*correlation)(i,j);

        hCov->SetBinContent(i+1, j+1, cov);
        hCovSq->SetBinContent(i+1, j+1, ((cov>0)-(cov<0))*sqrt(fabs(cov)));
        hCorr->SetBinContent(i+1, j+1, corr);
      }

    }

    // Take away the stat box
    gStyle->SetOptStat(0);
    // Make pretty correlation colors (red to blue)
    const int NRGBs = 5;
    TColor::InitializeColors();
    Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
    Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
    Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
    TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
    gStyle->SetNumberContours(255);

    c0->cd();
    c0->Clear();
    hCov->Draw("colz");
    c0->SetRightMargin(0.15);
    c0->Print(canvasname);

    c0->cd();
    c0->Clear();
    hCorr->Draw("colz");
    c0->SetRightMargin(0.15);
    c0->Print(canvasname);

    file->cd();
    hCov->Write("postfit_cov_plot");
    hCovSq->Write("postfit_covsq_plot");
    hCorr->Write("postfit_corr_plot");
  }

  // Then close the pdf file
  std::cout << "Closing pdf " << canvasname << std::endl;
  canvasname+="]";
  c0->Print(canvasname);

  // Write all the nice vectors
  file->cd();
  mean_vec->Write("postfit_params_arit");
  err_vec->Write("postfit_errors_arit");
  gaus_mean_vec->Write("postfit_params_gauss");
  gaus_err_vec->Write("postfit_errors_gauss");
  HPD_mean_vec->Write("postfit_params_HPD");
  HPD_err_vec->Write("postfit_errors_HPD");
  HPD_err_p_vec->Write("postfit_errors_HPD_pos");
  HPD_err_m_vec->Write("postfit_errors_HPD_neg");
  covariance->Write("postfit_cov");
  correlation->Write("postfit_corr");

  file->Close();
}

void DrawComp(std::string inputFile1, std::string title1, std::string inputFile2, std::string title2) {

  std::cout << "File for study:       " << inputFile1 << std::endl;
  std::cout << "File for nominal:     " << inputFile2 << std::endl;

  // Open the chain for file 1
  TChain* chain = new TChain("posteriors","");
  chain->Add(inputFile1.c_str());
  // Open the chain for file 2
  TChain* chain2 = new TChain("posteriors","");
  chain2->Add(inputFile2.c_str());

  int nEntries = chain->GetMaximum("step");
  int nEntries2 = chain2->GetMaximum("step");
  // Use first 1/7 as burn-in
  int cut = nEntries/BurnInCut;
  int cut2 = nEntries2/BurnInCut;
  std::stringstream ss;
  ss << "step > " << cut;
  std::stringstream ss2;
  ss2 << "step > " << cut2;

  // What entries are we interested in?
  // Here specifically 20k steps in, and events that have a logL < 5000
  std::string stepcut = ss.str();
  std::string stepcut2 = ss2.str();

  // Loop over interesting params; only need to do this for one of the files granted we have same params
  // Get the list of branches
  TObjArray* brlis = (TObjArray*)chain->GetListOfBranches();
  // Get the number of branches
  int nbr = brlis->GetEntries();
  std::cout << "# of branches: " << nbr << std::endl;
  // Make an array of TStrings
  TString bnames[nbr];
  // Array of doubles containing the nominal values (first step for xsec)
  std::vector<double> nom;
  nom.resize(nbr);

  // Have a counter for how many of each systematic we have
  int ndraw = 0;
  int nflux = 0;
  int nxsec = 0;

  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  chain->SetBranchStatus("*", false);
  chain->SetBranchStatus("step", true);
  for(int i=0; i < nbr; i++) {

    // Get the TBranch and its name
    TBranch* br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();

    // Get first entry of this branch (i.e. systematic); this is the first reconfigure so should be the fake-data which it was generated for (for xsec this is true...!)
    if(bname.BeginsWith("LogL") || bname.BeginsWith("ndd_")) continue;

    // If we're on beam systematics
    if(bname.BeginsWith("xsec_")) {
      //plotXsec = true;
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      chain->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      ndraw++;
      //nxsec++;

    }
    /*
    else if(bname.BeginsWith("b_")) {

      plotFlux = true;
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      nflux++;
      ndraw++;

    }
    */
  }
  TFile *TempFile = new TFile(inputFile1.c_str(), "open");
  TTree* Settings = (TTree*)(TempFile->Get("Settings"));
  // Get the xsec covariance matrix
  std::string *XsecInput = 0;
  Settings->SetBranchAddress("XsecCov", &XsecInput);
  Settings->GetEntry(0);
  delete Settings;
  TempFile->Close();
  delete TempFile;

  XsecCov = "../"+*XsecInput;

  chain->GetEntry(0);
  nom.resize(nbr);

  for(int i = 0; i < ndraw; ++i) {
     std::string tempString = GetXsecName(i);
      if (tempString.rfind("b_", 0) == 0) {
        nflux++;
      }
      else nxsec++;
  }
  
  std::cout << "# useful entries (flux, xsec, fsi): " << ndraw << std::endl;
  std::cout << "# useful entries (flux): " << nflux << std::endl;
  std::cout << "# useful entries (xsec): " << nxsec << std::endl;
  std::cout << "************************************************" << std::endl;

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  // Open a TCanvas to write the posterior onto
  TCanvas* c0 = new TCanvas("c0", "c0", 0, 0, 1024, 1024);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c0->SetGrid();
  c0->SetBottomMargin(0.1);
  c0->SetTopMargin(0.05);
  c0->SetRightMargin(0.03);
  c0->SetLeftMargin(0.10);

  // Write to a PDF file
  // Strip out the initial ../ if can find
  // Strip out the last .root
  inputFile1 = inputFile1.substr(0, inputFile1.find(".root")-1);
  while (inputFile2.find("/") != std::string::npos) {
    inputFile2 = inputFile2.substr(inputFile2.find("/")+1, inputFile2.find(".root")-1);
  }

  TString canvasname = inputFile1+"_"+inputFile2+".pdf[";
  c0->Print(canvasname);
  // Once the pdf file is open no longer need to bracket
  canvasname.ReplaceAll("[","");

  // We fit with this gaussian
  TF1 *gauss = new TF1("gauss","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",-5,5);
  gauss->SetLineColor(kBlue-1);
  gauss->SetLineStyle(kDashed);

  TF1 *gauss2 = new TF1("gauss2","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",-5,5);
  gauss2->SetLineColor(kRed);
  gauss2->SetLineStyle(kDashed);

  TVectorD* mean_vec = new TVectorD(ndraw);
  TVectorD* err_vec = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec = new TVectorD(ndraw);
  TVectorD* gaus_err_vec = new TVectorD(ndraw); 
  TVectorD* HPD_mean_vec  = new TVectorD(ndraw);
  TVectorD* HPD_err_vec  = new TVectorD(ndraw);
  TVectorD* HPD_err_p_vec   = new TVectorD(ndraw); 
  TVectorD* HPD_err_m_vec   = new TVectorD(ndraw); 

  TVectorD* mean_vec2 = new TVectorD(ndraw);
  TVectorD* err_vec2 = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec2 = new TVectorD(ndraw);
  TVectorD* gaus_err_vec2 = new TVectorD(ndraw); 
  TVectorD* HPD_mean_vec2  = new TVectorD(ndraw);
  TVectorD* HPD_err_vec2  = new TVectorD(ndraw);
  TVectorD* HPD_err_p_vec2   = new TVectorD(ndraw); 
  TVectorD* HPD_err_m_vec2   = new TVectorD(ndraw); 

  TVectorD* xsec_nom  = new TVectorD(nxsec);

  // Output file to write to
  TString rootfilename = inputFile1+"_"+inputFile2+".root";
  TFile* file = new TFile(rootfilename, "RECREATE");
  TDirectory *params = file->mkdir("hpost");

  int nbins = 70;
  // ndraw is number of draws we want to do
  for(int i=0; i<ndraw; ++i) {

    bool isXsec = false;
    file->cd();

    double asimovLine = 1.0;

    std::string tempString;
    if (i<ndraw-nflux) { //first go xsec then flux
      isXsec = true;

        tempString = GetXsecName(i);

        double central, prior, down, up;
        GetXsecLimits(i, central, prior, down, up);
        (*xsec_nom)(i) = central;
        asimovLine = prior;
      }
      else //now flux parameters
      {
        TString TempText= Form("b_%i", i-nxsec);
        tempString = TempText.Data();  
      }

    file->cd();
    std::cout << bnames[i] << std::endl;

    double maxi = chain->GetMaximum(bnames[i]);
    double mini = chain->GetMinimum(bnames[i]);
    double maxi2 = chain2->GetMaximum(bnames[i]);
    double mini2 = chain2->GetMinimum(bnames[i]);

    // This holds the posterior density
    TH1D *hpost = new TH1D(bnames[i], bnames[i], nbins, mini, maxi);
    hpost->SetMinimum(0);
    hpost->SetTitle(bnames[i]);
    TH1D *hpost2 = new TH1D(bnames[i]+"_2", bnames[i]+"_2", nbins, mini2, maxi2);
    hpost2->SetMinimum(0);
    hpost2->SetTitle(bnames[i]+"_2");
    hpost->GetYaxis()->SetTitle("Steps");
    hpost->GetYaxis()->SetTitleOffset(1.1);
    hpost->GetYaxis()->SetNoExponent(false);

    // Project bnames[i] onto hpost, applying stepcut
    chain->Project(bnames[i], bnames[i], stepcut.c_str());
    chain2->Project(bnames[i]+"_2", bnames[i], stepcut2.c_str());

    // Area normalise the distributions
    hpost->Scale(1./hpost->Integral(), "width");
    hpost2->Scale(1./hpost2->Integral(), "width");

    hpost->Smooth();
    hpost2->Smooth();

    // Get the characteristics of the hpost
    double mean, rms;
    GetArithmetic(hpost, mean, rms);
    double peakval, sigma_p, sigma_m, sigma_hpd;
    GetHPD(hpost, peakval, sigma_hpd, sigma_p, sigma_m);
    double gauss_mean, gauss_rms;
    GetGaussian(hpost, gauss, gauss_mean, gauss_rms);

    // Get the characteristics of the hpost
    double mean2, rms2;
    GetArithmetic(hpost2, mean2, rms2);
    double peakval2, sigma_p2, sigma_m2, sigma_hpd2;
    GetHPD(hpost2, peakval2, sigma_hpd2, sigma_p2, sigma_m2);
    double gauss_mean2, gauss_rms2;
    GetGaussian(hpost2, gauss2, gauss_mean2, gauss_rms2);

    std::cout << mean << " +/- " << rms << " (" << peakval << " + " << sigma_p << " - " << sigma_m << ")" << std::endl;
    std::cout << mean2 << " +/- " << rms2 << " (" << peakval2 << " + " << sigma_p2 << " - " << sigma_m2 << ")" << std::endl;

    // Set some nice colours
    hpost->SetLineColor(kBlue-1);
    hpost->SetLineWidth(2);

    hpost2->SetLineColor(kRed);
    hpost2->SetLineWidth(2);

    // Now make the TLine for the asimov
    TLine *asimov = new TLine(asimovLine, hpost->GetMinimum(), asimovLine, hpost->GetMaximum());
    asimov->SetLineColor(kBlack);
    asimov->SetLineWidth(2);
    asimov->SetLineStyle(kDashed);

    // Make a nice little TLegend
    TLegend *leg = new TLegend(0.12, 0.7, 0.6, 0.97);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    TString rebinLeg = Form("#splitline{%s}{#mu = %.2f, #sigma = %.2f}", title1.c_str(), mean, rms);
    TString nombinLeg = Form("#splitline{%s}{#mu = %.2f, #sigma = %.2f}", title2.c_str(), mean2, rms2);
    TString asimovLeg = Form("Central, #mu = %.2f", asimovLine);
    leg->AddEntry(asimov, asimovLeg, "l");
    leg->AddEntry(hpost,  rebinLeg, "l");
    leg->AddEntry(hpost2, nombinLeg, "l");

    TLine *hpd = new TLine(peakval, hpost->GetMinimum(), peakval, hpost->GetMaximum());
    hpd->SetLineColor(hpost->GetLineColor());
    hpd->SetLineWidth(2);
    hpd->SetLineStyle(kSolid);

    TLine *hpd2 = new TLine(peakval2, hpost2->GetMinimum(), peakval2, hpost2->GetMaximum());
    hpd2->SetLineColor(hpost2->GetLineColor());
    hpd2->SetLineWidth(2);
    hpd2->SetLineStyle(kSolid);

    if (isXsec && ((*xsec_nom)(i)) != 0) {
      mean = mean / ((*xsec_nom)(i));
      mean2 = mean2 / ((*xsec_nom)(i));
      rms = rms / ((*xsec_nom)(i));
      rms2 = rms2 / ((*xsec_nom)(i));
      gauss_mean = gauss_mean / ((*xsec_nom)(i));
      gauss_mean2 = gauss_mean2 / ((*xsec_nom)(i));
      gauss_rms = gauss_rms / ((*xsec_nom)(i));
      gauss_rms2 = gauss_rms2 / ((*xsec_nom)(i));
      peakval = peakval / ((*xsec_nom)(i));
      sigma_hpd = sigma_hpd / ((*xsec_nom)(i));
      sigma_p = sigma_p / ((*xsec_nom)(i));
      sigma_m = sigma_m / ((*xsec_nom)(i));
      peakval2 = peakval2 / ((*xsec_nom)(i));
      sigma_hpd2 = sigma_hpd2 / ((*xsec_nom)(i));
      sigma_p2 = sigma_p2 / ((*xsec_nom)(i));
      sigma_m2 = sigma_m2 / ((*xsec_nom)(i));
    } else if (isXsec && ((*xsec_nom)(i)) == 0) {
      mean = mean + 1.0;
      mean2 = mean2 + 1.0;
      gauss_mean = gauss_mean + 1.0;
      gauss_mean2 = gauss_mean2 + 1.0;
      peakval = peakval + 1.0;
      peakval2 = peakval2 + 1.0;
    }

    (*mean_vec)(i)      = mean;
    (*mean_vec2)(i)     = mean2;
    (*err_vec)(i)       = rms;
    (*err_vec2)(i)      = rms2;
    (*HPD_mean_vec)(i)  = peakval;
    (*HPD_err_vec)(i)  = sigma_hpd;
    (*HPD_err_p_vec)(i)  = sigma_p;
    (*HPD_err_m_vec)(i)  = sigma_m;

    (*gaus_mean_vec)(i) = gauss_mean;
    (*gaus_mean_vec2)(i)= gauss_rms;
    (*gaus_err_vec)(i)  = gauss_mean2;
    (*gaus_err_vec2)(i) = gauss_rms2;
    (*HPD_mean_vec2)(i)  = peakval2;
    (*HPD_err_vec2)(i)  = sigma_hpd2;
    (*HPD_err_p_vec2)(i)  = sigma_p2;
    (*HPD_err_m_vec2)(i)  = sigma_m2;

    gauss->SetLineColor(kBlue-1);
    gauss->SetLineStyle(kDashed);
    gauss2->SetLineColor(kRed);
    gauss2->SetLineStyle(kDashed);

    // Set the maximum properly
    double max = std::max(hpost->GetMaximum(), hpost2->GetMaximum());
    hpost->SetMaximum(1.5*max);
    hpost2->SetMaximum(hpost->GetMaximum());
    hpost->SetTitle(tempString.c_str());
    hpost->GetXaxis()->SetTitle((std::string(hpost->GetTitle())+" rel. nom").c_str());

    hpost->Draw();
    hpost2->Draw("same");
    asimov->Draw("same");
    hpd->Draw("same");
    hpd2->Draw("same");
    leg->Draw("same");
    c0->Print(canvasname);

    // cd into params directory in root file
    params->cd();
    hpost->Write();
    hpost2->Write();

    delete hpost;
    delete hpost2;
    delete asimov;
    delete leg;
    delete hpd;
    delete hpd2;

  } // End the for ndraw loop

  file->cd();

  // vector of the mean values
  std::vector<double> val;
  std::vector<double> val2;
  // vector of the error values
  std::vector<double> vale;
  std::vector<double> vale2;
  // vector of the index
  std::vector<double> ind;
  std::vector<double> ind2;
  // No idea...
  std::vector<double> inde;
  std::vector<double> inde2;

  for (int i = 0; i < mean_vec->GetNrows(); i++) {
    val.push_back((*mean_vec)(i));
    val2.push_back((*mean_vec2)(i));
    vale.push_back((*err_vec)(i));
    vale2.push_back((*err_vec2)(i));
    ind.push_back(i);
    ind2.push_back(i);
    inde.push_back(0.0);
    inde2.push_back(0.0);
  }

  TH1D *prefit = MakePrefit(ndraw ,nflux);
  file->cd();

  // Make a TH1D of the central values and the errors
  TH1D *paramPlot = new TH1D("paramPlot", "paramPlot", ndraw, 0, ndraw);
  paramPlot->SetName("mach3params");
  paramPlot->SetTitle((stepcut+"_"+stepcut2).c_str());
  paramPlot->SetFillStyle(3144);
  paramPlot->SetFillColor(kYellow-3);
  paramPlot->SetMarkerColor(kMagenta-3);
  paramPlot->SetMarkerStyle(22);
  paramPlot->SetLineColor(paramPlot->GetFillColor());
  paramPlot->SetMarkerSize(prefit->GetMarkerSize());

  // Make a TH1D of the central values and the errors
  TH1D *paramPlot2 = (TH1D*)paramPlot->Clone();
  paramPlot2->SetName("mach3params_2");
  paramPlot2->SetFillStyle(3344);
  paramPlot2->SetFillColor(kBlue-1);
  paramPlot2->SetMarkerColor(paramPlot2->GetFillColor());
  paramPlot2->SetMarkerStyle(33);
  paramPlot2->SetLineColor(paramPlot2->GetFillColor());
  paramPlot2->SetMarkerSize(prefit->GetMarkerSize());

  TLegend *graphLeg = new TLegend(0.33, 0.10, 0.80, 0.40);
  graphLeg->AddEntry(prefit, "Prefit", "fp");
  graphLeg->AddEntry(paramPlot,   title1.c_str(), "fp");
  graphLeg->AddEntry(paramPlot2,  title2.c_str(), "fp");
  graphLeg->SetFillColor(0);
  graphLeg->SetFillStyle(0);
  graphLeg->SetLineColor(0);
  graphLeg->SetLineWidth(0);
  graphLeg->SetLineStyle(0);
  graphLeg->SetBorderSize(0);

  // Set labels and data
  for (int i = 0; i < ndraw; ++i) {

    paramPlot->SetBinContent(i+1, (*mean_vec)(i));
    paramPlot->SetBinError(i+1, (*err_vec)(i));

    paramPlot2->SetBinContent(i+1, (*mean_vec2)(i));
    paramPlot2->SetBinError(i+1, (*err_vec2)(i));

    // Don't care about labelling beam for now
    if (i < nflux) {
      paramPlot->GetXaxis()->SetBinLabel(i+1, GetXsecName(i).c_str());
      paramPlot2->GetXaxis()->SetBinLabel(i+1, GetXsecName(i).c_str());
    }
  }
  file->cd();

  // Plot the flux parameters (0 to 100) if enabled
  // Have already looked through the branches earlier
  if (plotFlux == true) {
    prefit->GetXaxis()->SetRangeUser(0, nxsec);
    prefit->GetXaxis()->SetTitle("");

    prefit->GetYaxis()->SetRangeUser(0.7, 1.3);
    prefit->GetYaxis()->SetTitle("Variation");

    prefit->SetTitle((stepcut+"_"+stepcut2).c_str());
    int nplots = 4;
    for (int i = 0; i < nplots; ++i) {
      int bot = (double(i)/double(nplots))*nflux + nxsec;
      int top = (double(i+1)/double(nplots))*nflux + nxsec;

      prefit->GetXaxis()->SetRangeUser(bot, top);
      paramPlot->GetXaxis()->SetRangeUser(bot, top);
      paramPlot2->GetXaxis()->SetRangeUser(bot, top);

      prefit->Write(Form("param_flux_prefit_%i", i));
      paramPlot->Write(Form("param_flux_%i", i));
      paramPlot2->Write(Form("param_flux_HPD_%i", i));

      prefit->Draw("e2");
      paramPlot->Draw("e2, same");
      paramPlot2->Draw("e2, same");

      graphLeg->Draw("same");
      c0->Write("param_flux_canv");
      c0->Print(canvasname);
      c0->Clear();
    }
  }

  file->cd();
  // Plot the xsec parameters (100 to ~125)
  // Have already looked through the branches earlier
  if (plotXsec == true) {
    prefit->GetXaxis()->SetTitle("");
    prefit->GetYaxis()->SetRangeUser(-0.5, 2.0);
    prefit->GetXaxis()->LabelsOption("v");

    prefit->GetXaxis()->SetRangeUser(0, nxsec);
    paramPlot->GetXaxis()->SetRangeUser(0, nxsec);
    paramPlot2->GetXaxis()->SetRangeUser(0, nxsec);

    prefit->Draw("e2");
    paramPlot->Draw("e2,same");
    paramPlot2->Draw("e2,same");
    graphLeg->Draw("same");
    c0->Print(canvasname);
    c0->Write("param_xsec_canv");
    c0->Clear();

    paramPlot->Write("param_xsec");
    paramPlot2->Write("param_xsec_2");
  }

  // Finally draw the parameter plot onto the PDF
  // Close the .pdf file with all the posteriors
  c0->cd();
  c0->Clear();


  file->cd();
  // Close the pdf file
  std::cout << "Closing pdf " << canvasname << std::endl;
  canvasname+="]";
  c0->Print(canvasname);

  // Write all the nice vectors
  mean_vec->Write("FitParameters");
  mean_vec2->Write("FitParameters_2");

  err_vec->Write("FitErrors");
  err_vec2->Write("FitErrors_2");

  gaus_mean_vec->Write("FitParameters_gauss");
  gaus_mean_vec2->Write("FitParameters_gauss_2");

  gaus_err_vec->Write("FitErrors_gauss");
  gaus_err_vec2->Write("FitErrors_gauss_2");

  HPD_mean_vec->Write("FitParameters_HPD");
  HPD_mean_vec2->Write("FitParameters_HPD_2");

  HPD_err_vec->Write("FitErrors_HPD");
  HPD_err_vec2->Write("FitErrors_HPD_2");

  HPD_err_p_vec->Write("FitErrors_HPD_p");
  HPD_err_p_vec2->Write("FitErrors_HPD_p_2");

  HPD_err_m_vec->Write("FitErrors_HPD_m");
  HPD_err_m_vec2->Write("FitErrors_HPD_m_2");

  paramPlot->Write("postfit_params_plot");
  paramPlot2->Write("postfit_params_plot_2");

  file->Close();
}

void DrawComp(std::string inputFile1, std::string title1, std::string inputFile2, std::string title2, std::string inputFile3, std::string title3) {

  std::cout << "File for study:       " << inputFile1 << std::endl;
  std::cout << "File for study2:      " << inputFile2 << std::endl;
  std::cout << "File for study3:      " << inputFile3 << std::endl;

  // Open the chain for file 1
  TChain* chain = new TChain("posteriors","");
  chain->Add(inputFile1.c_str());
  // Open the chain for file 2
  TChain* chain2 = new TChain("posteriors","");
  chain2->Add(inputFile2.c_str());
  // Open the chain for file 3
  TChain* chain3 = new TChain("posteriors","");
  chain3->Add(inputFile3.c_str());

  int nEntries = chain->GetMaximum("step");
  int nEntries2 = chain2->GetMaximum("step");
  int nEntries3 = chain3->GetMaximum("step");

  // Use first 1/4 as burn-in
  int cut = nEntries/BurnInCut;
  int cut2 = nEntries2/BurnInCut;
  int cut3 = nEntries3/BurnInCut;
  std::stringstream ss;
  ss << "step > " << cut;
  std::stringstream ss2;
  ss2 << "step > " << cut2;
  std::stringstream ss3;
  ss3 << "step > " << cut3;

  // What entries are we interested in?
  // Here specifically 20k steps in, and events that have a logL < 5000
  std::string stepcut = ss.str();
  std::string stepcut2 = ss2.str();
  std::string stepcut3 = ss3.str();

  // Loop over interesting params; only need to do this for one of the files granted we have same params
  //
  // Get the list of branches
  TObjArray* brlis = (TObjArray*)chain->GetListOfBranches();
  // Get the number of branches
  int nbr = brlis->GetEntries();
  std::cout << "# of branches: " << nbr << std::endl;
  // Make an array of TStrings
  TString bnames[nbr];
  // Array of doubles containing the nominal values (first step for xsec)
  std::vector<double> nom;
  nom.resize(nbr);

  // Have a counter for how many of each systematic we have
  int ndraw = 0;
  int nflux = 0;
  int nxsec = 0;

  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  chain->SetBranchStatus("*", false);
  chain->SetBranchStatus("step", true);
  chain2->SetBranchStatus("*", false);
  chain2->SetBranchStatus("step", true);
  chain3->SetBranchStatus("*", false);
  chain3->SetBranchStatus("step", true);
  for(int i=0; i < nbr; i++) {

    // Get the TBranch and its name
    TBranch* br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();

    // Get first entry of this branch (i.e. systematic); this is the first reconfigure so should be the fake-data which it was generated for (for xsec this is true...!)
    if(bname.BeginsWith("LogL") || bname.BeginsWith("ndd_")) continue;

    // If we're on beam systematics
    if(bname.BeginsWith("xsec_")) {

      plotXsec = true;
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      chain2->SetBranchStatus(bname, true);
      chain2->SetBranchAddress(bname, &(nom[i]));
      chain3->SetBranchStatus(bname, true);
      chain3->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      ndraw++;
      //nxsec++;

    }
    /*   
    else if(bname.BeginsWith("b_")) {

      plotFlux = true;
      chain->SetBranchStatus(bname, true);
      chain->SetBranchAddress(bname, &(nom[i]));
      chain2->SetBranchStatus(bname, true);
      chain2->SetBranchAddress(bname, &(nom[i]));
      chain3->SetBranchStatus(bname, true);
      chain3->SetBranchAddress(bname, &(nom[i]));
      bnames[ndraw]=bname;
      nflux++;
      ndraw++;
    }
    */
  }
  TFile *TempFile = new TFile(inputFile1.c_str(), "open");
  TTree* Settings = (TTree*)(TempFile->Get("Settings"));
  // Get the xsec covariance matrix
  std::string *XsecInput = 0;
  Settings->SetBranchAddress("XsecCov", &XsecInput);
  Settings->GetEntry(0);
  delete Settings;
  TempFile->Close();
  delete TempFile;

  XsecCov = "../"+*XsecInput;
  chain->GetEntry(0);
  nom.resize(nbr);

 for(int i = 0; i < ndraw; ++i) {
     std::string tempString = GetXsecName(i);
      if (tempString.rfind("b_", 0) == 0) {
        nflux++;
      }
      else nxsec++;
  }
  
  std::cout << "# useful entries (flux, xsec, fsi): " << ndraw << std::endl;
  std::cout << "# useful entries (flux): " << nflux << std::endl;
  std::cout << "# useful entries (xsec): " << nxsec << std::endl;
  std::cout << "************************************************" << std::endl;

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  // Open a TCanvas to write the posterior onto
  TCanvas* c0 = new TCanvas("c0", "c0", 0, 0, 1024, 1024);
  c0->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c0->SetBottomMargin(0.1);
  c0->SetTopMargin(0.02);
  c0->SetRightMargin(0.03);
  c0->SetLeftMargin(0.10);


  // Write to a PDF file
  // Strip out the initial ../ if can find
  // Then strip out any additional /
  // Strip out the last .root
  inputFile1 = inputFile1.substr(0, inputFile1.find(".root"));
  std::cout << inputFile1 << std::endl;
  while (inputFile2.find("/") != std::string::npos) {
    inputFile2 = inputFile2.substr(inputFile2.find("/")+1, inputFile2.find(".root")-1);
  }
  std::cout << inputFile2 << std::endl;
  while (inputFile3.find("/") != std::string::npos) {
    inputFile3 = inputFile3.substr(inputFile3.find("/")+1, inputFile3.find(".root")-1);
  }
  std::cout << inputFile3 << std::endl;
  TString canvasname = inputFile1+"_"+inputFile2+"_"+inputFile3+".pdf[";
  c0->Print(canvasname);
  // Once the pdf file is open no longer need to bracket
  canvasname.ReplaceAll("[","");

  // We fit with this gaussian
  TF1 *gauss = new TF1("gauss","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",   -5, 5);
  TF1 *gauss2 = new TF1("gauss2","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])", -5, 5);
  TF1 *gauss3 = new TF1("gauss3","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])", -5, 5);

  TVectorD* mean_vec = new TVectorD(ndraw);
  TVectorD* err_vec = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec = new TVectorD(ndraw);
  TVectorD* gaus_err_vec = new TVectorD(ndraw); 

  TVectorD* mean_vec2 = new TVectorD(ndraw);
  TVectorD* err_vec2 = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec2 = new TVectorD(ndraw);
  TVectorD* gaus_err_vec2 = new TVectorD(ndraw); 

  TVectorD* mean_vec3 = new TVectorD(ndraw);
  TVectorD* err_vec3 = new TVectorD(ndraw); 
  TVectorD* gaus_mean_vec3 = new TVectorD(ndraw);
  TVectorD* gaus_err_vec3 = new TVectorD(ndraw); 

  TVectorD* xsec_nom  = new TVectorD(nxsec);

  // Output file to write to
  TString rootfilename = inputFile1+"_"+inputFile2+"_"+inputFile3+".root";
  TFile* file = new TFile(rootfilename, "RECREATE");
  file->cd();
  TDirectory *params = file->mkdir("hpost");

  int nbins = 70;
  bool isXsec = false;

  // ndraw is number of draws we want to do
  for(int i=0; i<ndraw; ++i) {

    file->cd();
    isXsec = false;
    double asimovLine = 0.0;

    std::string tempString;
    if (i<ndraw-nflux) { //first go xsec then flux
      isXsec = true;

        tempString = GetXsecName(i);

        double central, prior, down, up;
        GetXsecLimits(i, central, prior, down, up);
        (*xsec_nom)(i) = central;
        asimovLine = prior;
      }
      else //now flux parameters
      {
        TString TempText= Form("b_%i", i-nxsec);
        tempString = TempText.Data();  
      }
    file->cd();
    std::cout << bnames[i] << std::endl;

    double maxi = chain->GetMaximum(bnames[i]);
    double mini = chain->GetMinimum(bnames[i]);
    double maxi2 = chain2->GetMaximum(bnames[i]);
    double mini2 = chain2->GetMinimum(bnames[i]);
    double maxi3 = chain3->GetMaximum(bnames[i]);
    double mini3 = chain3->GetMinimum(bnames[i]);

    // This holds the posterior density
    TH1D *hpost = new TH1D(bnames[i], bnames[i], nbins, mini, maxi);
    hpost->SetMinimum(0);
    hpost->GetYaxis()->SetTitle("Steps (area norm.)");
    hpost->GetYaxis()->SetTitleOffset(1.1);
    hpost->GetYaxis()->SetNoExponent(false);
    TH1D *hpost2 = new TH1D(bnames[i]+"_2", bnames[i]+"_2", nbins, mini2, maxi2);
    hpost2->SetMinimum(0);
    TH1D *hpost3 = new TH1D(bnames[i]+"_3", bnames[i]+"_3", nbins, mini3, maxi3);
    hpost3->SetMinimum(0);

    hpost->SetTitle(bnames[i]);
    hpost2->SetTitle(bnames[i]+"_2");
    hpost3->SetTitle(bnames[i]+"_3");

    // Project bnames[i] onto hpost, applying stepcut
    chain->Project(bnames[i], bnames[i], stepcut.c_str());
    chain2->Project(bnames[i]+"_2", bnames[i], stepcut2.c_str());
    chain3->Project(bnames[i]+"_3", bnames[i], stepcut3.c_str());

    // Area normalise the distributions
    hpost->Scale(1./hpost->Integral(), "width");
    hpost2->Scale(1./hpost2->Integral(), "width");
    hpost3->Scale(1./hpost3->Integral(), "width");

    hpost->Smooth();
    hpost2->Smooth();
    hpost3->Smooth();

    hpost->SetTitle(tempString.c_str());
    hpost->GetXaxis()->SetTitle(tempString.c_str());
    hpost2->SetTitle(tempString.c_str());
    hpost2->GetXaxis()->SetTitle(tempString.c_str());
    hpost3->SetTitle(tempString.c_str());
    hpost3->GetXaxis()->SetTitle(tempString.c_str());

    // Get the characteristics of the hpost
    double mean = hpost->GetMean();
    double rms = hpost->GetRMS();
    double peakval = hpost->GetBinCenter(hpost->GetMaximumBin());
    double gauss_mean = gauss->GetParameter(1);
    double gauss_rms = gauss->GetParameter(2);

    double mean2 = hpost2->GetMean();
    double rms2 = hpost2->GetRMS();
    double peakval2 = hpost2->GetBinCenter(hpost2->GetMaximumBin());
    double gauss_mean2 = gauss2->GetParameter(1);
    double gauss_rms2 = gauss2->GetParameter(2);

    double mean3 = hpost3->GetMean();
    double rms3 = hpost3->GetRMS();
    double peakval3 = hpost3->GetBinCenter(hpost3->GetMaximumBin());
    double gauss_mean3 = gauss3->GetParameter(1);
    double gauss_rms3 = gauss3->GetParameter(2);

    // Set the range for the Gaussian fit
    gauss->SetRange(mean - 1.5*rms , mean + 1.5*rms);
    gauss2->SetRange(mean2 - 1.5*rms2 , mean2 + 1.5*rms2);
    gauss3->SetRange(mean3 - 1.5*rms3 , mean3 + 1.5*rms3);

    // Set the starting parameters close to RMS and peaks of the histograms
    gauss->SetParameters(hpost->GetMaximum()*rms*sqrt(2*3.14), peakval, rms);
    gauss2->SetParameters(hpost2->GetMaximum()*rms2*sqrt(2*3.14), peakval2, rms2);
    gauss3->SetParameters(hpost3->GetMaximum()*rms3*sqrt(2*3.14), peakval3, rms3);

    // Set some nice colours
    hpost->SetLineColor(kBlue-1);
    hpost->SetLineWidth(2);
    gauss->SetLineColor(kBlue-1);
    gauss->SetLineStyle(kDashed);

    hpost2->SetLineColor(kRed);
    hpost2->SetLineWidth(2);
    gauss2->SetLineColor(kRed);
    gauss2->SetLineStyle(kDashed);

    hpost3->SetLineColor(kGreen+2);
    hpost3->SetLineWidth(2);
    gauss3->SetLineColor(kGreen+2);
    gauss3->SetLineStyle(kDashed);

    // Perform the fit
    hpost->Fit("gauss","Rq");
    hpost2->Fit("gauss2","Rq");
    hpost3->Fit("gauss3","Rq");

    // Now make the TLine for the asimov
    TLine *asimov = new TLine(asimovLine, hpost->GetMinimum(), asimovLine, hpost->GetMaximum());
    asimov->SetLineColor(kBlack);
    asimov->SetLineWidth(2);
    asimov->SetLineStyle(kDashed);

    // Make a nice little TLegend
    TLegend *leg = new TLegend(0.12, 0.7, 0.8, 0.97);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    TString nombinLeg = Form("#splitline{%s}{#mu = %.2f, #sigma = %.2f}", title1.c_str(), mean, rms);
    TString rebinLeg = Form("#splitline{%s}{#mu = %.2f, #sigma = %.2f}", title2.c_str(), mean2, rms2);
    TString rebinLeg2 = Form("#splitline{%s}{#mu = %.2f, #sigma = %.2f}", title3.c_str(), mean3, rms3);
    TString asimovLeg = Form("Central, #mu = %.2f", asimovLine);
    leg->AddEntry(asimov, asimovLeg, "l");
    leg->AddEntry(hpost, nombinLeg, "l");
    leg->AddEntry(hpost2, rebinLeg, "l");
    leg->AddEntry(hpost3, rebinLeg2, "l");

    TLine *hpd = new TLine(peakval, hpost->GetMinimum(), peakval, hpost->GetMaximum());
    hpd->SetLineColor(hpost->GetLineColor());
    hpd->SetLineWidth(2);
    hpd->SetLineStyle(kSolid);

    TLine *hpd2 = new TLine(peakval2, hpost2->GetMinimum(), peakval2, hpost2->GetMaximum());
    hpd2->SetLineColor(hpost2->GetLineColor());
    hpd2->SetLineWidth(2);
    hpd2->SetLineStyle(kSolid);

    TLine *hpd3 = new TLine(peakval3, hpost3->GetMinimum(), peakval3, hpost3->GetMaximum());
    hpd3->SetLineColor(hpost3->GetLineColor());
    hpd3->SetLineWidth(2);
    hpd3->SetLineStyle(kSolid);
    
    if (isXsec && ((*xsec_nom)(i)) != 0) {
      mean = mean / ((*xsec_nom)(i));
      rms = rms / ((*xsec_nom)(i));
      gauss_mean = gauss_mean / ((*xsec_nom)(i));
      gauss_rms = gauss_rms / ((*xsec_nom)(i));

      mean2 = mean2 / ((*xsec_nom)(i));
      rms2 = rms2 / ((*xsec_nom)(i));
      gauss_mean2 = gauss_mean2 / ((*xsec_nom)(i));
      gauss_rms2 = gauss_rms2 / ((*xsec_nom)(i));

      mean3 = mean3 / ((*xsec_nom)(i));
      rms3 = rms3 / ((*xsec_nom)(i));
      gauss_mean3 = gauss_mean3 / ((*xsec_nom)(i));
      gauss_rms3 = gauss_rms3 / ((*xsec_nom)(i));
    } else if (isXsec && ((*xsec_nom)(i)) == 0) {
      mean = mean + 1.0;
      gauss_mean = gauss_mean + 1.0;
      mean2 = mean2 + 1.0;
      gauss_mean2 = gauss_mean2 + 1.0;
      mean3 = mean3 + 1.0;
      gauss_mean3 = gauss_mean3 + 1.0;
    }

    // Save the mean and RMS for histo and Gaussian fit
    (*mean_vec)(i)      = mean;
    (*err_vec)(i)       = rms;
    (*gaus_mean_vec)(i) = gauss_mean;
    (*gaus_err_vec)(i)  = gauss_rms;

    (*mean_vec2)(i)     = mean2;
    (*err_vec2)(i)      = rms2;
    (*gaus_mean_vec2)(i)= gauss_mean2;
    (*gaus_err_vec2)(i) = gauss_rms2;

    (*mean_vec3)(i)     = mean3;
    (*err_vec3)(i)      = rms3;
    (*gaus_mean_vec3)(i)= gauss_mean3;
    (*gaus_err_vec3)(i) = gauss_rms3;

    file->cd();

    // Do the difference between 1 and 2
    double maximum = 0;
    maximum = std::max(hpost->GetMaximum(), hpost3->GetMaximum());
    maximum = std::max(maximum, hpost3->GetMaximum());

    hpost->SetMaximum(1.3*maximum);
    hpost2->SetMaximum(hpost->GetMaximum());
    hpost3->SetMaximum(hpost->GetMaximum());

    file->cd();
    hpost->Draw();
    hpost2->Draw("same");
    hpost3->Draw("same");
    asimov->Draw("same");
    leg->Draw("same");
    c0->Print(canvasname);

    // cd into params directory in root file
    params->cd();
    hpost->Write();
    hpost2->Write();
    hpost3->Write();

    delete hpost;
    delete asimov;
    delete leg;
  } // End the for ndraw loop


  file->cd();

  // vector of the mean values
  std::vector<double> val;
  // vector of the error values
  std::vector<double> vale;
  // vector of the index
  std::vector<double> ind;
  std::vector<double> inde;

  std::vector<double> val2;
  std::vector<double> vale2;
  std::vector<double> ind2;
  std::vector<double> inde2;

  std::vector<double> val3;
  std::vector<double> vale3;
  std::vector<double> ind3;
  std::vector<double> inde3;

  for(int i=0; i<mean_vec->GetNrows(); i++) {

    val.push_back((*mean_vec)(i));
    vale.push_back((*err_vec)(i));
    ind.push_back(i);
    inde.push_back(0.0);

    val2.push_back((*mean_vec2)(i));
    vale2.push_back((*err_vec2)(i));
    ind2.push_back(i);
    inde2.push_back(0.0);

    val3.push_back((*mean_vec3)(i));
    vale3.push_back((*err_vec3)(i));
    ind3.push_back(i);
    inde3.push_back(0.0);
  }
  TH1D *prefit = MakePrefit(ndraw ,nflux);
  file->cd();

  // Make a TH1D of the central values and the errors
  TH1D *paramPlot = new TH1D("paramPlot", "paramPlot", ndraw, 0, ndraw);
  paramPlot->SetName("mach3params");
  paramPlot->SetTitle(stepcut.c_str());
  paramPlot->SetFillStyle(3001);
  paramPlot->SetFillColor(kBlue-1);
  paramPlot->SetMarkerColor(kBlue+10);
  paramPlot->SetLineColor(paramPlot->GetFillColor());
  paramPlot->SetMarkerStyle(20);
  paramPlot->SetMarkerSize(prefit->GetMarkerSize());

  // Make a TH1D of the central values and the errors
  TH1D *paramPlot2 = (TH1D*)paramPlot->Clone();
  paramPlot2->SetName("mach3params_2");
  paramPlot2->SetTitle(paramPlot2->GetName());
  paramPlot2->SetFillColor(kOrange-5);
  paramPlot2->SetMarkerStyle(23);
  paramPlot2->SetLineWidth(2);
  paramPlot2->SetFillStyle(3244);
  paramPlot2->SetMarkerColor(kOrange+8);
  paramPlot2->SetLineColor(paramPlot2->GetFillColor());
  paramPlot2->SetMarkerSize(prefit->GetMarkerSize()*0.75);

  // Make a TH1D of the central values and the errors
  TH1D *paramPlot3 = (TH1D*)paramPlot->Clone();
  paramPlot3->SetName("mach3params_3");
  paramPlot3->SetTitle(paramPlot3->GetName());
  paramPlot3->SetMarkerColor(kBlack);
  paramPlot3->SetLineColor(paramPlot3->GetMarkerColor());
  paramPlot3->SetMarkerStyle(25);
  paramPlot3->SetLineWidth(2);
  paramPlot3->SetMarkerSize(prefit->GetMarkerSize()*0.5);
  paramPlot3->SetFillColor(0);
  paramPlot3->SetFillStyle(0);

  TLegend *graphLeg = new TLegend(0.33, 0.12, 0.80, 0.40);
  graphLeg->SetFillColor(0);
  graphLeg->SetFillStyle(0);
  graphLeg->SetLineColor(0);
  graphLeg->SetLineStyle(0);
  graphLeg->AddEntry(prefit, "Prefit", "fp");
  graphLeg->AddEntry(paramPlot,   title1.c_str(), "fp");
  graphLeg->AddEntry(paramPlot2,  title2.c_str(), "fp");
  graphLeg->AddEntry(paramPlot3,  title3.c_str(), "lfep");
  graphLeg->Draw("same");

  // Set labels and data
  for (int i = 0; i < ndraw; ++i) {

    paramPlot->SetBinContent(i+1, (*mean_vec)(i));
    paramPlot->SetBinError(i+1, (*err_vec)(i));

    paramPlot2->SetBinContent(i+1, (*mean_vec2)(i));
    paramPlot2->SetBinError(i+1, (*err_vec2)(i));

    paramPlot3->SetBinContent(i+1, (*mean_vec3)(i));
    paramPlot3->SetBinError(i+1, (*err_vec3)(i));

    file->cd();
    // Don't care about labelling beam for now
    if (i < nflux) {
      paramPlot->GetXaxis()->SetBinLabel(i+1, GetXsecName(i ).c_str());
      paramPlot2->GetXaxis()->SetBinLabel(i+1, GetXsecName(i ).c_str());
      paramPlot3->GetXaxis()->SetBinLabel(i+1, GetXsecName(i ).c_str());
    }
  }
  file->cd();

  // Plot the flux parameters (0 to 100) if enabled
  // Have already looked through the branches earlier
  if (plotFlux == true) {
    prefit->GetXaxis()->SetRangeUser(0, nflux);
    prefit->GetXaxis()->SetTitle("");

    prefit->GetYaxis()->SetRangeUser(0.7, 1.3);
    prefit->GetYaxis()->SetTitle("Variation");

    prefit->SetTitle((stepcut+"_"+stepcut2).c_str());
    int nplots = 4;
    for (int i = 0; i < nplots; ++i) {
      int bot = (double(i)/double(nplots))*nflux + nxsec;
      int top = (double(i+1)/double(nplots))*nflux + nxsec;

      prefit->GetXaxis()->SetRangeUser(bot, top);
      paramPlot->GetXaxis()->SetRangeUser(bot, top);
      paramPlot2->GetXaxis()->SetRangeUser(bot, top);
      paramPlot3->GetXaxis()->SetRangeUser(bot, top);

      prefit->Write(Form("param_flux_prefit_%i", i));
      paramPlot->Write(Form("param_flux_%i", i));
      paramPlot2->Write(Form("param_flux_gaus_%i", i));
      paramPlot3->Write(Form("param_flux_HPD_%i", i));

      prefit->Draw("e2");
      paramPlot->Draw("e2, same");
      paramPlot2->Draw("e2, same");
      paramPlot3->Draw("e1, same");

      graphLeg->Draw("same");
      c0->Write("param_flux_canv");
      c0->Print(canvasname);
      c0->Clear();
    }
  }

  file->cd();
  // Plot the xsec parameters (100 to ~125)
  // Have already looked through the branches earlier
  if (plotXsec == true) {
    prefit->GetYaxis()->SetRangeUser(-0.5, 2.0);
    prefit->GetXaxis()->SetTitle("");
    prefit->GetXaxis()->SetRangeUser(0, nxsec);
    prefit->GetXaxis()->LabelsOption("v");

    prefit->GetXaxis()->SetRangeUser(0, nxsec);
    paramPlot->GetXaxis()->SetRangeUser(0, nxsec);
    paramPlot2->GetXaxis()->SetRangeUser(0, nxsec);
    paramPlot3->GetXaxis()->SetRangeUser(0, nxsec);

    prefit->Draw("e2");
    paramPlot->Draw("e2,same");
    paramPlot2->Draw("e2,same");
    paramPlot3->Draw("e1,same");
    graphLeg->Draw("same");

    paramPlot->Write("param_xsec");
    paramPlot2->Write("param_xsec_2");
    paramPlot3->Write("param_xsec_3");
    c0->Print(canvasname);
    c0->Clear();
  }

  // Then close the pdf file
  std::cout << "Closing pdf " << canvasname << std::endl;
  canvasname+="]";
  c0->Print(canvasname);

  // Write all the nice vectors
  file->cd();
  mean_vec->Write("FitParameters");
  mean_vec2->Write("FitParameters_2");
  mean_vec3->Write("FitParameters_3");
  err_vec->Write("FitErrors");
  err_vec2->Write("FitErrors_2");
  err_vec3->Write("FitErrors_3");
  gaus_mean_vec->Write("FitParameters_gauss");
  gaus_mean_vec2->Write("FitParameters_gauss_2");
  gaus_mean_vec3->Write("FitParameters_gauss_3");
  gaus_err_vec->Write("FitErrors_gauss");
  gaus_err_vec2->Write("FitErrors_gauss_2");
  gaus_err_vec3->Write("FitErrors_gauss_3");
  paramPlot->Write("postfit_params_plot");
  paramPlot2->Write("postfit_params_plot_2");
  paramPlot3->Write("postfit_params_plot_3");

  file->Close();
}


// *****************************
// Make the prefit plots
TH1D* MakePrefit(int ndraw, int nFlux) {
  // *****************************

  // Get Xsec matrix (xsec+flux params)
  TFile *XsecFile = new TFile(XsecCov.c_str(), "open");
  XsecFile->cd();
  // Get the matrix
  TMatrixDSym *XsecMatrix = (TMatrixDSym*)(XsecFile->Get("xsec_cov"));
  // And the central priors
  TVectorD *XsecPrior = (TVectorD*)(XsecFile->Get("xsec_param_prior"));
  // And the nominal values
  TVectorD *XsecNominal = (TVectorD*)(XsecFile->Get("xsec_param_nom"));
  TMatrixT<double> *XsecID = (TMatrixD*)(XsecFile->Get("xsec_param_id"));

  // Now make a TH1D of it
  int nXsec = XsecPrior->GetNrows() - nFlux;

  // Now drive in the actually fill
  double* XsecPrefitArray = new double[nXsec];
  // And the errors from the covariance matrix
  double* XsecErrorArray = new double[nXsec];

  for (int i = 0; i < nXsec; ++i) {
    // Normalise the prior relative the nominal: this is just by choice for BANFF comparisons
    if (((*XsecNominal)(i)) != 0) {
      XsecPrefitArray[i] = ((*XsecPrior)(i)) / ((*XsecNominal)(i));
      XsecErrorArray[i] = sqrt((*XsecMatrix)(i,i))/((*XsecNominal)(i));
    } else {
      XsecPrefitArray[i] = ((*XsecPrior)(i)) + 1.0;
      XsecErrorArray[i] = sqrt((*XsecMatrix)(i,i));
    }
  }
  
  double *FluxPrefitArray = new double[nFlux];
  double *FluxErrorArray = new double[nFlux];
  for (int i = 0; i < nFlux; ++i) {
    // The nominals are equal to the priors are equal to 1 for flux params
    FluxPrefitArray[i] = 1.0;
    FluxErrorArray[i] = sqrt((*XsecMatrix)(i+nXsec,i+nXsec));
  }
  
  if (ndraw != nFlux+nXsec) {
    std::cerr << "Number of requested drawn parameters are not equal to flux and xsec input" << std::endl;
    std::cerr << "ndraw = " << ndraw << std::endl;
    std::cerr << "nFlux = " << nFlux << "   " << "nXsec = " << nXsec << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  TH1D *PreFitPlot = new TH1D("Prefit", "Prefit", ndraw, 0, ndraw);
  for (int i = 0; i < PreFitPlot->GetNbinsX() + 1; ++i) {
    PreFitPlot->SetBinContent(i+1, 0);
    PreFitPlot->SetBinError(i+1, 0);
  }
  for (int i = 0; i < nXsec; ++i) {
    PreFitPlot->SetBinContent(i+1, XsecPrefitArray[i]);
    PreFitPlot->SetBinError(i+1, XsecErrorArray[i]);

    PreFitPlot->GetXaxis()->SetBinLabel(i+1, GetXsecName(i).c_str());
    // Also scale xsec entries relative the error
  }
  // Set the flux values
  for (int i = nXsec; i < nFlux+nXsec; ++i) {
    PreFitPlot->SetBinContent(i+1, FluxPrefitArray[i-nXsec]);
    PreFitPlot->SetBinError(i+1, FluxErrorArray[i-nXsec]);

    if (i % 5 == 0) PreFitPlot->GetXaxis()->SetBinLabel(i+1, Form("Flux %i", i-nXsec));
  }


  PreFitPlot->GetYaxis()->SetTitle("Variation rel. nominal");

  PreFitPlot->SetDirectory(0);

  PreFitPlot->SetFillStyle(1001);
  PreFitPlot->SetFillColor(kRed-3);
  PreFitPlot->SetMarkerStyle(21);
  PreFitPlot->SetMarkerSize(2.7);
  PreFitPlot->SetMarkerColor(kWhite);
  PreFitPlot->SetLineColor(PreFitPlot->GetFillColor());

  PreFitPlot->GetXaxis()->LabelsOption("v");

  delete XsecMatrix;
  delete XsecPrior;
  delete XsecNominal;
  delete XsecID;
  delete[] XsecPrefitArray;
  delete[] XsecErrorArray;

  XsecFile->Close();
  delete XsecFile;

  return PreFitPlot;
}


// **************************
// Function to get limits for a parameter from the input
// Hard-code to +/-2 sigma
void GetXsecLimits(const int param, double &central, double &prior, double &down_error, double &up_error) {
  // **************************

  central = 0.0;
  double error = 0.0;
  double sigmas = 1.0;

  // Open the input covariance file to get the pre-fit error
  TFile *covFile = new TFile(XsecCov.c_str(), "OPEN");
  covFile->cd();
  TMatrixDSym *covMatrix = (TMatrixDSym*)(covFile->Get("xsec_cov"));

  TVectorD* prior_a = (TVectorD*)(covFile->Get("xsec_param_prior"));
  TVectorD* nominal = (TVectorD*)(covFile->Get("xsec_param_nom"));
  TVectorD* lower   = (TVectorD*)(covFile->Get("xsec_param_lb"));
  TVectorD* upper   = (TVectorD*)(covFile->Get("xsec_param_ub"));
  TMatrixT<double> *id = (TMatrixD*)(covFile->Get("xsec_param_id"));

  if ( ((*nominal)(param)) != 0 ) {
    // Set the central to the nominal
    central = ((*nominal)(param));
    prior = ((*prior_a)(param));
    error = sigmas * sqrt(((*covMatrix)(param,param))) / central;
  } else {
    central = ((*nominal)(param));
    prior = ((*prior_a)(param));
    error = sigmas * sqrt(((*covMatrix)(param,param)));
  }

  // We might be passed the valid range of parameter
  // Do a check to see if this is true
  if (central - error <  ((*lower)(param))) {
    down_error = (*lower)(param);
  } else {
    down_error = central - error;
  }

  if (central + error > ((*upper)(param))) {
    up_error = (*upper)(param);
  } else {
    up_error = central + error;
  }

  // Normalisation parameter are unlikely of going above 2 sigma
  if ((*id)(param,0) == -1) {
    up_error = central + 2*sqrt((*covMatrix)(param, param));
  }

  delete covMatrix;
  delete prior_a;
  delete nominal;
  delete id;

  covFile->Close();
  delete covFile;
}



// **************************
// Just a converter from xsec_i to a real parameter name
// Use this when the conversion wasn't written to file :(
std::string GetXsecName(int xsecParamNo) {
  // **************************

  CheckXsecCov();
  // Open the input covariance file to get the pre-fit error
  TFile *covFile = new TFile(XsecCov.c_str(), "OPEN");
  covFile->cd();

  // Get the array of cross-section parameter names
  TObjArray *xsec_param_names = (TObjArray*)(covFile->Get("xsec_param_names"));

  std::string returnString;
  if (xsecParamNo >= xsec_param_names->GetEntries()) {
    std::cerr << "You gave me a cross-section parameter number which was greater than what the xsec covariance matrix has dimensions" << std::endl;
    std::cerr << "Please revise!" << std::endl;
    std::cerr << "xsecParamNo = " << xsecParamNo << std::endl;
    std::cerr << "GetEntries  = " << xsec_param_names->GetEntries() << std::endl;
    returnString = "INVALID";
  } else {
    returnString = std::string(((TObjString*)xsec_param_names->At(xsecParamNo))->GetString());
  }

  covFile->Close();

  return returnString;

} // End function


// **************************
// Function to check the cross-section covariance matrix
void CheckXsecCov() {
  // **************************

  if (XsecCov.empty()) {
    std::cerr << "Have not got a name for the cross-section covariance matrix" << std::endl;
    std::cerr << "Please specify as a global variable at the top of " << __FILE__ << std::endl;
    throw;
  }

}

// Get the highest posterior density from a TH1D
void GetHPD(TH1D * const hpost, double &central, double &error, double &error_pos, double &error_neg) {

  // Get the bin which has the largest posterior density
  int MaxBin = hpost->GetMaximumBin();
  // And it's value
  double peakval = hpost->GetBinCenter(MaxBin);

  // The total integral of the posterior
  double integral = hpost->Integral();

  // Keep count of how much area we're covering
  double sum = 0.0;

  // Counter for current bin
  int CurrBin = MaxBin;
  while (sum/integral < 0.6827/2.0 && CurrBin < hpost->GetNbinsX()+1) {
    sum += hpost->GetBinContent(CurrBin);
    CurrBin++;
  }
  double sigma_p = fabs(hpost->GetBinCenter(MaxBin)-hpost->GetBinCenter(CurrBin));
  // Reset the sum
  sum = 0.0;

  // Reset the bin counter
  CurrBin = MaxBin;
  // Counter for current bin
  while (sum/integral < 0.6827/2.0 && CurrBin >= 0) {
    sum += hpost->GetBinContent(CurrBin);
    CurrBin--;
  }
  double sigma_m = fabs(hpost->GetBinCenter(CurrBin)-hpost->GetBinCenter(MaxBin));

  // Now do the double sided HPD
  sum = 0.0;
  int LowBin = MaxBin-1;
  int HighBin = MaxBin+1;
  double LowCon = 0.0;
  double HighCon = 0.0;
  while (sum/integral < 0.6827 && LowBin >= 0 && HighBin < hpost->GetNbinsX()+1) {

    // Get the slice
    LowCon = hpost->GetBinContent(LowBin);
    HighCon = hpost->GetBinContent(HighBin);

    // If we're on the last slice and the lower contour is larger than the upper
    if ((sum+LowCon+HighCon)/integral > 0.6827 && LowCon > HighCon) {
      sum += LowCon;
      break;
      // If we're on the last slice and the upper contour is larger than the lower
    } else if ((sum+LowCon+HighCon)/integral > 0.6827 && HighCon >= LowCon) {
      sum += HighCon;
      break;
    } else {
      sum += LowCon + HighCon;
    }

    LowBin--;
    HighBin++;
  }

  double sigma_hpd = 0.0;
  if (LowCon > HighCon) {
    sigma_hpd = fabs(hpost->GetBinCenter(LowBin)-hpost->GetBinCenter(MaxBin));
  } else {
    sigma_hpd = fabs(hpost->GetBinCenter(HighBin)-hpost->GetBinCenter(MaxBin));
  }

  central = peakval;
  error = sigma_hpd;
  error_pos = sigma_p;
  error_neg = sigma_m;

}

// **************************
// Get the mean and RMS of a 1D posterior
void GetArithmetic(TH1D * const hpost, double &mean, double &error) {
  // **************************
  mean = hpost->GetMean();
  error = hpost->GetRMS();
}

// **************************
// Get Gaussian characteristics
void GetGaussian(TH1D *& hpost, TF1 *& gauss, double &central, double &error) {
  // **************************

  double mean = hpost->GetMean();
  double err = hpost->GetRMS();
  double peakval = hpost->GetBinCenter(hpost->GetMaximumBin());

  // Set the range for the Gaussian fit
  gauss->SetRange(mean - 1.5*err , mean + 1.5*err);
  // Set the starting parameters close to RMS and peaks of the histograms
  gauss->SetParameters(hpost->GetMaximum()*err*sqrt(2*3.14), peakval, err);

  // Perform the fit
  hpost->Fit(gauss->GetName(),"Rq");
  hpost->SetStats(0);

  central = gauss->GetParameter(1);
  error = gauss->GetParameter(2);

}
