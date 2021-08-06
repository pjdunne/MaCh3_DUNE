#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"

#include "TFile.h"
#include "TH1D.h"

#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>

using namespace std;
int main(int argc, char * argv[] ) {

  double total_prob_e=0.0;
  double total_prob_mu=0.0;
  double energy;
  double e_start, e_end, e_step;
  int i, j ;

  // Binning	
  int NBinsEnergy = 10000;

  double BasePath = 1300;
  double Density = 2.84;

  // Energy Range
  double EnergyBins[NBinsEnergy+1];
  e_start = 0.01;
  e_end = 100.;
  //e_step = log10(e_end/e_start)/double(NBinsEnergy);
  e_step = (e_end-e_start)/double(NBinsEnergy);
	
  // Oscillation Parameters
  bool kSquared = true;   // using sin^2(x) variables?
  double DM2 = 2.451e-3;
  double dm2 = 7.39e-5;

  //th23 = 0.866 --> sin^2(0.866)=0.58025
  double Theta23 = 0.58025; // sin^2 (theta23)
  //th13 = 0.150 --> sin^2(0.150)=0.02233
  double Theta13 = 0.02233; // sin^2 (theta13)
  //th12 = 0.5903 --> sin^2(0.5903)=0.30982
  double Theta12 = 0.30982; // sin^2 (theta12)

  std::cout << "Using          " << std::endl
            << "      DM2      " <<  DM2      << std::endl
            << "      Theta23  " <<  Theta23  << std::endl
            << "      Theta13  " <<  Theta13  << std::endl
            << "      dm2      " <<  dm2      << std::endl
            << "      Theta12  " <<  Theta12  << std::endl;

  BargerPropagator   * bNu; 

  bNu = new BargerPropagator( );
  bNu->UseMassEigenstates( false );

  double Entry = e_start;
  for(i=0; i<NBinsEnergy; i++ ) {
    //Entry = e_start*pow( 10.0 , double(i)*e_step );
    Entry = e_start+double(i)*e_step;
    EnergyBins[i] = Entry;
  }

  EnergyBins[NBinsEnergy] = EnergyBins[NBinsEnergy-1]*1.001;

  stringstream ssE;
  TH1D * nu_histos[2][3];

  ///////////////////////////
  /// e to e
  ssE.str(""); ssE <<  "P(#nu_{e} #rightarrow #nu_{e})" << " L = " << BasePath ; 
  TH1D * le2eE = new TH1D("nue_x_nue_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// e to mu 
  ssE.str(""); ssE <<  "P(#nu_{e} #rightarrow #nu_{#mu})" << " L = " << BasePath ; 
  TH1D * le2muE = new TH1D("nue_x_numu_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// e to tau 
  ssE.str(""); ssE <<  "P(#nu_{e} #rightarrow #nu_{#tau})" << " L = " << BasePath ; 
  TH1D * le2tauE = new TH1D("nue_x_nutau_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// mu to e
  ssE.str(""); ssE <<  "P(#nu_{#mu} #rightarrow #nu_{e})" << " L = " << BasePath ; 
  TH1D * lmu2eE = new TH1D("numu_x_nue_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// mu to mu 
  ssE.str(""); ssE <<  "P(#nu_{#mu} #rightarrow #nu_{#mu})" << " L = " << BasePath ; 
  TH1D * lmu2muE = new TH1D("numu_x_numu_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// mu to tau 
  ssE.str(""); ssE <<  "P(#nu_{#mu} #rightarrow #nu_{#tau})" << " L = " << BasePath ; 
  TH1D * lmu2tauE = new TH1D("numu_x_nutau_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  nu_histos[0][0] = le2eE; 
  nu_histos[0][1] = le2muE; 
  nu_histos[0][2] = le2tauE; 

  nu_histos[1][0] = lmu2eE;
  nu_histos[1][1] = lmu2muE;
  nu_histos[1][2] = lmu2tauE;

  for ( i = 0 ; i <= NBinsEnergy ; i ++ ) {
    //energy = e_start*pow(10.0, double(i)*e_step);
    energy = e_start + double(i)*e_step;
    bNu->SetMNS(Theta12, Theta13, Theta23, dm2, DM2, 0.0 , energy, kSquared); // Compare with 3Flavor05.ps//
    bNu->propagateLinear( 1, BasePath, Density );

    total_prob_e = 0.0;
    total_prob_mu = 0.0;
    for(int m=1; m<=3; m++) {
      total_prob_e += bNu->GetProb(1, m); // Normalize the Probabilities //
      total_prob_mu += bNu->GetProb(2, m); // Normalize the Probabilities //
    }

    if ( total_prob_e>1.00001 || total_prob_e<0.99998 || total_prob_mu>1.00001 || total_prob_mu<0.99998) {
      cerr << "ERROR Prob:" << "Energy: "<< energy << " " << endl;   abort();
    }

    for( j = 0 ; j < 3 ; j++ ) {
      //histos[j][1]->SetBinContent( i+1 , bNu->GetProb(2,j+1) );	
      nu_histos[0][j]->Fill( energy, bNu->GetProb(1,j+1) );	
      nu_histos[1][j]->Fill( energy, bNu->GetProb(2,j+1) );	
    }
  } // End Energy Loop 

  // Now do antineutrino:

  TH1D * nubar_histos[2][3];

  ///////////////////////////
  /// e to e
  ssE.str(""); ssE <<  "P(#bar{#nu}_{e} #rightarrow #bar{#nu}_{e})" << " L = " << BasePath ; 
  TH1D * le2eE_bar = new TH1D("nuebar_x_nuebar_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// e to mu 
  ssE.str(""); ssE <<  "P(#bar{#nu}_{e} #rightarrow #bar{#nu}_{#mu})" << " L = " << BasePath ; 
  TH1D * le2muE_bar = new TH1D("nuebar_x_numubar_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// e to tau 
  ssE.str(""); ssE <<  "P(#bar{#nu}_{e} #rightarrow #bar{#nu}_{#tau})" << " L = " << BasePath ; 
  TH1D * le2tauE_bar = new TH1D("nuebar_x_nutaubar_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// mu to e
  ssE.str(""); ssE <<  "P(#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e})" << " L = " << BasePath ; 
  TH1D * lmu2eE_bar = new TH1D("numubar_x_nuebar_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// mu to mu 
  ssE.str(""); ssE <<  "P(#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{#mu})" << " L = " << BasePath ; 
  TH1D * lmu2muE_bar = new TH1D("numubar_x_numubar_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  ///////////////////////////
  /// mu to tau 
  ssE.str(""); ssE <<  "P(#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{#tau})" << " L = " << BasePath ; 
  TH1D * lmu2tauE_bar = new TH1D("numubar_x_nutaubar_oscProbVsE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );

  nubar_histos[0][0] = le2eE_bar;
  nubar_histos[0][1] = le2muE_bar;
  nubar_histos[0][2] = le2tauE_bar;

  nubar_histos[1][0] = lmu2eE_bar;
  nubar_histos[1][1] = lmu2muE_bar;
  nubar_histos[1][2] = lmu2tauE_bar;

  for ( i = 0 ; i <= NBinsEnergy ; i ++ ) {
    //energy = e_start*pow(10.0, double(i)*e_step);
    energy = e_start + double(i)*e_step;
    bNu->SetMNS(Theta12, Theta13, Theta23, dm2, DM2, 0.0 , energy, kSquared, -1); // Compare with 3Flavor05.ps//
    bNu->propagateLinear( -1, BasePath, Density );

    total_prob_e = 0.0;
    total_prob_mu = 0.0;
    for(int m=1; m<=3; m++) {
      total_prob_e += bNu->GetProb(1, m); // Normalize the Probabilities //
      total_prob_mu += bNu->GetProb(2, m); // Normalize the Probabilities //
    }

    if ( total_prob_e>1.00001 || total_prob_e<0.99998 || total_prob_mu>1.00001 || total_prob_mu<0.99998)
    { cerr << "ERROR Prob:" << "Energy: "<< energy << " " << endl;   abort();   }

    for( j = 0 ; j < 3 ; j++ ) {
      //histos[j][1]->SetBinContent( i+1 , bNu->GetProb(2,j+1) );	
      nubar_histos[0][j]->Fill( energy, bNu->GetProb(1,j+1) );	
      nubar_histos[1][j]->Fill( energy, bNu->GetProb(2,j+1) );	
    }

   } // End Energy Loop 

   TFile *tmp = new TFile("DUNE_oscWeights.root", "recreate");
   tmp->cd();

   for( j = 0 ; j < 3 ; j++ ){
     nu_histos[0][j]->Write();	
     nubar_histos[0][j]->Write();	
     nu_histos[1][j]->Write();	
     nubar_histos[1][j]->Write();	
   }

   tmp->Write();
   tmp->Close();

   cout << endl<<"Done Cowboy!" << endl;
   return 0;
}

