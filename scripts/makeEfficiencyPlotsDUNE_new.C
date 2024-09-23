#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>

void makeEfficiencyPlotsDUNE_new(TString inputfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
   TFile* file = new TFile(inputfile);
   TList* list = file->GetListOfKeys();
   

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   TCanvas* c0 = new TCanvas("c0","c0",0,0,700,900);
   c0->Divide(1,2);
   c0->Print("efficiencyvars.ps[");
   std::vector<std::vector<TString>> kinematicvariables;

  std::vector<std::string> plotvariables = {"TrueNeutrinoEnergy"};

//  std::vector<std::string> plotvariables = {"RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "RecoXPos", "TrueYPos", "RecoYPos", "TrueZPos", "RecoZPos", "TrueRad", "RecoRad", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy"};

  std::vector<std::string> x_axis;
//   std::vector<std::string> plotvariables = {"TrueNeutrinoEnergy", "RecoNeutrinoEnergy", "TrueMinusRecoEnergy", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "TrueYPos", "TrueZPos", "NMuons", "TrueMinusRecoEnergyRatio", "RecoLepEnergy"};
   int n_vars = plotvariables.size();
   for(int iplots =0; iplots<n_vars; iplots++){
      std::vector<TString> keynames;
      kinematicvariables.push_back(keynames);
      std::string xaxisname= plotvariables[iplots];
      if(plotvariables[iplots].find("Ratio")!= std::string::npos){xaxisname += "";}
      else if(plotvariables[iplots].find("Energy")!= std::string::npos){xaxisname += " (GeV)";}
      else if(plotvariables[iplots].find("Pos")!= std::string::npos){xaxisname += " (cm)";}
      else if(plotvariables[iplots].find("Rad")!= std::string::npos){xaxisname += " (cm)";}
      else if(plotvariables[iplots].find("Muons")!= std::string::npos){xaxisname += "";}
      x_axis.push_back(xaxisname);
   }

  for(int i=0; i<list->GetEntries(); i++)
   {
     bool boolcontinue = false;
     //std::cout << "entry " << i << std::endl;
     TString keyname = list->At(i)->GetName();
        if(keyname.Contains("_NDGAr_")){
        for(int iplots=0; iplots<n_vars; iplots++){
          TString keynametest = plotvariables[iplots]+"_NDGAr_";
          if(keyname.Contains(keynametest)){
          kinematicvariables[iplots].push_back(keyname);
          std::cout<<"iplots: "<<iplots<<" keyname: "<<keyname<<std::endl;
          boolcontinue = true;
          break;
          }
        }

        if(!boolcontinue){continue;}
      }
   }
   std::vector<int> ccqe = {0};
   std::vector<int> ccdis = {2};
   std::vector<int> ccres ={3};
   std::vector<int> ccmec ={9};
   std::vector<int> ccother ={1, 4, 5, 6, 7, 8, 10, 11, 12, 13};
   std::vector<int> nc = {14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
   
   std::vector<std::vector<int>> modebreakdown;
   modebreakdown.push_back(ccqe);
   modebreakdown.push_back(ccdis);
   modebreakdown.push_back(ccres);
   modebreakdown.push_back(ccmec);
   modebreakdown.push_back(ccother);
   modebreakdown.push_back(nc);

   std::vector<std::string> modenames = {"CCQE", "CCDIS", "CCRES", "CCMEC", "CCOTHER", "NC"};
   std::vector<int> modecolours = {900-1, 616-6, 880-1, 860-5, 432-7, kGreen+2};


   TH1D* trueselected = (TH1D*)file->Get("TrueNeutrinoEnergy_NDGAr_all_TrueSelected");
   TH1D* alltrue = (TH1D*)file->Get("TrueNeutrinoEnergy_NDGAr_all_AllTrue");;
//   TH1D* efficiency = (TH1D*)file->Get("TrueNeutrinoEnergy_NDGAr_all_TrueSelected");
   std::vector<double> x1, y1, ex1, ey1;

   int nbins = trueselected->GetNbinsX();
   double truselec =0;
   double alltru = 0;
   double eff, energy, error, error2, error3;
   for(int bin =1; bin<nbins+1; bin++){
     error = pow(alltrue->GetBinContent(bin), 0.5);
     error2 = pow(trueselected->GetBinContent(bin), 0.5);
     trueselected->SetBinError(bin,error2);
     alltrue->SetBinError(bin,error);
     if(error !=0 && error2 !=0){
     error3 = pow((pow(1/error, 2)+(pow(1/error2, 2))), 0.5);
     }
     else{error3 = 0;}
     energy = trueselected->GetXaxis()->GetBinCenter(bin);

     if(alltrue->GetBinContent(bin) != 0) {
       eff = trueselected->GetBinContent(bin)/alltrue->GetBinContent(bin);
       y1.push_back(eff);
       x1.push_back(energy);
       ex1.push_back(0.1);
       ey1.push_back(error3);
     }
   }
   TGraphErrors* efficiency = new TGraphErrors(x1.size(), &x1[0], &y1[0], &ex1[0], &ey1[0]);

   auto legend = new TLegend(0.8,0.8,0.9,0.9);
   
   c0->cd(1);
   alltrue->SetMarkerStyle(kPlus);
   alltrue->SetMarkerColor(kBlue);
   alltrue->SetLineColor(kBlue);
   alltrue->Draw("HIST");
//   alltrue->Sumw2();
//   leg2->Draw("SAME");
//   c0->Update();
//   c0->cd(2);
    trueselected->SetMarkerStyle(kPlus);
    trueselected->SetMarkerColor(kRed);
    trueselected->SetLineColor(kRed);
    trueselected->Draw("EP SAME");
    legend->AddEntry(alltrue, "All True CC1mu", "l");
    legend->AddEntry(trueselected, "True Selected CC1mu", "lp");
    legend->Draw();
   c0->cd(2);
//   efficiency->Divide(trueselected, alltrue);
   efficiency->SetMarkerColor(kRed);
   efficiency->SetMarkerStyle(kPlus);
   efficiency->SetTitle("Efficiency");
//   efficiency->Sumw2();
//   trueselected->Divide(alltrue);
//   trueselected->SetTitle("Efficiency");
   efficiency->Draw("AP");
//   c0->Update();
   c0->Print("efficiencyvars.ps");
   std::cout<<"HERE"<<std::endl;
   c0->Print("efficiencyvars.ps]");

}
