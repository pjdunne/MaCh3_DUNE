#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>

void makeEfficiencyvsRadiusPlots(TString inputfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
   TFile* file = new TFile(inputfile);
   TList* list = file->GetListOfKeys();
   

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   TCanvas* c0 = new TCanvas("c0","c0",0,0,700,900);
   c0->Divide(1,4);
   c0->Print("efficiencyvsradius.ps[");
   std::vector<std::vector<TString>> kinematicvariables;

  std::vector<std::string> plotvariables = {"RecoRad"};

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

   TH1D* trueselected = (TH1D*)file->Get("RecoRad_NDGAr_all_TrueSelected");
   TH1D* alltrue = (TH1D*)file->Get("RecoRad_NDGAr_all_AllTrue");;
//   TH1D* efficiency = (TH1D*)trueselected->Clone();
   std::vector<double> x, y, ex, ey, x1, y1, ex1, ey1, x2, y2, ex2, ey2;

//   TGraphErrors* efficiency = new TGraphErrors();
//   TGraphErrors* efficiency_binbybin = new TGraphErrors();
//   TGraphErrors* effvsrad2 = new TGraphErrors();

//   TH1D* efficiency_binbybin = (TH1D*)trueselected->Clone();
//   efficiency->Reset();
//   efficiency_binbybin->Reset();

   int nbins = trueselected->GetNbinsX();
   double truselec =0;
   double alltru = 0;
   double eff, rad, absolute_eff, eff_per5cm, trusel_5, alltru_5, error, error2, error3, error4, error5;
   for(int bin =1; bin<nbins+1; bin++){
     truselec = truselec + trueselected->GetBinContent(bin);
     alltru = alltru + alltrue->GetBinContent(bin);
     trusel_5 = trusel_5 + trueselected->GetBinContent(bin);
     alltru_5 = alltru_5 + alltrue->GetBinContent(bin);
     eff = truselec/alltru;
     error = pow(alltrue->GetBinContent(bin), 0.5);
     error2 = pow(trueselected->GetBinContent(bin), 0.5);
     trueselected->SetBinError(bin,error2);
     alltrue->SetBinError(bin,error);
     if(error !=0 && error2 !=0){
     error3 = pow((pow(1/error, 2)+(pow(1/error2, 2))), 0.5);
     }
     else{error3 = 0;}
     if(truselec !=0 && alltru != 0){
     error4 = pow((1/truselec + 1/alltru), 0.5);
     }
     else{error4 = 0;}
     rad = trueselected->GetXaxis()->GetBinCenter(bin);
     std::cout<<"error3: "<<error3<<std::endl;
     std::cout<<"error4: "<<error4<<std::endl;
     std::cout<<"rad: "<<rad<<" truselec_5: "<<trusel_5<<" alltru_5: "<<alltru_5<<" eff: "<<eff<<std::endl;
     if(std::isnan(eff)){eff =0;}
//     efficiency->SetBinContent(bin, eff);
     //efficiency->AddPoint(rad, eff);
     y.push_back(eff);
     x.push_back(rad);
     ex.push_back(2.5);
     ey.push_back(error4);
     std::cout<<"here"<<std::endl;
     //efficiency->SetPointError(efficiency->GetN(), 1, error4);
//     efficiency->SetBinError(bin, 1/pow(bin, 0.5));
     if(alltrue->GetBinContent(bin) != 0) {
       absolute_eff = trueselected->GetBinContent(bin)/alltrue->GetBinContent(bin);
       y1.push_back(absolute_eff);
       x1.push_back(rad);
       ex1.push_back(2.5);
       ey1.push_back(error3);
//       efficiency_binbybin->AddPoint(rad, absolute_eff);
//       std::cout<<"point binbybin: "<<efficiency_binbybin->GetN()<<std::endl;
//       efficiency_binbybin->SetPointError((int)(efficiency_binbybin->GetN()),1, error3);
//       std::cout<<"geterror: "<<efficiency_binbybin->GetErrorX((int)(efficiency_binbybin->GetN()))<<" "<<efficiency_binbybin->GetErrorX(bin)<<std::endl;
       //effvsrad2->AddPoint(pow(rad, 2), absolute_eff);
       y2.push_back(absolute_eff);
       x2.push_back(pow(rad, 2));
       ex2.push_back(pow(2.5, 2));
       ey2.push_back(error3);
    
     }
/*
     if((bin+1)%5==0){
       if(alltru_5 !=0){
         eff_per5cm = trusel_5/alltru_5;
         if(trusel_5 != 0){
         error5 = pow((1/trusel_5 + 1/alltru_5), 0.5);
         }
         else{error5 = 0;}
         y2.push_back(eff_per5cm);
         x2.push_back(pow(rad-2.5, 2));
         ex2.push_back(pow(2.5, 2));
         ey2.push_back(error5);
         std::cout<<"y2: "<<eff_per5cm<<" x2: "<<(pow(rad-2.5, 2))<<" ex2: "<<pow(2.5, 2)<<" ey2: "<<error5<<std::endl;
         //effvsrad2->AddPoint(pow(rad-2.5, 2), eff_per5cm);
         //effvsrad2->SetPointError(effvsrad2->GetN(), pow(2.5, 2), error5);
       }
       trusel_5 =0;
       alltru_5 =0;
     }*/
   }
   auto legend = new TLegend(0.8,0.8,0.9,0.9);
   TGraphErrors* efficiency = new TGraphErrors(x.size(), &x[0], &y[0], &ex[0], &ey[0]);
   TGraphErrors* efficiency_binbybin = new TGraphErrors(x1.size(), &x1[0], &y1[0], &ex1[0], &ey1[0]);
   TGraphErrors* effvsrad2 = new TGraphErrors(x2.size(), &x2[0], &y2[0], &ex2[0], &ey2[0]);

   std::cout<<" points: "<<efficiency->GetN()<<" "<<efficiency_binbybin->GetN()<<" "<<effvsrad2->GetN()<<std::endl;
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
//   trueselected->Sumw2();
//   leg3->Draw("SAME");
   c0->Update();
//   c0->Print("efficiencyvsradius.ps");
   c0->cd(2);
//   efficiency->Divide(trueselected, alltrue);
//   efficiency->SetTitle("Efficiency");
//   efficiency->Sumw2();
//   trueselected->Divide(alltrue);
//   trueselected->SetTitle("Efficiency");
   efficiency->SetTitle("Cumulative Efficiency");
   efficiency->GetXaxis()->SetTitle("Radius /cm");
   efficiency->GetYaxis()->SetTitle("Efficiency");
   efficiency->SetMarkerStyle(kPlus);
   efficiency->SetMarkerSize(0.1);
   efficiency->Draw("AP");
   c0->cd(3);
//   efficiency_binbybin->Divide(trueselected, alltrue);
//   efficiency_binbybin->Sumw2();
   efficiency_binbybin->SetTitle("Absolute Efficiency");
   efficiency_binbybin->GetXaxis()->SetTitle("Radius /cm");
   efficiency_binbybin->GetYaxis()->SetTitle("Efficiency");
   efficiency_binbybin->SetMarkerStyle(kPlus);
   efficiency_binbybin->SetMarkerSize(0.1);
   efficiency_binbybin->Draw("AP");
   c0->Update();
   c0->cd(4);
   gPad->SetLogx();
   effvsrad2->SetTitle("Efficiency vs Radius Squared");
   effvsrad2->GetXaxis()->SetTitle("Radius Squared /cm^2");
   effvsrad2->GetYaxis()->SetTitle("Efficiency");
   effvsrad2->SetMarkerStyle(kPlus);
   effvsrad2->SetMarkerSize(0.3);   
   effvsrad2->Draw("AP");
   c0->Update();
//   c0->SetLogx(0);
   c0->Print("efficiencyvsradius.ps");
   std::cout<<"HERE"<<std::endl;
   c0->Print("efficiencyvsradius.ps]");

}
