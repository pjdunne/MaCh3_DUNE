#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>

void makeEfficiencyvsEnergyPlotsDUNE(TString inputfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
   TFile* file = new TFile(inputfile);
   TList* list = file->GetListOfKeys();
   

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   TCanvas* c0 = new TCanvas("c0","c0",0,0,700,900);
   c0->Divide(1,2);
   c0->Print("efficiencyenergyvars.ps[");
   std::vector<std::vector<TString>> kinematicvariables;

  std::vector<std::string> plotvariables = {"NTrueMuons", "NRecoMuons"};

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
   std::cout<<"here"<<std::endl;
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
   THStack* recostack = new THStack("recostack","NRecoMuons with Reco Cuts");    
   TLegend* leg2 = new TLegend(0.65,0.55,0.8,0.8);
   THStack* truestack = new THStack("truestack","NTrueMuons with Truth Cuts"); 
   TLegend* leg3 = new TLegend(0.65,0.55,0.8,0.8);
   TH1D* truenmuons;
   TH1D* reconmuons;
   TH1D* effvsenergy = new TH1D("effvsenergy", "Efficiency vs True Neutrino Energy", 10, 0, 10);
   std::cout<<"kinematicvariables.size(): "<<kinematicvariables.size()<<std::endl; 
   int energyval =0;
   bool onefinish = false;
   bool twofinish = false;
   for(int ienergy =0; ienergy<10; ienergy++){
   for(int ivars=0; ivars<kinematicvariables.size(); ivars++){
    std::cout<<"kinematicvariables[ivars].size(): "<<kinematicvariables[ivars].size()<<std::endl;  
    for(int ikeynames = 0; ikeynames<kinematicvariables[ivars].size(); ikeynames++){
      //std::cout<<"name: "<<kinematicvariables[ivars][ikeynames]<<std::endl;
      //std::cout<<"ikeyname: "<<ikeynames<<std::endl;
      std::string energy = "_TRUTH_" + std::to_string(energyval);
      std::string energy2 = "_RECO_" + std::to_string(energyval);
      TH1D* hist3;
      TH1D* hist4;
//     int counter =0;
      if(kinematicvariables[ivars][0].Contains("NTrueMuons")){
        if(kinematicvariables[ivars][ikeynames].Contains(energy.c_str())){
         hist3 =(TH1D*)file->Get(kinematicvariables[ivars][ikeynames]);
         std::cout<<"true name: "<<kinematicvariables[ivars][ikeynames]<<std::endl;
         onefinish = true;
//         if(counter ==0){truenmuons = (TH1D*)hist3->Clone("NTrueMuons"); truenmuons->Reset();counter++;}
//         truenmuons->Add(hist3);
         }
       }
//     counter =0;
     if(kinematicvariables[ivars][0].Contains("NRecoMuons")){
       if(kinematicvariables[ivars][ikeynames].Contains(energy2.c_str())){
         hist4 =(TH1D*)file->Get(kinematicvariables[ivars][ikeynames]);
         std::cout<<"reco name: "<<kinematicvariables[ivars][ikeynames]<<std::endl;
         twofinish = true;
//         if(counter ==0){reconmuons = (TH1D*)hist4->Clone("NRecoMuons"); reconmuons->Reset();counter++;}
//           reconmuons->Add(hist4);
         }
       }
   if(onefinish && twofinish){
   std::cout<<"hist4: "<<hist4->GetNbinsX()<<std::endl;
   std::cout<<"hist4 content: "<<hist4->GetMaximumBin()<<" content: "<<hist4->GetBinContent(2)<<std::endl;

//   c0->cd(0);
   hist4->Divide(hist3);
   std::cout<<"hist4 content: "<<hist4->GetMaximumBin()<<" content: "<<hist4->GetBinContent(2)<<std::endl;
   double efficiencyval = hist4->GetBinContent(2);
   std::cout<<"efficiencyval: "<<efficiencyval<<std::endl;
   effvsenergy->SetBinContent(energyval+1, efficiencyval);
   onefinish = false;
   twofinish = false; 
   energyval++;
   }
   }
   c0->cd(1);
   effvsenergy->GetXaxis()->SetTitle("True Neutrino Energy /GeV");
   effvsenergy->GetYaxis()->SetTitle("Efficiency for reconstructing muons");
   effvsenergy->Draw("HIST");
   }
   }
/*
//   leg2->Draw("SAME");
   //c0->Update();
   c0->cd(2);
   truestack->Draw("HIST");
   leg3->Draw("SAME");
   //c0->Update();
   c0->cd(3);
   hist4->Divide(hist3);
   hist4->SetTitle("Number of Reco Muons with Reco Cuts/Number of True Muons with truth cuts");
   hist4->Draw("HIST");*/
   c0->Update();
   c0->Print("efficiencyenergyvars.ps");
   std::cout<<"HERE"<<std::endl;
   c0->Print("efficiencyenergyvars.ps]");

}
