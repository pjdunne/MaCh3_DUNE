#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>

void makeSelectionVarsPlotsDUNE(TString inputfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
   TFile* file = new TFile(inputfile);
   TList* list = file->GetListOfKeys();
   

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   TCanvas* c0 = new TCanvas("c0","c0",0,0,700,900);
   c0->Divide(1,2);
   c0->Print("selectionvars.ps[");
   std::vector<std::vector<TString>> kinematicvariables;

   std::vector<std::string> plotvariables = {"NMuonsRecoOverTruth"};
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

   for(int ivars=0; ivars<kinematicvariables.size(); ivars++){
     THStack* stack = new THStack("stack","");
     TLegend* leg1 = new TLegend(0.65,0.55,0.8,0.8);
     std::vector<TH1D*> Hists(modebreakdown.size());
     for(int ikeynames = 0; ikeynames<kinematicvariables[ivars].size(); ikeynames++){
       int index = kinematicvariables[ivars][ikeynames].Last('_');
       int length = kinematicvariables[ivars][ikeynames].Length();
       TString fullname = kinematicvariables[ivars][ikeynames];
       TString name  = fullname.Remove(index, length);
       std::cout<<"name: "<<name<<" index: "<<index<<" length: "<<length<<std::endl;
       stack->SetTitle(name);
       TH1D* hist =(TH1D*)file->Get(kinematicvariables[ivars][ikeynames]);
       std::cout<<kinematicvariables[ivars][ikeynames]<<" mode: "<<fullname.Remove(0, 10)<<std::endl;
       fullname = kinematicvariables[ivars][ikeynames];
       int mode = fullname.Remove(0, index+1).Atoi();
       std::cout<<"mode: "<<mode<<std::endl;
         for(int modes =0; modes<modebreakdown.size(); modes++){
         if(std::find(modebreakdown[modes].begin(), modebreakdown[modes].end(), mode)!=modebreakdown[modes].end()){
           if(Hists[modes] ==NULL){ Hists[modes] = hist;}
           else{Hists[modes]->Add(hist);}
         }
       }
     }
     for(int ihists=0; ihists<Hists.size(); ihists++){
       TH1D* histogram = Hists[ihists];
       if(histogram != NULL){
         histogram->SetFillColor(modecolours[ihists]);
         histogram->SetLineColor(modecolours[ihists]);
         stack->Add(histogram);
         leg1->AddEntry(histogram,modenames[ihists].c_str()); 
       }
     }
     c0->cd(1);
     stack->Draw("HIST");
     leg1->Draw("SAME"); 
     c0->Update();
     c0->Print("selectionvars.ps");
     delete stack;
     delete leg1;
     Hists.clear();
   }
   std::cout<<"HERE"<<std::endl;
   c0->Print("selectionvars.ps]");

}
