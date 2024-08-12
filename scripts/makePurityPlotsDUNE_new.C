#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>

void makePurityPlotsDUNE_new(TString inputfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
   TFile* file = new TFile(inputfile);
   TList* list = file->GetListOfKeys();
   

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   TCanvas* c0 = new TCanvas("c0","c0",0,0,700,900);
   c0->Divide(2,2);
   c0->Print("purityvars.ps[");
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
  /*
   for(int ivars=0; ivars<kinematicvariables.size(); ivars++){
     THStack* stack = new THStack("stack","");
     TLegend* leg1 = new TLegend(0.65,0.55,0.8,0.8);
     std::vector<TH1D*> HistsReco(modebreakdown.size());
     std::vector<TH1D*> HistsTruth(modebreakdown.size());
     for(int ikeynames = 0; ikeynames<kinematicvariables[ivars].size(); ikeynames++){
       if(kinematicvariables[ivars][ikeynames].Contains("_RECO")){
         int index = kinematicvariables[ivars][ikeynames].First('_');
         int index2 = kinematicvariables[ivars][ikeynames].Last('_');
         int length = kinematicvariables[ivars][ikeynames].Length();
         TString fullname = kinematicvariables[ivars][ikeynames];
         TString name  = fullname.Remove(index+6, length)+"_RecoCuts";
         std::cout<<"name: "<<name<<" index: "<<index<<" length: "<<length<<std::endl;
         stack->SetTitle(name);
         std::cout<<"name: "<<kinematicvariables[ivars][ikeynames]<<std::endl;
         TH1D* hist =(TH1D*)file->Get(kinematicvariables[ivars][ikeynames]);
         //std::cout<<kinematicvariables[ivars][ikeynames]<<" mode: "<<fullname.Remove(0, 10)<<std::endl;
         fullname = kinematicvariables[ivars][ikeynames];
         fullname.Remove(index2, length);
         int mode = fullname.Remove(0, index+7).Atoi();
         std::cout<<"mode: "<<mode<<std::endl;
           for(int modes =0; modes<modebreakdown.size(); modes++){
           if(std::find(modebreakdown[modes].begin(), modebreakdown[modes].end(), mode)!=modebreakdown[modes].end()){
             if(HistsReco[modes] ==NULL){ HistsReco[modes] = hist;}
             else{HistsReco[modes]->Add(hist);}
           }
         }
       }
     }
     for(int ihists=0; ihists<HistsReco.size(); ihists++){
       TH1D* histogram = HistsReco[ihists];
       if(histogram != NULL){
         histogram->SetFillColor(modecolours[ihists]);
         histogram->SetLineColor(modecolours[ihists]);
         stack->Add(histogram);
         leg1->AddEntry(histogram,modenames[ihists].c_str());
         if(kinematicvariables[ivars][0].Contains("NRecoMuons")){recostack->Add(histogram); leg2->AddEntry(histogram,modenames[ihists].c_str());}
       }
     }
     c0->cd(1);
     stack->Draw("HIST");
     leg1->Draw("SAME"); 
//     c0->Update();
//     c0->Print("efficiencyvars.ps");
//     delete stack;
//     leg1->Clear();
//     Hists.clear();

     THStack* stack2 = new THStack("stack2","");
     TLegend* leg0 = new TLegend(0.65,0.55,0.8,0.8);
     for(int ikeynames = 0; ikeynames<kinematicvariables[ivars].size(); ikeynames++){
       if(kinematicvariables[ivars][ikeynames].Contains("_TRUTH")){
         int index = kinematicvariables[ivars][ikeynames].First('_');
         int index2 = kinematicvariables[ivars][ikeynames].Last('_');
         int length = kinematicvariables[ivars][ikeynames].Length();
         TString fullname = kinematicvariables[ivars][ikeynames];
         TString name  = fullname.Remove(index+6, length)+"_TruthCuts";
         std::cout<<"name: "<<name<<" index: "<<index<<" length: "<<length<<std::endl;
         stack2->SetTitle(name);
         TH1D* hist1 =(TH1D*)file->Get(kinematicvariables[ivars][ikeynames]);
         std::cout<<kinematicvariables[ivars][ikeynames]<<" mode: "<<fullname.Remove(0, 10)<<std::endl;
         fullname = kinematicvariables[ivars][ikeynames];
         fullname.Remove(index2, length);
         int mode = fullname.Remove(0, index+7).Atoi();
         std::cout<<"mode: "<<mode<<" name: "<<name<<std::endl;
           for(int modes =0; modes<modebreakdown.size(); modes++){
           if(std::find(modebreakdown[modes].begin(), modebreakdown[modes].end(), mode)!=modebreakdown[modes].end()){
             if(HistsTruth[modes] ==NULL){ HistsTruth[modes] = hist1;}
             else{HistsTruth[modes]->Add(hist1);}
           }
         }
       }
     }
     for(int ihists=0; ihists<HistsTruth.size(); ihists++){
       TH1D* histogram1 = HistsTruth[ihists];
       if(histogram1 != NULL){
         histogram1->SetFillColor(modecolours[ihists]);
         histogram1->SetLineColor(modecolours[ihists]);
         stack2->Add(histogram1);
         leg0->AddEntry(histogram1,modenames[ihists].c_str());
         if(kinematicvariables[ivars][0].Contains("NTrueMuons")){truestack->Add(histogram1); leg3->AddEntry(histogram1,modenames[ihists].c_str());}
       }
     }
     c0->cd(2);
     stack2->Draw("HIST");
     leg0->Draw("SAME"); 
     c0->Update();
     c0->Print("efficiencyvars.ps");
     delete stack;
     delete stack2;
     delete leg1;
     delete leg0;
     HistsReco.clear();
     HistsTruth.clear();
     //}
     
     int counter =0;
     if(kinematicvariables[ivars][0].Contains("NTrueMuons")){
       for(int ikeynames = 0; ikeynames<kinematicvariables[ivars].size(); ikeynames++){
         if(kinematicvariables[ivars][ikeynames].Contains("_TRUTH")){
           TH1D* hist3 =(TH1D*)file->Get(kinematicvariables[ivars][ikeynames]);
           if(counter ==0){truenmuons = (TH1D*)hist3->Clone("NTrueMuons"); truenmuons->Reset();counter++;}
           truenmuons->Add(hist3);
         }
       }
     }
     counter =0;
     if(kinematicvariables[ivars][0].Contains("NRecoMuons")){
       for(int ikeynames = 0; ikeynames<kinematicvariables[ivars].size(); ikeynames++){
         if(kinematicvariables[ivars][ikeynames].Contains("_RECO")){
           TH1D* hist4 =(TH1D*)file->Get(kinematicvariables[ivars][ikeynames]);
           if(counter ==0){reconmuons = (TH1D*)hist4->Clone("NRecoMuons"); reconmuons->Reset();counter++;}
           reconmuons->Add(hist4);
         }
       }
     }
   }
*/
   TH1D* trueselected = (TH1D*)file->Get("TrueNeutrinoEnergy_NDGAr_all_TrueSelected");
   TH1D* allselected = (TH1D*)file->Get("TrueNeutrinoEnergy_NDGAr_all_AllSelected");;
   TH1D* purity = (TH1D*)file->Get("TrueNeutrinoEnergy_NDGAr_all_TrueSelected");

   c0->cd(1);
   allselected->Draw("HIST");
   allselected->Sumw2();
//   leg2->Draw("SAME");
//   c0->Update();
   c0->cd(2);
   trueselected->Draw("HIST");
   trueselected->Sumw2();
//   leg3->Draw("SAME");
   c0->Update();
   c0->Print("purityvars.ps");
   c0->cd(3);
   purity->Divide(trueselected, allselected);
   purity->SetTitle("Purity");
   purity->Sumw2();
//   trueselected->Divide(alltrue);
//   trueselected->SetTitle("Efficiency");
   purity->Draw("E P*");
//   c0->Update();
   c0->Print("purityvars.ps");
   std::cout<<"HERE"<<std::endl;
   c0->Print("purityvars.ps]");

}
