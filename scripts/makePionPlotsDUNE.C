#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>

void makePionPlotsDUNE(TString inputfile, bool gif, bool allnpions, int npions = 0)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
  TCanvas* c0 = new TCanvas("c0","c0",0,0,700,900);
  c0->Divide(1,2);
  TCanvas* c1 = new TCanvas("c1","c1",0,0,600,600);
//  c1->Divide(1,2);
  c0->Print("pionvars.ps[");
  int ipionmax;
  int yaxismax;
  int yaxismin;
  
  if(allnpions){ipionmax = 4; npions = 0;}
  else{ipionmax = npions+1;}
  std::cout<<"ipionmax: "<<ipionmax<<" npions: "<<npions<<std::endl; 
  for(int ipion = npions; ipion<ipionmax; ipion++){
   for(int ithreshold = 0; ithreshold<710; ithreshold=ithreshold+10){
   if(ipion == 0 && ithreshold > 0){continue;}
   if((ithreshold > 150) && (ithreshold % 50 != 0)){continue;}
   std::string filenamestart(inputfile);
   std::string fullfile = filenamestart+"_"+std::to_string(ipion)+"_pions_"+std::to_string(ithreshold)+"_MeV.root";
   TFile* file = new TFile(fullfile.c_str());
   TList* list = file->GetListOfKeys();
   

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   std::vector<std::vector<TString>> kinematicvariables;

//   std::vector<std::string> plotvariables = {"NMuonsRecoOverTruth"};
//  std::vector<std::string> plotvariables = {"RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "RecoXPos", "TrueYPos", "RecoYPos", "TrueZPos", "RecoZPos", "TrueRad", "RecoRad", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy"};
//  std::vector<std::string> plotvariables = {"IdealNeutrinoRecoEnergy", "TrueNeutrinoEnergy", "TrueMinusIdealRecoEnergy", "TrueMinusIdealRecoEnergyRatio", "RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "ChargedPionMultiplicity", "NRecoPions", "NRecoParticles", "InFDV", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy", "LepPT", "LepPZ", "LepRecoPT", "LepRecoPZ", "PiRecoEnergy", "PiTrueEnergy", "PiRecoMomentum", "PiTrueMomentum", "MuonPiRecoAngle", "MuonPiAngle", "PiZRecoAngle", "PiZAngle"};

//  std::vector<std::string> plotvariables = {"TrueMinusIdealRecoEnergy", "TrueMinusIdealRecoEnergyRatio", "TrueQ2", "TrueW"};
  std::vector<std::string> plotvariables = {"TrueMinusIdealRecoEnergy"};
//  std::vector<std::string> plotvariables = {"IdealNeutrinoRecoEnergy", "RecoNeutrinoEnergy"};


//{"RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "NRecoPions", "NRecoParticles", "InFDV", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy", "LepPT", "LepPZ", "LepRecoPT", "LepRecoPZ", "PiRecoEnergy", "PiTrueEnergy", "PiRecoMomentum", "PiTrueMomentum", "MuonPiRecoAngle", "MuonPiAngle", "PiZRecoAngle", "PiZAngle"};

  std::vector<std::string> x_axis;
//   std::vector<std::string> plotvariables = {"TrueNeutrinoEnergy", "RecoNeutrinoEnergy", "TrueMinusRecoEnergy", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "TrueYPos", "TrueZPos", "NMuons", "TrueMinusRecoEnergyRatio", "RecoLepEnergy"};
   int n_vars = plotvariables.size();
   for(int iplots =0; iplots<n_vars; iplots++){
      std::vector<TString> keynames;
      kinematicvariables.push_back(keynames);
      std::string xaxisname= plotvariables[iplots];
      if(plotvariables[iplots].find("TrueMinus")!= std::string::npos){
      if(plotvariables[iplots].find("Ratio")!= std::string::npos){xaxisname = "(E_{#nu} - E_{rec})/E_{#nu}";}
      else{xaxisname = "E_{#nu} - E_{rec} (GeV)";}
      }
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
   std::vector<int> modecolours = {900-1, 616-6,  kGreen+2, 860-5, 432-7, 880-1};
   std::string histnamefull = "";
   for(int ivars=0; ivars<kinematicvariables.size(); ivars++){
     THStack* stack = new THStack("stack","");
     TLegend* leg1 = new TLegend(0.65,0.55,0.8,0.8);
     std::vector<TH1D*> Hists(modebreakdown.size());
     for(int ikeynames = 0; ikeynames<kinematicvariables[ivars].size(); ikeynames++){
       int index = kinematicvariables[ivars][ikeynames].Last('_');
       int length = kinematicvariables[ivars][ikeynames].Length();
       TString fullname = kinematicvariables[ivars][ikeynames];
       TString name  = fullname.Remove(index, length);
       std::string histnamestart(name);
       histnamefull = histnamestart+"_"+std::to_string(ipion)+"_pions_"+std::to_string(ithreshold)+"_MeV";
       std::cout<<"name: "<<name<<" index: "<<index<<" length: "<<length<<std::endl;
       stack->SetTitle(histnamefull.c_str());
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
//         TH1D* histrebin = (TH1D*)(histogram->Rebin(2, "histrebin"));
//         stack->Add(histrebin);
//         leg1->AddEntry(histrebin,modenames[ihists].c_str()); 
         stack->Add(histogram);
         leg1->AddEntry(histogram,modenames[ihists].c_str()); 
       }
     }
     if(ithreshold == 0){yaxismax = stack->GetMaximum()+100;}
     c0->cd(1);
   //  if(ithreshold == 0){stack->Draw("HIST");}
  //   else{ stack->Draw("HIST SAME");}
     stack->Draw("HIST");
     std::cout<<"yaxismax = "<<yaxismax<<std::endl;
//     stack->GetYaxis()->SetRangeUser(0, yaxismax);
     stack->SetTitle(histnamefull.c_str());
     stack->SetMaximum(yaxismax);
     if(plotvariables[0] == "TrueMinusIdealRecoEnergy"){
       stack->GetXaxis()->SetRangeUser(-0.2, 1.0);
     }
//     gPad->RedrawAxis();
//     gPad->GetTickx();
//     gPad->SetTicks(0,0);
//     gPad->GetTicky();
//     gPad->SetTicks(0,0);
     std::cout<<"yaxismax = "<<yaxismax<<std::endl;
     stack->GetXaxis()->SetTitle(x_axis[ivars].c_str());
     leg1->Draw("SAME"); 
//     gPad->Modified();
//     gPad->Update();
     c0->Update();
     c0->Print("pionvars.ps");
     if(gif){
     c1->cd();
     stack->Draw("HIST");
     std::cout<<"yaxismax = "<<yaxismax<<std::endl;
     stack->SetTitle(histnamefull.c_str());
     stack->SetMaximum(yaxismax);
     if(plotvariables[0] == "TrueMinusIdealRecoEnergy"){
       std::cout<<"here:"<<std::endl;
       stack->GetXaxis()->SetRangeUser(-0.2, 1.0);
     }
     std::cout<<"yaxismax = "<<yaxismax<<std::endl;
     stack->GetXaxis()->SetTitle(x_axis[ivars].c_str());
     leg1->Draw("SAME");     
     c1->Update();
     c1->Print("pion.gif+25");}
     delete stack;
     delete leg1;
     Hists.clear();
   }
   }
   }
   std::cout<<"HERE"<<std::endl;
   c0->Print("pionvars.ps]");

}
