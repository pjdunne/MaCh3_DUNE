#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>

void makeAcceptanceCorrectionPlotsDUNE_2D(TString inputfile, TString inputfile2)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
  TCanvas* c0 = new TCanvas("c0","c0",0,0,600,600);
  c0->Divide(1,2);
  TCanvas* c1 = new TCanvas("c1","c1",0,0,600,600);
//  c1->Divide(1,2);
  c0->Print("acceptancecorrectionvars_2d.ps[");
  int ipionmax;
  int yaxismax;
  int yaxismin;
  
//  if(allnpions){ipionmax = 4; npions = 0;}
//  else{ipionmax = npions+1;}
/*  std::cout<<"ipionmax: "<<ipionmax<<" npions: "<<npions<<std::endl; 
  for(int ipion = npions; ipion<ipionmax; ipion++){
   for(int ithreshold = 0; ithreshold<710; ithreshold=ithreshold+10){
   if(ipion == 0 && ithreshold > 0){continue;}
   if((ithreshold > 150) && (ithreshold % 50 != 0)){continue;}
*/
   std::vector<TString> inputfilelist;
   inputfilelist.push_back(inputfile);
   std::vector<TH2D*> histlist;
   std::vector<std::string> histnames;
   std::vector<std::string> x_axis;

   if(inputfile2){inputfilelist.push_back(inputfile2);}
     
   for(int i_files =0; i_files<inputfilelist.size(); i_files++){
   std::string filenamestart(inputfilelist[i_files]);
   std::string fullfile = filenamestart+".root";
   TFile* file = new TFile(fullfile.c_str());
   TList* list = file->GetListOfKeys();
   

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   std::vector<std::vector<TString>> kinematicvariables;

//   std::vector<std::string> plotvariables = {"NMuonsRecoOverTruth"};
//  std::vector<std::string> plotvariables = {"RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "RecoXPos", "TrueYPos", "RecoYPos", "TrueZPos", "RecoZPos", "TrueRad", "RecoRad", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy"};
//  std::vector<std::string> plotvariables = {"IdealNeutrinoRecoEnergy", "TrueNeutrinoEnergy", "TrueMinusIdealRecoEnergy", "TrueMinusIdealRecoEnergyRatio", "RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "ChargedPionMultiplicity", "NRecoPions", "NRecoParticles", "InFDV", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy", "LepPT", "LepPZ", "LepRecoPT", "LepRecoPZ", "PiRecoEnergy", "PiTrueEnergy", "PiRecoMomentum", "PiTrueMomentum", "MuonPiRecoAngle", "MuonPiAngle", "PiZRecoAngle", "PiZAngle"};

//  std::vector<std::string> plotvariables = {"TrueMinusIdealRecoEnergy", "TrueMinusIdealRecoEnergyRatio", "TrueQ2", "TrueW"};
  std::vector<std::string> plotvariables = {"TrueQ3", "TrueQ0"};
//  std::vector<std::string> plotvariables = {"IdealNeutrinoRecoEnergy", "RecoNeutrinoEnergy"};


//{"RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "NRecoPions", "NRecoParticles", "InFDV", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy", "LepPT", "LepPZ", "LepRecoPT", "LepRecoPZ", "PiRecoEnergy", "PiTrueEnergy", "PiRecoMomentum", "PiTrueMomentum", "MuonPiRecoAngle", "MuonPiAngle", "PiZRecoAngle", "PiZAngle"};

//   std::vector<std::string> plotvariables = {"TrueNeutrinoEnergy", "RecoNeutrinoEnergy", "TrueMinusRecoEnergy", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "TrueYPos", "TrueZPos", "NMuons", "TrueMinusRecoEnergyRatio", "RecoLepEnergy"};
   int n_vars = plotvariables.size();
   for(int iplots =0; iplots<n_vars; iplots++){
      std::vector<TString> keynames;
      kinematicvariables.push_back(keynames);
      std::string xaxisname= plotvariables[iplots];
      if(plotvariables[iplots].find("TrueMinus")!= std::string::npos){
      if(plotvariables[iplots].find("Ratio")!= std::string::npos){xaxisname = "(E_nu - E_rec)/E_nu";}
      else{xaxisname = "E_nu - E_rec (GeV)";}
      }
      else if(plotvariables[iplots].find("Energy")!= std::string::npos){xaxisname += " (GeV)";}
      else if(plotvariables[iplots].find("Pos")!= std::string::npos){xaxisname += " (cm)";}
      else if(plotvariables[iplots].find("Rad")!= std::string::npos){xaxisname += " (cm)";}
      else if(plotvariables[iplots].find("Muons")!= std::string::npos){xaxisname += "";}
      else if(plotvariables[iplots].find("TrueQ2")!= std::string::npos){xaxisname += " (GeV^{2})";}
      else if(plotvariables[iplots].find("TrueQ0")!= std::string::npos){xaxisname += " (GeV)";}
      else if(plotvariables[iplots].find("TrueQ3")!= std::string::npos){xaxisname += " (GeV)";}
      else if(plotvariables[iplots].find("TrueW")!= std::string::npos){xaxisname += " (GeV^{2})";}
      x_axis.push_back(xaxisname);
   }

  for(int i=0; i<list->GetEntries(); i++)
   {
     bool boolcontinue = false;
     //std::cout << "entry " << i << std::endl;
     TString keyname = list->At(i)->GetName();
        if(keyname.Contains("_NDGAr_")){
        std::cout<<"here ndgar"<<std::endl;
        for(int iplots=0; iplots<n_vars; iplots = iplots+2){
          TString keynametest = plotvariables[iplots]+"_"+plotvariables[iplots+1]+"_NDGAr_";
          std::cout<<"keynametest: "<<keynametest<<std::endl;
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
     for(int ikeynames = 0; ikeynames<kinematicvariables[ivars].size(); ikeynames++){
     std::cout<<"here2 "<<kinematicvariables[ivars][ikeynames]<<std::endl;
       int index = kinematicvariables[ivars][ikeynames].Last('_');
       int length = kinematicvariables[ivars][ikeynames].Length();
       TString fullname = kinematicvariables[ivars][ikeynames];
       TString name  = fullname.Remove(index, length);
//       std::string histnamestart(name);
       histnamefull = fullname;

       TH2D* hist = (TH2D*)file->Get(kinematicvariables[ivars][ikeynames]);
       histlist.push_back(hist);
       histnames.push_back(histnamefull);
//       TH2D* histrebin = (TH2D*)(hist->Rebin2D(1, 1, "histrebin"));
  

     }
   }
   }
   for(int i_hists = 0; i_hists<histlist.size()+1; i_hists++){
     if(i_hists<2){
       std::string finalname = histnames[i_hists];
       if(i_hists==0){finalname=finalname+"_AcceptedEvents";}
       if(i_hists==1){finalname=finalname+"_AllEvents";}
       c0->cd(i_hists);
       histlist[i_hists]->Draw("COLZ");
       std::cout<<"draw here"<<std::endl;
       histlist[i_hists]->SetTitle(finalname.c_str());
       histlist[i_hists]->GetXaxis()->SetTitle(x_axis[0].c_str());
       histlist[i_hists]->GetYaxis()->SetTitle(x_axis[1].c_str());
       c0->Print("acceptancecorrectionvars_2d.ps");
     }
     if(i_hists==2){
       TH2D* histratios = (TH2D*)histlist[0]->Clone();
       std::cout<<"cloned"<<std::endl;
       histratios->GetXaxis()->SetTitle(x_axis[0].c_str());
       histratios->GetYaxis()->SetTitle(x_axis[1].c_str());
       histratios->SetTitle("Acceptance vs Q3 vs Q0");
       histratios->Divide(histlist[1]);
       for(int i_x =0; i_x<histratios->GetNbinsX(); i_x++){
         for(int i_y = 0; i_y<histratios->GetNbinsY(); i_y++){
           int binnum = histratios->GetBin(i_x+1, i_y+1);
           double numerator = histlist[0]->GetBinContent(binnum);
           double denominator = histlist[1]->GetBinContent(binnum);
           if(numerator == 0 && denominator != 0){histratios->SetBinContent(binnum, 1e-6);}
         }
       }
       std::cout<<"divided"<<std::endl;
       c0->cd(i_hists);
       histratios->Draw("COLZ");
       c0->Print("acceptancecorrectionvars_2d.ps");
     }
/*       c0->cd(1);
       hist->Draw("COLZ");
       hist->SetTitle(histnamefull.c_str());
       hist->GetXaxis()->SetTitle(x_axis[0].c_str());
       hist->GetYaxis()->SetTitle(x_axis[1].c_str());
//     histrebin->SetMaximum(yaxismax);
//     c0->Update();
       c0->Print("acceptancecorrectionvars_2d.ps");
*/
     }
   
   //}
   std::cout<<"HERE"<<std::endl;
   c0->Print("acceptancecorrectionvars_2d.ps]");

}
