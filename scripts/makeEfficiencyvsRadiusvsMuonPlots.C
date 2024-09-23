#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>
//input name = first part of the file name 
//efficiencybool = 1 for efficiency and 0 for purity
//radius = 1 for plot against reco rad and 0 for plot against true neutrino energy

void makeEfficiencyvsRadiusvsMuonPlots(TString inputname, bool efficiencybool, bool onedhist, bool cumulative, bool acceptance, TString plotvariable, TString plotvariable2 = "",TString plotvariable3="")
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
  std::vector<TH1D*> onedhistlist;
  std::vector<std::string> plotvariables; 
  std::string xaxisname;
  std::vector<std::string> xaxis;
  std::vector<bool> efficiencyboolsvec;
  std::cout<<"onedhist: "<<onedhist<<std::endl;
  if(onedhist){
   TFile* file = new TFile(inputname);
   TList* list = file->GetListOfKeys();

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   std::vector<std::vector<TString>> kinematicvariables;

  std::string plotvariablestr(plotvariable);
  plotvariables.push_back(plotvariablestr);
  if(plotvariable2!=""){
    std::string plotvariable2str(plotvariable2);
    plotvariables.push_back(plotvariable2str);
  }


//  std::vector<std::string> plotvariables = {"RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "RecoXPos", "TrueYPos", "RecoYPos", "TrueZPos", "RecoZPos", "TrueRad", "RecoRad", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy"};

  std::vector<std::string> x_axis;
//   std::vector<std::string> plotvariables = {"TrueNeutrinoEnergy", "RecoNeutrinoEnergy", "TrueMinusRecoEnergy", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "TrueYPos", "TrueZPos", "NMuons", "TrueMinusRecoEnergyRatio", "RecoLepEnergy"};
   int n_vars = plotvariables.size();
   for(int iplots =0; iplots<n_vars; iplots++){
      std::vector<TString> keynames;
      kinematicvariables.push_back(keynames);
      xaxisname= plotvariables[iplots];
      if(plotvariables[iplots].find("Ratio")!= std::string::npos){xaxisname += "";}
      else if(plotvariables[iplots].find("Energy")!= std::string::npos){xaxisname += " (GeV)";}
      else if(plotvariables[iplots].find("Pos")!= std::string::npos){xaxisname += " (cm)";}
      else if(plotvariables[iplots].find("Rad")!= std::string::npos){xaxisname += " (cm)";}
      else if(plotvariables[iplots].find("Muons")!= std::string::npos){xaxisname += "";}
      xaxis.push_back(xaxisname);
      std::cout<<"pushed back: "<<xaxis.size()<<" "<<plotvariables[iplots]<<std::endl;
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

   if(acceptance){
     efficiencyboolsvec.push_back(0);
     efficiencyboolsvec.push_back(1);
   }
   else if(!acceptance){
     efficiencyboolsvec.push_back(efficiencybool);
   }
   for(int ihists = 0; ihists<efficiencyboolsvec.size(); ihists++){

     std::string denomname=plotvariablestr;
     std::string trueselecname = plotvariables[0]+"_NDGAr_all_TrueSelected";
     if(efficiencyboolsvec[ihists]){
       denomname = plotvariables.back()+"_NDGAr_all_AllTrue";
     }
     if(!efficiencyboolsvec[ihists]){
       denomname = plotvariables.back()+"_NDGAr_all_AllSelected";
     }

     TH1D* trueselectedfile = (TH1D*)file->Get(trueselecname.c_str());
     TH1D* denominatorfile = (TH1D*)file->Get(denomname.c_str());
     TH1D* trueselected = (TH1D*)(trueselectedfile);
     TH1D* denominator = (TH1D*)(denominatorfile);
     
     if(plotvariable == "TrueNeutrinoEnergy"){
     trueselected = (TH1D*)(trueselectedfile->Rebin(10, "trueselected"));
     denominator = (TH1D*)(denominatorfile->Rebin(10, "denominator"));
     }
     TH1D* efficiency = (TH1D*)trueselected->Clone();
     TGraph* effvsrad2 = new TGraph();
     TH1D* efficiency_binbybin = (TH1D*)trueselected->Clone();
     efficiency->Reset();
     efficiency_binbybin->Reset();

     int nbins = trueselected->GetNbinsX();
     double truselec =0;
     double alldenom = 0;
     double result, rad, absolute_eff, error, error2, error3, eff_per5cm, trusel_5, alltru_5;
     for(int bin =1; bin<nbins+1; bin++){
       error = pow(denominator->GetBinContent(bin), 0.5);
       error2 = pow(trueselected->GetBinContent(bin), 0.5);
       trueselected->SetBinError(bin,error2);
       denominator->SetBinError(bin,error);
       if(error !=0 && error2 !=0){
         error3 = pow((pow(1/error, 2)+(pow(1/error2, 2))), 0.5);
       }
       else{error3 = 0;}
     
       if(cumulative){
         truselec = truselec + trueselected->GetBinContent(bin);
         alldenom = alldenom + denominator->GetBinContent(bin);
       }
       else if(!cumulative){
         truselec = trueselected->GetBinContent(bin);
         alldenom = denominator->GetBinContent(bin);
       }
       result = truselec/alldenom;
       std::cout<<"truselec: "<<truselec<<" denom: "<<alldenom<<" result: "<<result<<std::endl;
       if(std::isnan(result)){result =0;}
       efficiency->SetBinContent(bin, result);
       efficiency->SetBinError(bin, error3);
       }
       std::cout<<"HERE NK"<<std::endl; 
       onedhistlist.push_back(efficiency);
     }
  }
  else if(!onedhist){
  if(acceptance){
     efficiencyboolsvec.push_back(0);
     efficiencyboolsvec.push_back(1);
   }
   else if(!acceptance){
     efficiencyboolsvec.push_back(efficiencybool);
   }
 
  for(int imuon =0; imuon < 10; imuon++){
   std::string inputnamestr(inputname);
   std::string inputfile;
   inputfile = inputnamestr + "_"+imuon+".root";
   TFile* file = new TFile(inputfile.c_str());
   TList* list = file->GetListOfKeys();
   

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   std::vector<std::vector<TString>> kinematicvariables;

  std::string plotvariablestr(plotvariable);
  std::string plotvariable2str(plotvariable2);
  plotvariables.push_back(plotvariablestr);
  plotvariables.push_back(plotvariable2str);

//  std::vector<std::string> plotvariables = {"RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "RecoXPos", "TrueYPos", "RecoYPos", "TrueZPos", "RecoZPos", "TrueRad", "RecoRad", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy"};

  std::vector<std::string> x_axis;
//   std::vector<std::string> plotvariables = {"TrueNeutrinoEnergy", "RecoNeutrinoEnergy", "TrueMinusRecoEnergy", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "TrueYPos", "TrueZPos", "NMuons", "TrueMinusRecoEnergyRatio", "RecoLepEnergy"};
   int n_vars = plotvariables.size();
   for(int iplots =0; iplots<n_vars; iplots++){
      std::vector<TString> keynames;
      kinematicvariables.push_back(keynames);
      xaxisname= plotvariables[iplots];
      if(plotvariables[iplots].find("Ratio")!= std::string::npos){xaxisname += "";}
      else if(plotvariables[iplots].find("Energy")!= std::string::npos){xaxisname += " (GeV)";}
      else if(plotvariables[iplots].find("Pos")!= std::string::npos){xaxisname += " (cm)";}
      else if(plotvariables[iplots].find("Rad")!= std::string::npos){xaxisname += " (cm)";}
      else if(plotvariables[iplots].find("Muons")!= std::string::npos){xaxisname += "";}
      xaxis.push_back(xaxisname);
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
  
   for(int ihists = 0; ihists<efficiencyboolsvec.size(); ihists++){

   std::string denomname=plotvariablestr;
   std::string trueselecname = plotvariables[0]+"_NDGAr_all_TrueSelected";

   for(int nplotvars =0; nplotvars<plotvariables.size(); nplotvars++){ std::cout<<"plot vars: "<<plotvariables[nplotvars]<<std::endl;}
   if(efficiencyboolsvec[ihists]){
     denomname = plotvariables.back()+"_NDGAr_all_AllTrue";
   }
   if(!efficiencyboolsvec[ihists]){
     denomname = plotvariables[0]+"_NDGAr_all_AllSelected";
   }
  
   std::cout<<"denom name: "<<denomname<<std::endl;
   TH1D* trueselectedfile = (TH1D*)file->Get(trueselecname.c_str());
   TH1D* denominatorfile = (TH1D*)file->Get(denomname.c_str());
   std::cout<<"trueselectedfile: "<<trueselectedfile->GetNbinsX()<<std::endl;
   TH1D* trueselected = (TH1D*)(trueselectedfile->Rebin(3, "trueselected"));
   TH1D* denominator = (TH1D*)(denominatorfile->Rebin(3, "denominator"));
   std::cout<<"trueselected: "<<trueselected->GetNbinsX()<<std::endl;
//   TH1D* trueselected = (TH1D*)file->Get(trueselecname.c_str());
//   TH1D* denominator = (TH1D*)file->Get(denomname.c_str());
   TH1D* efficiency = (TH1D*)trueselected->Clone();
   TGraph* effvsrad2 = new TGraph();
   TH1D* efficiency_binbybin = (TH1D*)trueselected->Clone();
   efficiency->Reset();
   efficiency_binbybin->Reset();

   int nbins = trueselected->GetNbinsX();
   double truselec =0;
   double alldenom = 0;
   double result, rad, absolute_eff, error, error2, error3, eff_per5cm, trusel_5, alltru_5;
   for(int bin =1; bin<nbins+1; bin++){
     std::cout<<"bin: "<<bin<<std::endl;
     error = pow(denominator->GetBinContent(bin), 0.5);
     error2 = pow(trueselected->GetBinContent(bin), 0.5);
     trueselected->SetBinError(bin,error2);
     denominator->SetBinError(bin,error);
     if(error !=0 && error2 !=0){
     error3 = pow((pow(1/error, 2)+(pow(1/error2, 2))), 0.5);
     }
     else{error3 = 0;}
     
     if(cumulative){
       truselec = truselec + trueselected->GetBinContent(bin);
       alldenom = alldenom + denominator->GetBinContent(bin);
     }
     else if(!cumulative){
       truselec = trueselected->GetBinContent(bin);
       alldenom = denominator->GetBinContent(bin);
     }
     result = truselec/alldenom;
     std::cout<<"truselec: "<<truselec<<" denom: "<<alldenom<<" result: "<<result<<std::endl;
     if(std::isnan(result)){result =0;}
     efficiency->SetBinContent(bin, result);
     efficiency->SetBinError(bin, error3);
     }
   onedhistlist.push_back(efficiency);
   }
   }
   }
   if(!onedhist){
   std::cout<<"number muon hists"<< onedhistlist.size()<<std::endl;
   int nxbins = onedhistlist[0]->GetNbinsX();
   std::vector<double> muonscore = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
   int nybins = muonscore.size()-1;
   TH2D* twodhist = new TH2D("twodhist", "", nxbins, onedhistlist[0]->GetXaxis()->GetBinLowEdge(1), onedhistlist[0]->GetXaxis()->GetBinLowEdge(nxbins)+onedhistlist[0]->GetXaxis()->GetBinWidth(nxbins), nybins, &muonscore[0]);
   int nhists = efficiencyboolsvec.size(); 
   std::string outfilename ="";
   if(acceptance){nhists = efficiencyboolsvec.size()+1; outfilename = "acceptancevs"+plotvariables[0]+"vs"+plotvariable3;}
     else if(!acceptance){ 
       if(efficiencybool){
         outfilename = "efficiencyvs"+plotvariables[0]+"vs"+plotvariable3;
       }
       else if(!efficiencybool){
         outfilename = "purityvs"+plotvariables[0]+"vs"+plotvariable3;
       }
     }
     TCanvas* c2 = new TCanvas("c2","c2",0,0,510,510);
     auto legend = new TLegend(0.7,0.7,0.8,0.8);
     c2->SetRightMargin(0.18);
//     c2->Divide(2,2);
     std::string openfile = outfilename+".ps[";
     std::string closefile = outfilename+".ps]";
     std::string finalfile = outfilename+".ps";
     std::cout<<"HERE NK 5"<<std::endl;
     c2->Print(openfile.c_str());
   
   for(int ihists = 0; ihists<nhists; ihists++){ 
   std::cout<<"start loop: "<<std::endl;
   int nstart;
   if(ihists<efficiencyboolsvec.size()){
   if(efficiencyboolsvec[ihists]){
     nstart = 2;
   }
   else if(!efficiencyboolsvec[ihists]){
     nstart = 1;
   }
   for(int imuon =nstart; imuon<nybins*2+1; imuon=imuon+2){
     for(int ienergy = 1; ienergy<nxbins+1; ienergy++){
//      std::cout<<"ienergy: "<<ienergy<<" imuon: "<<imuon<<" content: "<<onedhistlist[imuon-1]->GetBinContent(ienergy)<<std::endl;
         float muonbin;
         if(efficiencyboolsvec[ihists]){
          muonbin= imuon/2;
         }
         else if(!efficiencyboolsvec[ihists]){
          muonbin =(imuon+1)/2;
         } 
         std::cout<<"HERE NK 6"<<" imuon: "<<imuon<<"muonbin: "<<muonbin<<" ienergy: "<<ienergy<<" content: "<<onedhistlist[imuon-1]->GetBinContent(ienergy)<<std::endl;
         twodhist->SetBinContent(ienergy, muonbin, onedhistlist[imuon-1]->GetBinContent(ienergy));
     }
   }
   std::cout<<"HERE NK 8"<<std::endl;
   }
   else if(ihists==efficiencyboolsvec.size()){
   for(int imuon =1; imuon<nybins*2+1; imuon=imuon+2){
     for(int ienergy = 1; ienergy<nxbins+1; ienergy++){
//      std::cout<<"ienergy: "<<ienergy<<" imuon: "<<imuon<<" content: "<<onedhistlist[imuon-1]->GetBinContent(ienergy)<<std::endl;
         int muonbin;
         if(efficiencyboolsvec[ihists]){
          muonbin= imuon/2;
         }
         else if(!efficiencyboolsvec[ihists]){
          muonbin =(imuon+1)/2;
         }              
//         std::cout<<"HERE NK 6 new"<<std::endl;
         twodhist->SetBinContent(ienergy, (imuon+1)/2, (onedhistlist[imuon-1]->GetBinContent(ienergy)*onedhistlist[imuon]->GetBinContent(ienergy)));
     }
   }
   }
   std::cout<<"HERE NK 9"<<std::endl;
   std::string histname="";
   std::string xaxisname = "";
   if(ihists<efficiencyboolsvec.size()){
     if(efficiencyboolsvec[ihists]){
       histname = "Efficiency";
     }
     else if(!efficiencyboolsvec[ihists]){
       histname = "Purity";
     }
   }
   else{histname = "Acceptance";}
 
   std::cout<<xaxis.back()<<std::endl;
   xaxisname = xaxis.back();
   std::cout<<"HERE NK 7"<<std::endl;
   c2->cd(ihists+1);
   std::cout<<"ihists: "<<ihists<<"histname: "<<histname<<std::endl;
   
   
   twodhist->SetTitle(histname.c_str());
   twodhist->GetXaxis()->SetTitle(xaxisname.c_str());
   twodhist->GetYaxis()->SetTitle(plotvariable3);
   twodhist->Draw("COLZ");
   c2->Update();
   //twodhist->Reset();
   c2->Print(finalfile.c_str());

   }
   c2->Print(finalfile.c_str());
   std::cout<<"HERE"<<std::endl;
   c2->Print(closefile.c_str());
   }
   else if(onedhist){
     int nhists = efficiencyboolsvec.size(); 
     std::string outfilename ="";
     if(acceptance){nhists = efficiencyboolsvec.size()+2; outfilename = "acceptancevs"+plotvariables[0];}
     else if(!acceptance){ 
       if(efficiencybool){
         outfilename = "efficiencyvs"+plotvariables[0];
       }
       else if(!efficiencybool){
         outfilename = "purityvs"+plotvariables[0];
       }
     }
     TCanvas* c1 = new TCanvas("c1","c1",0,0,700,900);
     auto legend = new TLegend(0.8,0.8,0.9,0.9);

     c1->Divide(2,2);
     std::string openfile = outfilename+".ps[";
     std::string closefile = outfilename+".ps]";
     std::string finalfile = outfilename+".ps";
     std::cout<<"HERE NK 5"<<std::endl;
     c1->Print(openfile.c_str());
     TH1D* acceptance = (TH1D*)onedhistlist.back()->Clone();
     acceptance->Reset();
     if(onedhistlist.size()>=2){
       acceptance->Multiply(onedhistlist[0], onedhistlist[1]);
     }   
     for(int ihists = 0; ihists<nhists; ihists++){ 
       std::cout<<"HERE NK"<<" onehistlist.size(): "<<onedhistlist.size()<<std::endl;
       std::string histname="";
       std::cout<<"HERE NK 2 new "<<xaxis.size()<<std::endl;
       xaxisname = xaxis.back();
       std::cout<<"HERE NK 2 new"<<std::endl;
       TH1D* acceptance = (TH1D*)onedhistlist[0]->Clone();
       if(ihists<efficiencyboolsvec.size()){
         if(efficiencyboolsvec[ihists]){
           histname = "Efficiency";
         }
         else if(!efficiencyboolsvec[ihists]){
           histname = "Purity";
         }
       }
       else{histname = "Acceptance";}
       std::cout<<"HERE NK 2"<<std::endl;

       std::cout<<"HERE NK 3"<<std::endl;
       std::cout<<"HERE NK 4"<<std::endl;

       std::cout<<"hist name: "<<histname<<std::endl;
       c1->cd(ihists+1);
       if(ihists<efficiencyboolsvec.size()){
         onedhistlist[ihists]->SetTitle(histname.c_str());
         onedhistlist[ihists]->GetXaxis()->SetTitle(xaxisname.c_str());
         onedhistlist[ihists]->Draw("P");
//         c1->Update();
       }
       else if(ihists==efficiencyboolsvec.size()){
//         acceptance = onedhistlist[ihists-2]*onedhistlist[ihists-1];
         acceptance->Multiply(onedhistlist[ihists-2], onedhistlist[ihists-1]);
         acceptance->SetTitle(histname.c_str());
         acceptance->GetXaxis()->SetTitle(xaxisname.c_str());
         acceptance->Draw("P");
//         c1->Update();
       }
       else if(ihists==efficiencyboolsvec.size()+1){
         acceptance->Multiply(onedhistlist[ihists-3], onedhistlist[ihists-2]);
         acceptance->SetMarkerStyle(kPlus);
         acceptance->SetMarkerColor(kBlack);
         acceptance->SetLineColor(kBlack);
         acceptance->GetYaxis()->SetRangeUser(0, 1);
         acceptance->Draw("P");
         onedhistlist[ihists-3]->SetMarkerStyle(kPlus);
         onedhistlist[ihists-3]->SetMarkerColor(kRed);
         onedhistlist[ihists-3]->SetLineColor(kRed);
         onedhistlist[ihists-3]->Draw("P SAME");
         onedhistlist[ihists-2]->SetMarkerStyle(kPlus);
         onedhistlist[ihists-2]->SetMarkerColor(kBlue);
         onedhistlist[ihists-2]->SetLineColor(kBlue);
         onedhistlist[ihists-2]->Draw("P SAME");
         legend->AddEntry(onedhistlist[ihists-3], "Purity", "lp");
         legend->AddEntry(onedhistlist[ihists-2], "Efficiency", "lp");
         legend->AddEntry(acceptance, "Acceptance", "lp");
         legend->Draw();

//         c1->Update();
       }
    }       
   c1->Print(finalfile.c_str());
   std::cout<<"HERE"<<std::endl;
   c1->Print(closefile.c_str()); 
   }
}
