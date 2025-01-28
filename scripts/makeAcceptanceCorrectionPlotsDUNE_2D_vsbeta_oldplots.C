#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include <iostream>

void makeAcceptanceCorrectionPlotsDUNE_2D_vsbeta_oldplots(TString inputfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlueRedYellow);
  gStyle->SetPaintTextFormat("1.2f");
  TCanvas* c0 = new TCanvas("c0","c0",0,0,600,600);
  c0->SetMargin(0.15, 0.15, 0.10, 0.10);
//  c0->Divide(1,2);
  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,800);
//  c1->Divide(1,2);
  c0->Print("acceptancecorrectionvars_2d.ps[");
  int ipionmax;
  int yaxismax;
  int yaxismin;
  int ncd =0;
//  if(allnpions){ipionmax = 4; npions = 0;}
//  else{ipionmax = npions+1;}
/*  std::cout<<"ipionmax: "<<ipionmax<<" npions: "<<npions<<std::endl; 
  for(int ipion = npions; ipion<ipionmax; ipion++){
   for(int ithreshold = 0; ithreshold<710; ithreshold=ithreshold+10){
   if(ipion == 0 && ithreshold > 0){continue;}
   if((ithreshold > 150) && (ithreshold % 50 != 0)){continue;}
*/
//  for(int i_radius = 170; i_radius<240; i_radius=i_radius+20){
//   for(int i_pi0eff = 50; i_pi0eff<110; i_pi0eff=i_pi0eff+10){ 
//   for(int i_gammeff = 50; i_gammaeff<110; i_gammaeff=i_gammaeff+10){
   std::vector<double> beta = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
   std::vector<double> Enu = {0.00, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 4.00, 4.50, 5.00};

  std::vector<TString> fidrad = {"160_0", "170_0", "180_0", "190_0", "200_0", "210_0", "220_0", "230_0", "240_0", "249_45"};
  std::vector<TString> ecalcontainment = {"noecalcontainment", "withecalcontainment"};

  std::vector<TH1D*> histvector;
   std::vector<TString> histnamesfracenu;
  for(int i_fidrad =0; i_fidrad<fidrad.size(); i_fidrad++){
  for(int i_ecal =0; i_ecal<ecalcontainment.size(); i_ecal++){
  TH2D* histfracbeta = new TH2D("histfracbeta", "histfracbeta", Enu.size()-1, &Enu[0], beta.size()-1, &beta[0]);

  TH1D* histfracenu = new TH1D("histfracenu", "histfracenu", Enu.size()-1, &Enu[0]);

  for(int i_enu = 0; i_enu<Enu.size()-1; i_enu++){
   TString fidradpath = fidrad[i_fidrad].ReplaceAll(".", "_");
//   TString fidradpath = filepathfidrad.c_str();
//   fidradpath.ReplaceAll(".", "_");
   TString filepath = "fidrad"+fidradpath+"_Bfield0_5_"+ecalcontainment[i_ecal]+"/";
   std::cout<<"filepath: "<<filepath<<std::endl;
   std::vector<TString> inputfilelist;
   std::string enubin = "Enubin"+std::to_string(i_enu);
   TString inputfilename = inputfile+"_fidrad"+fidradpath+"_Bfield0_5_"+enubin+"_"+ecalcontainment[i_ecal];
//   inputfilename.ReplaceAll("Enubin0", enubin.c_str());
   inputfilelist.push_back(filepath+inputfilename);
   std::cout<<"inputfile: "<<inputfilename<<std::endl;
   inputfilename.ReplaceAll("AcceptedEvents", "AllEvents");
   inputfilelist.push_back(filepath+inputfilename);
   std::vector<TH2D*> histlist;
   std::vector<std::string> histnames;
   std::vector<std::string> x_axis;
   int histlistsize;
   bool highestpt = true;
   std::vector<std::string> plotvariables = {"TrueQ3", "TrueQ0"};

   histlistsize = plotvariables.size() + plotvariables.size()/2;
   for(int i_files =0; i_files<inputfilelist.size(); i_files++){
   std::string filenamestart(inputfilelist[i_files]);
   std::string fullfile = filenamestart+".root";
   TFile* file = new TFile(fullfile.c_str());
   TList* list = file->GetListOfKeys();
   if(filenamestart.find("NotAccepted")!= std::string::npos){highestpt = false;}

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   std::vector<std::vector<TString>> kinematicvariables;

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
      else if(plotvariables[iplots].find("HighestpT")!= std::string::npos && !highestpt){xaxisname.erase(0,9); xaxisname.insert(0, "Rejected");}
      if(plotvariables[iplots].find("Momentum")!= std::string::npos){xaxisname += " (GeV)";}
      if(plotvariables[iplots].find("Length")!= std::string::npos){xaxisname += " (m)";}
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
       histnamefull = fullname+", fidrad = "+fidradpath+" B Field = 0.5, "+(std::to_string(Enu[i_enu])).c_str()+"< Enu < "+(std::to_string(Enu[i_enu+1])).c_str()+" "+ecalcontainment[i_ecal];

       TH2D* hist = (TH2D*)file->Get(kinematicvariables[ivars][ikeynames]);
       histlist.push_back(hist);
       histnames.push_back(histnamefull);
//       TH2D* histrebin = (TH2D*)(hist->Rebin2D(1, 1, "histrebin"));
  

     }
   }
   }

   if(inputfilelist.size() == 1){
   for(int i_hists = 0; i_hists<histlistsize; i_hists++){
       std::string finalname = histnames[i_hists];
       int rebinnum = 2;
       if((std::string(histlist[i_hists]->GetName()).find("HighestpTLength") != std::string::npos)){rebinnum = 5;}
       if((std::string(histlist[i_hists]->GetName()).find("TrueSquaredRad") != std::string::npos) || (std::string(histlist[i_hists]->GetName()).find("TrueXPos") != std::string::npos)){rebinnum = 10;}
       TH2D* histrebin = (TH2D*)(histlist[i_hists]->Rebin2D(rebinnum, rebinnum, "histrebin"));
       c0->cd(i_hists);
       histrebin->Draw("COLZ");
       std::cout<<"draw here"<<std::endl;
       histrebin->SetTitle(finalname.c_str());
       histrebin->GetXaxis()->SetTitle(x_axis[i_hists*2].c_str());
       histrebin->GetYaxis()->SetTitle(x_axis[i_hists*2+1].c_str());
       c0->Print("acceptancecorrectionvars_2d.ps");
   }
   }
   else{
   std::cout<<"histlistsize: "<<histlistsize<<" num hists: "<<histlist.size()<<std::endl;
   for(int i_hists = 0; i_hists<histlistsize; i_hists++){
     if(i_hists<histlistsize-plotvariables.size()/2){
       std::cout<<"ihists: "<<i_hists<<" histlistsize-plotvariables.size()/2: "<<histlistsize-plotvariables.size()/2<<std::endl;
       std::string finalname = histnames[i_hists];
       if(i_hists==0){finalname=finalname+" AcceptedEvents";}
       if(i_hists==1){finalname=finalname+" AllEvents";}
       int rebinnum = 2;
       std::cout<<finalname.c_str()<<std::endl;
       if(std::string(histlist[i_hists]->GetName()).find("TrueRad") != std::string::npos){rebinnum = 4;}
       if((std::string(histlist[i_hists]->GetName()).find("TrueSquaredRad") != std::string::npos) || (std::string(histlist[i_hists]->GetName()).find("TrueXPos") != std::string::npos)){rebinnum = 10;}

       TH2D* histrebin = (TH2D*)(histlist[i_hists]->Rebin2D(rebinnum, rebinnum, "histrebin"));
       c0->cd(i_hists);
       histrebin->SetContour(100);
//       histrebin->Draw("COLZ");
       std::cout<<"draw here"<<std::endl;
       histrebin->SetTitle(finalname.c_str());
       std::cout<<"here"<<std::endl;
       histrebin->GetXaxis()->SetTitle(x_axis[i_hists*2].c_str());
       histrebin->GetYaxis()->SetTitle(x_axis[i_hists*2+1].c_str());
//       c0->Print("acceptancecorrectionvars_2d.ps");
     }
     else{
       TH2D* histratios = (TH2D*)histlist[i_hists-plotvariables.size()]->Clone();
       int rebinnum = 1;
       if((std::string(histlist[i_hists-plotvariables.size()]->GetName()).find("TrueRad") != std::string::npos)){rebinnum =4;}
       if((std::string(histlist[i_hists-plotvariables.size()]->GetName()).find("TrueQ0") != std::string::npos)){rebinnum =2;}
       if((std::string(histlist[i_hists-plotvariables.size()]->GetName()).find("TrueXPos") != std::string::npos) || (std::string(histlist[i_hists-plotvariables.size()]->GetName()).find("TrueSquaredRad") != std::string::npos)){rebinnum = 4;}
 //      if(Enu[i_enu]==4.0 || Enu[i_enu]==3.5){rebinnum = 2;}
       if(Enu[i_enu]==4.5){rebinnum = 3;}
       TH2D* histdenom = (TH2D*)(histlist[i_hists-plotvariables.size()/2]->Rebin2D(rebinnum, rebinnum, "histrebin"));
       TH2D* histrebin = (TH2D*)(histratios->Rebin2D(rebinnum, rebinnum, "histrebin"));
       std::cout<<"cloned"<<std::endl;
       histrebin->GetXaxis()->SetTitle(x_axis[(i_hists-plotvariables.size())*2].c_str());
       histrebin->GetYaxis()->SetTitle(x_axis[(i_hists-plotvariables.size())*2+1].c_str());
       TString titleacceptance = "Acceptance vs Q3 vs Q0, fidrad = "+fidradpath.ReplaceAll("_", ".")+", Bfield = 0.5, "+(std::to_string(Enu[i_enu])).c_str()+"    < Enu < "+(std::to_string(Enu[i_enu+1])).c_str()+" "+ecalcontainment[i_ecal];
       histrebin->SetTitle(titleacceptance);
       histrebin->Divide(histdenom);

       for(int i_x =0; i_x<histratios->GetNbinsX(); i_x++){
         for(int i_y = 0; i_y<histratios->GetNbinsY(); i_y++){
           int binnum = histrebin->GetBin(i_x+1, i_y+1);
           double value = histrebin->GetBinContent(binnum);
           double denominator = histdenom->GetBinContent(binnum);
           if(value == 0 && denominator != 0){histrebin->SetBinContent(binnum, 1e-6);}
         }
       }
       for(int i_beta = 0; i_beta<beta.size(); i_beta++){
         int count = 0;
         int total = 0;
         for(int i_x =0; i_x<histrebin->GetNbinsX(); i_x++){
           for(int i_y = 0; i_y<histrebin->GetNbinsY(); i_y++){
             if(histrebin->GetBinContent(i_x+1, i_y+1)!=0){
               total++;
//               std::cout<<"histrebin->GetBinContent(i_x+1, i_y+1): "<<histrebin->GetBinContent(i_x+1, i_y+1)<<std::endl;
               if(histrebin->GetBinContent(i_x+1, i_y+1)<=beta[i_beta+1]){
               count++;
               }
             }
           }
         }
         std::cout<<"count: "<<count<<" total: "<<total<<std::endl;
         double fraction = (double)count/(double)total;
         histfracbeta->SetBinContent(i_enu+1, i_beta+1, fraction);
         if(beta[i_beta+1] == 0.1){histfracenu->SetBinContent(i_enu+1, fraction);};
         std::cout<<"i_enu+1: "<<i_enu+1<<" i_beta+1: "<<i_beta+1<<" fraction: "<<fraction<<"Enu: "<<Enu[i_enu]<<"beta: "<<beta[i_beta+1]<<std::endl;
       }
       std::cout<<"divided"<<std::endl;
       c0->cd(i_hists);
       histrebin->SetContour(100);
//       histrebin->Draw("COLZ");
       histrebin->SetMarkerSize(0.2);
       histrebin->Draw("TEXT SAME");
//       histratios->Draw("COLZ");
//       c0->Print("acceptancecorrectionvars_2d.ps");
       ncd =i_hists;
     }
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
   }
//   }
//   } 
   //}
   ncd++;
   TString title = "fiducial rad = "+fidrad[i_fidrad].ReplaceAll("_", ".")+" cm";
   c0->cd(ncd);
   histfracbeta->SetContour(100);
   histfracbeta->Draw("COLZ0");
   histfracbeta->GetXaxis()->SetTitle("E#nu (GeV)");
   histfracbeta->GetYaxis()->SetTitle("#beta threshold");
   histfracbeta->SetTitle(title);
   c0->Print("acceptancecorrectionvars_2d.ps");

   histvector.push_back(histfracenu);
   histnamesfracenu.push_back(title);
     }
   }

   std::vector<double> xfracline = {0, 2.5, 4.5};
   std::vector<double> yfracline = {0.03,0.03, 0.07};
   TGraph* fracline = new TGraph(xfracline.size(), &xfracline[0], &yfracline[0]); 
   TLegend *legend = new TLegend(0.151,0.7,0.75,0.88);
//   TLegend *legend = new TLegend(0.65,0.7,0.85,0.9);
   c0->cd(ncd+1);
//   histfracenu->SetContour(100);
   histvector[0]->GetYaxis()->SetRangeUser(0.0, 0.15);  
   histvector[0]->SetLineStyle(2);
   histvector[0]->SetLineColor(1);
   histvector[0]->Draw("L");
   histvector[0]->GetXaxis()->SetTitle("E_{#nu} (GeV)");
   histvector[0]->GetYaxis()->SetTitle("Fraction");
   histvector[0]->SetTitle("Fraction of Phase Space with less than 10% Acceptance");
//   legend->AddEntry(histvector[0],histnamesfracenu[0]);
     
   for(int i_histenu = 1; i_histenu<histvector.size(); i_histenu++){
     int linecol = int(i_histenu/2)+1;
     if(linecol == 10){linecol = 12;}
     histvector[i_histenu]->SetLineColor(linecol);
     if(i_histenu%2==0){histvector[i_histenu]->SetLineStyle(2);}
     histvector[i_histenu]->Draw("L SAME");
     if(i_histenu%2==1){
     legend->AddEntry(histvector[i_histenu],histnamesfracenu[i_histenu]);
     }
   }
    TLine *dummy_line = new TLine(0, 0, 1, 1); //Arbitrary points
    dummy_line->SetLineStyle(2); // Dashed
    dummy_line->SetLineColor(kGray);
  // Add the dummy line to the legend
    legend->AddEntry(dummy_line, "Without ECal Containment", "l");

   fracline->SetLineWidth(3);
   fracline->Draw("L SAME");
   legend->SetNColumns(2);
   legend->SetBorderSize(0);
   legend->SetTextSize(0.018);
   legend->SetFillColor(0);
   legend->SetFillStyle(1001);
   legend->Draw();
   c0->Print("acceptancecorrectionvars_2d.ps");
   std::cout<<"HERE"<<std::endl;
   c0->Print("acceptancecorrectionvars_2d.ps]");
}
