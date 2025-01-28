#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include <iostream>
void rebin_2d_based_on_events(TH2D* h, int min_events, std::vector<double> &new_bins_x, std::vector<double> &new_bins_y) {

    // Rebinning in x
    new_bins_x.push_back(h->GetXaxis()->GetXmin());
/*    double current_bin_content_sum_x = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            current_bin_content_sum_x += h->GetBinContent(i, j);
        }
        if (current_bin_content_sum_x >= min_events) {
            new_bins_x.push_back(h->GetXaxis()->GetBinUpEdge(i));
            std::cout<<"new bin x: "<<h->GetXaxis()->GetBinUpEdge(i)<<std::endl;
            current_bin_content_sum_x = 0;
        }
    }*/
    if (new_bins_x.back() < h->GetXaxis()->GetXmax()){
        new_bins_x.push_back(h->GetXaxis()->GetXmax());
    }

    // Rebinning in y
    new_bins_y.push_back(h->GetYaxis()->GetXmin());
/*    double current_bin_content_sum_y = 0;
    for (int j = 1; j <= h->GetNbinsY(); ++j) {
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            current_bin_content_sum_y += h->GetBinContent(i, j);
        }
        if (current_bin_content_sum_y >= min_events) {
            new_bins_y.push_back(h->GetYaxis()->GetBinUpEdge(j));
            current_bin_content_sum_y = 0;
            std::cout<<"new bin y: "<<h->GetYaxis()->GetBinUpEdge(j)<<std::endl;
        }
    }*/
     for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
          if (h->GetBinContent(i,j) >= min_events){
            std::cout<<" i, j: "<<i<<" , "<<j<<std::endl;
            bool x_exists = false;
            bool y_exists = false;
            for (double x : new_bins_x){
              if (h->GetXaxis()->GetBinUpEdge(i) == x) x_exists = true;
            }
            if (!x_exists) new_bins_x.push_back(h->GetXaxis()->GetBinUpEdge(i));

            for (double y : new_bins_y){
              if (h->GetYaxis()->GetBinUpEdge(j) == y) y_exists = true;
            }
            if (!y_exists) new_bins_y.push_back(h->GetYaxis()->GetBinUpEdge(j));
          }
        }
    }
    if (new_bins_y.back() < h->GetYaxis()->GetXmax()){
        new_bins_y.push_back(h->GetYaxis()->GetXmax());
    }


    if (new_bins_x.size() <= 1 || new_bins_y.size() <= 1) {
        std::cout << "Warning: No rebinning necessary or possible. Returning original histogram." << std::endl;
    }
}

void makeAcceptanceCorrectionPlotsDUNE_2D_vsbeta(TString inputfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlueRedYellow);
  gStyle->SetPaintTextFormat("1.2f");
  TCanvas* c0 = new TCanvas("c0","c0",0,0,600,600);
  c0->SetMargin(0.15, 0.15, 0.10, 0.10);
//  c0->Divide(1,2);
  TCanvas* c1 = new TCanvas("c1","c1",0,0,600,600);
  c1->SetMargin(0.15, 0.15, 0.10, 0.10);
//  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,800);
//  c1->Divide(1,2);
  c0->Print("acceptancecorrectionvars_2d.ps[");
  int ipionmax;
  int yaxismax;
  int yaxismin;
  int ncd =0;
  std::vector<double> new_bins_x;
  std::vector<double> new_bins_y;
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
   std::vector<double> Enuerrors = {0.75, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,0.125, 0.125, 0.25, 0.25, 0.25};
   std::vector<double> deltaR;
//  std::vector<TString> fidrad = {"70.0", "90.0"};
//  std::vector<TString> fidrad = {"30.0", "50.0", "70.0", "90.0", "110.0"};
  std::vector<TString> fidrad = {"90.0"};
  std::vector<TString> samplingfreq = {"20.0"};
//  std::vector<TString> samplingfreq = {"2.0", "20.0"};

//  std::vector<TString> samplingfreq = {"2.0", "10.0", "20.0"};
//  std::vector<TString> spatialres = {"2.5"};
  std::vector<TString> spatialres = {"1.0"};
//  std::vector<TString> spatialres = {"1.0", "2.0", "2.5", "3.0"};
//  std::vector<TString> bfield = {"0.5", "0.4", "0.3", "0.2", "0.1"};
  std::vector<TString> bfield = {"0.5", "0.4", "0.3"};
//  std::vector<TString> bfield = {"0.5"};

//  std::vector<TString> instrumentedrad = {"199.45", "249.45"};  
  std::vector<TString> instrumentedrad = {"249.45"};  
  std::vector<std::vector<TString>> fidradnum;
  std::vector<TString> fidradnumvect;
  for(int i_str0=0; i_str0<instrumentedrad.size(); i_str0++){
    fidradnum.push_back(fidradnumvect);
    for(int i_str=0; i_str<fidrad.size(); i_str++){
       double fidraddouble = std::stod(std::string(instrumentedrad[i_str0])) - std::stod(std::string(fidrad[i_str]));
       std::stringstream ss1;
       ss1 << std::fixed << std::setprecision(2) << fidraddouble;
       fidradnum[i_str0].push_back((ss1.str()).c_str());
    }
  }
//  std::vector<TString> fidrad = {"160_0", "170_0", "180_0", "190_0", "200_0", "210_0", "220_0", "230_0", "240_0", "249_45"};
//  std::vector<TString> ecalcontainment = {"noecalcontainment", "withecalcontainment"};
  std::vector<TString> ecalcontainment = {"withecalcontainment"};

  std::vector<TH1D*> histvector;
  std::vector<TString> histnamesfracenu;
  std::vector<double> bin_i;
  std::vector<double> weight_i;
  std::vector<std::vector<double>> surroundingbins;
  std::vector<std::vector<double>> surroundingbins_weights;

  double threshold = 0.1;
  for(int i_bfield = 0; i_bfield<bfield.size(); i_bfield++){
  for(int i_freq = 0; i_freq<samplingfreq.size(); i_freq++){
  for(int i_res = 0; i_res<spatialres.size(); i_res++){
  for(int i_rad = 0; i_rad<instrumentedrad.size(); i_rad++){
  for(int i_fidrad =0; i_fidrad<fidrad.size(); i_fidrad++){
  for(int i_ecal =0; i_ecal<ecalcontainment.size(); i_ecal++){
  TH2D* histfracbeta = new TH2D("histfracbeta", "histfracbeta", Enu.size()-1, &Enu[0], beta.size()-1, &beta[0]);
  TH2D* histfracbeta2 = new TH2D("histfracbeta2", "histfracbeta2", Enu.size()-1, &Enu[0], beta.size()-1, &beta[0]);
  TH1D* histfracenu = new TH1D("histfracenu", "histfracenu", Enu.size()-1, &Enu[0]);

  for(int i_enu = 0; i_enu<Enu.size()-1; i_enu++){
  bool morethan75percent = true;
   TString fidradpath = fidradnum[i_rad][i_fidrad].ReplaceAll(".", "_");
//   TString fidradpath = filepathfidrad.c_str();
//   fidradpath.ReplaceAll(".", "_");
//   TString filepath = "fidrad"+fidradpath+"_Bfield0_5_"+ecalcontainment[i_ecal]+"/";
   std::cout<<"fidradpath: "<<fidradpath<<std::endl;
   std::vector<TString> inputfilelist;
   std::string enubin = "Enubin"+std::to_string(i_enu);
   TString filepath = inputfile+spatialres[i_res].ReplaceAll(".", "_")+"mmpointres_"+samplingfreq[i_freq].ReplaceAll(".", "_")+"hz_fidrad"+fidradpath+"_Bfield"+bfield[i_bfield].ReplaceAll(".", "_")+"_instrumentedrad"+instrumentedrad[i_rad].ReplaceAll(".", "_")+"_"+ecalcontainment[i_ecal]+"/";
   TString inputfilename = "Selections_CC_AcceptanceCorrection_AcceptedEvents_"+enubin;
//   TString inputfilename = inputfile+"_fidrad"+fidradpath+"_Bfield0_5_"+enubin+"_"+ecalcontainment[i_ecal];
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
   std::cout<<fullfile<<std::endl;
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
//        std::cout<<"here ndgar"<<std::endl;
        for(int iplots=0; iplots<n_vars; iplots = iplots+2){
          TString keynametest = plotvariables[iplots]+"_"+plotvariables[iplots+1]+"_NDGAr_";
//          std::cout<<"keynametest: "<<keynametest<<std::endl;
          if(keyname.Contains(keynametest)){
          kinematicvariables[iplots].push_back(keyname);
//          std::cout<<"iplots: "<<iplots<<" keyname: "<<keyname<<std::endl;
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
       histnamefull = fullname+", fidrad = "+fidradpath+" B Field ="+bfield[i_bfield].ReplaceAll("_", ".")+" "+(std::to_string(Enu[i_enu])).c_str()+"< Enu < "+(std::to_string(Enu[i_enu+1])).c_str()+" "+ecalcontainment[i_ecal];

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
       int rebinnum = 1;
       TH2D* histrebin;
/*       if(i_hists==1){
       new_bins_x.clear();
       new_bins_y.clear();
       rebin_2d_based_on_events(histlist[i_hists], 10, new_bins_x, new_bins_y);
       std::sort(new_bins_x.begin(), new_bins_x.end());
       std::sort(new_bins_y.begin(), new_bins_y.end());
       if (new_bins_x.size() <= 1 || new_bins_y.size() <= 1){
       histrebin=(TH2D*)(histlist[i_hists]->Clone());
       }
       int n_new_bins_x = new_bins_x.size() - 1;
       std::cout<<"n_new_bins_x: "<<n_new_bins_x<<std::endl;
       double* bins_array_x = new double[n_new_bins_x + 1];
       std::copy(new_bins_x.begin(), new_bins_x.end(), bins_array_x);
       for(int i_x; i_x<n_new_bins_x + 1; i_x++){std::cout<<"x bin: "<<bins_array_x[i_x]<<std::endl;}
       int n_new_bins_y = new_bins_y.size() - 1;
       double* bins_array_y = new double[n_new_bins_y + 1];
       std::copy(new_bins_y.begin(), new_bins_y.end(), bins_array_y);
       for(int i_y; i_y<n_new_bins_y + 1; i_y++){std::cout<<"y bin: "<<bins_array_y[i_y]<<std::endl;}
       std::cout<<"here!!!"<<std::endl;
       TH2D* histrebin= new TH2D("histrebin", histlist[i_hists]->GetTitle(), n_new_bins_x, bins_array_x, n_new_bins_y, bins_array_y);
       std::cout<<"here!!!!"<<std::endl;
       for (int i = 1; i <= histlist[i_hists]->GetNbinsX(); ++i) {
        for (int j = 1; j <= histlist[i_hists]->GetNbinsY(); ++j) {
            histrebin->Fill(histlist[i_hists]->GetXaxis()->GetBinCenter(i), histlist[i_hists]->GetYaxis()->GetBinCenter(j), histlist[i_hists]->GetBinContent(i,j));
        }
    }
//       delete[] bins_array_x;
//       delete[] bins_array_y;
       }*/
//       else{
       if(std::string(histlist[i_hists]->GetName()).find("TrueRad") != std::string::npos){rebinnum = 4;}
       if((std::string(histlist[i_hists]->GetName()).find("TrueSquaredRad") != std::string::npos) || (std::string(histlist[i_hists]->GetName()).find("TrueXPos") != std::string::npos)){rebinnum = 10;}
       histrebin = (TH2D*)(histlist[i_hists]->Rebin2D(rebinnum, rebinnum, "histrebin"));
       std::cout<<"here   !!!"<<std::endl;
//       }
       std::cout<<finalname.c_str()<<std::endl;
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
/*       double binlessthan10 =0;
       double binstot=0;
       for(int i_x =0; i_x<histlist[i_hists-plotvariables.size()/2]->GetNbinsX(); i_x++){
         for(int i_y = 0; i_y<histlist[i_hists-plotvariables.size()/2]->GetNbinsY(); i_y++){
           int binnum = histlist[i_hists-plotvariables.size()/2]->GetBin(i_x+1, i_y+1);
           double value = histlist[i_hists-plotvariables.size()/2]->GetBinContent(binnum);
           if(value < 10 && value != 0){binlessthan10++;}
           if(value != 0){binstot++;}
         }
       }      
       int rebinnum = 1;
       if((double)(binlessthan10)/binstot > 0.85){rebinnum=2;}
       TH2D* histdenomtry = (TH2D*)(histlist[i_hists-plotvariables.size()/2]->Rebin2D(rebinnum, rebinnum, "histrebin"));
       binlessthan10 =0;
       binstot=0;
       for(int i_x =0; i_x<histdenomtry->GetNbinsX(); i_x++){
         for(int i_y = 0; i_y<histdenomtry->GetNbinsY(); i_y++){
           int binnum = histdenomtry->GetBin(i_x+1, i_y+1);
           double value = histdenomtry->GetBinContent(binnum);
           if(value < 10 && value != 0){binlessthan10++;}
           if(value != 0){binstot++;}
         }
       }      
       if((double)(binlessthan10)/binstot > 0.85){rebinnum=3;}
*/
       int rebinnum = 2;
//       if(Enu[i_enu]==0.0){rebinnum = 1;}
       if(Enu[i_enu]==4.5){rebinnum = 3;}

       if((std::string(histlist[i_hists-plotvariables.size()]->GetName()).find("TrueRad") != std::string::npos)){rebinnum =4;}
//       if((std::string(histlist[i_hists-plotvariables.size()]->GetName()).find("TrueQ0") != std::string::npos)){rebinnum =2;}
       if((std::string(histlist[i_hists-plotvariables.size()]->GetName()).find("TrueXPos") != std::string::npos) || (std::string(histlist[i_hists-plotvariables.size()]->GetName()).find("TrueSquaredRad") != std::string::npos)){rebinnum = 4;}
 //      if(Enu[i_enu]==4.0 || Enu[i_enu]==3.5){rebinnum = 2;}
 //      if(Enu[i_enu]==4.5){rebinnum = 4;}
/*       new_bins_x.clear();
       new_bins_y.clear();
       rebin_2d_based_on_events(histlist[i_hists-plotvariables.size()/2], 10, new_bins_x, new_bins_y);
       std::sort(new_bins_x.begin(), new_bins_x.end());
       std::sort(new_bins_y.begin(), new_bins_y.end());
       if (new_bins_x.size() <= 1 || new_bins_y.size() <= 1){
       TH2D* histrebin=(TH2D*)(histlist[i_hists]->Clone());
       }
       int n_new_bins_x = new_bins_x.size() - 1;
       int n_new_bins_y = new_bins_y.size() - 1;
       std::cout<<"n_new_bins_x: "<<n_new_bins_x<<std::endl;
       std::cout<<"n_new_bins_y: "<<n_new_bins_y<<std::endl;
       double* bins_array_x = new double[n_new_bins_x + 1];
       std::copy(new_bins_x.begin(), new_bins_x.end(), bins_array_x);
       for(int i_x; i_x<n_new_bins_x + 1; i_x++){std::cout<<"x bin: "<<bins_array_x[i_x]<<std::endl;}
       double* bins_array_y = new double[n_new_bins_y + 1];
       std::copy(new_bins_y.begin(), new_bins_y.end(), bins_array_y);
       for(int i_y; i_y<n_new_bins_y + 1; i_y++){std::cout<<"y bin: "<<bins_array_y[i_y]<<std::endl;}
       std::cout<<"here!!!"<<std::endl;
       TH2D* histrebin= new TH2D("histrebin", histratios->GetTitle(), n_new_bins_x, bins_array_x, n_new_bins_y, bins_array_y);
       TH2D* histdenom= new TH2D("histdenom", histlist[i_hists-plotvariables.size()/2]->GetTitle(), n_new_bins_x, bins_array_x, n_new_bins_y, bins_array_y);
       std::cout<<"here!!!!"<<std::endl;
       for (int i = 1; i <= histratios->GetNbinsX(); ++i) {
        for (int j = 1; j <= histratios->GetNbinsY(); ++j) {
            histrebin->Fill(histratios->GetXaxis()->GetBinCenter(i), histratios->GetYaxis()->GetBinCenter(j), histratios->GetBinContent(i,j));
            histdenom->Fill(histlist[i_hists-plotvariables.size()/2]->GetXaxis()->GetBinCenter(i), histlist[i_hists-plotvariables.size()/2]->GetYaxis()->GetBinCenter(j), histlist[i_hists-plotvariables.size()/2]->GetBinContent(i,j));
        }
    }*/
//       delete[] bins_array_x;
//       delete[] bins_array_y;
       TH2D* histdenom = (TH2D*)(histlist[i_hists-plotvariables.size()/2]->Rebin2D(rebinnum, rebinnum, "histrebin"));
       TH2D* histrebin = (TH2D*)(histratios->Rebin2D(rebinnum, rebinnum, "histrebin"));
       std::cout<<"cloned"<<std::endl;
       histrebin->GetXaxis()->SetTitle(x_axis[(i_hists-plotvariables.size())*2].c_str());
       histrebin->GetYaxis()->SetTitle(x_axis[(i_hists-plotvariables.size())*2+1].c_str());
       TString titleacceptance = "Acceptance vs Q3 vs Q0, fidrad = "+fidradpath.ReplaceAll("_", ".")+", B Field = "+bfield[i_bfield].ReplaceAll("_", ".")+"T, "+(std::to_string(Enu[i_enu])).c_str()+" < Enu < "+(std::to_string(Enu[i_enu+1])).c_str()+" "+ecalcontainment[i_ecal];
       histrebin->SetTitle(titleacceptance);
       histrebin->Sumw2();
       histdenom->Sumw2();
       histrebin->Divide(histdenom);
       double sum = 0;
       double sum_weights =0;
       int n_valid_bins = 0;
       for(int i_x =0; i_x<histrebin->GetNbinsX(); i_x++){
         for(int i_y = 0; i_y<histrebin->GetNbinsY(); i_y++){
           int binnum = histrebin->GetBin(i_x+1, i_y+1);
           double value = histrebin->GetBinContent(binnum);
           double denominator = histdenom->GetBinContent(binnum);
           if(value == 0 && denominator != 0){histrebin->SetBinContent(binnum, 1e-6);}
           double bin_content = histrebin->GetBinContent(binnum);
           if(bin_content>0){
             double weight = pow(histdenom->GetBinContent(binnum), 0.5);
             weight_i.push_back(weight);
             sum += bin_content*weight;
             sum_weights += weight;
             bin_i.push_back(bin_content);
             n_valid_bins++;
             std::vector<double> nearbins;
             std::vector<double> nearbins_weights;
/*             int binnum1 = histrebin->GetBin(i_x, i_y+1);
             if(histrebin->GetBinContent(binnum1) !=0){nearbins.push_back(histrebin->GetBinContent(binnum1));nearbins_weights.push_back(pow(histdenom->GetBinContent(binnum1), 0.5));}
             int binnum2 = histrebin->GetBin(i_x, i_y);
             if(histrebin->GetBinContent(binnum2) !=0){nearbins.push_back(histrebin->GetBinContent(binnum2));nearbins_weights.push_back(pow(histdenom->GetBinContent(binnum2), 0.5));}
             int binnum3 = histrebin->GetBin(i_x+1, i_y);
             if(histrebin->GetBinContent(binnum3) !=0){nearbins.push_back(histrebin->GetBinContent(binnum3));nearbins_weights.push_back(pow(histdenom->GetBinContent(binnum3), 0.5));} 
             int binnum4 = histrebin->GetBin(i_x+2, i_y);
             if(histrebin->GetBinContent(binnum4) !=0){nearbins.push_back(histrebin->GetBinContent(binnum4));nearbins_weights.push_back(pow(histdenom->GetBinContent(binnum4), 0.5));}
             int binnum5 = histrebin->GetBin(i_x+2, i_y+1);
             if(histrebin->GetBinContent(binnum5) !=0){nearbins.push_back(histrebin->GetBinContent(binnum5));nearbins_weights.push_back(pow(histdenom->GetBinContent(binnum5), 0.5));}
             int binnum6 = histrebin->GetBin(i_x+2, i_y+2);
             if(histrebin->GetBinContent(binnum6) !=0){nearbins.push_back(histrebin->GetBinContent(binnum6));nearbins_weights.push_back(pow(histdenom->GetBinContent(binnum6), 0.5));}
             int binnum7 = histrebin->GetBin(i_x+1, i_y+2);
             if(histrebin->GetBinContent(binnum7) !=0){nearbins.push_back(histrebin->GetBinContent(binnum7));nearbins_weights.push_back(pow(histdenom->GetBinContent(binnum7), 0.5));}
             int binnum8 = histrebin->GetBin(i_x, i_y+2);
             if(histrebin->GetBinContent(binnum8) !=0){nearbins.push_back(histrebin->GetBinContent(binnum8));nearbins_weights.push_back(pow(histdenom->GetBinContent(binnum8), 0.5));}
             surroundingbins.push_back(nearbins);
             surroundingbins_weights.push_back(nearbins_weights);
             nearbins.clear();
             nearbins_weights.clear();*/
           }
         }
       }
//       double mean = sum/sum_weights;
/*       std::vector<double> stderrors;
       std::vector<double> stddevs;
       std::vector<double> means;
       //std::cout<<"surroundingbins.size(): "<<surroundingbins.size()<<std::endl;
       if(surroundingbins.size()>0){
       //std::cout<<"HERE"<<std::endl;
       for(int i_err=0; i_err<surroundingbins.size(); i_err++){
         //std::cout<<"HERE"<<std::endl;
         double summean = bin_i[i_err];
         double sumweights = weight_i[i_err];
         //std::cout<<"summean: "<<summean<<" sumweights: "<<sumweights<<std::endl;
         for(int i_bins; i_bins<surroundingbins[i_err].size(); i_bins++){
           summean += surroundingbins[i_err][i_bins]*surroundingbins_weights[i_err][i_bins];
           sumweights += surroundingbins_weights[i_err][i_bins];
         }
         double meanweighted = summean/sumweights;
         double sumnumerator =weight_i[i_err]*pow((bin_i[i_err]-meanweighted), 2);
         for(int i_bins; i_bins<surroundingbins[i_err].size(); i_bins++){
           sumnumerator += surroundingbins_weights[i_err][i_bins]*pow((surroundingbins[i_err][i_bins]-meanweighted), 2);
         }
         double M = surroundingbins_weights[i_err].size();
         double stddev_i = pow(sumnumerator/(sumweights*(M/(M+1))), 0.5);
         //std::cout<<"surroundingbins_weights[i_err].size(): "<<surroundingbins_weights[i_err].size()<<std::endl;
         //std::cout<<"(sumweights*(surroundingbins_weights[i_err].size()/(surroundingbins_weights[i_err].size()+1))): "<<(sumweights*(surroundingbins_weights[i_err].size()/(surroundingbins_weights[i_err].size()+1)))<<std::endl;
         //std::cout<<"weights fraction: "<<(surroundingbins_weights[i_err].size()/(surroundingbins_weights[i_err].size()+1))<<std::endl;

         //std::cout<<"stddev_i: "<<stddev_i<<" meanweighted: "<<meanweighted<<" sumnumerator/: "<<sumnumerator<<" sumweights: "<<sumweights<<std::endl;
         double stderror_i = stddev_i/pow(M, 0.5);
         stderrors.push_back(stderror_i);
         stddevs.push_back(stddev_i);
         means.push_back(meanweighted);
         }
        }
        else{stderrors.push_back(1.0);}
*/        
//       double stddev = std::sqrt(sum_sq / n_valid_bins - mean * mean);
//       std::cout<<"stddev: "<<stddev<<std::endl;
//       histrebin->Sumw2();
       for(int i_beta = 0; i_beta<beta.size(); i_beta++){
         int count = 0;
         int total = 0;
         int count2 =0;
         int count3 =0;
         double probsum = 0;
         for(int i_x =0; i_x<histrebin->GetNbinsX(); i_x++){
           for(int i_y = 0; i_y<histrebin->GetNbinsY(); i_y++){
             if(histrebin->GetBinContent(i_x+1, i_y+1)!=0){
               total++;
//               std::cout<<"histrebin->GetBinContent(i_x+1, i_y+1): "<<histrebin->GetBinContent(i_x+1, i_y+1)<<" histdenom->GetBinContent(i_x+1, i_y+1): "<<histdenom->GetBinContent(i_x+1, i_y+1)<<std::endl;
//               std::cout<<"error bin "<<i_x+1<<", "<<i_y+1<<": "<<histrebin->GetBinError(i_x+1, i_y+1)<<std::endl;
//               histrebin->SetBinError(i_x+1, i_y+1, histrebin->GetBinContent(i_x+1, i_y+1)*(1/(pow(histdenom->GetBinContent(i_x+1, i_y+1), 0.5))));
               if(histrebin->GetBinContent(i_x+1, i_y+1)<=beta[i_beta+1]){
               count++;
               }
               if(histrebin->GetBinContent(i_x+1, i_y+1)<=beta[i_beta+1] && histrebin->GetBinContent(i_x+1, i_y+1)>beta[i_beta]){
               count2++;
               }
               //std::cout<<"hist bin content: "<<histrebin->GetBinContent(i_x+1, i_y+1)<<" bin error: "<<histrebin->GetBinError(i_x+1, i_y+1)<<std::endl;
//               if(histrebin->GetBinContent(i_x+1, i_y+1)-histrebin->GetBinError(i_x+1, i_y+1)<=beta[i_beta+1]){
//               count3++;
//               }
               //std::cout<<" total: "<<total<<std::endl;
               //std::cout<<"stderrors size: "<<stderrors.size()<<std::endl;
               //std::cout<<"stderrors[total-1]: "<<stderrors[total-1]<<std::endl;
/*               if(histrebin->GetBinContent(i_x+1, i_y+1)-stderrors[total-1]<=beta[i_beta+1]){
               std::cout<<" bin content: "<<histrebin->GetBinContent(i_x+1, i_y+1)<<std::endl;
               count3++;
               }*/
/*               double mean = means[total-1];
               double stddev = stddevs[total-1];
               if(stddev == 0 || std::isnan(stddev)){stddev = 0.01;}
               if(std::isinf(stddev)){stddev = 1;}
               double prob = (std::erf((0.1-mean)/(pow(2, 0.5)*stddev))-std::erf((0.0-mean)/(pow(2, 0.5)*stddev)))/(std::erf((1.0-mean)/(pow(2, 0.5)*stddev))-std::erf((0.0-mean)/(pow(2, 0.5)*stddev)));
               //std::cout<<"stddev: "<<stddev<<" mean: "<<mean<<" std::erf((0.1-mean)/(pow(2, 0.5)*stddev): "<<std::erf((0.1-mean)/(pow(2, 0.5)*stddev))<<" std::erf((0.0-mean)/(pow(2, 0.5)*stddev)): "<<std::erf((0.0-mean)/(pow(2, 0.5)*stddev))<<" std::erf((1.0-mean)/(pow(2, 0.5)*stddev)): "<<std::erf((1.0-mean)/(pow(2, 0.5)*stddev))<<std::endl;
               //std::cout<<"prob: "<<prob<<std::endl;
               probsum += prob*(1-prob);*/
             }
           }
         }
/*         bin_i.clear();
         weight_i.clear();
         surroundingbins.clear();
         surroundingbins_weights.clear();*/

//         std::cout<<"count: "<<count<<" total: "<<total<<" probsum: "<<probsum<<std::endl;
         double fraction = (double)count/(double)total;
         double fraction2 = (double)(count2)/(double)total;
//         double fraction3 = (double)(count3)/(double)total;
         
//         double sigma_frac = pow(probsum, 0.5)/(double)(total);
//         std::cout<<"sigma_frac: "<<sigma_frac<<std::endl;
         if(fraction2 == 0){fraction2 = 1e-6;}
         count2 =0;
         histfracbeta->SetBinContent(i_enu+1, i_beta+1, fraction);
         histfracbeta2->SetBinContent(i_enu+1, i_beta+1, fraction2);
         if(beta[i_beta+1] == threshold){
         std::vector<double> fracbelowthreshold;
         double sum =0;
         std::default_random_engine generator;
         for(int i_trials =0; i_trials<100; i_trials++){
           int countbelow = 0;
           for(int i_x =0; i_x<histrebin->GetNbinsX(); i_x++){
             for(int i_y = 0; i_y<histrebin->GetNbinsY(); i_y++){
               double acceptanceprob = histrebin->GetBinContent(i_x+1, i_y+1);
//               std::cout<<"acceptanceprob: "<<acceptanceprob<<std::endl;
               if(acceptanceprob!=0){
                 int ntrials = histdenom->GetBinContent(i_x+1, i_y+1);
                 std::binomial_distribution<int> distribution(ntrials,acceptanceprob);
                 int naccepted = distribution(generator);
//                 std::cout<<"naccepted: "<<naccepted<<" ntrials: "<<ntrials<<std::endl;
                 double frac = (double)(naccepted)/(double)(ntrials);
                 if(frac<threshold){countbelow++;}
               }
             }
           }
           double fracbelow = (double)(countbelow)/(double)(total);
//           std::cout<<"fracbelow: "<<fracbelow<<std::endl;
           fracbelowthreshold.push_back(fracbelow);
           sum += fracbelow;
         }
         double meanfrac = sum/fracbelowthreshold.size();
//         std::cout<<"meanfrac: "<<meanfrac<<std::endl;
         double sumsquare = 0;
         for(int i_frac = 0; i_frac<fracbelowthreshold.size(); i_frac++){
           std::cout<<"fracbelowthreshold[i_frac]: "<<fracbelowthreshold[i_frac]<<std::endl;
           sumsquare += pow(fracbelowthreshold[i_frac] - meanfrac, 2);
         }
         double stddev_new = pow(sumsquare/(fracbelowthreshold.size()-1), 0.5);
         std::cout<<"stddev_new: "<<stddev_new<<std::endl;
         
         histfracenu->SetBinContent(i_enu+1, fraction);
         histfracenu->SetBinError(i_enu+1, stddev_new);
//          std::cout<<"count: "<<count<<" total: "<<total<<std::endl; 
//          std::cout<<"i_enu+1: "<<i_enu+1<<" i_beta+1: "<<i_beta+1<<" fraction: "<<fraction<<"Enu: "<<Enu[i_enu]<<"beta: "<<beta[i_beta+1]<<std::endl;};
       }
       std::cout<<"divided"<<std::endl;
       c0->cd(i_hists);
       histrebin->SetContour(100);
       histrebin->Draw("COLZ");
       histrebin->SetMarkerSize(0.2);
       histrebin->Draw("TEXT SAME");
//       histratios->Draw("COLZ");
       c0->Print("acceptancecorrectionvars_2d.ps");
       ncd =i_hists;
     }
     }
/*      c0->cd(1);
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
   } 
//   }
//   } 
   //}
   ncd++;
   TString titlelegend = "Bfield = "+bfield[i_bfield].ReplaceAll("_", ".")+ " T";
//   TString titlelegend = "Res = "+spatialres[i_res].ReplaceAll("_", ".")+" mm";
//   TString titlelegend = "Freq = "+samplingfreq[i_freq].ReplaceAll("_", ".")+" MHz, Fid rad = "+fidrad[i_fidrad].ReplaceAll("_", ".")+" cm from edge";

//   TString titlelegend = "Fid rad = "+fidrad[i_fidrad].ReplaceAll("_", ".")+" cm from edge, Instrumented rad = "+instrumentedrad[i_rad].ReplaceAll("_", ".");

//   TString titlelegend = "Fid rad = "+fidrad[i_fidrad].ReplaceAll("_", ".")+" cm from edge";
   TString title = "Bfield = "+bfield[i_bfield].ReplaceAll("_", ".")+ "Res = "+spatialres[i_res].ReplaceAll("_", ".")+" mm, Freq = "+samplingfreq[i_freq].ReplaceAll("_", ".")+" MHz, Fid rad = "+fidrad[i_fidrad].ReplaceAll("_", ".")+" cm from edge";


//   TString title = "Resolution = "+spatialres[i_res].ReplaceAll("_", ".")+" mm Sampling Freq = "+samplingfreq[i_freq].ReplaceAll("_", ".")+" MHz, Fid rad = "+fidradnum[i_rad][i_fidrad].ReplaceAll("_", ".")+" cm Bfield = "+bfield[i_bfield].ReplaceAll("_", ".")+" Instrumented rad = "+instrumentedrad[i_rad].ReplaceAll("_", ".");
//   TString title = "fidrad = "+fidrad[i_fidrad].ReplaceAll("_", ".")+", Bfield = 0.5, "+ecalcontainment[i_ecal];
   c0->cd(ncd);
   histfracbeta->SetContour(200);
   histfracbeta->Draw("COLZ0");
   histfracbeta->GetXaxis()->SetTitle("E#nu (GeV)");
   histfracbeta->GetYaxis()->SetTitle("#beta threshold");
   histfracbeta->SetTitle(title);
   c0->Print("acceptancecorrectionvars_2d.ps");
   ncd++;
   c0->cd(ncd);
   histfracbeta2->SetContour(200);
   histfracbeta2->Draw("COLZ0");
   histfracbeta2->GetXaxis()->SetTitle("E#nu (GeV)");
   histfracbeta2->GetYaxis()->SetTitle("#beta range");
   histfracbeta2->SetTitle(title);
   c0->Print("acceptancecorrectionvars_2d.ps");

   histvector.push_back(histfracenu);
   histnamesfracenu.push_back(titlelegend);
     }
   }
   }
   } 
   }
   }
   std::vector<double> xfracline = {0, 2.5, 4.5};
   std::vector<double> yfracline = {0.03,0.03, 0.07};
   TGraph* fracline = new TGraph(xfracline.size(), &xfracline[0], &yfracline[0]);
   std::vector<Color_t> colours = {kRed+2, kRed+1, kRed,kRed-7, kOrange+2, kOrange+1, kOrange, kOrange-7, kBlue+2, kBlue+1, kBlue, kBlue-7, kGreen+2, kGreen+1, kGreen, kGreen-7, kMagenta+2, kMagenta+1, kMagenta, kMagenta-7, kCyan+2, kCyan+1, kCyan, kCyan-7};
   TLegend *legend = new TLegend(0.151,0.65,0.35,0.88);
//  TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
   c0->cd(ncd+1);
   gStyle->SetPalette(kRainbow);
   c0->Update();
//   histfracenu->SetContour(100);
   histvector[0]->GetYaxis()->SetRangeUser(0.0, 0.14);  
//   histvector[0]->SetLineStyle(2);
   histvector[0]->SetLineColor(colours[0]);
   histvector[0]->SetLineWidth(1);
//   histvector[0]->Draw("P E");
   histvector[0]->Draw("L");
   histvector[0]->GetXaxis()->SetTitle("E_{#nu} (GeV)");
   histvector[0]->GetYaxis()->SetTitle("Fraction");
   histvector[0]->SetTitle("Fraction of Phase Space with less than 10% Acceptance");
   legend->AddEntry(histvector[0],histnamesfracenu[0]);
     
   for(int i_histenu = 1; i_histenu<histvector.size(); i_histenu++){
     int linecol = int(i_histenu+1);
//     int linecol = int(i_histenu/2)+1;
     if(linecol == 10){linecol = 12;}
     histvector[i_histenu]->SetLineColor(colours[i_histenu*4+1]);
     histvector[i_histenu]->SetLineWidth(1);
//     if(i_histenu%2==0){histvector[i_histenu]->SetLineStyle(2);}
//     histvector[i_histenu]->Draw("P E SAME");
     histvector[i_histenu]->Draw("L SAME");
     legend->AddEntry(histvector[i_histenu],histnamesfracenu[i_histenu]);
   }
   fracline->SetLineWidth(3);
   fracline->Draw("L SAME");
   legend->SetBorderSize(0);
   legend->SetTextSize(0.025);
   legend->SetFillColor(0);
   legend->SetFillStyle(1001);
   legend->Draw();
   c0->Print("acceptancecorrectionvars_2d.ps");
   c0->Print("acceptancecorrectionvars_2d.png");
   std::cout<<"HERE"<<std::endl;
   c0->Print("acceptancecorrectionvars_2d.ps]");

}
