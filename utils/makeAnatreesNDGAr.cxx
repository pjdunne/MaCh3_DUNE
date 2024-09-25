#include <TFile.h>
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <stdio.h>
#include <iostream>
#include <cmath>
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TDirectory.h"
//#include "Framework/Ntuple/NtpMCEventRecord.h"
//#include "EVGCore/EventRecord.h"
//#include "nusystematics/artless/response_helper.hh"
//#include "duneanaobj/StandardRecord/StandardRecord.h"
#include <filesystem>

/*
void Split(std::string infile, int initial_nuPDG, int final_nuPDG, bool numuselec) 
{  

   const char* path  = "/vols/dune/nk3717/data/NDGAr_10kCAFs";
   const char* outpath  = "/vols/dune/nk3717/data/NDGAr_10kCAFs/Outputs";
   
   caf::StandardRecord * sr = new caf::StandardRecord();
   genie::NtpMCEventRecord * gr = new genie::NtpMCEventRecord();
   caf::StandardRecord * srnew = new caf::StandardRecord();
   genie::NtpMCEventRecord * grnew = new genie::NtpMCEventRecord();
//   caf::StandardRecord * sr;
//  sr = &sr1;
//   std::cout<<sizeof(sr1)<<std::endl;
   int nutype;
   int nutypeunosc;
   double cvn;
   float cvn_threshold;

//   std::vector<float> muonparticle;  //mc.nu.prim[0].p.Mag() thats the magnitude of the primary muon energy. plot this as histogram
   std::string initial_nu_name; 
   std::string final_nu_name; 
   std::string selec_name;
   std::string outfile_name; 
   std::string cvn_name;

   if(numuselec){ selec_name = "numuselec";}
   else{ selec_name = "nueselec";}
   if (initial_nuPDG == 14) {initial_nu_name = "numu_x_";}    
   if (initial_nuPDG == -14) {initial_nu_name = "numubar_x_";}    
   if (initial_nuPDG == 12) {initial_nu_name = "nue_x_";}    
   if (initial_nuPDG == -12) {initial_nu_name = "nuebar_x_";}    
   if (initial_nuPDG == 16) {initial_nu_name = "nutau_x_";}    
   if (initial_nuPDG == -16) {initial_nu_name = "nutaubar_x_";}    

   if (final_nuPDG == 14) {final_nu_name = "numu_";}    
   if (final_nuPDG == -14) {final_nu_name = "numubar_";}    
   if (final_nuPDG == 12) {final_nu_name = "nue_";}    
   if (final_nuPDG == -12) {final_nu_name = "nuebar_";}    
   if (final_nuPDG == 16) {final_nu_name = "nutau_";}    
   if (final_nuPDG == -16) {final_nu_name = "nutaubar_";}    

   outfile_name = infile.substr(0, 13) + "_" + initial_nu_name + final_nu_name + selec_name + ".root";
   
   char filepath[200];
   char outfilepath[200];
   sprintf(outfilepath, "%s/%s", outpath, outfile_name.c_str());
   TFile * outf = new TFile(outfilepath, "RECREATE");
   int neventstot =0;

   TFile * tfstructure = new TFile("/vols/dune/nk3717/data/NDGAr_10kCAFs/gar_caf_13229676_1_20240615.root", "READ");

   TTree * treestructure = NULL;
   tfstructure->GetObject("cafTree", treestructure);
   TTree * gtreestructure = NULL;
   tfstructure->GetObject("genieEvt", gtreestructure);
   TTree * newtree = treestructure->CloneTree(0);
   newtree->SetDirectory(outf);
   TTree * gnewtree = gtreestructure->CloneTree(0);
   gnewtree->SetDirectory(outf);
   std::cout<<"cloned"<<std::endl;
   outf->cd();
   for (auto const& dir_entry : std::filesystem::directory_iterator{path}){
     if(dir_entry.is_regular_file()){
       sprintf(filepath, "%s/%s", path, dir_entry.path().filename().c_str());
       std::cout<<"Opening CAF: "<<filepath<<std::endl;
       TFile * tf = new TFile(filepath, "READ");
       TTree * tree = NULL;
       tf->GetObject("cafTree", tree);
       TTree * gtree = NULL;
       tf->GetObject("genieEvt", gtree);
       TBranch * branch = tree ->GetBranch("rec");
       branch ->SetAddress(&sr); 
       gtree->SetBranchStatus("*",1);
       TBranch * gbranch = gtree->GetBranch("genie_record");
       gbranch->SetAddress(&gr);
       gtree->SetAlias("genie_record", "gmcrec");
       std::cout<<"nevents new tree"<<newtree->GetEntries()<<std::endl;
       std::cout<<"nevents gnew tree"<<gnewtree->GetEntries()<<std::endl;

       int N = tree->GetEntries();
       int GN = gtree ->GetEntries();
       std::cout << "N = " << N <<" GN = "<<GN<< std::endl;
       std::cout<<"nevents new tree"<<newtree->GetEntries()<<std::endl;
       std::cout<<"nevents gnew tree"<<gnewtree->GetEntries()<<std::endl;

       for (int i = 0; i < N; i++) {
         std::cout<<"before get entry"<<std::endl;
         tree->GetEntry(i);
         std::cout<<"after get entry i: "<<i<<std::endl;
         if(numuselec) {
           cvn_name = "cvnnumu";
           cvn = (double)(sr->common.ixn.gsft[0].nuhyp.cvn.numu);
           cvn_threshold = -0.5; 
           selec_name = "numuselec";}
         else {
           cvn_name = "cvnnue";
           cvn = (double)(sr->common.ixn.gsft[0].nuhyp.cvn.nue);
           cvn_threshold = -0.85; 
           selec_name = "nueselec";}
         std::cout<<"pdgorig "<<(double)(sr->mc.nu[0].pdgorig)<<" pdg "<<(double)(sr->mc.nu[0].pdg)<<" cvn "<<cvn<< " cvn threshold " << cvn_threshold <<std::endl;
    
         if(std::isnan(cvn)){cvn =0;}
         std::cout<<"pdgorig "<<(double)(sr->mc.nu[0].pdgorig)<<" pdg "<<(double)(sr->mc.nu[0].pdg)<<" cvn "<<cvn<< " cvn threshold " << cvn_threshold <<std::endl;
         if((sr->mc.nu[0].pdgorig) == initial_nuPDG and (sr->mc.nu[0].pdg) == final_nuPDG and cvn >= cvn_threshold) {
           tree->GetEntry(i);
           newtree->Fill();
           gtree->GetEntry(i);
           gnewtree->Fill();
         }  
       }
    
       int n_selection =  newtree->GetEntries();
       neventstot = neventstot+n_selection;
       int gn_selection =  gnewtree->GetEntries();
       std::cout << "Check: number of events in " << outfile_name << " is: " << n_selection << std::endl;
       std::cout << "Check: number of events in GENIE TREE is: " << gn_selection << std::endl;
      delete tf;
    }
  }
  if(neventstot !=0){
    outf->cd();
    newtree->Write("cafTree");
    gnewtree->Write("gtree");
    std::cout << "Writing output ROOT File" << std::endl;
    outf->Write();
    outf->Close();
  }
  else {
    std::cout << "No events inside " << outfile_name << " , so no output ROOT file " << std::endl;
    remove(outfilepath);
  }
  
  delete sr;
  delete gr;
  delete tfstructure;
}
*/

void Split(std::string infile, int initial_nuPDG, int final_nuPDG, bool numuselec) 
{  

   const char* path  = "/vols/dune/nk3717/data/NDGAr_100kCAFs/AnaTrees/";
   const char* outpath  = "/vols/dune/nk3717/data/NDGAr_100kCAFs/AnaTreesOutputs";

   vector<int> *ntypes=0;
   int nutype;
   int nutypeunosc;
   double cvn;
   float cvn_threshold;

   std::string initial_nu_name; 
   std::string final_nu_name; 
   std::string selec_name;
   std::string outfile_name; 
   std::string cvn_name;

   if(numuselec){ selec_name = "numuselec";}
   else{ selec_name = "nueselec";}
   if (initial_nuPDG == 14) {initial_nu_name = "numu_x_";}    
   if (initial_nuPDG == -14) {initial_nu_name = "numubar_x_";}    
   if (initial_nuPDG == 12) {initial_nu_name = "nue_x_";}    
   if (initial_nuPDG == -12) {initial_nu_name = "nuebar_x_";}    
   if (initial_nuPDG == 16) {initial_nu_name = "nutau_x_";}    
   if (initial_nuPDG == -16) {initial_nu_name = "nutaubar_x_";}    

   if (final_nuPDG == 14) {final_nu_name = "numu_";}    
   if (final_nuPDG == -14) {final_nu_name = "numubar_";}    
   if (final_nuPDG == 12) {final_nu_name = "nue_";}    
   if (final_nuPDG == -12) {final_nu_name = "nuebar_";}    
   if (final_nuPDG == 16) {final_nu_name = "nutau_";}    
   if (final_nuPDG == -16) {final_nu_name = "nutaubar_";}    

   outfile_name = infile.substr(0, 13) + "_" + initial_nu_name + final_nu_name + selec_name + "_geant"+".root";
   
   char filepath[200];
   char outfilepath[200];
   sprintf(filepath, "%s/%s", path, infile.c_str());
   sprintf(outfilepath, "%s/%s", outpath, outfile_name.c_str());


   TFile * tf = new TFile(filepath, "READ");
   TFile * outf = new TFile(outfilepath, "RECREATE");     

   TDirectory* dir = tf->GetDirectory("anatree");
   TTree * tree = (TTree*)dir->Get("GArAnaTree");
   if(!tree){std::cout<<"Tree not found"<<std::endl;}
   TTree * newtree = tree->CloneTree(0);
   std::cout<<"cloned"<<std::endl;
   tree->SetBranchStatus("NType", 1);
   tree->SetBranchAddress("NType", &ntypes);
  
   int N = tree->GetEntries();
   std::cout << "N = " << N << std::endl; 
   for (int i = 0; i < N; i++) {

     tree->GetEntry(i);

/*     if(numuselec) {
       cvn_name = "cvnnumu";
       cvn = (double)(sr->common.ixn.gsft[0].nuhyp.cvn.numu);
       cvn_threshold = -0.5; 
       selec_name = "numuselec";}
     else {
       cvn_name = "cvnnue";
       cvn = (double)(sr->common.ixn.gsft[0].nuhyp.cvn.nue);
       cvn_threshold = -0.85; 
       selec_name = "nueselec";}
*/
//     if(std::isnan(cvn)){cvn =0;}
//     std::cout<<"pdgorig "<<(double)(sr->mc.nu[0].pdgorig)<<" pdg "<<(double)(sr->mc.nu[0].pdg)<<" cvn "<<cvn<< " cvn threshold " << cvn_threshold <<std::endl;
     if(!ntypes){std::cout<<"Ntypes not found"<<std::endl;}
     int NType = ntypes->at(0);
     //std::cout<<"NType: "<<NType<<std::endl;     
     if(NType == initial_nuPDG and NType == final_nuPDG) {
       tree->GetEntry(i);
       newtree->Fill();
     }
   }

   int n_selection =  newtree->GetEntries();
   std::cout << "Check: number of events in " << outfile_name << " is: " << n_selection << std::endl;
   
   if (n_selection != 0) {
     newtree->Write("GArAnaTree");
     std::cout << "Writing output ROOT File" << std::endl;
     outf->Write();
     outf->Close();
   }
   
   else {
     std::cout << "No events inside " << outfile_name << " , so no output ROOT file " << std::endl;
     remove(outfilepath);
   }
  tree->ResetBranchAddresses();
}

/*
void Plot(std::string infile) 
{  

   const char* path  = "/vols/dune/nk3717/data/NDGAr_100kCAFs";
//   const char* outpath  = "/vols/dune/nk3717/MaCh3_refactortry/MaCh3_DUNE/inputs/DUNE_NDGAr_CAF_files/";
   
//   caf::StandardRecord * sr = new caf::StandardRecord();
   std::vector<float> muonparticle;  //mc.nu.prim[0].p.Mag() thats the magnitude of the primary muon energy. plot this as histogram
   int n_int;
   char filepath[200];
   sprintf(filepath, "%s/%s", path, infile.c_str());

   TFile * tf = new TFile(filepath, "READ");

   TTree * tree = NULL;
   tf->GetObject("cafTree", tree);
   TBranch * branch = tree ->GetBranch("rec");
   branch ->SetAddress(&sr); 
   int N = tree->GetEntries();
   std::cout << "N = " << N << std::endl; 
   for (int i = 0; i < N; i++) {
     tree->GetEntry(i);
     //n_int = (int)(sr->mc.nu.size());
     n_int =1;
     std::cout<<"n_int "<<n_int<<std::endl;
     for(int j = 0; j <n_int; j++){
       muonparticle.push_back((float)(sr->mc.nu[0].pdgorig)); //sr->mc.nu[j].prim[0].p.Mag()));
       std::cout<<"muon energy "<<(float)(sr->mc.nu[0].pdgorig) << std::endl;
     }
   }
  TH1F *h1 = new TH1F("h1", "PDG Number", 40, -20, 20);
  int n_energies = muonparticle.size();
  for(int i_muon =0; i_muon < n_energies; i_muon++){
    h1->Fill(muonparticle[i_muon]);
  }
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  h1->Draw();
  h1->SetTitle("Original PDG; PDG Number; Counts");
  c1->SaveAs("pdgorig.png");
  tree->ResetBranchAddresses();
  delete tf;
  //delete tree;
  //delete branch;
  delete sr;
}

*/
int makeAnatreesNDGAr(std::string infile)
{
/*
  if( (argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1])) ) {
    std::cout << "Help yourself by looking at the source code to see what the options are." << std::endl;
    return 0;
  }

  std::string infile;
  infile = argv[1];*/
//  Plot(infile);
 
  // NUMU_X_NUMU
//  Split(infile, 14, 14, true); //numu_x_numu_numuselec
//  Split(infile, 14, 14, false); //numu_x_numu_nueselec
//  Split(infile, -14, -14, true); //numubar_x_numubar_numuselec
//  Split(infile, -14, -14, false); //numubar_x_numubar_nueselec
  
  // NUMU_X_NUE
//  Split(infile, 14, 12, true); //numu_x_nue_numuselec
//  Split(infile, 14, 12, false); //numu_x_nue_nueselec
//  Split(infile, -14, -12, true); //numubar_x_nuebar_numuselec
//  Split(infile, -14, -12, false); //numubar_x_nuebar_nueselec
  
  // NUMU_X_NUTAU
//  Split(infile, 14, 16, true); //numu_x_nutau_numuselec
//  Split(infile, 14, 16, false); //numu_x_nutau_nueselec
//  Split(infile, -14, -16, true); //numubar_x_nutaubar_numuselec
//  Split(infile, -14, -16, false); //numubar_x_nutaubar_nueselec
  
  // NUE_X_NUE
  Split(infile, 12, 12, true); //nue_x_nue_numuselec 
//  Split(infile, 12, 12, false); //nue_x_nue_nueselec
  Split(infile, -12, -12, true); //nuebar_x_nuebar_numuselec
//  Split(infile, -12, -12, false); //nuebar_x_nuebar_nueselec
  
  // NUE_X_NUMU
//  Split(infile, 12, 14, true); //nue_x_numu_numuselec
//  Split(infile, 12, 14, false); //nue_x_numu_nueselec
//  Split(infile, -12, -14, true); //nuebar_x_numubar_numuselec
//  Split(infile, -12, -14, false); //nuebar_x_numubar_nueselec
  
  //NUE_X_NUTAU
//  Split(infile, 12, 16, true); //nue_x_nutau_numuselec
//  Split(infile, 12, 16, false); //nue_x_nutau_nueselec
//  Split(infile, -12, -16, true); //nuebar_x_nutaubar_numuselec
//  Split(infile, -12, -16, false); //nuebar_x_nutaubar_nueselec

  return 0;
}




  
