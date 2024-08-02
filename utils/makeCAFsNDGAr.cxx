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
//#include "Ntuple/NtpMCEventRecord.h"
//#include "EVGCore/EventRecord.h"
//#include "nusystematics/artless/response_helper.hh"
#include "duneanaobj/StandardRecord/StandardRecord.h"


void Split(std::string infile, int initial_nuPDG, int final_nuPDG, bool numuselec) 
{  

   const char* path  = "/vols/dune/nk3717/data/NDGAr_newtestCAFs";
   const char* outpath  = "/vols/dune/nk3717/data/NDGAr_newtestCAFs";
   
   caf::StandardRecord * sr = new caf::StandardRecord();
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

 //  if(numuselec){ selec_name = "numuselec";}
 //  else{ selec_name = "nueselec";}
   if (initial_nuPDG == 14) {initial_nu_name = "numu_x_";}    
   if (initial_nuPDG == -14) {initial_nu_name = "numubar_x_";}    
   if (initial_nuPDG == 12) {initial_nu_name = "nue_x_";}    
   if (initial_nuPDG == -12) {initial_nu_name = "nuebar_x_";}    
   if (initial_nuPDG == 16) {initial_nu_name = "nutau_x_";}    
   if (initial_nuPDG == -16) {initial_nu_name = "nutaubar_x_";}    

   if (final_nuPDG == 14) {final_nu_name = "numu";}    
   if (final_nuPDG == -14) {final_nu_name = "numubar";}    
   if (final_nuPDG == 12) {final_nu_name = "nue";}    
   if (final_nuPDG == -12) {final_nu_name = "nuebar";}    
   if (final_nuPDG == 16) {final_nu_name = "nutau";}    
   if (final_nuPDG == -16) {final_nu_name = "nutaubar";}    

//   outfile_name = infile.substr(0, 13) + "_" + initial_nu_name + final_nu_name + selec_name + ".root";
   outfile_name = infile.substr(0, 13) + "_" + initial_nu_name + final_nu_name + ".root";  

   char filepath[200];
   char outfilepath[200];
   sprintf(filepath, "%s/%s", path, infile.c_str());
   sprintf(outfilepath, "%s/%s", outpath, outfile_name.c_str());


   TFile * tf = new TFile(filepath, "READ");
   TFile * outf = new TFile(outfilepath, "RECREATE");     

   TTree * tree = NULL;
   tf->GetObject("cafTree", tree);
//   TTree * tree = (TTree*) tf->GetObject("cafTree", tree);
//   TTree * gtree = (TTree*) tf->Get( "genieEvt" );
   TTree * gtree = NULL;
   tf->GetObject("genieEvt", gtree);
   TTree * newtree = tree->CloneTree(0);
   TTree * gnewtree = gtree->CloneTree(100);
   std::cout<<"cloned"<<std::endl;
   TBranch * branch = tree ->GetBranch("rec");
//   tree->SetBranchStatus("*", 1);
   branch ->SetAddress(&sr); 
//   tree->SetBranchStatus("rec", 1);
//   tree->SetBranchAddress("rec", &sr, &branch);
   gtree->SetBranchStatus("*",1);
   TBranch * gbranch = gtree->GetBranch("genie_record");
//   tree->SetBranchAddress("nuPDGunosc",&nutypeunosc);
//   tree->SetBranchAddress("nuPDG",&nutype);
//   tree->SetBranchAddress(cvn_name.c_str(),&cvn);


//   TLeaf* leafptr = tree->GetBranch("common.ixn.dlp.nuhyp.cvn.numu")->GetLeaf("common.ixn.dlp.nuhyp.cvn.numu");
//   TLeaf* leafptr2 = tree->GetBranch("common.ixn.dlp.nuhyp.cvn.nue")->GetLeaf("common.ixn.dlp.nuhyp.cvn.nue");
//   TLeaf* leafptr3 = tree->GetBranch("mc.nu.pdg")->GetLeaf("mc.nu.pdg");
//   TLeaf* leafptr4 = tree->GetBranch("mc.nu.pdgorig")->GetLeaf("mc.nu.pdgorig");
//
   int N = tree->GetEntries();
   int GN = gtree ->GetEntries();
   std::cout << "N = " << N <<" GN = "<<GN<< std::endl; 
   for (int i = 0; i < N; i++) {
//     tree->SetBranchStatus("*", 1);
//     muonparticle.push_back((float)(sr->mc.nu[0].prim[0].p.E));
//     std::cout<<"muon energy "<< muonparticle[i] << std::endl;
//   tree->SetBranchStatus("rec", 1);
//     tree->SetBranchAddress("rec", &sr);
//     std::cout<<"start loop"<<std::endl;

     tree->GetEntry(i);
//     gtree->GetEntry(i);
//     std::cout<<"get entries"<<std::endl;
     // if( i % 100 == 0 ) printf( "Event %d of %d...\n", i, N );
/*     if(numuselec) {
       cvn_name = "cvnnumu";
//       cvn = (double)(leafptr->GetValue(i));
       cvn = (double)(sr->common.ixn.dlp[0].nuhyp.cvn.numu);
//       std::cout<<"cvn "<<cvn<<std::endl;
       cvn_threshold = -0.5; 
       selec_name = "numuselec";}
     else {
       cvn_name = "cvnnue";
       cvn = (double)(sr->common.ixn.dlp[0].nuhyp.cvn.nue);
//       cvn = (double)(leafptr2->GetValue(i));
//       std::cout<<"cvn "<<cvn<<std::endl;
       cvn_threshold = -0.85; 
       selec_name = "nueselec";}
*/
     if(std::isnan(cvn)){cvn =0;}
     std::cout<<"pdgorig "<<(double)(sr->mc.nu[0].pdgorig)<<" pdg "<<(double)(sr->mc.nu[0].pdg)<<" cvn "<<cvn<< " cvn threshold " << cvn_threshold <<std::endl;
//     std::cout<<"pdgorig "<<(double)(leafptr4->GetValue(i))<<" pdg "<<(double)(leafptr3->GetValue(i))<<" cvn "<<cvn<< " cvn threshold " << cvn_threshold <<std::endl;
     if((sr->mc.nu[0].pdgorig) == initial_nuPDG and (sr->mc.nu[0].pdg) == final_nuPDG and cvn >= cvn_threshold) {
//     if((double)(leafptr4->GetValue(i)) == initial_nuPDG and (double)(leafptr3->GetValue(i))== final_nuPDG and cvn > cvn_threshold) {
       std::cout<<"here"<<std::endl;
       tree->GetEntry(i);
       std::cout<<"here1"<<std::endl;
       newtree->Fill();
       std::cout<<"here2"<<std::endl;
//       gbranch->GetEntry(i);
//       gtree->GetEntry(i);
//       std::cout<<"here3"<<std::endl;
//       gnewtree->Fill();
//       std::cout<<"here4"<<std::endl;
     }  
   }
   std::cout<<"here3"<<std::endl;
   gtree->GetEntry(0);
   gnewtree->Fill(); 
   std::cout<<"here4"<<std::endl;

   int n_selection =  newtree->GetEntries();
   int gn_selection =  gnewtree->GetEntries();
   std::cout << "Check: number of events in " << outfile_name << " is: " << n_selection << std::endl;
   std::cout << "Check: number of events in GENIE TREE is: " << gn_selection << std::endl;
   
   if (n_selection != 0) {
     std::cout<<"here2"<<std::endl;
     newtree->Write("cafTree");
     gnewtree->Write("genieEvt");
     std::cout << "Writing output ROOT File" << std::endl;
     outf->Write();
     outf->Close();
     std::cout<<"here2"<<std::endl;
   }
   
   else {
     std::cout << "No events inside " << outfile_name << " , so no output ROOT file " << std::endl;
     remove(outfilepath);
   }
  tree->ResetBranchAddresses();
  delete tf;
  //delete tree;
  //delete branch;
  delete sr;
}


void Plot(std::string infile) 
{  

   const char* path  = "/vols/dune/nk3717/data/NDGAr_newtestCAFs";
//   const char* outpath  = "/vols/dune/nk3717/MaCh3_refactortry/MaCh3_DUNE/inputs/DUNE_NDGAr_CAF_files/";
   
   caf::StandardRecord * sr = new caf::StandardRecord();
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


int main( int argc, char const *argv[] )
{

  if( (argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1])) ) {
    std::cout << "Help yourself by looking at the source code to see what the options are." << std::endl;
    return 0;
  }

  std::string infile;
  infile = argv[1];
//  Plot(infile);
 
  // NUMU_X_NUMU
  Split(infile, 14, 14, true); //numu_x_numu_numuselec
//  Split(infile, 14, 14, false); //numu_x_numu_nueselec
  Split(infile, -14, -14, true); //numubar_x_numubar_numuselec
//  Split(infile, -14, -14, false); //numubar_x_numubar_nueselec
  
  // NUMU_X_NUE
  Split(infile, 14, 12, true); //numu_x_nue_numuselec
//  Split(infile, 14, 12, false); //numu_x_nue_nueselec
  Split(infile, -14, -12, true); //numubar_x_nuebar_numuselec
//  Split(infile, -14, -12, false); //numubar_x_nuebar_nueselec
  
  // NUMU_X_NUTAU
  Split(infile, 14, 16, true); //numu_x_nutau_numuselec
//  Split(infile, 14, 16, false); //numu_x_nutau_nueselec
  Split(infile, -14, -16, true); //numubar_x_nutaubar_numuselec
//  Split(infile, -14, -16, false); //numubar_x_nutaubar_nueselec
  
  // NUE_X_NUE
  Split(infile, 12, 12, true); //nue_x_nue_numuselec 
//  Split(infile, 12, 12, false); //nue_x_nue_nueselec
  Split(infile, -12, -12, true); //nuebar_x_nuebar_numuselec
//  Split(infile, -12, -12, false); //nuebar_x_nuebar_nueselec
  
  // NUE_X_NUMU
  Split(infile, 12, 14, true); //nue_x_numu_numuselec
//  Split(infile, 12, 14, false); //nue_x_numu_nueselec
  Split(infile, -12, -14, true); //nuebar_x_numubar_numuselec
//  Split(infile, -12, -14, false); //nuebar_x_numubar_nueselec
  
  //NUE_X_NUTAU
  Split(infile, 12, 16, true); //nue_x_nutau_numuselec
//  Split(infile, 12, 16, false); //nue_x_nutau_nueselec
  Split(infile, -12, -16, true); //nuebar_x_nutaubar_numuselec
//  Split(infile, -12, -16, false); //nuebar_x_nutaubar_nueselec

}




  
