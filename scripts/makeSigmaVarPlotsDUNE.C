#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>

TString getXsecParName(TString keyname)
{
  if (keyname.EndsWith("_par_0"))
    keyname.ReplaceAll("_par_0","_0 MaQE");
  else if (keyname.EndsWith("_par_1"))
    keyname.ReplaceAll("_par_1","_1 VecFFCCQEshape");
  else if (keyname.EndsWith("_par_2"))
    keyname.ReplaceAll("_par_2","_2 CCQEPauliSupViaKF");
  else if (keyname.EndsWith("_par_3"))
    keyname.ReplaceAll("_par_3","_3 MaCCRES");
  else if (keyname.EndsWith("_par_4"))
    keyname.ReplaceAll("_par_4","_4 MvCCRES");
  else if (keyname.EndsWith("_par_5"))
    keyname.ReplaceAll("_par_5","_5 MaNCRES");
  else if (keyname.EndsWith("_par_6"))
    keyname.ReplaceAll("_par_6","_6 MvNCRES");
  else if (keyname.EndsWith("_par_7"))
    keyname.ReplaceAll("_par_7","_7 ThetaDelta2NPi");
  else if (keyname.EndsWith("_par_8"))
    keyname.ReplaceAll("_par_8","_8 AhtBY");
  else if (keyname.EndsWith("_par_9"))
    keyname.ReplaceAll("_par_9","_9 BhtBY");
  else if (keyname.EndsWith("_par_10"))
    keyname.ReplaceAll("_par_10","_10 CV1uBY");
  else if (keyname.EndsWith("_par_11"))
    keyname.ReplaceAll("_par_11","_11 CV2uBY");
  else if (keyname.EndsWith("_par_12"))
    keyname.ReplaceAll("_par_12 ","_12 FrCEx_pi");
  else if (keyname.EndsWith("_par_13"))
    keyname.ReplaceAll("_par_13","_13 FrElas_pi");
  else if (keyname.EndsWith("_par_14"))
    keyname.ReplaceAll("_par_14","_14 FrInel_pi");
  else if (keyname.EndsWith("_par_15"))
    keyname.ReplaceAll("_par_15","_15 FrAbs_pi");
  else if (keyname.EndsWith("_par_16"))
    keyname.ReplaceAll("_par_16","_16 FrPiProd_pi");
  else if (keyname.EndsWith("_par_17"))
    keyname.ReplaceAll("_par_17","_17 FrCEx_N");
  else if (keyname.EndsWith("_par_18"))
    keyname.ReplaceAll("_par_18","_18 FrElas_N");
  else if (keyname.EndsWith("_par_19"))
    keyname.ReplaceAll("_par_19","_19 FrInel_N");
  else if (keyname.EndsWith("_par_20"))
    keyname.ReplaceAll("_par_20","_20 FrAbs_N");
  else if (keyname.EndsWith("_par_21"))
    keyname.ReplaceAll("_par_21","_21 FrPiProd_N");
  else if (keyname.EndsWith("_par_22"))
    keyname.ReplaceAll("_par_22","_22 E2p2h_A_nu");
  else if (keyname.EndsWith("_par_23"))
    keyname.ReplaceAll("_par_23","_23 E2p2h_B_nu");
  else if (keyname.EndsWith("_par_24"))
    keyname.ReplaceAll("_par_24","_24 E2p2h_A_nubar");
  else if (keyname.EndsWith("_par_25"))
    keyname.ReplaceAll("_par_25","_25 E2p2h_B_nubar");
  else if (keyname.EndsWith("_par_26"))
    keyname.ReplaceAll("_par_26","_26 C12ToAr40_2p2hScaling_nu");
  else if (keyname.EndsWith("_par_27"))
    keyname.ReplaceAll("_par_27","_27 C12ToAr40_2p2hScaling_nubar");
  else if (keyname.EndsWith("_par_28"))
    keyname.ReplaceAll("_par_28","_28 NR_nu_n_CC_2pi");
  else if (keyname.EndsWith("_par_29"))
    keyname.ReplaceAll("_par_29","_29 NR_nu_n_CC_3pi");
  else if (keyname.EndsWith("_par_30"))
    keyname.ReplaceAll("_par_30","_30 NR_nu_p_CC_2pi");
  else if (keyname.EndsWith("_par_31"))
    keyname.ReplaceAll("_par_31","_31 NR_nu_p_CC_3pi");
  else if (keyname.EndsWith("_par_32"))
    keyname.ReplaceAll("_par_32","_32 NR_nu_np_CC_1pi");
  else if (keyname.EndsWith("_par_33"))
    keyname.ReplaceAll("_par_33","_33 NR_nu_n_NC_1pi");
  else if (keyname.EndsWith("_par_34"))
    keyname.ReplaceAll("_par_34","_34 NR_nu_n_NC_2pi");
  else if (keyname.EndsWith("_par_35"))
    keyname.ReplaceAll("_par_35","_35 NR_nu_n_NC_3pi");
  else if (keyname.EndsWith("_par_36"))
    keyname.ReplaceAll("_par_36","_36 NR_nu_p_NC_1pi");
  else if (keyname.EndsWith("_par_37"))
    keyname.ReplaceAll("_par_37","_37 NR_nu_p_NC_2pi");
  else if (keyname.EndsWith("_par_38"))
    keyname.ReplaceAll("_par_38","_38 NR_nu_p_NC_3pi");
  else if (keyname.EndsWith("_par_39"))
    keyname.ReplaceAll("_par_39","_39 NR_nubar_n_CC_1pi");
  else if (keyname.EndsWith("_par_40"))
    keyname.ReplaceAll("_par_40","_40 NR_nubar_n_CC_2pi");
  else if (keyname.EndsWith("_par_41"))
    keyname.ReplaceAll("_par_41","_41 NR_nubar_n_CC_3pi");
  else if (keyname.EndsWith("_par_42"))
    keyname.ReplaceAll("_par_42","_42 NR_nubar_p_CC_1pi");
  else if (keyname.EndsWith("_par_43"))
    keyname.ReplaceAll("_par_43","_43 NR_nubar_p_CC_2pi");
  else if (keyname.EndsWith("_par_44"))
    keyname.ReplaceAll("_par_44","_44 NR_nubar_p_CC_3pi");
  else if (keyname.EndsWith("_par_45"))
    keyname.ReplaceAll("_par_45","_45 NR_nubar_n_NC_1pi");
  else if (keyname.EndsWith("_par_46"))
    keyname.ReplaceAll("_par_46","_46 NR_nubar_n_NC_2pi");
  else if (keyname.EndsWith("_par_47"))
    keyname.ReplaceAll("_par_47","_47 NR_nubar_n_NC_3pi");
  else if (keyname.EndsWith("_par_48"))
    keyname.ReplaceAll("_par_48","_48 NR_nubar_p_NC_1pi");
  else if (keyname.EndsWith("_par_49"))
    keyname.ReplaceAll("_par_49","_49 NR_nubar_p_NC_2pi");
  else if (keyname.EndsWith("_par_50"))
    keyname.ReplaceAll("_par_50","_50 NR_nubar_p_NC_3pi");
  else if (keyname.EndsWith("_par_51"))
    keyname.ReplaceAll("_par_51","_51 BeRPA_A");
  else if (keyname.EndsWith("_par_52"))
    keyname.ReplaceAll("_par_52","_52 BeRPA_B");
  else if (keyname.EndsWith("_par_53"))
    keyname.ReplaceAll("_par_53","_53 BeRPA_D");
  else if (keyname.EndsWith("_par_54"))
    keyname.ReplaceAll("_par_54","_54 nuenuebar_xsec_ratio");
  else if (keyname.EndsWith("_par_55"))
    keyname.ReplaceAll("_par_55","_55 nuemumu_xsec_ratio");

  return keyname;
}

void makeSigmaVarPlotsDUNE(TString inputfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
   TFile* file = new TFile(inputfile);
   TList* list = file->GetListOfKeys();
   
   TH1D* fhc_nue_nom = (TH1D*)file->Get("fhc_nue_osc");
   fhc_nue_nom->SetLineColor(kBlack);
   fhc_nue_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
   TH1D* fhc_numu_nom = (TH1D*)file->Get("fhc_numu_osc");
   fhc_numu_nom->SetLineColor(kBlack);
   fhc_numu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
   TH1D* rhc_nue_nom = (TH1D*)file->Get("rhc_nue_osc");
   rhc_nue_nom->SetLineColor(kBlack);
   rhc_nue_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
   TH1D* rhc_numu_nom = (TH1D*)file->Get("rhc_numu_osc");
   rhc_numu_nom->SetLineColor(kBlack);
   rhc_numu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   TCanvas* c0 = new TCanvas("c0","c0",0,0,700,900);
   c0->Divide(1,2);
   c0->Print("sigmavar.ps[");
  
   for(int i=0; i<list->GetEntries(); i++)
   {
     //std::cout << "entry " << i << std::endl;
     TString keyname = list->At(i)->GetName();
      if(keyname=="fhc_nue_osc" || keyname=="fhc_numu_osc" 
	 || keyname=="rhc_nue_osc" || keyname=="rhc_numu_osc" || !keyname.Contains("n3"))
	 continue;

      TH1D* hVar[4];
      hVar[0]=(TH1D*)file->Get(keyname);
      keyname.ReplaceAll("n3","n1");
      hVar[1]=(TH1D*)file->Get(keyname);
      keyname.ReplaceAll("n1","p1");
      hVar[2]=(TH1D*)file->Get(keyname);
      keyname.ReplaceAll("p1","p3");
      hVar[3]=(TH1D*)file->Get(keyname);

      hVar[0]->SetLineColor(kRed);
      hVar[1]->SetLineColor(kMagenta);
      hVar[2]->SetLineColor(kCyan);
      hVar[3]->SetLineColor(kBlue);
      
      TH1D* divisor;
      if(keyname.Contains("fhc_numu"))
      {
	 divisor=fhc_numu_nom;
	 //divisor->GetXaxis()->SetRange(1,70);
	 //for(int j=0; j<4; j++)
	    //hVar[j]->GetXaxis()->SetRange(1,70);
      }
      if(keyname.Contains("fhc_nue"))
	 divisor=fhc_nue_nom;
      if(keyname.Contains("rhc_nue"))
	 divisor=rhc_nue_nom;
      if(keyname.Contains("rhc_numu"))
	 divisor=rhc_numu_nom;

      keyname.ReplaceAll("_sig_p3","");
      
      // if xsec, find the name of the parameter
      if(keyname.Contains("xsec"))
	keyname = getXsecParName(keyname);
      divisor->SetTitle(keyname);
      if (i == 840) {
      	std::cout << "Title: " << keyname << std::endl; }

      
      c0->cd(1);
      
      // Set the axis range to max value of all variations
      std::vector<double> maxlist;
      maxlist.push_back(divisor->GetMaximum());
      for(int j=0; j<4; j++) {
         maxlist.push_back(hVar[j]->GetMaximum()); }
      double maxmax = *std::max_element(maxlist.begin(), maxlist.end());
      divisor->GetYaxis()->SetRangeUser(0,1.05*maxmax);
      

     divisor->DrawCopy("HIST");
     for(int j=0; j<4; j++) 
	 hVar[j]->DrawCopy("HIST SAME");
     
      c0->Update();
//      getchar();

      double max=1.6, min=0.4;
      for(int j=0; j<4; j++)
      {
	 hVar[j]->Divide(divisor);
	 if(hVar[j]->GetMinimum(0.1)<min)
	    min = hVar[j]->GetMinimum();
	 if(hVar[j]->GetMaximum()>max)
	    max = hVar[j]->GetMaximum();
      }
      if(min==0)
	 min=(2-max);
      max*=1.02;
      min*=0.98;

      TH1D* divcopy = (TH1D*)divisor->Clone("divcopy");
      divcopy->Divide(divisor);
      divcopy->SetMaximum(max);
      divcopy->SetMinimum(min);
      
      c0->cd(2);
      divcopy->Draw("HIST");
      for(int j=0; j<4; j++)
	 hVar[j]->DrawCopy("HIST SAME");

      TLegend* leg = new TLegend(0.65,0.55,0.8,0.8);
      leg->AddEntry(divcopy,"nominal","L");
      leg->AddEntry(hVar[0],"-3#sigma","L");
      leg->AddEntry(hVar[1],"-1#sigma","L");
      leg->AddEntry(hVar[2],"+1#sigma","L");
      leg->AddEntry(hVar[3],"+3#sigma","L");
      leg->SetFillColor(0);
      c0->cd(1);
      leg->Draw("SAME");

      c0->Update();
//      getchar();
      c0->Print("sigmavar.ps");

   }

   c0->Print("sigmavar.ps]");

}
