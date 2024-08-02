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
    keyname.ReplaceAll("_par_0","_0 Total Energy Scale FD");
  else if (keyname.EndsWith("_par_1"))
    keyname.ReplaceAll("_par_1","_1 Total Energy Scale Sqrt FD");
  else if (keyname.EndsWith("_par_2"))
    keyname.ReplaceAll("_par_2","_2 Total Energy Scale Inverse Sqrt FD");
  else if (keyname.EndsWith("_par_3"))
    keyname.ReplaceAll("_par_3","_3 Charged Hadron Energy Scale FD");
  else if (keyname.EndsWith("_par_4"))
    keyname.ReplaceAll("_par_4","_4 Charged Hadron Energy Scale Sqrt FD");
  else if (keyname.EndsWith("_par_5"))
    keyname.ReplaceAll("_par_5","_5 Charged Hadron Energy Scale Inverse Sqrt FD");
  else if (keyname.EndsWith("_par_6"))
    keyname.ReplaceAll("_par_6","_6 Muon Energy Scale FD");
  else if (keyname.EndsWith("_par_7"))
    keyname.ReplaceAll("_par_7","_7 Muon Energy Scale Sqrt FD");
  else if (keyname.EndsWith("_par_8"))
    keyname.ReplaceAll("_par_8","_8 Muon Energy Scale Inverse Sqrt FD");
  else if (keyname.EndsWith("_par_9"))
    keyname.ReplaceAll("_par_9","_9 Neutron Energy Scale FD");
  else if (keyname.EndsWith("_par_10"))
    keyname.ReplaceAll("_par_10","_10 Neutron Energy Scale Sqrt FD");
  else if (keyname.EndsWith("_par_11"))
    keyname.ReplaceAll("_par_11","_11 Neutron Energy Scale Inverse Sqrt FD");
  else if (keyname.EndsWith("_par_12"))
    keyname.ReplaceAll("_par_12 ","_12 EM Shower Energy Scale FD");
  else if (keyname.EndsWith("_par_13"))
    keyname.ReplaceAll("_par_13","_13 EM Shower Energy Scale Sqrt FD");
  else if (keyname.EndsWith("_par_14"))
    keyname.ReplaceAll("_par_14","_14 EM Shower Energy Scale Inverse Sqrt FD");
  else if (keyname.EndsWith("_par_15"))
    keyname.ReplaceAll("_par_15","_15 Charged Hadron Resolution FD");
  else if (keyname.EndsWith("_par_16"))
    keyname.ReplaceAll("_par_16","_16 Muon Resolution FD");
  else if (keyname.EndsWith("_par_17"))
    keyname.ReplaceAll("_par_17","_17 Neutron Resolution FD");
  else if (keyname.EndsWith("_par_18"))
    keyname.ReplaceAll("_par_18","_18 EM Shower Resolution FD");
  else if (keyname.EndsWith("_par_19"))
    keyname.ReplaceAll("_par_19","_19 MaQE");
  else if (keyname.EndsWith("_par_20"))
    keyname.ReplaceAll("_par_20","_20 VecFFCCQEshape");
  else if (keyname.EndsWith("_par_21"))
    keyname.ReplaceAll("_par_21","_21 CCQEPauliSupViaKF");
  else if (keyname.EndsWith("_par_22"))
    keyname.ReplaceAll("_par_22","_22 MaCCRES");
  else if (keyname.EndsWith("_par_23"))
    keyname.ReplaceAll("_par_23","_23 MvCCRES");
  else if (keyname.EndsWith("_par_24"))
    keyname.ReplaceAll("_par_24","_24 MaNCRES");
  else if (keyname.EndsWith("_par_25"))
    keyname.ReplaceAll("_par_25","_25 MvNCRES");
  else if (keyname.EndsWith("_par_26"))
    keyname.ReplaceAll("_par_26","_26 ThetaDelta2NPi");
  else if (keyname.EndsWith("_par_27"))
    keyname.ReplaceAll("_par_27","_27 AhtBY");
  else if (keyname.EndsWith("_par_28"))
    keyname.ReplaceAll("_par_28","_28 BhtBY");
  else if (keyname.EndsWith("_par_29"))
    keyname.ReplaceAll("_par_29","_29 CV1uBY");
  else if (keyname.EndsWith("_par_30"))
    keyname.ReplaceAll("_par_30","_30 CV2uBY");
  else if (keyname.EndsWith("_par_31"))
    keyname.ReplaceAll("_par_31 ","_31 FrCEx_pi");
  else if (keyname.EndsWith("_par_32"))
    keyname.ReplaceAll("_par_32","_32 FrElas_pi");
  else if (keyname.EndsWith("_par_33"))
    keyname.ReplaceAll("_par_33","_33 FrInel_pi");
  else if (keyname.EndsWith("_par_34"))
    keyname.ReplaceAll("_par_34","_34 FrAbs_pi");
  else if (keyname.EndsWith("_par_35"))
    keyname.ReplaceAll("_par_35","_35 FrPiProd_pi");
  else if (keyname.EndsWith("_par_36"))
    keyname.ReplaceAll("_par_36","_36 FrCEx_N");
  else if (keyname.EndsWith("_par_37"))
    keyname.ReplaceAll("_par_37","_37 FrElas_N");
  else if (keyname.EndsWith("_par_38"))
    keyname.ReplaceAll("_par_38","_38 FrInel_N");
  else if (keyname.EndsWith("_par_39"))
    keyname.ReplaceAll("_par_39","_39 FrAbs_N");
  else if (keyname.EndsWith("_par_40"))
    keyname.ReplaceAll("_par_40","_40 FrPiProd_N");
  else if (keyname.EndsWith("_par_41"))
    keyname.ReplaceAll("_par_41","_41 E2p2h_A_nu");
  else if (keyname.EndsWith("_par_42"))
    keyname.ReplaceAll("_par_42","_42 E2p2h_B_nu");
  else if (keyname.EndsWith("_par_43"))
    keyname.ReplaceAll("_par_43","_43 E2p2h_A_nubar");
  else if (keyname.EndsWith("_par_44"))
    keyname.ReplaceAll("_par_44","_44 E2p2h_B_nubar");
  else if (keyname.EndsWith("_par_45"))
    keyname.ReplaceAll("_par_45","_45 C12ToAr40_2p2hScaling_nu");
  else if (keyname.EndsWith("_par_46"))
    keyname.ReplaceAll("_par_46","_46 C12ToAr40_2p2hScaling_nubar");
  else if (keyname.EndsWith("_par_47"))
    keyname.ReplaceAll("_par_47","_47 NR_nu_n_CC_2pi");
  else if (keyname.EndsWith("_par_48"))
    keyname.ReplaceAll("_par_48","_48 NR_nu_n_CC_3pi");
  else if (keyname.EndsWith("_par_49"))
    keyname.ReplaceAll("_par_49","_49 NR_nu_p_CC_2pi");
  else if (keyname.EndsWith("_par_50"))
    keyname.ReplaceAll("_par_50","_50 NR_nu_p_CC_3pi");
  else if (keyname.EndsWith("_par_51"))
    keyname.ReplaceAll("_par_51","_51 NR_nu_np_CC_1pi");
  else if (keyname.EndsWith("_par_52"))
    keyname.ReplaceAll("_par_52","_52 NR_nu_n_NC_1pi");
  else if (keyname.EndsWith("_par_53"))
    keyname.ReplaceAll("_par_53","_53 NR_nu_n_NC_2pi");
  else if (keyname.EndsWith("_par_54"))
    keyname.ReplaceAll("_par_54","_54 NR_nu_n_NC_3pi");
  else if (keyname.EndsWith("_par_55"))
    keyname.ReplaceAll("_par_55","_55 NR_nu_p_NC_1pi");
  else if (keyname.EndsWith("_par_56"))
    keyname.ReplaceAll("_par_56","_56 NR_nu_p_NC_2pi");
  else if (keyname.EndsWith("_par_57"))
    keyname.ReplaceAll("_par_57","_57 NR_nu_p_NC_3pi");
  else if (keyname.EndsWith("_par_58"))
    keyname.ReplaceAll("_par_58","_58 NR_nubar_n_CC_1pi");
  else if (keyname.EndsWith("_par_59"))
    keyname.ReplaceAll("_par_59","_59 NR_nubar_n_CC_2pi");
  else if (keyname.EndsWith("_par_60"))
    keyname.ReplaceAll("_par_60","_60 NR_nubar_n_CC_3pi");
  else if (keyname.EndsWith("_par_61"))
    keyname.ReplaceAll("_par_61","_61 NR_nubar_p_CC_1pi");
  else if (keyname.EndsWith("_par_62"))
    keyname.ReplaceAll("_par_62","_62 NR_nubar_p_CC_2pi");
  else if (keyname.EndsWith("_par_63"))
    keyname.ReplaceAll("_par_63","_63 NR_nubar_p_CC_3pi");
  else if (keyname.EndsWith("_par_64"))
    keyname.ReplaceAll("_par_64","_64 NR_nubar_n_NC_1pi");
  else if (keyname.EndsWith("_par_65"))
    keyname.ReplaceAll("_par_65","_65 NR_nubar_n_NC_2pi");
  else if (keyname.EndsWith("_par_66"))
    keyname.ReplaceAll("_par_66","_66 NR_nubar_n_NC_3pi");
  else if (keyname.EndsWith("_par_67"))
    keyname.ReplaceAll("_par_67","_67 NR_nubar_p_NC_1pi");
  else if (keyname.EndsWith("_par_68"))
    keyname.ReplaceAll("_par_68","_68 NR_nubar_p_NC_2pi");
  else if (keyname.EndsWith("_par_69"))
    keyname.ReplaceAll("_par_69","_69 NR_nubar_p_NC_3pi");
  else if (keyname.EndsWith("_par_70"))
    keyname.ReplaceAll("_par_70","_70 BeRPA_A");
  else if (keyname.EndsWith("_par_71"))
    keyname.ReplaceAll("_par_71","_71 BeRPA_B");
  else if (keyname.EndsWith("_par_72"))
    keyname.ReplaceAll("_par_72","_72 BeRPA_D");
  else if (keyname.EndsWith("_par_73"))
    keyname.ReplaceAll("_par_73","_73 nuenuebar_xsec_ratio");
  else if (keyname.EndsWith("_par_74"))
    keyname.ReplaceAll("_par_74","_74 nuemumu_xsec_ratio");

  return keyname;
}

void makeSigmaVarPlotsDUNE(TString inputfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
   TFile* file = new TFile(inputfile);
   TList* list = file->GetListOfKeys();
   

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   TCanvas* c0 = new TCanvas("c0","c0",0,0,700,900);
   c0->Divide(1,2);
   c0->Print("sigmavarnewreco.ps[");
   std::vector<std::vector<TString>> kinematicvariables;

  std::vector<std::string> plotvariables = {"RecoNeutrinoEnergy", "TrueNeutrinoEnergy", "TrueMinusRecoEnergy", "TrueMinusRecoEnergyRatio", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "RecoXPos", "TrueYPos", "RecoYPos", "TrueZPos", "RecoZPos", "TrueRad", "RecoRad", "NTrueMuons", "NRecoMuons", "RecoLepEnergy", "TrueLepEnergy"};

//   std::vector<std::string> plotvariables = {"TrueNeutrinoEnergy", "RecoNeutrinoEnergy", "TrueMinusRecoEnergy", "PionMultiplicity", "NRecoParticles", "InFDV", "TrueXPos", "TrueYPos", "TrueZPos", "NMuons", "TrueMinusRecoEnergyRatio", "RecoLepEnergy"};
   int n_vars = plotvariables.size();
   for(int iplots =0; iplots<n_vars; iplots++){
      std::vector<TString> keynames;
      kinematicvariables.push_back(keynames);
   }
//   std::vector<TString> trueEkeynames;
//   std::vector<TString> recoEkeynames;
//   std::vector<TString> truerecoEkeynames;
//   std::vector<TString> pionkeynames;
//   std::vector<TString> nrecoparticleskeynames;
//   std::vector<TString> infdvkeynames;
//   std::vector<TString> truexposkeynames;
//   std::vector<TString> trueyposkeynames;
//   std::vector<TString> truezposkeynames;
   for(int i=0; i<list->GetEntries(); i++)
   {
     bool boolcontinue = false;
     //std::cout << "entry " << i << std::endl;
     TString keyname = list->At(i)->GetName();
        if(keyname=="FHC_nue_unosc" || keyname=="FHC_numu_unosc" 
	 || keyname=="RHC_nue_unosc" || keyname=="RHC_numu_unosc" || !keyname.Contains("n3") || keyname.Contains("_NDGAr_")){
        for(int iplots=0; iplots<n_vars; iplots++){
          TString keynametest = plotvariables[iplots]+"_NDGAr_";
          if(keyname.Contains(keynametest)){
          kinematicvariables[iplots].push_back(keyname);
          std::cout<<"iplots: "<<iplots<<" keyname: "<<keyname<<std::endl;
          boolcontinue = true;
          break;
          }
        }
/*
        if(keyname.Contains("RecoNeutrinoEnergy_NDGAr")){
        recoEkeynames.push_back(keyname);
        }
        if(keyname.Contains("TrueMinusRecoEnergy_NDGAr")){
        truerecoEkeynames.push_back(keyname);
        }
        if(keyname.Contains("PionMultiplicity_NDGAr")){
        pionkeynames.push_back(keyname);
        }
        if(keyname.Contains("NRecoParticles_NDGAr")){
        nrecoparticleskeynames.push_back(keyname);
        }
        if(keyname.Contains("InFDV_NDGAr")){
        infdvkeynames.push_back(keyname);
        }
        if(keyname.Contains("TrueXPos_NDGAr")){
        truexposkeynames.push_back(keyname);
        }
        if(keyname.Contains("TrueYPos_NDGAr")){
        trueyposkeynames.push_back(keyname);
        }
        if(keyname.Contains("TrueZPos_NDGAr")){
        truezposkeynames.push_back(keyname);
        }
        */
        if(!boolcontinue){continue;}
        /*
        stack->SetTitle("TrueNeutrinoEnergy_NDGAr ;True Energy (GeV);");
        TH1D* hist =(TH1D*)file->Get(keyname);
        hist->SetFillColor(icolour);
        icolour++;
        c0->cd(1);
        //c0->Update();
        stack->Add(hist);
        std::cout<<"Added new hist to stack"<<std::endl;
        }
        else{continue;}
        stack->Draw();
        c0->Update();
        c0->Print("sigmavarnewreco.ps");
	 continue;*/
      }
      /*else{
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
      if(keyname.Contains("FHC_numu")) {
        TH1D* fhc_numu_nom = (TH1D*)file->Get("FHC_numu_unosc");
        fhc_numu_nom->SetLineColor(kBlack);
        fhc_numu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
	    divisor=fhc_numu_nom;
	 //divisor->GetXaxis()->SetRange(1,70);
	 //for(int j=0; j<4; j++)
	    //hVar[j]->GetXaxis()->SetRange(1,70);
      }
      if(keyname.Contains("FHC_nue")) {
        TH1D* fhc_nue_nom = (TH1D*)file->Get("FHC_nue_unosc");
        fhc_nue_nom->SetLineColor(kBlack);
        fhc_nue_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
	    divisor=fhc_nue_nom;
	  }
      if(keyname.Contains("RHC_nue")) {
        TH1D* rhc_nue_nom = (TH1D*)file->Get("RHC_nue_unosc");
        rhc_nue_nom->SetLineColor(kBlack);
        rhc_nue_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
	    divisor=rhc_nue_nom;
	  }
      if(keyname.Contains("RHC_numu")) {
       TH1D* rhc_numu_nom = (TH1D*)file->Get("RHC_numu_unosc");
       rhc_numu_nom->SetLineColor(kBlack);
       rhc_numu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
	   divisor=rhc_numu_nom;
	  }
      if(keyname.Contains("ND_FHC_CCnumu")) {
       TH1D* nd_fhc_ccnumu_nom = (TH1D*)file->Get("ND_FHC_CCnumu_unosc");
       nd_fhc_ccnumu_nom->SetLineColor(kBlack);
       nd_fhc_ccnumu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
	   divisor=nd_fhc_ccnumu_nom;
	  }
      if(keyname.Contains("ND_RHC_CCnumu")) {
       TH1D* nd_rhc_ccnumu_nom = (TH1D*)file->Get("ND_RHC_CCnumu_unosc");
       nd_rhc_ccnumu_nom->SetLineColor(kBlack);
       nd_rhc_ccnumu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
	   divisor=nd_rhc_ccnumu_nom;
	  }
      if(keyname.Contains("NDGAr_FHC_")) {
       TH1D* ndgar_fhc_ccnumu_nom = (TH1D*)file->Get("NDGAr_FHC_unosc");
       ndgar_fhc_ccnumu_nom->SetLineColor(kBlack);
       ndgar_fhc_ccnumu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
	   divisor=ndgar_fhc_ccnumu_nom;
	  }
      if(keyname.Contains("NDGAr_RHC_")) {
       TH1D* ndgar_rhc_ccnumu_nom = (TH1D*)file->Get("NDGAr_RHC_unosc");
       ndgar_rhc_ccnumu_nom->SetLineColor(kBlack);
       ndgar_rhc_ccnumu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
	   divisor=ndgar_rhc_ccnumu_nom;
	  }
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
//      maxlist.push_back(divisor->GetMaximum());
      for(int j=0; j<4; j++) {
         maxlist.push_back(hVar[j]->GetMaximum());
         std::cout<<"maxlist ["<<j<<"]:"<<maxlist[j]<<std::endl; }
      double maxmax = *std::max_element(maxlist.begin(), maxlist.end());
      std::cout<<"maxmax: "<<maxmax<<std::endl;
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
      }
   c0->Print("sigmavarnewreco.ps");*/
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

//   kinematicvariables.push_back(trueEkeynames);
//   kinematicvariables.push_back(recoEkeynames);
//   kinematicvariables.push_back(truerecoEkeynames);
//   kinematicvariables.push_back(pionkeynames);
//   kinematicvariables.push_back(nrecoparticleskeynames);
//   kinematicvariables.push_back(infdvkeynames);
//   kinematicvariables.push_back(truexposkeynames);
//   kinematicvariables.push_back(trueyposkeynames);
//   kinematicvariables.push_back(truezposkeynames);


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
     c0->Print("sigmavarnewreco.ps");
     delete stack;
     delete leg1;
     Hists.clear();
   }
/*
   THStack* stack1 = new THStack("stack1","");
   TLegend* leg2 = new TLegend(0.65,0.55,0.8,0.8);
   int icolour = 1;
   for(int ikeynames = 0; ikeynames<recoEkeynames.size(); ikeynames++){
     stack1->SetTitle("RecoNeutrinoEnergy_NDGAr ;Reco Energy (GeV);");
     TH1D* hist1 =(TH1D*)file->Get(recoEkeynames[ikeynames]);
     if(hist1->Integral() != 0){
       hist1->SetFillColor(icolour);
       icolour++;
       if(icolour ==10){icolour++;}
       stack1->Add(hist1);
       leg2->AddEntry(hist1,recoEkeynames[ikeynames].Remove(0, 25));
     }
   }
   c0->cd(1);
   stack1->Draw("HIST");
   leg2->Draw("SAME"); 
   c0->Update();
   c0->Print("sigmavarnewreco.ps");
   //delete stack;
   delete stack1;*/
   std::cout<<"HERE"<<std::endl;
   c0->Print("sigmavarnewreco.ps]");

}
