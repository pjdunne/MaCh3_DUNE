////////////////////////////////////////////////////////////////////////////
//
// Quick plotting script to show the affect of marginalising 2D distributions 
// Give this exe a reduced chain and it'll draw th23 vs. dm32 and show the 1D
// and 2D best-fit points and make panels either side of the main canvas
// showing the 1D distributions 
//  
/////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMarker.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TROOT.h"
#include "TColor.h"
#include "TImage.h"

void draw_LLH_scans(TString dir, TString infile){

  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPalette(51,0);
  gStyle->SetNumberContours(104);

  TFile *f = new TFile(dir+infile, "READ");

  TH2D* dis_llh = 0;
  TH2D* dis_llh_cont = 0;
  //TGraph* dis_true = 0;
  //TGraph* dis_bestfit = 0;

  //TH2D* app_llh = 0;
  //TH2D* app_llh_cont = 0;
  //TGraph* app_true = 0;
  //TGraph* app_bestfit = 0;

  f->GetObject("dCP Th23 LLH", dis_llh);
  //f->GetObject("true", dis_true);
  //f->GetObject("bestfit", dis_bestfit);

  //f->GetObject("llh_scan_dcpth13", app_llh);
  //f->GetObject("true_dcpth13", app_true);
  //f->GetObject("bestfit_dcpth13", app_bestfit);

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
  c1->cd();
  
  dis_llh->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  dis_llh->GetXaxis()->SetLabelSize(0.05);
  dis_llh->GetYaxis()->SetLabelSize(0.05);
  dis_llh->GetXaxis()->SetTitleSize(0.06);
  dis_llh->GetXaxis()->SetTitleOffset(0.9);
  dis_llh->GetYaxis()->SetTitle("dCP");
  dis_llh->GetYaxis()->SetTitleSize(0.06);
  dis_llh->GetYaxis()->SetTitleOffset(0.8);
  dis_llh->GetYaxis()->SetTickLength(0);
  dis_llh->GetXaxis()->SetTickLength(0);
  dis_llh->GetYaxis()->SetMaxDigits(3);

  dis_llh->GetZaxis()->SetTitle("-LLH");
  dis_llh->GetZaxis()->SetTitleSize(0.06);
  dis_llh->GetZaxis()->SetLabelSize(0.04);
  dis_llh->GetZaxis()->SetTitleOffset(0.8);
  dis_llh->Draw("colz");
  //gStyle->SetNumberContours(25);
  //Clone histograms so we can draw nice contours
  std::cout << "Maximum is " << dis_llh->GetMaximum() << std::endl;
  std::cout << "Minimum is " << dis_llh->GetMinimum() << std::endl;
  double llh_min = dis_llh->GetMinimum();
  double contours_dis[10] = {80, 160, 240, 320, 400, 480, 560, 640, 720, 800};
  //double contours_dis[10] = {llh_min, llh_min+1, llh_min+4, llh_min+9, llh_min+16, llh_min+25, llh_min+36, llh_min+49, llh_min+64, llh_min+81};
  dis_llh_cont = (TH2D*)dis_llh->Clone("dis_llh_cont");
  dis_llh_cont->SetContour(10, contours_dis);
  dis_llh_cont->SetLineColor(kGray+2);
  dis_llh_cont->SetLineWidth(2);
  dis_llh_cont->SetLineStyle(1);
  dis_llh_cont->Draw("CONT3 SAMES");

  //Now draw true and best-fit points
  //dis_true->SetMarkerStyle(29);
  //dis_true->SetMarkerColor(kGray+2);
  //dis_true->SetMarkerSize(2);
  //dis_true->Draw("P SAMES");

  //dis_bestfit->SetMarkerStyle(29);
  //dis_bestfit->SetMarkerSize(2);
  //dis_bestfit->SetMarkerColor(kRed+1);
  //dis_bestfit->Draw("P SAMES");

  //TLegend *leg_dis = new TLegend(0.48, 0.48, 0.85, 0.59);
  //leg_dis->AddEntry(dis_true,    Form("True values:      sin^{2}#theta_{23} = %1.3f, #Delta m^{2}_{32} = %2.3e", dis_true->GetPointX(0), dis_true->GetPointY(0)), "p");
  //leg_dis->AddEntry(dis_bestfit, Form("Best-fit values:  sin^{2}#theta_{23} = %1.3f, #Delta m^{2}_{32} = %2.3e", dis_bestfit->GetPointX(0), dis_bestfit->GetPointY(0)), "p");
  //leg_dis->Draw();

  c1->Print("Test_LLH_dis.pdf");
  c1->Print("Test_LLH_dis.png");

  //////////////////////////////////////// 
  /*

  app_llh->GetXaxis()->SetTitle("sin^{2}#theta_{13}");
  app_llh->GetXaxis()->SetLabelSize(0.05);
  app_llh->GetYaxis()->SetLabelSize(0.05);
  app_llh->GetXaxis()->SetTitleSize(0.06);
  app_llh->GetXaxis()->SetTitleOffset(0.9);
  app_llh->GetYaxis()->SetTitle("#delta_{CP} [radians]");
  app_llh->GetYaxis()->SetTitleSize(0.06);
  app_llh->GetYaxis()->SetTitleOffset(0.8);
  app_llh->GetYaxis()->SetTickLength(0);
  app_llh->GetXaxis()->SetTickLength(0);
  app_llh->GetYaxis()->SetMaxDigits(3);

  app_llh->GetZaxis()->SetTitle("-LLH");
  app_llh->GetZaxis()->SetTitleSize(0.06);
  app_llh->GetZaxis()->SetLabelSize(0.04);
  app_llh->GetZaxis()->SetTitleOffset(0.8);
  app_llh->Draw("colz");
  //gStyle->SetNumberContours(25);
  //Clone histograms so we can draw nice contours
  std::cout << "Maximum is " << app_llh->GetMaximum() << std::endl;
  double contours_app[10] = {3, 6, 9, 12, 15, 18, 21, 24, 27, 30};
  app_llh_cont = (TH2D*)app_llh->Clone("app_llh_cont");
  app_llh_cont->SetContour(10, contours_app);
  app_llh_cont->SetLineColor(kGray+2);
  app_llh_cont->SetLineWidth(2);
  app_llh_cont->SetLineStyle(0);
  app_llh_cont->Draw("CONT3 SAMES");

  //Now draw true and best-fit points
  app_true->SetMarkerStyle(29);
  app_true->SetMarkerColor(kGray+2);
  app_true->SetMarkerSize(2);
  app_true->Draw("P SAMES");

  app_bestfit->SetMarkerStyle(29);
  app_bestfit->SetMarkerSize(2);
  app_bestfit->SetMarkerColor(kRed+1);
  app_bestfit->Draw("P SAMES");

  TLegend *leg_app = new TLegend(0.45, 0.82, 0.85, 0.94);
  leg_app->AddEntry(app_true,    Form("True values:      sin^{2}#theta_{13} = %2.3e, #delta_{CP} = %1.3f", app_true->GetPointX(0), app_true->GetPointY(0)), "p");
  leg_app->AddEntry(app_bestfit, Form("Best-fit values:  sin^{2}#theta_{13} = %2.3e, #delta_{CP} = %1.3f", app_bestfit->GetPointX(0), app_bestfit->GetPointY(0)), "p");
  leg_app->Draw(); 

  c1->Print("Test_LLH_app.pdf"); */

  return;

}
