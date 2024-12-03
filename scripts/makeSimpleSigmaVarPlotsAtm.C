

const std::string outdir = "sigma_variation_plots";

std::vector<std::string> vars = {"delta_cp"};
// std::vector<std::string> vars = {"delta_cp", "sin2th_23", "sin2th_13", "delm2_23"};

std::vector<const char*> dcp_angles = {"0", "#pi/4", "#pi/2", "3#pi/4", "#pi",
        "-3#pi/4", "-#pi/2", "-#pi/4"};

void makeSimpleSigmaVarPlotsAtm(std::string filename) {
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas;
    TFile *f = TFile::Open(filename.c_str());
    for (auto var : vars) {
        for (std::string s : {"numuselec", "nueselec"}) {
            TLegend* leg = new TLegend(0.65,0.55,0.8,0.8);
            gPad->SetLogx();
            if (var == "delta_cp") {
                TH1D *h = (TH1D*) f->Get((var + "/" + s + "/" + "Variation_0").c_str());
                // h->Scale(1.0, "width");
                h->GetXaxis()->SetRangeUser(0, 10);
                h->SetTitle((var + " " + s).c_str());
                h->GetXaxis()->SetTitle("True Energy (GeV)");
                h->SetLineColor(kBlack);
                h->Draw("hist");

                leg->AddEntry(h,dcp_angles[0],"L");
                for (int i = 1; i < 8; i++) {
                    TH1D *ho = (TH1D*) f->Get((var + "/" + s + "/Variation_" + std::to_string(i)).c_str());
                    ho->SetLineColor(kBlack + i);
                    leg->AddEntry(ho,dcp_angles[i],"L");
                //     ho->Scale(1.0, "width");
                    ho->Draw("hist same");
                }
                leg->Draw("same");
                        
                c->SaveAs((outdir + "/" + var + s + ".png").c_str());
                
                h->GetXaxis()->SetRangeUser(0, 10);
                h->GetYaxis()->SetRangeUser(0.98, 1.02);
                h->SetTitle((var + " " + s + " ratios").c_str());
                h->GetXaxis()->SetTitle("True Energy (GeV)");
                h->SetLineColor(kBlack);
                TH1D *hn = (TH1D *) h->Clone(); 
                hn->Divide(h);
                hn->Draw("hist");

                leg = new TLegend(0.65,0.55,0.8,0.8);

+                 leg->AddEntry(h,dcp_angles[0],"L");
                for (int i = 1; i < 8; i++) {
                    TH1D *ho = (TH1D*) f->Get((var + "/" + s + "/Variation_" + std::to_string(i)).c_str());
                    ho->Divide(h);
                    ho->SetLineColor(kBlack + i);
                    leg->AddEntry(ho,dcp_angles[i],"L");
                //     ho->Scale(1.0, "width");
                    ho->Draw("hist same");
                }
                leg->Draw("same");
                        
                c->SaveAs((outdir + "/" + var + s + "_ratio.png").c_str());
            }
            else {
                TH1D *h = (TH1D*) f->Get((var + "/" + s + "/" + "Variation_2").c_str());
                // h->Scale(1.0, "width");
                h->GetXaxis()->SetRangeUser(0, 10);
                h->SetTitle((var + " " + s).c_str());
                h->GetXaxis()->SetTitle("True Energy (GeV)");
                h->SetLineColor(kBlack);
                h->Draw("hist");

                leg->AddEntry(h,"nominal","L");
                for (std::string variation: {"Variation_0", "Variation_4"}) {
                    TH1D *ho = (TH1D*) f->Get((var + "/" + s + "/" + variation).c_str());
                    if (variation == "Variation_4") {
                        ho->SetLineColor(kBlue);
                        leg->AddEntry(ho,"+3#sigma","L");
                    }
                    else {
                        ho->SetLineColor(kRed);
                        leg->AddEntry(ho,"-3#sigma","L");
                    }
                //     ho->Scale(1.0, "width");
                    ho->Draw("hist same");
                }
                leg->Draw("same");
                        
                c->SaveAs((outdir + "/" + var + s + ".png").c_str());


                h->GetXaxis()->SetRangeUser(0, 10);
                h->GetYaxis()->SetRangeUser(0.98, 1.02);
                h->SetTitle((var + " " + s + " ratios").c_str());
                h->GetXaxis()->SetTitle("True Energy (GeV)");
                h->SetLineColor(kBlack);
                TH1D *hn = (TH1D *) h->Clone(); 
                hn->Divide(h);
                hn->Draw("hist");

                leg = new TLegend(0.65,0.55,0.8,0.8);
                leg->AddEntry(h,"nominal","L");
                for (std::string variation: {"Variation_0", "Variation_4"}) {
                    TH1D *ho = (TH1D*) f->Get((var + "/" + s + "/" + variation).c_str());
                //     ho->Scale(1.0, "width");
                    ho->Divide(h);
                    if (variation == "Variation_4") {
                        ho->SetLineColor(kBlue);
                        leg->AddEntry(ho,"+3#sigma","L");
                    }
                    else {
                        ho->SetLineColor(kRed);
                        leg->AddEntry(ho,"-3#sigma","L");
                    }
                    ho->Draw("hist same");
                }
                leg->Draw("same");
                        
                c->SaveAs((outdir + "/" + var + s + "_ratio.png").c_str());
            }

        }
    }
}
