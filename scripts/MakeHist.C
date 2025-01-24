#include <vector>

void MakeHist() {

  std::vector<float> Binning;
  for (int i=0;i<10;i++) {
    Binning.push_back(i*0.05);
  }
  for (int i=0;i<55;i++) {
    Binning.push_back(0.5+i*0.1);
  }
  for (int i=0;i<10;i++) {
    Binning.push_back(6.0+i*0.2);
  }
  for (int i=0;i<10;i++) {
    Binning.push_back(8.0+i*0.5);
  }
  for (int i=0;i<7;i++) {
    Binning.push_back(13.0+i*1.0);
  }
  for (int i=0;i<5;i++) {
    Binning.push_back(20.0+i*5.0);
  }
  for (int i=0;i<3;i++) {
    Binning.push_back(60.0+i*20.0);
  }

  for (int i=0;i<Binning.size();i++) {
    std::cout << i << " " << Binning[i] << std::endl;
  }

  TFile* File = new TFile("OscBinning.root","RECREATE");
  TH1D* Hist = new TH1D("OscBinning","",Binning.size()-1,Binning.data());
  Hist->Write();
  File->Close();
}
