#include "samplePDFDUNE/samplePDFBlob.h"

#include "yaml-cpp/yaml.h"

#include "TFile.h"

#include <iostream>

int main(int argc, char *argv[]) {

  auto tf = std::unique_ptr<TFile>(TFile::Open(argv[1]));
  auto ttree = tf->Get<TTree>("caf");

  auto doc = YAML::LoadFile(argv[2]);

  auto exptblob = BuildExperimentBlob(ttree, doc);

  std::cout << exptblob.to_string() << std::endl;

  auto myerec = new_column_op<double>("xvar", [](ColumnarBlob const &exptblob) {
    auto eRecoP = exptblob.col<double>("eRecoP");
    auto eRecoN = exptblob.col<double>("eRecoN");
    auto eRecoPip = exptblob.col<double>("eRecoPip");
    auto eRecoPim = exptblob.col<double>("eRecoPim");
    auto eRecoPi0 = exptblob.col<double>("eRecoPi0");
    auto eRecoOther = exptblob.col<double>("eRecoOther");

    return eRecoP + eRecoN + eRecoPip + eRecoPim + eRecoPi0 + eRecoOther;
  });

  auto hasrecopi =
      new_column_op<int>("hasrecopi", [](ColumnarBlob const &exptblob) {
        auto eRecoPip = exptblob.col<double>("eRecoPip");
        auto eRecoPim = exptblob.col<double>("eRecoPim");
        auto eRecoPi0 = exptblob.col<double>("eRecoPi0");

        auto haspi = exptblob.new_column<int>();

        for (size_t i = 0; i < exptblob.num_rows; ++i) {
          haspi(i) =
              ((eRecoPip(i) > 0) || (eRecoPim(i) > 0) || (eRecoPi0(i) > 0));
        }

        return haspi;
      });

  auto isNC = new_column_op<int>("isNC", [](ColumnarBlob const &exptblob) {
    auto eRecoPip = exptblob.col<int>("isCC");

    auto nc = exptblob.new_column<int>();

    for (size_t i = 0; i < exptblob.num_rows; ++i) {
      nc(i) = !eRecoPip(i);
    }

    return nc;
  });

  auto evblob = Blob2Blob(
      exptblob, std::vector<ColumnOp>{myerec, hasrecopi,
                                      CopyColumn<double>("Ev", "etrue"),
                                      CopyColumn<int>("mode"), isNC});

  std::cout << evblob.to_string() << std::endl;

  std::vector<double> xbins(50), ybins(50);

  Eigen::VectorXd::Map(xbins.data(), 50) =
      Eigen::VectorXd::LinSpaced(50, 0, 10);
  Eigen::VectorXd::Map(ybins.data(), 50) =
      Eigen::VectorXd::LinSpaced(50, 0, 10);

  samplePDFDUNEBlob pdfb;
  pdfb.setupFDMC(evblob, 1, 1, true);
  pdfb.set1DBinning(xbins.size() - 1, xbins.data());
  pdfb.set2DBinning(xbins.size() - 1, xbins.data(), ybins.size() - 1,
                    ybins.data());

  covarianceOsc *osc = new covarianceOsc(
      "osc_cov", "../inputs/osc_covariance_DUNE_PDG2021_v1.root");

  osc->setParameters(std::vector<double>(6, 0));
  osc->acceptStep();

  pdfb.SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
  pdfb.reweight(osc->getPropPars());

  auto h = pdfb.get1DHist();
  TFile out("out.root", "RECREATE");
  out.WriteTObject(h, "myh");
  out.Write();
  out.Close();
}
