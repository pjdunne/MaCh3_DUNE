#include "samplePDFBlob.h"

#include <cassert>
#include <stdexcept>

samplePDFDUNEBlob::samplePDFDUNEBlob()
    : samplePDFBase(1) {

  _hPDF1D = new TH1D("hErec_nue", "Reconstructed Energy", 200, 0, 50.0);
  _hPDF2D = new TH2D("a", "b", 15, 0, 50.0 * 1000, 15, 0, 150);
  dathist = new TH1D("dat_nue", "", 200, 0, 50.0);
  dathist2d = new TH2D("dat2d_nue", "", 15, 0, 1500, 15, 0, 150);

  setXsecCov(nullptr);
}

samplePDFDUNEBlob::~samplePDFDUNEBlob() {}

void samplePDFDUNEBlob::setupFDMC(ColumnarBlob const &evblob, int nutype,
                                  int oscnutype, bool signal) {

  MCSamples.emplace_back();
  auto &core_sample = MCSamples.back();

  core_sample.nEvents = evblob.num_rows;
  core_sample.nutype = nutype;
  core_sample.oscnutype = oscnutype;
  core_sample.signal = signal;
  core_sample.SampleDetID = 0;

  core_sample.x_var = new double *[core_sample.nEvents];
  core_sample.y_var = new double *[core_sample.nEvents];
  core_sample.enu_s_bin = new unsigned int[core_sample.nEvents];
  core_sample.xvar_s_bin = new unsigned int[core_sample.nEvents];
  core_sample.yvar_s_bin = new unsigned int[core_sample.nEvents];
  core_sample.rw_etru = new double *[core_sample.nEvents];
  core_sample.XBin = new int[core_sample.nEvents];
  core_sample.YBin = new int[core_sample.nEvents];
  core_sample.NomXBin = new int[core_sample.nEvents];
  core_sample.NomYBin = new int[core_sample.nEvents];
  core_sample.XBin = new int[core_sample.nEvents];
  core_sample.YBin = new int[core_sample.nEvents];

  core_sample.rw_lower_xbinedge = new double[core_sample.nEvents];
  core_sample.rw_lower_lower_xbinedge = new double[core_sample.nEvents];
  core_sample.rw_upper_xbinedge = new double[core_sample.nEvents];
  core_sample.rw_upper_upper_xbinedge = new double[core_sample.nEvents];
  core_sample.mode = new int *[core_sample.nEvents];
  core_sample.xsec_w = new double[core_sample.nEvents];
  core_sample.flux_w = new double[core_sample.nEvents];
  core_sample.osc_w = new double[core_sample.nEvents];
  core_sample.isNC = new bool[core_sample.nEvents];
  core_sample.nxsec_norm_pointers = new int[core_sample.nEvents];
  core_sample.xsec_norm_pointers = new const double **[core_sample.nEvents];
  core_sample.xsec_norms_bins = new std::list<int>[core_sample.nEvents];
  core_sample.nxsec_spline_pointers = new int[core_sample.nEvents];
  core_sample.xsec_spline_pointers = new const double **[core_sample.nEvents];
  core_sample.ntotal_weight_pointers = new int[core_sample.nEvents];
  core_sample.total_weight_pointers = new double **[core_sample.nEvents];
  core_sample.Target = new int *[core_sample.nEvents];

  rw_etru = std::vector<double>(core_sample.nEvents);
  mode = std::vector<int>(core_sample.nEvents);
  Target = std::vector<int>(core_sample.nEvents, 40);
  x_var = std::vector<double>(core_sample.nEvents);

  Eigen::VectorXd::Map(rw_etru.data(), rw_etru.size()) =
      evblob.col<double>("etrue");
  Eigen::VectorXi::Map(mode.data(), mode.size()) = evblob.col<int>("mode");
  Eigen::VectorXd::Map(x_var.data(), x_var.size()) = evblob.col<double>("xvar");

  Eigen::VectorXi::Map(core_sample.nxsec_norm_pointers, core_sample.nEvents) =
      Eigen::VectorXi::Constant(core_sample.nEvents, 0);
  Eigen::VectorXi::Map(core_sample.nxsec_spline_pointers, core_sample.nEvents) =
      Eigen::VectorXi::Constant(core_sample.nEvents, 0);
  Eigen::VectorXi::Map(core_sample.ntotal_weight_pointers,
                       core_sample.nEvents) =
      Eigen::VectorXi::Constant(core_sample.nEvents, 0);

  Eigen::VectorXi::Map(core_sample.NomXBin, core_sample.nEvents) =
      Eigen::VectorXi::Constant(core_sample.nEvents, -1);
  Eigen::VectorXi::Map(core_sample.NomYBin, core_sample.nEvents) =
      Eigen::VectorXi::Constant(core_sample.nEvents, -1);

  Eigen::VectorXi::Map(core_sample.XBin, core_sample.nEvents) =
      Eigen::VectorXi::Constant(core_sample.nEvents, -1);
  Eigen::VectorXi::Map(core_sample.YBin, core_sample.nEvents) =
      Eigen::VectorXi::Constant(core_sample.nEvents, -1);

  Eigen::VectorXd::Map(core_sample.rw_lower_xbinedge, core_sample.nEvents) =
      Eigen::VectorXd::Constant(core_sample.nEvents, -1);

  Eigen::VectorXd::Map(core_sample.rw_lower_lower_xbinedge,
                       core_sample.nEvents) =
      Eigen::VectorXd::Constant(core_sample.nEvents, -1);

  Eigen::VectorXd::Map(core_sample.rw_upper_xbinedge, core_sample.nEvents) =
      Eigen::VectorXd::Constant(core_sample.nEvents, -1);

  Eigen::VectorXd::Map(core_sample.rw_upper_upper_xbinedge,
                       core_sample.nEvents) =
      Eigen::VectorXd::Constant(core_sample.nEvents, -1);

  Eigen::VectorXd::Map(core_sample.xsec_w, core_sample.nEvents) =
      Eigen::VectorXd::Constant(core_sample.nEvents, 1);
  Eigen::VectorXd::Map(core_sample.osc_w, core_sample.nEvents) =
      Eigen::VectorXd::Constant(core_sample.nEvents, 1);

  Eigen::VectorX<bool>::Map(core_sample.isNC, core_sample.nEvents) =
      evblob.col<int>("isNC").cast<bool>();
  Eigen::VectorXd::Map(core_sample.flux_w, core_sample.nEvents) =
      Eigen::VectorXd::Constant(core_sample.nEvents, 1);

  for (int iEvent = 0; iEvent < core_sample.nEvents; ++iEvent) {
    core_sample.rw_etru[iEvent] = &(rw_etru[iEvent]);
    core_sample.mode[iEvent] = &(mode[iEvent]);
    core_sample.Target[iEvent] = &(Target[iEvent]);

    core_sample.x_var[iEvent] = &(x_var[iEvent]);

    // only ever 1D in my samplePDF
    core_sample.y_var[iEvent] = &(x_var[iEvent]);
  }

  std::cout << "MCSamples.size(): " << MCSamples.size() << std::endl;
  std::cout << "MCSamples[0].nEvents: " << MCSamples[0].nEvents << std::endl;
}
