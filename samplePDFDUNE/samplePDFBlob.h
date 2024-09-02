#pragma once

#include "covariance/covarianceOsc.h"
#include "covariance/covarianceXsec.h"
#include "samplePDF/samplePDFFDBase.h"
#include "splines/splinesDUNE.h"

#include "StructsDUNE.h"

#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "Eigen/Core"

#include "yaml-cpp/yaml.h"

#include <sstream>
#include <stdexcept>
#include <string>

struct ColumnarBlob {

  std::vector<std::string> double_column_names;
  Eigen::ArrayXXd double_blob;

  std::vector<std::string> int_column_names;
  Eigen::ArrayXXi int_blob;

  int num_rows;

  ColumnarBlob(std::vector<std::string> dcolnames,
               std::vector<std::string> icolnames, int nrows)
      : double_column_names(dcolnames), int_column_names(icolnames),
        num_rows(nrows) {
    double_blob = Eigen::ArrayXXd::Zero(nrows, dcolnames.size());
    int_blob = Eigen::ArrayXXi::Zero(nrows, icolnames.size());
  }

  template <typename T> auto col(std::string const &cname) {
    if constexpr (std::is_same_v<T, double>) {
      for (size_t i = 0; i < double_column_names.size(); ++i) {
        if (cname == double_column_names[i]) {
          return double_blob.col(i);
        }
      }
      std::stringstream ss;
      ss << "Failed to find double column named: " << cname;
      throw std::runtime_error(ss.str());
    } else if constexpr (std::is_same_v<T, int>) {
      for (size_t i = 0; i < int_column_names.size(); ++i) {
        if (cname == int_column_names[i]) {
          return int_blob.col(i);
        }
      }
      std::stringstream ss;
      ss << "Failed to find int column named: " << cname;
      throw std::runtime_error(ss.str());
    }
    throw std::logic_error("invalid type");
  }

  template <typename T> auto col(std::string const &cname) const {
    if constexpr (std::is_same_v<T, double>) {
      for (size_t i = 0; i < double_column_names.size(); ++i) {
        if (cname == double_column_names[i]) {
          return double_blob.col(i);
        }
      }
      std::stringstream ss;
      ss << "Failed to find double column named: " << cname;
      throw std::runtime_error(ss.str());
    } else if constexpr (std::is_same_v<T, int>) {
      for (size_t i = 0; i < int_column_names.size(); ++i) {
        if (cname == int_column_names[i]) {
          return int_blob.col(i);
        }
      }
      std::stringstream ss;
      ss << "Failed to find int column named: " << cname;
      throw std::runtime_error(ss.str());
    }
    throw std::logic_error("invalid type");
  }

  template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> new_column() const {
    if constexpr (std::is_same_v<T, double>) {
      return Eigen::VectorXd::Zero(num_rows);
    } else if constexpr (std::is_same_v<T, int>) {
      return Eigen::VectorXi::Zero(num_rows);
    }
  }

  std::string to_string() {
    std::stringstream ss;

    ss << "double columns: ";
    for (auto const &dcn : double_column_names) {
      ss << dcn << ", ";
    }
    ss << "\nint columns: ";
    for (auto const &icn : int_column_names) {
      ss << icn << ", ";
    }
    ss << "\n";

    for (int i = 0; i < std::min(20, num_rows); ++i) {
      ss << "\ndouble.col[" << i << "]: ";
      for (int j = 0; j < double_blob.cols(); ++j) {
        ss << double_blob(i, j) << ", ";
      }

      ss << "\nints.col[" << i << "]: ";
      for (int j = 0; j < int_blob.cols(); ++j) {
        ss << int_blob(i, j) << ", ";
      }
    }
    ss << "\n";
    return ss.str();
  }
};

inline ColumnarBlob BuildExperimentBlob(TTree *ttree, YAML::Node const &yhf) {
  TTreeReader ttreader(ttree);

  std::vector<TTreeReaderValue<double>> doublevars;
  std::vector<TTreeReaderValue<int>> intvars;

  std::vector<std::string> doublecnames;
  std::vector<std::string> intcnames;

  for (auto branch : yhf["branches"]) {

    auto bname = branch[0].as<std::string>();
    auto btype = branch[1].as<char>();

    if (btype == 'd') {
      doublevars.emplace_back(ttreader, bname.c_str());
      doublecnames.push_back(bname);
      std::cout << "added double branch: " << bname << std::endl;
    } else if (btype == 'i') {
      intvars.emplace_back(ttreader, bname.c_str());
      intcnames.push_back(bname);
      std::cout << "added int branch: " << bname << std::endl;
    } else {
      throw "invalid branch type.";
    }
  }

  Long64_t numrows = ttree->GetEntries();
  std::cout << "Processing " << numrows << " entries." << std::endl;

  ColumnarBlob blob(doublecnames, intcnames, numrows);

  Long64_t i = 0;
  while (ttreader.Next()) {

    for (size_t j = 0; j < doublevars.size(); ++j) {
      // if (!i) {
      //   if (!doublevars[j].IsValid()) {
      //     std::stringstream ss;
      //     ss << "invalid branch name: " << doublecnames[j]
      //        << ", or type: double";
      //     throw std::runtime_error(ss.str());
      //   }
      // }
      blob.double_blob(i, j) = *doublevars[j];
    }
    for (size_t j = 0; j < intvars.size(); ++j) {
      // if (!i) {

      //   if (!intvars[j].IsValid()) {
      //     std::stringstream ss;
      //     ss << "invalid branch name: " << intcnames[j] << ", or type: int";
      //     throw std::runtime_error(ss.str());
      //   }
      // }
      blob.int_blob(i, j) = *intvars[j];
    }

    i++;
  }

  return blob;
}

struct ColumnOp {

  std::string column_name_double;
  std::function<Eigen::VectorXd(ColumnarBlob const &)> Op_double;

  std::string column_name_int;
  std::function<Eigen::VectorXi(ColumnarBlob const &)> Op_int;
};

template <typename T>
inline ColumnOp new_column_op(
    std::string const &ncolname,
    std::function<Eigen::Vector<T, Eigen::Dynamic>(ColumnarBlob const &)>
        operation) {
  if constexpr (std::is_same_v<T, double>) {
    return ColumnOp{ncolname, operation, "", nullptr};
  } else if constexpr (std::is_same_v<T, int>) {
    return ColumnOp{"", nullptr, ncolname, operation};
  } else {
    throw std::runtime_error("invalid new column op type");
  }
}

template <typename T>
inline auto CopyColumn(std::string const &column_name,
                       std::string const &new_column_name = "") {
  return new_column_op<T>(new_column_name.size() ? new_column_name
                                                 : column_name,
                          [=](ColumnarBlob const &exptblob) {
                            return exptblob.col<T>(column_name);
                          });
}

inline ColumnarBlob Blob2Blob(ColumnarBlob const &exptblob,
                              std::vector<ColumnOp> const &ops) {

  std::vector<std::string> double_column_names;
  std::vector<std::string> int_column_names;

  for (auto const &op : ops) {
    if (op.column_name_double.size()) {
      double_column_names.push_back(op.column_name_double);
    } else if (op.column_name_int.size()) {
      int_column_names.push_back(op.column_name_int);
    } else {
      throw std::runtime_error("operation doesn't specify output column name");
    }
  }

  ColumnarBlob evblob(double_column_names, int_column_names, exptblob.num_rows);

  for (auto const &op : ops) {
    if (op.column_name_double.size()) {
      evblob.col<double>(op.column_name_double) = op.Op_double(exptblob);
    } else {
      evblob.col<int>(op.column_name_int) = op.Op_int(exptblob);
    }
  }

  return evblob;
}

class samplePDFDUNEBlob : public samplePDFFDBase {
public:
  samplePDFDUNEBlob();
  ~samplePDFDUNEBlob();

  void setupFDMC(ColumnarBlob const &evblob, int nutype, int oscnutype,
                 bool signal);

  std::vector<double> rw_etru;
  std::vector<int> mode;
  std::vector<int> Target;
  std::vector<double> x_var;
  std::vector<double> y_var;

  void printPosteriors() {}
  double getCovLikelihood() { return 1; }
  double getDiscVar(int, int, int) { return 1; }
  int getNMCSamples() { return 1; }
  int getNEventsInSample(int) { return rw_etru.size(); }
  void setupWeightPointers() {}
  double CalcXsecWeightFunc(int, int) {}
  double ReturnKinematicParameter(std::string, int, int) {}
  std::vector<double> ReturnKinematicParameterBinning(std::string) {}
};
