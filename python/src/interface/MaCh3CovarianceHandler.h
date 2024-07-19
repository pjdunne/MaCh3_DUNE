#pragma once

// MaCh3 Includes
#include "covariance/covarianceBase.h"
#include "covariance/covarianceOsc.h"
#include "manager/manager.h"
#include "pyMaCh3.h"

// C++ Includes
#include <string>
#include <iostream>
#include <unordered_map>

// pybind
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;


/// @brief class that wraps around usual set of covariance matrices
class MaCh3CovarianceHandler{
 public:
    /// @brief Constructor
    MaCh3CovarianceHandler(){}
    virtual ~MaCh3CovarianceHander(){};

    /// @brief Add new covariance to handler object
    /// @param systematics_manager Node containing ["General"]["Systematics"]
    /// @param cov_label Assumes covariance is stored as <cov_label>FILE and <cov_label>Name
    add_covariance(YAML::Node& systematics_manager, std::string cov_label);

    /// @brief Lets covariance know which object is the oscillator
    /// @param cov_label name used in covariance_map
    set_oscillator(std::string cov_label);

 protected:
    /// Data structure for storing covariance objects
    std::unordered_map<std::string, covarianceBase*> covariance_map;

};