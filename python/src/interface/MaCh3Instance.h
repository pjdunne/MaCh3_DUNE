#pragma once
/*
HW : Simple demonstrator of a MaCh3 instance creator
*/

// C++ includes
#include <string>
#include <vector>

// MaCh3 Includes
#include "samplePDFDUNE/samplePDFDUNEBase.h"
#include "samplePDFDUNE/samplePDFDUNEBaseND.h"

// pybind
#include "pybind11/pybind11.h"

namespace py = pybind11;

class MaCh3Instance{
 public:
    /// @brief MaCh3Instance constructor
    /// @param yaml_config config file name
    MaCh3Instance(const std::string yaml_config);
    double propose_step(std::vector<double> new_step);

 protected:
    template<typename t>
    void setup_samples(std::vector<std::string> sample_names, double pot);

    // Important getters/setters/
    void set_parameter_values(std::vector<double> new_pars);
    double get_likelihood();
    // Hacky hack
    std::vector<int> parameter_indices;

    std::vector<samplePDFFDBase*> sample_vector;
    covarianceXsec *xsec;
    covarianceOsc *osc;

}; // MaCh3 Instance

// PyBind wrapping
PYBIND11_MODULE(MaCh3, m){
    py::class_<MaCh3Instance>(m, "MaCh3")
        .def(py::init<const std::string yaml_config>)
        .def("propose_step", &MaCh3Instance::propose_step);
}

