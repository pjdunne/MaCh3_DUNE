#pragma once
/*
HW : Simple demonstrator of a MaCh3 instance creator
*/

// C++ includes
#include <string>
#include <vector>
#include <type_traits>

// ROOT Includes
#include "TH1D.h"

// MaCh3 Includes
#include "samplePDFDUNE/samplePDFDUNEBase.h"
#include "samplePDFDUNE/samplePDFDUNEBaseND.h"

// pybind
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

// pyMaCh3
#include "pyMaCh3.h"

namespace py = pybind11;

class MaCh3Instance{
 public:
    /// @brief MaCh3Instance constructor
    /// @param yaml_config config file name
    MaCh3Instance(std::string yaml_config);

    std::vector<double> get_parameter_values();
    double propose_step(std::vector<double> new_step);

    void set_nominal_values(std::vector<double> new_vals);


    std::vector<double> get_nominal_values();
    std::vector<std::vector<std::vector<double>>>  get_event_hists();

 protected:
    template<class T>
    void setup_samples(std::vector<std::string> sample_names, double pot){
        for(std::string sample : sample_names){
            static_assert(std::is_base_of<samplePDFFDBase, T>::value, "If Error: samplePDFFDBase is not base of T");
            
            T* sample_pdf = new T(pot, sample.c_str(), xsec);

            sample_pdf->useNonDoubledAngles(true);
            sample_pdf->SetupOscCalc(osc->GetPathLength(), osc->GetDensity());

            // Add Data
            sample_pdf->reweight(osc->getPropPars());

            // Get name
            TString samp_name(sample_pdf->GetSampleName().c_str());

            switch(sample_pdf->GetBinningOpt()){
                case 1:
                    sample_pdf->addData((TH1D*)sample_pdf->get1DHist()->Clone(samp_name+"_asimov"));
                    break;
                case 2:
                    sample_pdf->addData((TH2D*)sample_pdf->get2DHist()->Clone(samp_name+"_asimov"));
                    break;
                default:
                    std::cerr<<"ERROR::Unknown sample binning opt, please pick 1 (1D) or 2 (2D)"<<std::endl;
                    std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
                    throw;
            }

            sample_vector.push_back(sample_pdf);
        }
    }

    // Important getters/setters/
    void set_parameter_values(std::vector<double> new_pars);
    double get_likelihood();
    std::vector<std::vector<double>> convert_hist_to_vec(std::shared_ptr<TH1D> input_hist);
    // Hacky hack
    std::vector<int> parameter_indices;

    std::vector<samplePDFFDBase*> sample_vector;
    covarianceXsec *xsec;
    covarianceOsc *osc;

}; // MaCh3 Instance


