#pragma once

/*
Base class instance of a MaCh3 Systematic
*/

// MaCh3 includes
#include "covariance/covarianceBase.h"
#include "MaCh3Object.h"

// C++ Includes
#include <iostream>
#include <vector>
#include <string>
#include <memory>


class MaCh3Systematic : public MaCh3Object<covarianceBase>{
   public:
    /// @brief MaCh3Systematic constructor
    MaCh3Systematic(std::string config_opt);

    /// @brief destructor
    virtual ~MaCh3Systematic(){}

    /// Likelihood getter
    double get_likelihood(){
        return _stored_object->GetLikelihood();
    }

    /// @brief Set the model to a new set of values
    void set_parameter_values(std::vector<double> parameter_values){
        /// Mute the couts
        std::cout.setstate(std::ios_base::failbit);

        if((int)parameter_values.size()!=_stored_object->getNpars()){
            std::cerr<<"ERROR::Trying to set wrong number of values"<<std::endl;
            std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
            throw;
        }

        for(int i=0; i<_stored_object->getNpars()){
            _stored_object->setParProp(i, parameter_values[i]);
        }

        /// Unmute them
        std::cout.clear();
    }

   protected:

    /// @brief Constructs dummy covarianceBase, since this is never actually likely to be used I'm giving it a dummy arg!
    std::unique_ptr<covariancebase> construct_default() override{
        return std::unique_ptr<samplePDFBase>(covarianceBase());
    }

};