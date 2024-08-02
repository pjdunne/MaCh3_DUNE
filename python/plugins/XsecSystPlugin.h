#pragma once

// C++ Includes
#include <memory>
#include <vector>
#include <string>

// MaCh3 includes
#include "../src/interface/MaCh3Plugin.h"
#include "covariance/covarianceXsec.h"

class XsecSystPlugin : public MaCh3Plugin<covarianceXsec>{
    XsecSystPlugin() : MaCh3Plugin() {
        override_default_constructor = true;
    }

    std::unique_ptr<covarianceXsec> construct_obj(std::string matrix_name, std::string file_name){
        auto xsec = std::make_unique<covarianceXsec>(matrix_name, file_name);
        return xsec;
    }


    void process(MaCh3Systematic * mach3_object){
        mach3_object->get_stored_object()->setParameters();
    }

};