#pragma once

// Include the input manager
#include "manager/manager.h"

// Include the samplePDFs
#include "samplePDF/samplePDFFDBase.h"

void MakeMaCh3DuneBeamInstance(manager *fitMan, std::vector<samplePDFFDBase*> &sample_vec,  covarianceXsec *&xsec, covarianceOsc *&osc);
