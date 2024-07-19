#include "MaCh3CovarianceHandler.h"

using py=pybind11


add_covariance(YAML::Node& systematics_manager, std::string cov_label){
    std::string file_name = GetFromManager<std::string>(systematics_manager[cov_label+"File"])
    std::string matrix_name = GetFromManager<std::string>(systematics_manager[cov_label+"Name"])


}