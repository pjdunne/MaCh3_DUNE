#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "pyMaCh3.h"

namespace py = pybind11;

void pyMaCh3Instance(py::module &);

PYBIND11_MODULE(_pyMaCh3, m) {

    py::bind_vector<std::vector<bool>>(m, "Vector_bool");
    py::bind_vector<std::vector<int>>(m, "Vector_int");
    py::bind_vector<std::vector<double>>(m, "Vector_double");
    py::bind_vector<std::vector<uint32_t>>(m, "Vector_uint32_t");

    py::bind_vector<std::vector<std::vector<double>> >(m, "Vector_Vector_double_t");

    py::implicitly_convertible<py::list, std::vector<bool>>();
    py::implicitly_convertible<py::list, std::vector<int>>();
    py::implicitly_convertible<py::list, std::vector<double>>();
    py::implicitly_convertible<py::list, std::vector<uint32_t>>();
    py::implicitly_convertible<py::list, std::vector<std::vector<double>>>();

    m.doc() = "MaCh3 implementation in python";

    pyMaCh3Instance(m);
}
