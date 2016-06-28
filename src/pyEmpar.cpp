#include <pybind11/pybind11.h>

#include <Empar/src/Empar.h>

namespace py = pybind11;


PYBIND11_PLUGIN(empar) {
  py::module m("empar", "sample c++ bindings module");
	m.def("run", &run, "main call");

	return m.ptr();
}
