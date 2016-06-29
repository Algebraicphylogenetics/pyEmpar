#include <pybind11/pybind11.h>

#include <Empar/src/Empar.h>
#include <Empar/src/model.h>
#include <Empar/src/parameters.h>
#include <Empar/src/tree.h>
#include <Empar/src/alignment.h>


namespace py = pybind11;


PYBIND11_PLUGIN(empar) {
  py::module m("empar", "sample c++ bindings module");
	m.def("run", &run, "main call");

  py::class_<Model>(m, "Model");
  py::class_<Parameters>(m, "Parameters");
  py::class_<Tree>(m, "Tree");
  py::class_<Counts>(m, "Counts");

  m.def("create_model", &create_model, "Create the model");
  m.def("read_tree", &read_tree, "Read the tree");
  m.def("create_parameters", &create_parameters, "Create the parameters");
  m.def("EMalgorithm", &EMalgorithm, "Run the EM algorithm");
	return m.ptr();
}
