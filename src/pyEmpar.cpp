#include <pybind11/pybind11.h>

#include <Empar/src/Empar.h>
#include <Empar/src/model.h>
#include <Empar/src/parameters.h>
#include <Empar/src/tree.h>
#include <Empar/src/alignment.h>
#include <Empar/src/em.h>
#include <Empar/src/state_list.h>
#include <Empar/src/matrix.h>


namespace py = pybind11;


PYBIND11_PLUGIN(empar) {
  py::module m("empar", "sample c++ bindings module");
	m.def("run", &run, "main call");

  py::class_<Model>(m, "Model");
  py::class_<Parameters>(m, "Parameters");
  py::class_<Tree>(m, "Tree");
  py::class_<Counts>(m, "Counts");
  py::class_<Counts>(m, "Counts");
  py::class_<StateList>(m, "StateList");
  py::class_<Matrix>(m, "Matrix");

  m.def("create_model", &create_model, "Create the model");
  m.def("read_tree", &read_tree, "Read the tree");
  m.def("create_parameters", &create_parameters, "Create the parameters");
  m.def("EMalgorithm", &EMalgorithm, "Run the EM algorithm");
  m.def("create_state_list", &create_state_list, "Create a list of states");

  m.def("create_matrix", &create_matrix, "Create matrix");
  m.def("log_likelihood_fast", &log_likelihood_fast, "Calculate loglik- fast version");
  m.def("joint_prob_fast", &joint_prob_fast, "Calculate joint probability");
	return m.ptr();
}
