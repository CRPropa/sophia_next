#include <pybind11/pybind11.h>

#include "sophia_interface.h"

namespace py = pybind11;

PYBIND11_MODULE(pysophia, m) {
	m.doc() = "SophiaNext python binding";
	py::class_<sophiaevent_output, std::shared_ptr<sophiaevent_output>>(m, "SophiaEventOutput")
		.def(py::init<>())
		.def_readonly("Nout", &sophiaevent_output::Nout)
		.def("getPartP", &sophiaevent_output::getPartP)
		.def("getPartID", &sophiaevent_output::getPartID);
	py::class_<sophia_interface, std::shared_ptr<sophia_interface>>(m, "SophiaInterface")
		.def(py::init<>())
		.def("sophiaevent", &sophia_interface::sophiaevent,
			py::arg("onProton"),
			py::arg("Ein"),
			py::arg("eps"),
			py::arg("declareChargedPionsStable") = false);
}

