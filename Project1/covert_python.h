#pragma once

void init_python(wchar_t* python_url,char* env_path= NULL) {
	Py_SetPythonHome(python_url);
	Py_Initialize();
	PyRun_SimpleString("from __future__ import print_function");
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.path.append('./')");
	PyRun_SimpleString(env_path);//"sys.path.append('F:/SWMM_lisflood/xiashenshiyan/anuga_bmi')"
	PyRun_SimpleString("import numpy as np");
	PyRun_SimpleString("from anuga_bmi import BmiAnuga");
}
int init_python_var(PyObject*& pModule, PyObject*& pFunc, PyObject*& pArgs, PyObject*& sww,
	PyObject*& pDict, PyObject*& pInstancesww, PyObject*& sw) {
	pModule = PyImport_ImportModule("anuga_bmi.anugaBMI");
	if (!pModule) {
		printf("Cant open python file!/n");
		return -1;
	}
	pDict = PyModule_GetDict(pModule);
	if (!pDict) {
		printf("Cant find dictionary./n");
		return -1;
	}
	sww = PyDict_GetItemString(pDict, "BmiAnuga");
	if (!sww) {
		printf("Cant find second class./n");
		return -1;
	}
	pInstancesww = PyInstanceMethod_New(sww);
	sw = PyObject_CallObject(pInstancesww, sw);
	if (!sw) {
		return -1;
	}
	pArgs = PyTuple_New(1);
	if (!pArgs) {
		return -1;
	}
	//PyObject_GetAttrString(sw, "_anuga");
	//if (anuga == NULL) {
	//	printf("can't get bmi_anuga().anuga");
	//	return -1;
	//}
	return 0;
}