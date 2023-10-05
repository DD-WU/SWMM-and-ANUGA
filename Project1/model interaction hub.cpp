#include <stdio.h>
#include<iostream>
#include <conio.h>
#include <math.h>
#include <map>
#include "bmi_swmm.h"
#include"consts.h"
#include"space.h"
#include"swmm5.h"
#include"Tools.h"
#include <sstream>
#include<Python.h>
#include<pybind11/pybind11.h>
#include<pybind11/embed.h>
#include <omp.h>
using namespace std;
namespace py = pybind11;
using namespace py::literals;

int main() {


//[Loading models]
	py::scoped_interpreter python;
	//查看系统路径
	py::module sys = py::module::import("sys");
	py::module anuga = py::module::import("anuga");
	sys.attr("path").attr("append")("F:\\SWMM_lisflood\\xiashenshiyan");
	py::object anugaBMI = py::module::import("anugaBMI");

	//INIT_SWMM
	// set path and read file
	string swmm_inp = "F:\\SWMM_lisflood\\xiashenshiyan\\Project1\\guihua.inp";
	string swmm_out = "F:\\SWMM_lisflood\\xiashenshiyan\\Project1\\guihua.rpt";
	string swmm_rpt = "F:\\SWMM_lisflood\\xiashenshiyan\\Project1\\guihua.out";
	vector<string> argv_swmm;
	double ElapsedTime = 0;
	double swmm_t = 0.0, anuga_time = 0.0;
	char* PS_name = NULL;//bci文件里的点QVAR后的名称
	bool start_flag = true;
	argv_swmm.push_back(swmm_inp);
	argv_swmm.push_back(swmm_out);
	argv_swmm.push_back(swmm_rpt);
	initialize_swmm(argv_swmm);
	py::object jiangbei = anugaBMI.attr("BmiAnuga")();
	jiangbei.attr("initialize_ANUGA")("jiangbei.yaml");
	//[Spatial and temporal matching]
		//read anuga time
	anuga_time = jiangbei.attr("get_current_time_ANUGA")().cast<double>();
	//read swmm time
	get_current_time_swmm(&ElapsedTime);
	swmm_t = ElapsedTime * MSECperDAY / 1000.0;

	//read anuga pipe region
	py::dict pipeKeys = jiangbei.attr("regions");
	//read swmm pipe node
	int p_number = Tools::readSWMMPoint(
		"F:\\SWMM_lisflood\\xiashenshiyan\\Project1\\guihua.inp",
		PS_name);
	// inflow and overflow array
	double* p_val = new double[p_number] {0};
	double* q_in = new double[p_number] {0};

	//[LOOP]
	while (swmm_t > 0 || start_flag)
	{
		//[Data exchange]
		start_flag = false;
		if (swmm_t > anuga_time) {
			jiangbei.attr("update_until")(anuga_time +1);
			anuga_time = jiangbei.attr("get_current_time_ANUGA")().cast<double>();
		}
		else if (swmm_t < anuga_time)
		{
			update_swmm(ElapsedTime);
			get_current_time_swmm(&ElapsedTime);
			swmm_t = ElapsedTime * MSECperDAY / 1000.0;
		}
		else
		{

			for (int i = 0; i < p_number; i++) p_val[i] = 0;
			set_inflow_swmm(q_in);
			get_overflow_swmm(p_val);

			for (int i = 0; i < p_number; i++)
			{
				bool flag = pipeKeys.attr("__contains__")(PS_name + i * 80).cast<bool>();
				// overflow  into anuga
				if (flag) {
					double area = jiangbei.attr("area_map")[PS_name + i * 80].cast<double>();
					double v = p_val[i];
					jiangbei.attr("set_rate_ANUGA")("name"_a = PS_name + i * 80, "value"_a = v);
				}
			}

			double stage = 0, elevation = 0;
			for (int i = 0; i < p_number; i++) q_in[i] = 0;
			for (int i = 0; i < p_number; i++) {
				bool flag = pipeKeys.attr("__contains__")(PS_name + i * 80).cast<bool>();
				// get swmm inflow
				if (flag) {

					stage = jiangbei.attr("get_stage_ANUGA")("name"_a = PS_name + i * 80).cast<double>();
					elevation = jiangbei.attr("get_elevation_ANUGA")("name"_a = PS_name + i * 80).cast<double>();
					double q_j = 0;
					if (stage - elevation>0)
					{
						q_j = Tools::waterDepthToInflow(stage - elevation, 0.6, 0.75, p_val[i]);
					}
					q_in[i] += q_j;
					jiangbei.attr("set_rate_ANUGA")("name"_a = PS_name + i * 80, "value"_a = -q_j);
				}
			}
			//[Update the time]
			update_swmm(ElapsedTime);
			jiangbei.attr("update_until")(anuga_time + 1);
			get_current_time_swmm(&ElapsedTime);
			swmm_t = ElapsedTime * MSECperDAY / 1000.0;
			anuga_time = jiangbei.attr("get_current_time_ANUGA")().cast<double>();
		}
	}

//[Export result]
finalize_swmm();
jiangbei.attr("finalize_ANUGA")();

return 0;

}