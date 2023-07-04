#include <stdio.h>
#include<iostream>
#include <conio.h>
#include <math.h>
#include <map>
#include "bmi_swmm.h"
#include "bmi_lisflood.h"
#include"consts.h"
#include"space.h"
#include"swmm5.h"
#include"Tools.h"
#include <sstream>
//int main() {
//	// set path and read file
//	string swmm_inp = "F:\\SWMM_lisflood\\xiashenshiyan\\Project1\\guihua.inp";
//	string swmm_out = "F:\\SWMM_lisflood\\xiashenshiyan\\Project1\\guihua.rpt";
//	string swmm_rpt = "F:\\SWMM_lisflood\\xiashenshiyan\\Project1\\guihua.out";
//	string lisflood_conf = "F:\\SWMM_lisflood\\xiashenshiyan\\Project1\\jiangbei.par";
//	vector<string> argv_lisflood;
//	argv_lisflood.push_back(lisflood_conf);
//	vector<string> argv_swmm;
//	argv_swmm.push_back(swmm_inp);
//	argv_swmm.push_back(swmm_out);
//	argv_swmm.push_back(swmm_rpt);
//
//	//INIT_MODEL
//	initialize_swmm(argv_swmm);//这个是zaijian全曼宁系数相等
//	initialize_lisflood(argv_lisflood);
//
//	//read lisflood time
//	double endtime = 0, starttime = 0, initstep = 0, t = 0;
//	get_end_time_lisflood(&endtime);
//	get_start_time_lisflood(&starttime);
//	get_time_step_lisflood(&initstep);
//	get_current_time_lisflood(&t);
//
//	//read swmm time
//	double ElapsedTime = 0, time = 0;
//	get_time_step_swmm(&ElapsedTime);
//
//	//read parameter
//	int i, j, k, ix, iy;//i,j,k 循环变量，ix,iy（bci文件里的点坐标）
//	char* PS_name = NULL;//bci文件里的点QVAR后的名称
//	int p_number = Tools::readCommonPoint("F:\\SWMM_lisflood\\xiashenshiyan\\Project1\\guihua.inp", PS_name); //读inp里的点数，并给ps_name赋值
//	int num = get_PSNum();//获取bci文件里的点数
//	double* dx_ = get_dx();//获取分辨率
//	double* p_val = new double[p_number] {0};//溢流数组
//	double* q_in = new double[p_number] {0};//入流数组	
//	int* c_count = new int[p_number] {0};//inp文件中一个点对应的bci文件点数	
//	c_count=Tools::calcCommonPoint(c_count,PS_name,p_number,num);
//	// 各走一步
//	update_swmm(ElapsedTime);
//	update_lisflood(initstep);
//	//更新时间
//	get_time_step_lisflood(&initstep);
//	get_current_time_lisflood(&t);
//	get_time_step_swmm(&ElapsedTime);
//	time = ElapsedTime * MSECperDAY / 1000;
//
//	//loop
//	while (ElapsedTime > 0 && t < endtime)
//	{
//		if (abs(time - t)<0.1) {//时间相等时
//			std::cout << "time:"<<t<<"\n";
//			for (i = 0; i < p_number; i++) p_val[i] = 0;
//			swmm_inflow(q_in, p_val);//给swmm赋入流q_in，并获取溢流p_val
//			
//			for (j = 0; j < num; j++)
//			{
//				for (i = 0; i < p_number; i++)
//				{
//					if (strcmp((get_PSName() + j * 80), PS_name + i * 80) == 0)//找到同名点
//					{
//						double* data = get_BCVar(j);   // 获取bdy文件所对应的时间序列的引用
//						k = 0;
//						while (data[k] != -1)//改序列里的流速值
//						{
//							if (k % 2 == 0) {
//								data[k] = p_val[i] / *dx_/c_count[i];    // overflow  into LISFLOOD 单位为m^2/s
//								if (c_count[i] == 0) std::cout << "报错" << "\n";
//							}
//							k++;
//						}
//					}
//				}
//			}
//			double q_j = 0;//溢流量
//			for (i = 0; i < p_number; i++) q_in[i] = 0;
//			for (j = 0; j < num; j++)
//				for (i = 0; i < p_number; i++)
//					if (strcmp((get_PSName() + j * 80), PS_name + i * 80) == 0)
//					{
//						//获取点位坐标
//						ix = get_BCxpi(j);
//						iy = get_BCypi(j);
//						if (c_count[i] > 1) {
//							if (0.04 < *dx_ * *dx_) {
//								q_j = Tools::waterDepthToInflow(throwWaterDepth(ix, iy), 0.04) / c_count[i]; //  leaving LISFLOOD  给溢流量赋值
//							}
//							else
//								q_j = Tools::waterDepthToInflow(throwWaterDepth(ix, iy), (*dx_) * (*dx_)) / c_count[i]; // leaving LISFLOOD
//						}
//						else
//						{
//							if (0.04 < *dx_ * *dx_) {
//								q_j = Tools::waterDepthToInflow(throwWaterDepth(ix, iy), 0.04) / c_count[i]; //  leaving LISFLOOD  给溢流量赋值
//							}
//							else
//								q_j = Tools::waterDepthToInflow(throwWaterDepth(ix, iy), (*dx_) * (*dx_)) / c_count[i]; // leaving LISFLOOD
//						}
//					
//						  //   dx*dx  A   smaller?    lisflood many to  swmm one source
//						q_in[i] = q_in[i] + q_j; // into SWMM 计算要给swmm多少水
//						if (p_val[i] > 0)//如果swmm有溢流就什么都不做
//						{
//						}
//						else//否则，就把刚刚循环改过的bdy里的流量改为负的刚刚算出来的值
//						{
//							double* data = get_BCVar(j);    //获取bdy里的第j个时间序列
//							k = 0;
//							while (data[k] != -1)
//							{
//								if (k % 2 == 0) {
//									data[k] = -q_j / *dx_;  //  leaving LISFLOOD 
//								}
//								k++;
//							}
//						}
//					}
//		}
//		else if (time > t)//swmm比lisflood快 
//		{
//			std::cout << "不统一了\n";
//
//			update_lisflood(initstep);
//			get_current_time_lisflood(&t);
//			continue;
//		}
//		else//lisflood比swmm快 
//		{
//			std::cout << "不统一了\n";
//			update_swmm(ElapsedTime);
//			get_time_step_swmm(&ElapsedTime);
//			time = ElapsedTime * MSECperDAY / 1000;
//			continue;
//		}
//		// 各走一步
//		update_swmm(ElapsedTime);
//		update_lisflood(initstep);
//		//更新时间
//		get_time_step_lisflood(&initstep);
//		get_current_time_lisflood(&t);
//		get_time_step_swmm(&ElapsedTime);
//		time = ElapsedTime * MSECperDAY / 1000;
//
//	};
//	while (ElapsedTime > 0)
//	{
//		update_swmm(ElapsedTime);
//		get_time_step_swmm(&ElapsedTime);
//		time = ElapsedTime * MSECperDAY / 1000;
//	}
//	while (t < endtime)
//	{
//		update_lisflood(initstep);
//		get_current_time_lisflood(&t);
//	}
//	//结束
//	finalize_swmm();
//	finalize_lisflood();
//	printf("");
//}