#pragma once
#include <string>
class Tools
{
public:
	static double waterDepthToInflow(double waterDepth, double A, double C1, double flag);//���ܼ���������С��λ��m^3/s��,����ֵΪ����ֵ
	static int readSWMMPoint(char* inputFile, char*& PS_name);//����ֵΪ��ϵ����
	static int* calcCommonPoint(int* c_count, char* PS_name, int p_number, int num);//ȷ����Ϲ�ϵ��һ��һ����һ�Զ�
};

