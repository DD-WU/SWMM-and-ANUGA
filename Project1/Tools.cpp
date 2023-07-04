#include "Tools.h"
#include"bmi_swmm.h"
#define g 9.8
#define C 0.5
#define PI 3.1415926
double Tools::waterDepthToInflow(double waterDepth, double A, double C1) {
	double r = sqrt(A / PI);
	if (r / waterDepth > C1) {
		return C * 2 * PI * r * sqrt(2 * g * waterDepth) * waterDepth; //m³/s
	}
	else
	{
		return C * A * sqrt(2 * g * waterDepth); //m³/s
	}

}

int Tools::readSWMMPoint(char* inputFile, char*& PS_name) {
	FILE* fp_swmm;
	char line[100] = "";
	PS_name = new char[5000 * 80];
	float  xc, yc;
	int p_number = 0;
	int i = 0;
	fp_swmm = fopen(inputFile, "r");
	if (fp_swmm == NULL) { printf("no inp files for swmm5\n"); exit(0); }
	while (fgets(line, 100, fp_swmm))
	{
		line[strlen(line) - 1] = 0;

		if (strcmp("[COORDINATES]", line) == 0)                    //[COORDINATES]  nodes  outfalls 
		{
			fgets(line, 100, fp_swmm);
			fgets(line, 100, fp_swmm);
			//for (i = 0; i < 10000; i++)     //  NOT p_number
			//{
			//	fscanf(fp_swmm, " %s %f  %f ", PS_name + i * 80, &xc, &yc);
			//	if (strcmp("[VERTICES]", PS_name + i * 80) == 0)  break;          // may changed according to inp file         
			//}
			//p_number = i;
			get_pipe_name_swmm(fp_swmm, PS_name, &p_number);
		}
	}
	return p_number;
};