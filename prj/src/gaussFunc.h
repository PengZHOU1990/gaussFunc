#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <Eigen/Eigen>
#include <iomanip>

using namespace std;
using namespace Eigen;

#define SQR2(x) ((x)*(x))
#define SQR3(x) (SQR2(x)*(x))
#define SQR4(x) (SQR3(x)*(x))
#define SQR6(x) (SQR4(x)*SQR2(x))
#define SQR8(x) (SQR4(x)*SQR4(x))
#define SQR10(x) (SQR8(x)*SQR2(x))

#define PI (3.1415926535898)
#define D2R (PI/180.)
#define R2D (180./PI)

// 可选的椭球类型
enum  enEarType
{
	SYSBJ54 = 0,
	SYSGDZ80,
	SYSWGS84,
	SYSCGCS2000,
};

// 椭球参数，建议通过两个参数（长半轴、扁率）构造椭球，即自定义的带2个参数的构造函数
struct stEarPara
{
	double a;
	double f;
	double e;
	double b;
	double e2;
	stEarPara(){}
	stEarPara(double da, double df)
	{
		a = da, f = df;
		e = sqrt(2 * f - f * f);
		b = a - a * f;
		e2 = sqrt(a*a - b* b) / b;
	}
};

// 坐标转换参数，可用于平面四参数&空间直角坐标7参数
struct stCoorTransPara
{
	double dX;
	double dY;
	double dZ;
	double rX;
	double rY;
	double rZ;
	double m;
};

// 几种常见椭球的集合参数，通过长半轴、扁率即可获取其他几何参数
// ref--孔祥元，P98； 潘正风，P8； ICD BDS，2013
static const stEarPara gBJ54Ear(6378245., 1. / 298.3);

static const stEarPara gGDZ80Ear(6378140., 1. / 298.257);

static const stEarPara gWGS84Ear(6378137., 1. / 298.257223563);

static const stEarPara gCGCS2000Ear(6378137., 1. / 298.257222101);


typedef struct stRawData
{
	string m_sN;
	double m_dX1;
	double m_dY1;
	double m_dZ1;
	double m_dX2;
	double m_dY2;
	double m_dZ2;
};

void printMat(const string &str, const int &iWid, const int &iFixed, ofstream &os, const Eigen::MatrixXd &mat);

void parseStr(const string &sSrc, const string &sParser, vector<string> &vsParsed);

string fTrimAll(string &sSrc);



// 高斯投影坐标正算不分带: 指定中央子午线,计算的dX dY不加500km，不加带号
// dB -- [-90, 90], deg; dL -- [-180, 180], deg; dL0 -- [0., 360.], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int _declspec(dllexport) gaussFw(enEarType earType, double dB, double dL, double dL0, double *dX, double *dY);

// 高斯投影坐标正算不分带: 指定中央子午线,计算的dX dY加500km，不加带号
// dB -- [-90, 90], deg; dL -- [-180, 180], deg; dL0 -- [0., 360.], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int _declspec(dllexport) gaussFwWithAdd(enEarType earType, double dB, double dL, double dL0, double *dX, double *dY);

// 高斯投影坐标正算6度带: 指定中央子午线,计算的dX dY加500km，加带号
// dB -- [-90, 90], deg; dL -- [-180, 180], deg; dL0 -- [0., 360.], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int _declspec(dllexport) gaussFwWithBelt(enEarType earType, double dB, double dL, double dL0, double *dX, double *dY);

// 高斯投影坐标反算: 指定中央子午线,输入的dx dy不加500km，不带带号
// dL0 -- [0., 360.], deg; dB -- [-90, 90], deg; dL -- [0, 360], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int _declspec(dllexport) gaussBw(enEarType earType, double dx, double dy, double dL0, double *dB, double *dL);

// 高斯投影坐标反算: 指定中央子午线,输入的dx dy加500km，不带带号
// dL0 -- [0., 360.], deg; dB -- [-90, 90], deg; dL -- [0, 360], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int _declspec(dllexport) gaussBwWithAdd(enEarType earType, double dx, double dy, double dL0, double *dB, double *dL);

// 求底点纬度
double calcuBf(enEarType earType, double dx);

// 高斯投影坐标反算: *另一种算法*指定中央子午线,输入的dx dy加500km，不带带号，与gaussBw1精度相当
// dL0 -- [0., 360.], deg; dB -- [-90, 90], deg; dL -- [0, 360], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int gaussBw1(enEarType earType, double dx, double dy, double dL0, double *dB, double *dL);


// 批量高斯投影坐标正反算不分带: 指定中央子午线,计算的dX dY加500km，不加带号
// imode -- 0正算， 1反算
// dB -- [-90, 90], deg; dL -- [-180, 180], deg; dL0 -- [0., 360.], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int _declspec(dllexport) gaussBatch1(enEarType earType, int imode, const char *sFnSrc, const char *sFnDst);

