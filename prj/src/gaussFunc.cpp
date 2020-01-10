#include "gaussFunc.h"

void printMat(const string &str, const int &iWid, const int &iFixed, ofstream &os, const Eigen::MatrixXd &mat)
{
	os << str << ":" << endl;
	for (int i = 0; i < mat.rows(); i++)
	{
		for (int j = 0; j < mat.cols(); j++)
		{
			os << fixed << setw(iWid) << setprecision(iFixed) << mat(i, j);
		}

		os << endl;
	}
	os << endl;

	return;
}

void parseStr(const string &sSrc, const string &sParser, vector<string> &vsParsed)
{
	int iPos = 0;
	string sTmp = "", s1(sSrc);
	if (sParser.length() < 1)
	{
		return;
	}

	if (sParser == " ")
	{
		stringstream ss;
		ss << s1;
		while (ss >> sTmp)
		{
			vsParsed.push_back(sTmp);
		}
	}
	else
	{
		while ((iPos = s1.find(sParser)) != string::npos)
		{
			sTmp = s1.substr(0, iPos);
			s1.erase(0, iPos + sParser.length());
			sTmp.length() > 0 ? vsParsed.push_back(sTmp) : 0;
		}
		if (vsParsed.size() > 0)
		{
			s1.length() > 0 ? vsParsed.push_back(s1) : 0;
		}
	}

	return;
}

string fTrimAll(string &sSrc)
{
	string sTmp(sSrc);
	if (sSrc.length() < 0)
	{
		return "";
	}
	if (string::npos != sTmp.find_first_not_of(" "))
	{
		sTmp.erase(0, sTmp.find_first_not_of(" "));
	}
	if (string::npos != sTmp.find_last_not_of(" "))
	{
		sTmp.erase(sTmp.find_last_not_of(" ") + 1, sTmp.length());
	}
	return sTmp;
}

// 高斯投影坐标正算不分带: 指定中央子午线,计算的dY（横坐标）不加500km，不加带号
// dB -- [-90, 90], deg; dL -- [-180, 180], deg; dL0 -- [0., 360.], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int gaussFw(enEarType earType, double dB, double dL, double dL0, double *dX, double *dY)
{
    stEarPara earPa;
    // 坐标检核
    if (dB < -90. || dB > 90. || dL < -180. || dL > 180. || dL0 < 0. || dL0 > 360.
        || (earType != SYSBJ54 && earType != SYSGDZ80 && earType != SYSWGS84 && earType != SYSCGCS2000))
    {
        return -1;
    }
    if (dL < 0.)
    {
        dL += 360.;
    }

    switch (earType)
    {
    case SYSBJ54:
        earPa = gBJ54Ear;
        break;
    case SYSGDZ80:
        earPa = gGDZ80Ear;
        break;
    case SYSWGS84:
        earPa = gWGS84Ear;
        break;
    case SYSCGCS2000:
        earPa = gCGCS2000Ear;
        break;
    default:
        break;
    }

    double dB1 = dB * D2R, dL1 = dL * D2R, dL00 = dL0 * D2R; // 转换为弧度
    double de2 = 2 * earPa.f - earPa.f * earPa.f;
    double dee = de2 * (1. - de2);
    double dN = earPa.a / sqrt(1. - de2 * sin(dB1)*sin(dB1));
    double dT = tan(dB1) * tan(dB1);
    double dC = dee * cos(dB1) * cos(dB1);
    double dA = (dL1 - dL00) * cos(dB1);
    double dM = earPa.a * ((1. - de2 / 4 - 3 * de2 * de2 / 64 - 5 * SQR3(de2) / 256) * dB1
        - (3 * de2 / 8 + 3 * de2 * de2 / 32 + 45 * SQR3(de2) / 1024) * sin(2 * dB1)
        + (15 * de2 * de2 / 256 + 45 * SQR3(de2) / 1024) * sin(4 * dB1)
        - (35 * SQR3(de2) / 3072) * sin(6 * dB1));

    *dX = dM + dN * tan(dB1) * (dA * dA / 2
        + (5 - dT + 9 * dC + 4 * dC * dC) * SQR4(dA) / 24
        + (61 - 58 * dT + dT * dT + 600 * dC - 330 * dee) * SQR6(dA) / 720);
    *dY = dN * (dA + (1. - dT + dC) * dA * dA * dA / 6
        + (5. - 18 * SQR3(dT) + 72 * dC - 58 * dee) * SQR4(dA)*dA / 120);

    return 1;
}

// 高斯投影坐标正算不分带: 指定中央子午线,计算的dY（横坐标）加500km，不加带号
// dB -- [-90, 90], deg; dL -- [-180, 180], deg; dL0 -- [0., 360.], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int gaussFwWithAdd(enEarType earType, double dB, double dL, double dL0, double *dX, double *dY)
{
	// 高斯-克吕格3°带投影，投影结果横坐标尚未进行偏移，中央子午线左侧坐标横坐标为负值
    int iRtn = gaussFw(earType, dB, dL, dL0, dX, dY);    
    if (1 != iRtn)
    {
        return iRtn;
    }

	// 遵照我国惯例，投影后的平面坐标东向偏移500km，规避横坐标为负值的情况
    *dY += 500000L;

    return 1;
}

// 高斯投影坐标正算6度带: 指定中央子午线,计算的dY（横坐标）加500km，加带号
// dB -- [-90, 90], deg; dL -- [-180, 180], deg; dL0 -- [0., 360.], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int gaussFwWithBelt(enEarType earType, double dB, double dL, double dL0, double *dX, double *dY)
{
    int iRtn = gaussFw(earType, dB, dL, dL0, dX, dY);
    if (1 != iRtn)
    {
        return iRtn;
    }

    int iBeltNum = int((dL0 + 3) / 6.);
    *dY += 1000000L * iBeltNum + 500000L;

    return 1;
}

// 高斯投影坐标反算: 指定中央子午线,输入的dy（横坐标）不加500km，不带带号
// dL0 -- [0., 360.], deg; dB -- [-90, 90], deg; dL -- [0, 360], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int gaussBw(enEarType earType, double dx, double dy, double dL0, double *dB, double *dL)
{
    stEarPara earPa;
    // 坐标检核
    if (dL0 < 0. || dL0 > 360.
        || (earType != SYSBJ54 && earType != SYSGDZ80 && earType != SYSWGS84 && earType != SYSCGCS2000))
    {
        return -1;
    }

    switch (earType)
    {
    case SYSBJ54:
        earPa = gBJ54Ear;
        break;
    case SYSGDZ80:
        earPa = gGDZ80Ear;
        break;
    case SYSWGS84:
        earPa = gWGS84Ear;
        break;
    case SYSCGCS2000:
        earPa = gCGCS2000Ear;
        break;
    default:
        break;
    }

    double de2 = 2 * earPa.f - earPa.f * earPa.f;
    double de1 = (1.0 - sqrt(1 - de2)) / (1.0 + sqrt(1 - de2));
    double dee = de2 / (1 - de2);
    double du = dx / (earPa.a * (1 - de2 / 4 - 3 * de2 * de2 / 64 - 5 * de2 * de2 * de2 / 256));
    double fai = du + (3 * de1 / 2 - 27 * de1 * de1 * de1 / 32) * sin(2 * du)
        + (21 * de1 * de1 / 16 - 55 * de1 * de1 * de1 * de1 / 32) * sin(4 * du)
        + (151 * de1 * de1 * de1 / 96) * sin(6 * du)
        + (1097 * de1 * de1 * de1 * de1 / 512) * sin(8 * du);
    double dC = dee * cos(fai) * cos(fai);
    double dT = tan(fai) * tan(fai);
    double dN = earPa.a / sqrt(1.0 - de2 * sin(fai) * sin(fai));
    double dR = earPa.a * (1 - de2) / sqrt((1 - de2 * sin(fai) * sin(fai)) * (1 - de2 * sin(fai) * sin(fai)) * (1 - de2 * sin(fai) * sin(fai)));
    double dD = dy / dN;

    //计算经度(Longitude) 纬度(Latitude)，单位：rad
    *dB = fai - (dN * tan(fai) / dR) * (dD * dD / 2 - (5 + 3 * dT + 10 * dC - 4 * dC * dC - 9 * dee) * dD * dD * dD * dD / 24
        + (61 + 90 * dT + 298 * dC + 45 * dT * dT - 256 * dee - 3 * dC * dC) * dD * dD * dD * dD * dD * dD / 720);
    *dL = (dD - (1 + 2 * dT + dC) * dD * dD * dD / 6
        + (5 - 2 * dC + 28 * dT - 3 * dC * dC + 8 * dee + 24 * dT * dT) * dD * dD * dD * dD * dD / 120) / cos(fai);

	// 单位换算为：度
    *dB = *dB * R2D;
    *dL = dL0 + *dL * R2D;
    if (*dL < 0.)
    {
        *dL += 360.;
    }
    return 1;
}

// 高斯投影坐标反算: 指定中央子午线,输入的dy（横坐标）加500km，不带带号
// dL0 -- [0., 360.], deg; dB -- [-90, 90], deg; dL -- [0, 360], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int gaussBwWithAdd(enEarType earType, double dx, double dy, double dL0, double *dB, double *dL)
{
	// 根据我国惯例，为避免横坐标为负值，高斯-克吕格平面直角坐标横坐标在投影后一般东向偏移了500km，
	// 因此，平面坐标反算经纬度时，应先去掉偏移的影响
    dy -= 500000.;
    return gaussBw(earType, dx, dy, dL0, dB, dL);
}

// 求底点纬度
double calcuBf(enEarType earType, double dx)
{
    stEarPara earPa;
    switch (earType)
    {
    case SYSBJ54:
        earPa = gBJ54Ear;
        break;
    case SYSGDZ80:
        earPa = gGDZ80Ear;
        break;
    case SYSWGS84:
        earPa = gWGS84Ear;
        break;
    case SYSCGCS2000:
        earPa = gCGCS2000Ear;
        break;
    default:
        break;
    }

    double c0 = 1 + earPa.e * earPa.e / 4 + 7 * SQR4(earPa.e) / 64 + 15 * SQR6(earPa.e) / 256
        + 579 * SQR8(earPa.e) / 16384 + 1515 * SQR10(earPa.e) / 65536
        + 16837 * SQR6(earPa.e) * SQR6(earPa.e) / 1048576
        + 48997 * SQR8(earPa.e) * SQR6(earPa.e) / 4194304
        + 9467419 * SQR8(earPa.e) * SQR8(earPa.e) / 1073741824;
    c0 = earPa.a / c0;

    double b0 = dx / c0;
    double d1 = 3 * earPa.e * earPa.e / 8 + 45 * SQR4(earPa.e) / 128 + 175 * SQR6(earPa.e) / 512
        + 11025 * SQR8(earPa.e) / 32768 + 43659 * SQR10(earPa.e) / 131072
        + 693693 * SQR6(earPa.e) * SQR6(earPa.e) / 2097152 + 10863435 * SQR8(earPa.e) * SQR6(earPa.e) / 33554432;
    double d2 = -21 * SQR4(earPa.e) / 64 - 277 * SQR6(earPa.e) / 384 - 19413 * SQR8(earPa.e) / 16384
        - 56331 * SQR10(earPa.e) / 32768 - 2436477 * SQR6(earPa.e) * SQR6(earPa.e) / 1048576
        - 196473 * SQR8(earPa.e) * SQR6(earPa.e) / 65536;
    double d3 = 151 * SQR6(earPa.e) / 384 + 5707 * SQR8(earPa.e) / 4096
        + 53189 * SQR10(earPa.e) / 163840 + 4599609 * SQR6(earPa.e) * SQR6(earPa.e) / 655360
        + 15842375 * SQR8(earPa.e) * SQR6(earPa.e) / 1048576;
    double d4 = -1097 * SQR8(earPa.e) / 2048 - 1687 * SQR10(earPa.e) / 640
        - 3650333 * SQR6(earPa.e) * SQR6(earPa.e) / 327680 - 114459079 * SQR8(earPa.e) * SQR6(earPa.e) / 27525120;
    double d5 = 8011 * SQR10(earPa.e) / 1024 + 874457 * SQR6(earPa.e) * SQR6(earPa.e) / 98304
        + 216344925 * SQR8(earPa.e) * SQR6(earPa.e) / 3670016;
    double d6 = -682193 * SQR6(earPa.e) * SQR6(earPa.e) / 245760
        - 46492223 * SQR8(earPa.e) * SQR6(earPa.e) / 1146880;
    double d7 = 36941521 * SQR8(earPa.e) * SQR6(earPa.e) / 3440640;

    double	bf;
    bf = b0 + sin(2 * b0) *
        (d1 + sin(b0) * sin(b0) *
        (d2 + sin(b0) * sin(b0) *
        (d3 + sin(b0) * sin(b0) *
        (d4 + sin(b0) * sin(b0) *
        (d5 + sin(b0) * sin(b0) *
        (d6 + d7 * sin(b0) * sin(b0)))))));
    return bf;
}

// 高斯投影坐标反算: *另一种算法*指定中央子午线,输入的dy（横坐标）加500km，不带带号，与gaussBw1精度相当
// dL0 -- [0., 360.], deg; dB -- [-90, 90], deg; dL -- [0, 360], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int gaussBw1(enEarType earType, double dx, double dy, double dL0, double *dB, double *dL)
{
    dy -= 500000.;
    double Bf = calcuBf(earType, dx);			 // 求底点纬度

    stEarPara earPa;
    switch (earType)
    {
    case SYSBJ54:
        earPa = gBJ54Ear;
        break;
    case SYSGDZ80:
        earPa = gGDZ80Ear;
        break;
    case SYSWGS84:
        earPa = gWGS84Ear;
        break;
    case SYSCGCS2000:
        earPa = gCGCS2000Ear;
        break;
    default:
        break;
    }

    double Mf = (earPa.a * (1 - earPa.e * earPa.e)) / sqrt((1 - earPa.e * earPa.e * sin(Bf) * sin(Bf))
        * (1 - earPa.e * earPa.e * sin(Bf) * sin(Bf)) * (1 - earPa.e * earPa.e * sin(Bf) * sin(Bf)));
    double Nf = earPa.a / sqrt(1 - earPa.e * earPa.e * sin(Bf) * sin(Bf));
    double tf = tan(Bf);
    double etaf_2 = earPa.e2 * earPa.e2 * cos(Bf) * cos(Bf);

    *dB = Bf - tf * SQR2(dy) / (2 * Mf * Nf) + tf * (5 + 3 * SQR2(tf)
        + etaf_2 - 9 * etaf_2 * SQR2(tf)) * SQR4(dy) / (24 * Mf * SQR3(Nf))
        - tf * (61 + 90 * SQR2(tf) + 45 * SQR4(tf)) * SQR6(dy) / (720 * Mf * SQR4(Nf)*Nf);
    *dL = dy / (Nf * cos(Bf))
        - SQR3(dy) * (1 + 2 * SQR2(tf) + etaf_2) / (6 * SQR3(Nf) * cos(Bf))
        + (5 + 28 * SQR2(tf) + 24 * SQR4(tf) + 6 * etaf_2
        + 8 * etaf_2 * SQR2(tf)) * SQR4(dy)*(dy) / (120 * SQR4(Nf)*Nf * cos(Bf));

    *dB *= R2D;
    *dL = *dL * R2D + dL0;

    return 1;
}


// 批量高斯投影坐标正反算不分带: 指定中央子午线,计算的dY（横坐标）加500km，不加带号
// imode -- 0正算， 1反算
// dB -- [-90, 90], deg; dL -- [-180, 180], deg; dL0 -- [0., 360.], deg;
// return: -1 -- 输入参数有误, 1 -- ok
// note: 输入的BL与输出的xy均属于同一椭球
int gaussBatch1(enEarType earType, int imode, const char *sFnSrc, const char *sFnDst)
{
	if (imode != 0 && imode != 1)
	{
		return -1;
	}
	if (earType != SYSBJ54 && earType != SYSGDZ80 && earType != SYSWGS84 && earType != SYSCGCS2000)
	{
		return -2;
	}
	FILE *fSrc;
	fSrc = fopen(sFnSrc, "r");
	if (NULL == fSrc)
	{
		return -3;
	}

	FILE *fDst;
	fDst = fopen(sFnDst, "w");
	if (NULL == fDst)
	{
		fclose(fSrc);
		return -4;
	}

	fprintf(fDst, "%% EarType     : %d, (0-BJ54, 1-GDZ80, 2-WGS84, 3-CGCS2000) \n", earType);
	fprintf(fDst, "%% ConvertType : %d, (0-高斯正算, 1-高斯反算) \n", imode);
	fprintf(fDst, "%% SourceFile  : %s\n", sFnSrc);
	fprintf(fDst, "%% DestFile    : %s\n", sFnDst);
	fprintf(fDst, "/********************************************************/\n");

	char sbuf[1024] = "";
	double d1, d2, d3, d4;
	int iL0, iRtn = 0;
	while (!feof(fSrc))
	{
		fgets(sbuf, 1024, fSrc);
		sscanf(sbuf, "%d,%lf,%lf", &iL0, &d1, &d2);
		if (imode == 0)
		{
			iRtn = gaussFwWithAdd(earType, d1, d2, double(iL0), &d3, &d4);
		}
		else if (imode == 1)
		{
			iRtn = gaussBwWithAdd(earType, d1, d2, double(iL0), &d3, &d4);
		}

		if (iRtn != 1)
		{
			d3 = d4 = 0.;
		}

		if (imode == 0)
		{
			fprintf(fDst, "%4d, %.8f, %.8f, %.4f, %.4f\n", iL0, d1, d2, d3, d4);
		}
		else if (imode == 1)
		{
			fprintf(fDst, "%4d, %.4f, %.4f, %.8f, %.8f\n", iL0, d1, d2, d3, d4);
		}

		d1 = d2 = d3 = d4 = iL0 = iRtn = 0;
	}


	fclose(fSrc);
	fclose(fDst);

	return 1;
}
