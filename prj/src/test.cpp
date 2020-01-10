#include "gaussFunc.h"


int gaussBatchTest();


int main(int argc, char *argv[])
{
	gaussBatchTest();
#
	return 0;
}


/*
高斯正反算批量处理测试
输入文件（参考值）格式：

点名，B，L，x，y

其中：
B-经度（十进制，度）；
L-纬度（十进制，度）；
x-高斯平面坐标系纵坐标（十进制，米）；
y-高斯平面坐标系横坐标（十进制，米）
*/
int gaussBatchTest()
{
	double dB, dL, x, y;

	int					iType, iCount;
	double				dL0;
	stCoorTransPara		coorTranspara;
	string				sOfn, sIfn, sLine;
	ifstream			ifs;
	ofstream			ofs;
	stRawData			rawDataTmp;
	vector<stRawData>   vRawData;
	vector<string>      vsTmp;

	iType = 3;
	dL0 = 113.;				// 中央子午线

	sOfn = "out.txt";		// 结果文件
	sIfn = "dataGauss.txt"; // 输入文件（参考值）

	// 打开文件
	ofs.open(sOfn, ios::out);
	if (!ofs)
	{
		cout << "# error: fail to open file: " << sOfn << endl;
		return 0;
	}

	ifs.open(sIfn, ios::in);
	if (!ifs)
	{
		cout << "# error: fail to open file: " << sIfn << endl;
		ofs.close();
		return 0;
	}


	// 读取输入文件（参考值）
	while (!ifs.eof())
	{
		getline(ifs, sLine);

		// 使用逗号分离字段，提取数据
		vsTmp.clear();
		parseStr(sLine, ",", vsTmp);

		// 去掉两头空格
		for (vector<string>::iterator it = vsTmp.begin(); it != vsTmp.end(); it++)
			*it = fTrimAll(*it);

		
		if (vsTmp.size() == 2)   // 预留给“使用七参数进行坐标转换”
		{
			if (vsTmp[0] == "dX")
				coorTranspara.dX = atof(vsTmp[1].c_str());
			else if (vsTmp[0] == "dY")
				coorTranspara.dY = atof(vsTmp[1].c_str());
			else if (vsTmp[0] == "dZ")
				coorTranspara.dZ = atof(vsTmp[1].c_str());
			else if (vsTmp[0] == "rX" || vsTmp[0] == "m*cos")
				coorTranspara.rX = atof(vsTmp[1].c_str());
			else if (vsTmp[0] == "rY")
				coorTranspara.rY = atof(vsTmp[1].c_str());
			else if (vsTmp[0] == "rZ")
				coorTranspara.rZ = atof(vsTmp[1].c_str());
			else if (vsTmp[0] == "m" || vsTmp[0] == "m*sin")
				coorTranspara.m = atof(vsTmp[1].c_str());
		}
		else if (vsTmp.size() == 5)  // 高斯正反算
		{
			rawDataTmp.m_sN = vsTmp[0];
			rawDataTmp.m_dX1 = atof(vsTmp[1].c_str());
			rawDataTmp.m_dY1 = atof(vsTmp[2].c_str());
			rawDataTmp.m_dX2 = atof(vsTmp[3].c_str());
			rawDataTmp.m_dY2 = atof(vsTmp[4].c_str());
			vRawData.push_back(rawDataTmp);
		}
		else if (vsTmp.size() == 7)  // 预留给“坐标转换参数求取”
		{
			rawDataTmp.m_sN = vsTmp[0];
			rawDataTmp.m_dX1 = atof(vsTmp[1].c_str());
			rawDataTmp.m_dY1 = atof(vsTmp[2].c_str());
			rawDataTmp.m_dZ1 = atof(vsTmp[3].c_str());
			rawDataTmp.m_dX2 = atof(vsTmp[4].c_str());
			rawDataTmp.m_dY2 = atof(vsTmp[5].c_str());
			rawDataTmp.m_dZ2 = atof(vsTmp[6].c_str());
			vRawData.push_back(rawDataTmp);
		}
	}

	// 正算算例（经纬度转平面坐标）
	// 结果文件的文件头
	ofs << "*--------------------BL2xy--------------------*" << endl
		<< "  pt      x(before)      y(before)       x(after)       y(after)    deltaX    deltaY" << endl;
	for (vector<stRawData>::const_iterator it = vRawData.begin(); it != vRawData.end(); it++)
	{
		gaussFwWithAdd(SYSCGCS2000, it->m_dX1, it->m_dY1, dL0, &x, &y);

		ofs << fixed << it->m_sN << ", "
			<< fixed << setw(13) << setprecision(5) << it->m_dX2 << ", "
			<< fixed << setw(13) << setprecision(5) << it->m_dY2 << ", "
			<< fixed << setw(13) << setprecision(5) << x << ", "
			<< fixed << setw(13) << setprecision(5) << y << ", "
			<< fixed << setw(8) << setprecision(5) << x - it->m_dX2 << ", "
			<< fixed << setw(8) << setprecision(5) << y - it->m_dY2
			<< endl;
	}

	// 反算算例（平面坐标转经纬度）
	// 结果文件的文件头
	ofs << endl << endl 
		<< "*--------------------xy2BL--------------------*" << endl
		<< "  pt      B(before)       L(before)       B(after)        L(after)         deltaB         deltaL" << endl;
	for (vector<stRawData>::const_iterator it = vRawData.begin(); it != vRawData.end(); it++)
	{
		gaussBwWithAdd(SYSCGCS2000, it->m_dX2, it->m_dY2, dL0, &dB, &dL);

		ofs << fixed << it->m_sN << ", "
			<< fixed << setw(13) << setprecision(10) << it->m_dX1 << ", "
			<< fixed << setw(14) << setprecision(10) << it->m_dY1 << ", "
			<< fixed << setw(13) << setprecision(10) << dB << ", "
			<< fixed << setw(14) << setprecision(10) << dL << ", "
			<< fixed << setw(13) << setprecision(10) << dB - it->m_dX1 << ", "
			<< fixed << setw(13) << setprecision(10) << dL - it->m_dY1
			<< endl;
	}

	// 清理工作
	ifs.close(), ofs.close();

	return 1;
}
