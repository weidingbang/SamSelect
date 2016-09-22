#include "stdafx.h"
#include <memory.h>
#include <vector>
#include <string>
class GetData
{
public:
	void Input(vector<subStrContext> str);
	int DataNum;	//数据个数
	int Dimension;		//数据维数
	double *s;		//相似度矩阵
	vector<subStrContext> DataSet;	//输入数据集C
	void GetS();		//相似度矩阵
	//double *DataSet;
	double p;		//聚类参考值
	int m_l;
	double GetP();	//获取参考值
	double ComputeSim(string a, string b);
	int get_charlen(char* a);
};
struct exemplar //聚类中心类
{
	int num;
	bool count_flag;
	int count;
};
class AP_Cluster :public GetData
{
public:
	AP_Cluster(vector<subStrContext> str, double rd, int l);
	int MaxIter;	//迭代次数
	double *r;		//吸引度矩阵
	double *a;		//归属度矩阵
	exemplar *exe;  //聚类中心
	double *first_max;
	double *second_max;
	int *first_max_k1;
	double *add_max;
	double x;		//平滑系数
	void Initial();
	void Get_max();
	void Get_add_max();
	void UpdateR_A();
	void Update_R();
	void Update_A();
	void Display();
	vector<vector<subStrContext> > CutClusterSet(vector<int> &center);
};