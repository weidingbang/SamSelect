#include "stdafx.h"
#include <memory.h>
#include <vector>
#include <string>
class GetData
{
public:
	void Input(vector<subStrContext> str);
	int DataNum;	//���ݸ���
	int Dimension;		//����ά��
	double *s;		//���ƶȾ���
	vector<subStrContext> DataSet;	//�������ݼ�C
	void GetS();		//���ƶȾ���
	//double *DataSet;
	double p;		//����ο�ֵ
	int m_l;
	double GetP();	//��ȡ�ο�ֵ
	double ComputeSim(string a, string b);
	int get_charlen(char* a);
};
struct exemplar //����������
{
	int num;
	bool count_flag;
	int count;
};
class AP_Cluster :public GetData
{
public:
	AP_Cluster(vector<subStrContext> str, double rd, int l);
	int MaxIter;	//��������
	double *r;		//�����Ⱦ���
	double *a;		//�����Ⱦ���
	exemplar *exe;  //��������
	double *first_max;
	double *second_max;
	int *first_max_k1;
	double *add_max;
	double x;		//ƽ��ϵ��
	void Initial();
	void Get_max();
	void Get_add_max();
	void UpdateR_A();
	void Update_R();
	void Update_A();
	void Display();
	vector<vector<subStrContext> > CutClusterSet(vector<int> &center);
};