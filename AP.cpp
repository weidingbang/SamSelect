#include "BaseClass.h"
#include "stdafx.h"
#include "AP.h"
#include <map>

//sequence filtering  threshold value
#define CLUFILVAL 5

void GetData::Input(vector<subStrContext> str)
{
	int i, j;
	Dimension = 1;
	DataNum = str.size();
	/*infile.open("irisdataset.txt");
	DataSet = new double[Dimension*DataNum];
	for (i = 0; i < DataNum; i++)
		for (j = 0; j < Dimension; j++)
			infile >> DataSet[i*Dimension + j];
	infile.close();*/
	DataSet = str;
}
void GetData::GetS()
{
	int i, j, k;
	double sum;
	s = new double[DataNum*DataNum];
	for (i = 0; i < DataNum; i++)
		for (j = 0; j < DataNum; j++)
		{
			sum = 0;
			/*
			sum = ComputeSim(DataSet[i].subStr, DataSet[j].subStr);
			if (sum == 1)
			{
				s[i*DataNum + j] = INF;
			}
			else
			{
				s[i*DataNum + j] = log(sum);
			}*/
			sum = ComputeSim(DataSet[i].subStr, DataSet[j].subStr);
			//sum *= sum;
			s[i*DataNum + j] = -sum;
		}
	double key = GetP();
	for (i = 0; i < DataNum; i++)
		s[i*DataNum + i] = key;
}
double GetData::GetP()  //pÖµÈ¡ÏàËÆ¶È¾ØÕó×îÐ¡Öµ
{
	int i, j;
	
/*
	p = 0;
	for (i = 0; i < DataNum; i++)
		for (j = 0; j < DataNum; j++)
		{
			if (i != j&&s[i*DataNum + j] < p)
				p = s[i*DataNum + j];
		}*/
	
	p = 0;
	for (i = 0; i < DataNum; i++)
		for (j = 0; j < DataNum; j++)
		{
			p += s[i*DataNum + j];
		}
	p = p / (DataNum * DataNum);

/*
	
	p = -1000;
	for (i = 0; i < DataNum; i++)
		for (j = 0; j < DataNum; j++)
		{
			if (i != j&&s[i*DataNum + j] > p)
				p = s[i*DataNum + j];
		}*/
	return p;
}

int GetData::get_charlen(char* a)
{
	int i = 0;
	while (a[i] != '\0')
	{
		i++;
	}
	return i;

}

double GetData::ComputeSim(string str_a, string str_b)
{
	int minhamdis = m_l + 1;
	for (int i = 0; i < str_a.size() - m_l + 1; ++i)
	{
		for (int j = 0; j < str_b.size() - m_l + 1; ++j)
		{
			int hamdis = 0;
			for (int l = 0; l < m_l; ++l)
			{
				if (str_a[i + l] != str_b[j + l])
				{
					hamdis++;
				}
			}
			if (hamdis < minhamdis)
			{
				minhamdis = hamdis;
			}
		}
	}
	return minhamdis;
}


void AP_Cluster::Initial()
{
	int i, j;
	exe = new exemplar[DataNum];
	r = new double[DataNum*DataNum];
	a = new double[DataNum*DataNum];
	first_max = new double[DataNum];
	second_max = new double[DataNum];
	first_max_k1 = new int[DataNum];
	add_max = new double[DataNum];
	MaxIter = 100;
	for (i = 0; i < DataNum; i++)
		for (j = 0; j < DataNum; j++)
		{
			r[i*DataNum + j] = 0;
			a[i*DataNum + j] = 0;
		}
	for (i = 0; i < DataNum; i++)
	{
		exe[i].count_flag = false;
		exe[i].count = 0;
		exe[i].num = 0;
	}
}
void AP_Cluster::Get_max()
{
	int i, k1;
	for (i = 0; i < DataNum; i++)
	{
		first_max[i] = s[i*DataNum + 0] + a[i*DataNum + 0];
		first_max_k1[i] = 0;
		for (k1 = 1; k1 < DataNum; k1++)
		{
			if (s[i*DataNum + k1] + a[i*DataNum + k1] > first_max[i])
			{
				first_max[i] = s[i*DataNum + k1] + a[i*DataNum + k1];
				first_max_k1[i] = k1;
			}
		}
		for (k1 = 0; k1 < DataNum; k1++)
			if (k1 != first_max_k1[i])
			{
				second_max[i] = s[i*DataNum + k1] + a[i*DataNum + k1];
				break;
			}
		for (k1 = 1; k1<DataNum; k1++)
			if (s[i*DataNum + k1] + a[i*DataNum + k1]>second_max[i] && k1 != first_max_k1[i])
				second_max[i] = s[i*DataNum + k1] + a[i*DataNum + k1];
	}
}
void AP_Cluster::Update_R()
{
	int i, k;
	Get_max();
	for (i = 0; i < DataNum; i++)
	{
		for (k = 0; k < DataNum; k++)
		{
			if (first_max_k1[i] == k)
				r[i*DataNum + k] = x*(s[i*DataNum + k] - second_max[i]) + (1 - x)*r[i*DataNum + k];
			else r[i*DataNum + k] = x*(s[i*DataNum + k] - first_max[i]) + (1 - x)*r[i*DataNum + k];
		}
	}
}
void AP_Cluster::Get_add_max()
{
	int k, i1;
	for (k = 0; k < DataNum; k++)
	{
		add_max[k] = 0;
		for (i1 = 0; i1 < DataNum; i1++)
		{
			if (i1 != k&&r[i1*DataNum + k]>0)
				add_max[k] += r[i1*DataNum + k];
		}
	}
}
void AP_Cluster::Update_A()
{
	int i, k;
	double sum;
	Get_add_max();
	for (i = 0; i < DataNum; i++)
	{
		for (k = 0; k < DataNum; k++)
		{
			if (i == k)
				a[i*DataNum + k] = x*add_max[k] + (1 - x)*a[i*DataNum + k];
			else
			{
				if (r[i*DataNum + k]>0)
					sum = add_max[k] - r[i*DataNum + k];
				else sum = add_max[k];
				if (r[k*DataNum + k] + sum < 0)
					a[i*DataNum + k] = x*(r[k*DataNum + k] + sum) + (1 - x)*a[i*DataNum + k];
				else a[i*DataNum + k] = (1 - x)*a[i*DataNum + k];
			}
		}
	}
}
void AP_Cluster::UpdateR_A()
{
	//int i, k;
	double max;
	int t = 1;
	int m;
	while (t < MaxIter)
	{
		Update_R();
		Update_A();
		for (int i = 0; i < DataNum; i++)
		{
			if (!exe[i].count_flag)
			{
				for (int k = 0; k < DataNum; k++)
				{
					if (k == 0)
					{
						max = a[i*DataNum + k] + r[i*DataNum + k];
						m = 0;
					}
					if (k != 0 && a[i*DataNum + k] + r[i*DataNum + k] > max)
					{
						max = a[i*DataNum + k] + r[i*DataNum + k];
						m = k;
					}
				}
				if (t == 1)
					exe[i].num = m;
				if (t != 1 && exe[i].num == m)
				{
					exe[i].num = m;
					exe[i].count++;
				}
				if (t != 1 && exe[i].num != m)
				{
					exe[i].num = m;
					exe[i].count = 0;
				}
				if (exe[i].count == 50)
					exe[i].count_flag = true;
			}
		}
		t++;
	}
	//Display();
}

vector<vector<subStrContext> > AP_Cluster::CutClusterSet(vector<int> &center)
{
	vector<vector<subStrContext> > result;
	vector<vector<int> > result_number;

	for (int d = 0; d < DataNum; ++d)
	{
		int data_row = -1;
		int data_flag = 0;

		int exe_data_row = -1;
		int exe_data_flag = 0;

		for (int i = 0; i < result_number.size(); ++i)
		{
			for (int j = 0; j < result_number[i].size(); ++j)
			{
				if (d == result_number[i][j])
				{
					data_row = i;
					data_flag = 1;
				}
				if (exe[d].num == result_number[i][j])
				{

					exe_data_row =  i;
					exe_data_flag = 1;
				}
			}
		}

		//cout<<d<<" "<<data_flag<<" "<<exe_data_flag<<" "<<data_row<<" "<<exe_data_row<<endl;
		if (data_flag == 0 && exe_data_flag == 0)
		{
			vector<int> temp;
			temp.push_back(d);

			if (d != exe[d].num)
			{
				temp.push_back(exe[d].num);
			}
			result_number.push_back(temp);
		}

		if (data_flag == 1 && exe_data_flag == 0)
		{
			result_number[data_row].push_back(exe[d].num);
		}

		if (data_flag == 0 && exe_data_flag == 1)
		{
			result_number[exe_data_row].push_back(d);
		}

		if (data_flag == 1 && exe_data_flag == 1)
		{
			if (data_row != exe_data_row)
			{
				for (int i = 0; i < result_number[exe_data_row].size(); ++i)
				{
					result_number[data_row].push_back(result_number[exe_data_row][i]);
				}

				//for(it = a.begin();it != a.end(); it++)
				//{
					//if (1 == *it)
					//a.earse(it);
				//}
				int count = -1;
				for(vector<vector<int> >::iterator iter=result_number.begin(); iter!=result_number.end(); )
				{
					count++;
     				if(count == exe_data_row)
          				iter = result_number.erase(iter);
      				else
            			iter++;
				}
			}

			
		}
	}


	vector<subStrContext> temp;
	for (int i = 0; i < result_number.size(); ++i)
	{
		result.push_back(temp);
	}

	for (int i = 0; i < result_number.size(); ++i)
	{

		if (result_number[i].size() > CLUFILVAL)
		{
			map<int, int> mp;
			int maxfreqNum = 0;

			for (int j = 0; j < result_number[i].size(); ++j)
			{
				if (++mp[exe[result_number[i][j]].num] >= mp[0])
				{
					maxfreqNum = exe[result_number[i][j]].num;
				}
			}

			center.push_back(maxfreqNum);

			for (int j = 0; j < result_number[i].size(); ++j)
			{
				result[i].push_back(DataSet[result_number[i][j]]);
			}
		}	
	}

	vector<vector<subStrContext> > final;
	for (int i = 0; i < result.size(); ++i)
	{
		if (result[i].size() > CLUFILVAL)
		{
			final.push_back(result[i]);
		}
	}
	return final;
}

void AP_Cluster::Display()
{
	int i, j;
	cout << "ID" << "    " << "exemplar" << "        " << "data" << endl;
	for (i = 0; i < DataNum; i++)
	{
		cout << i << "        " << exe[i].num;
		/*for (j = 0; j < Dimension; j++)
			cout << ' ' << setw(8) << DataSet[i*Dimension + j];*/
		cout << ' ' << "        " << DataSet[i].subStr << endl;
		cout << endl;
	}
}

AP_Cluster::AP_Cluster(vector<subStrContext> str, double rd, int l)
{
	x = rd;
	m_l = l;
	Input(str);
	GetS();
	Initial();
	UpdateR_A();
	//CutClusterSet();
}
