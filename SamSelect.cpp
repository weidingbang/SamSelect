/*============================================
# Filename: SamSelect.cpp
# Ver 1.0 2016-08-26
# Copyright (C) 2016 weidingbang (weidingbang1992@163.com)
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include"SamSelect.h"
#include <stack>
#include <fstream>
#include <time.h>
#include <math.h>

SamSelect::SamSelect()
{

}

int SamSelect::get_charlen(const char* a)
{
	int i = 0;
	while (a[i] != '\0')
	{
		i++;
	}
	return i;

}

void SamSelect::SeqFiltering(MotifConfig mc)
{
	m_fileName = mc.inputFile.c_str();
	m_l = mc.l;
	m_d = mc.d;
	m_q = mc.q_seq_num;
	
	clock_t be = clock();
	FM *fm = NULL;
	fm=new FM(m_fileName);

	//存储第一个查询表的具体表项
	ifstream inFile(m_fileName);
	if (!inFile)//打开文件失败
		cout<<"file open failed"<<endl;
	//将文件内容读入到string类型的vector容器
	//每一行存储为该容器的一个元素
	string s;
	int count = 0;
	while (!inFile.eof())
	{
		count++;
		getline(inFile, s);
		if(count % 2 == 1)
		{
			continue;
		}
		filestring.push_back(s);
	}
	filestring[filestring.size() - 1] += " ";

	m_t = filestring.size();
	m_n = filestring[0].size() - 1;

	int index_table1_length = 8;
	alllmer_1mismatch_count = fm->CountingAllOneMisLmer(index_table1_length);


	GenSubString();
	
	//测试聚类之前正确性
	//CorrentnessBeforeAP();
	
	//聚类过程
	vector<int> temp_center;

	AP_Cluster AP(vssc, 0.00001, m_l);
	temp_result = AP.CutClusterSet(temp_center);

	AP_Result temp;
	for (int i = 0; i < temp_result.size(); ++i)
	{
		for (int j = 0; j < temp_result[i].size(); ++j)
		{
			temp.sub_result.push_back(temp_result[i][j]);
		}
		temp.center = temp_center[i];
		ap_result.push_back(temp);
		temp.sub_result.clear();
	}

	AP_Result p;
	for (int i = 0; i < ap_result.size() - 1; i++)
	{
		for (int j = 0; j < ap_result.size() - 1 - i; j++)
		{
			if (ap_result[j].sub_result.size() < ap_result[j + 1].sub_result.size())
			{
				p = ap_result[j];
				ap_result[j] = ap_result[j + 1];
				ap_result[j + 1] = p;
			}
		}
	}

	//combining after AP
	for (int i = 0; i < ap_result.size(); ++i)
	{
		int cluster_flag = 0;
		string center_str = vssc[ap_result[i].center].subStr;
		while(cluster_flag == 0)
		{
			vector<int> temp_rate;
			for (int j = i + 1; j < ap_result.size(); ++j)
			{
				float rate = HamDisOneDNum(center_str, ap_result[j].sub_result);
				if (rate > 0.8)
				{
					temp_rate.push_back(1);
				}
				else
					temp_rate.push_back(0);
			}

			int combine_flag = 0;
			for (int i = 0; i < temp_rate.size(); ++i)
			{
				if (temp_rate[i] == 1)
				{
					combine_flag = 1;
				}
			}

			int del = 0;
			if (combine_flag == 0)
			{
				cluster_flag = 1;
				break;
			}
			else
			{
				for (int j = 0; j < temp_rate.size(); ++j)
				{
					if (temp_rate[j] == 1)
					{
						for (int k = 0; k < ap_result[i + j -del + 1].sub_result.size(); ++k)
						{
							ap_result[i].sub_result.push_back(ap_result[i + j -del + 1].sub_result[k]);
						}

						int count = -1;
						for(vector<AP_Result>::iterator iter=ap_result.begin(); iter!=ap_result.end(); )
						{
							count++;
     						if(count == i + j - del + 1)
     						{
          						iter = ap_result.erase(iter);
          						del++;
     						}
      						else
            					iter ++ ;
						}
					}
					
				}
			}
			if (ap_result.size() == 1)
			{
				cluster_flag = 1;
			}


		}
	}

	//for (int i = 0; i < ap_result.size(); ++i)
	//{
		//cout<<ap_result[i].sub_result.size()<<endl;
		//for (int j = 0; j < ap_result[i].sub_result.size(); ++j)
		//{
			//cout<<ap_result[i].sub_result[j].subStr<<" ";
		//}
		//cout<<endl;
	//}

	m_top = 100;
	Rank(m_top);

	//测试聚类之后正确性
	//CorrentnessAfterAP();
	

	clock_t en = clock();
	cout<<"build time is "<<(en-be)/CLOCKS_PER_SEC<<"s"<<endl<<endl;
}

void SamSelect::Rank(int top)
{
	
}

float SamSelect::HamDisOneDNum(string str, vector<subStrContext> vs)
{
	int count = 0;
	for (int i = 0; i < vs.size(); ++i)
	{
		int minhamdis = m_l + 1;
		for (int m = 0; m < str.size() - m_l + 1; ++m)
		{
			for (int n = 0; n < vs[i].subStr.size() - m_l + 1; ++n)
			{
				int hamdis = 0;
				for (int l = 0; l < m_l; ++l)
				{
					if (str[m + l] != vs[i].subStr[n + l])
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

		//判断条件，海明距离小于等于m_d时，认为两个子串是属于同一个模体的两个模体实例
		if (minhamdis <= m_d)
		{
			count++;
		}
	}
	float rate  = (float) count / vs.size();
	return rate;
}

int SamSelect::HamDisNum(string str_a, string str_b)
{

}

int SamSelect::HamDisTwoDNum(vector<subStrContext> a, vector<subStrContext> b)
{
	int count = 0;
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < b.size(); ++j)
		{
			int minhamdis = m_l + 1;
			for (int m = 0; m < a[i].subStr.size() - m_l + 1; ++m)
			{
				for (int n = 0; n < b[j].subStr.size() - m_l + 1; ++n)
				{
					int hamdis = 0;
					for (int l = 0; l < m_l; ++l)
					{
						if (a[i].subStr[m + l] != b[j].subStr[n + l])
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

			//判断条件，海明距离小于等于m_d时，认为两个子串是属于同一个模体的两个模体实例
			if (minhamdis <= m_d)
			{
				count++;
			}
		}
	}
	return count;
}

// void SamSelect::CorrentnessAfterAP()
// {
// 	vector<int> cluster_after_MotIns_line;
// 	vector<string> string_true_MotIns_line;
// 	vector<int> int_true_MotIns_line;
// 	int temp = -1;
// 	for (int i = 0; i < m_t; ++i)
// 	{
// 		cluster_after_MotIns_line.push_back(temp);
// 	}
// 	for (int i = 0; i < ap_result.size(); ++i)
// 	{
// 		if (ap_result[i].size() >= 10)
// 		{
// 			for (int j = 0; j < ap_result[i].size(); ++j)
// 			{
// 				cluster_after_MotIns_line[ap_result[i][j].row] = 1;
// 			}
// 		}
// 	}
	
// 	string motif_fileName = "";
// 	string data_fileName(m_fileName);
// 	string temp_fileName = "RanGenData/motif/motif";

// 	for (int i = 20; i < data_fileName.size(); ++i)
// 	{
// 		motif_fileName += data_fileName[i];
// 	}
// 	motif_fileName = temp_fileName + motif_fileName;

// 	ifstream inFile_motif(motif_fileName.c_str());
// 	if (!inFile_motif)//打开文件失败
// 		cout<<"file open failed"<<endl;
// 	//将文件内容读入到string类型的vector容器
// 	//每一行存储为该容器的一个元素
// 	string motif_s;
// 	int motif_count = 0;
// 	while (!inFile_motif.eof())
// 	{
// 		motif_count++;
// 		getline(inFile_motif, motif_s);

// 		if (motif_count > 4)
// 		{
// 			if(motif_count % 2 == 1)
// 			{
// 				string_true_MotIns_line.push_back(motif_s);
// 			}
// 		}
		
// 	}

// 	for (int i = 0; i < string_true_MotIns_line.size(); ++i)
// 	{
// 		int_true_MotIns_line.push_back(atoi(string_true_MotIns_line[i].c_str()));
// 	}

// 	int true_instance_num = 0, corrent_num = 0, wrong_num = 0;
// 	for (int i = 0; i < cluster_after_MotIns_line.size(); ++i)
// 	{
// 		if ((cluster_after_MotIns_line[i] >= 0) && (int_true_MotIns_line[i] >= 0))
// 		{
// 			corrent_num++;
// 		}
// 		if ((cluster_after_MotIns_line[i] >= 0) && (int_true_MotIns_line[i] < 0))
// 		{
// 			wrong_num++;
// 		}
// 		if (int_true_MotIns_line[i] >= 0)
// 		{
// 			true_instance_num++;
// 		}
// 	}

// 	float corrent_rate = (float)corrent_num / (corrent_num + wrong_num);
// 	cout<<"file: "<<m_fileName<<endl;
// 	cout<<"High frequency substring number after AP: "<<vssc.size()<<endl;
// 	cout<<"total sequence: "<<m_t<<endl;
// 	cout<<"true motif instance num "<<true_instance_num<<endl;
// 	cout<<"counting corrent num: "<<corrent_num<<endl;
// 	cout<<"counting wrong num "<<wrong_num<<endl;
// 	cout<<"rate "<<corrent_rate<<endl;
// }

void SamSelect::CorrentnessBeforeAP()
{
	vector<int> cluster_before_MotIns_line;
	vector<string> string_true_MotIns_line;
	vector<int> int_true_MotIns_line;
	int temp = -1;
	for (int i = 0; i < m_t; ++i)
	{
		cluster_before_MotIns_line.push_back(temp);
	}
	for (int i = 0; i < vssc.size(); ++i)
	{
		cluster_before_MotIns_line[vssc[i].row] = 1;
	}
	
	string motif_fileName = "";
	string data_fileName(m_fileName);
	string temp_fileName = "RanGenData/motif/motif";

	for (int i = 20; i < data_fileName.size(); ++i)
	{
		motif_fileName += data_fileName[i];
	}
	motif_fileName = temp_fileName + motif_fileName;

	ifstream inFile_motif(motif_fileName.c_str());
	if (!inFile_motif)//打开文件失败
		cout<<"file open failed"<<endl;
	//将文件内容读入到string类型的vector容器
	//每一行存储为该容器的一个元素
	string motif_s;
	int motif_count = 0;
	while (!inFile_motif.eof())
	{
		motif_count++;
		getline(inFile_motif, motif_s);

		if (motif_count > 4)
		{
			if(motif_count % 2 == 1)
			{
				string_true_MotIns_line.push_back(motif_s);
			}
		}
		
	}

	for (int i = 0; i < string_true_MotIns_line.size(); ++i)
	{
		int_true_MotIns_line.push_back(atoi(string_true_MotIns_line[i].c_str()));
	}

	int true_instance_num = 0, corrent_num = 0, wrong_num = 0;
	for (int i = 0; i < cluster_before_MotIns_line.size(); ++i)
	{
		if ((cluster_before_MotIns_line[i] >= 0) && (int_true_MotIns_line[i] >= 0))
		{
			corrent_num++;
		}
		if ((cluster_before_MotIns_line[i] >= 0) && (int_true_MotIns_line[i] < 0))
		{
			wrong_num++;
		}
		if (int_true_MotIns_line[i] >= 0)
		{
			true_instance_num++;
		}
	}

	float corrent_rate = (float)corrent_num / (corrent_num + wrong_num);
	cout<<"file: "<<m_fileName<<endl;
	cout<<"High frequency substring number before AP: "<<vssc.size()<<endl;
	cout<<"total sequence: "<<m_t<<endl;
	cout<<"true motif instance num "<<true_instance_num<<endl;
	cout<<"counting corrent num: "<<corrent_num<<endl;
	cout<<"counting wrong num "<<wrong_num<<endl;
	cout<<"rate "<<corrent_rate<<endl;
}

void SamSelect::GenSubString()
{
	
	int occx = OccRx() + OccMx() + 1;

	//对阈值大于occx的相邻子串合并(当且仅当只有一个子串阈值大于occx，不需要合并)
	//如果子串长度大于等于l，不用处理
	//如果子串长度小于l，在左右两端分别加上l/2
	int temp = m_l/2;

	for (int i = 1; i < alllmer_1mismatch_count.size(); ++i)
	{
		for (int j = 0; j < alllmer_1mismatch_count[i].size(); ++j)
		{
			int num = 0;
			while (alllmer_1mismatch_count[i][j] > occx)
			{
				if (j < alllmer_1mismatch_count[i].size())
				{
					num++;
					j++;
				}
				else
					break;
				
			}

			if (num > 0)
			{
				//存储对应的高频子串
				ssc.row = i;
				for (int l = 0; l < 12 + num - 1; ++l)
				{
					ssc.subStr += filestring[i][j - num + l];
				}

				if (ssc.subStr.size() >= m_l)
				{
					vssc.push_back(ssc);
					ssc.subStr.clear();
				}

				else
				{

					if ((j - num + 1 < temp)||(j + 12 - 1 + temp > filestring[i].size() - 2))
					{
						if (j - num + 1 < temp)
						{
							for (int k = 0; k < j - num + 1; ++k)
							{
								ssc.subStr = filestring[i][j - num - 1 - k] + ssc.subStr;
							}
							for (int k = 0; k < temp + temp - (j - num + 1); ++k)
							{
								ssc.subStr = ssc.subStr + filestring[i][j + 12 - 1 + k];
							}
						}

						if (j + 12 - 1 + temp > filestring[i].size() - 2)
						{
							for (int k = 0; k < (filestring[i].size() - 2) - (j + 12 - 1) + 1; ++k)
							{
								ssc.subStr = ssc.subStr + filestring[i][j + 12 - 1 + k];
							}
							for (int k = 0; k < temp + temp - ((filestring[i].size() - 2) - (j + 12 - 1) + 1); ++k)
							{
								ssc.subStr = filestring[i][j - num - 1 - k] + ssc.subStr;
							}
						}
					}
					else
					{
						for (int k = 0; k < temp; ++k)
						{
							ssc.subStr = filestring[i][j - num -1 - k] + ssc.subStr;
							ssc.subStr = ssc.subStr + filestring[i][j + 12 - 1 + k];
						}
					}

					vssc.push_back(ssc);
					ssc.subStr.clear();
				}
			}
			
		}
	}
}

float SamSelect::OccRx()
{
	float pk = 0, result;

	for (int i = 0; i <= 1; ++i)
	 {
	 	float temp;
	 	temp = (Combine(12, i) * pow(3, i))/ pow(4, 12);
	 	pk += temp;
	 }

	 return m_t * (m_n - 12 + 1) * pk;
}

int SamSelect::Combine(int x, int y)
{
	int result = 1;
	for (int i = 0; i < y; ++i)
	{
		result = result * (x - i);
	}
	for (int i = 0; i < y; ++i)
	{
		result = result / (y - i);
	}

	return result;
}

float SamSelect::OccMx()
{
	//变量pk'
	float pk = 0.01;
	m_qpercent = ((float)m_q)/m_t;
	float result = m_t * m_qpercent * pk;

	return result;
}