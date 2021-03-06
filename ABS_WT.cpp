/*============================================
# Filename: ABS_WT.cpp
# Ver 1.0 2014-06-08
# Copyright (C) 2014 ChenLonggang (chenlonggang.love@163.com)
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include"ABS_WT.h"
#include<string>
#include <stack>
#include<map>
#include <fstream>
#include <time.h>
#include <math.h>

u64 GetBits(u64 * buff,int &index,int bits)
{
	if((index & 0x3f) + bits < 65)
		return (buff[index>>6] >>( 64 -((index&0x3f) + bits))) & ((0x01ull<<bits)- 1);
	int first = 64 - (index &0x3f);
	int second = bits - first;
	u64 high = (buff[index>>6] & ((0x01ull<<first)-1)) << second;
	return high + (buff[(index>>6)+1]>>(64-second));
}

int Zeros(u16 x,ABS_FM *t)
{
	if(t->Z[x>>8]==8)
		return t->Z[x>>8]+t->Z[(uchar)x];
	else
		return t->Z[x>>8];
}

int GammaDecode(u64 * buff,int & index,ABS_FM * t)
{
	u32 x = GetBits(buff,index,32);
	int runs = Zeros(x>>16,t);
	int bits = (runs<<1)+1;
	index = index + bits;
	return x>>(32-bits);
}

ABS_FM::ABS_FM(const char * filename,int block_size,int D)
{
	this->block_size = block_size;
	this->D =D;
	this->T=NULL;
	T = Getfile(filename);
	this->fileName = filename;
	Inittable();
	
	//protein
//	for (int i = 0; i < 20; i++)
//	{
//		char c = 'A' + i;
//		stringstream stream;
//		stream << c;
//		string str;
//		str = stream.str();
//		vprotein.push_back(str);
//	}

	//DNA
	string str = "A";
	vDNA.push_back(str);
	str = "C";
	vDNA.push_back(str);
	str = "G";
	vDNA.push_back(str);
	str = "T";
	vDNA.push_back(str);
}

ABS_FM::~ABS_FM()
{
	DestroyWaveletTree();
	if(T)
		delete [] T;
	if(bwt)
		delete [] bwt;
	if(SAL)
		delete SAL;
	if(RankL)
		delete RankL;
	if(C)
		delete [] C;
	if(code)
		delete [] code;
	if(Z)
		delete [] Z;
	if(R)
		delete [] R;
}

int ABS_FM::SizeInByte()
{
	return TreeSizeInByte(root) + SAL->GetMemorySize() + RankL->GetMemorySize();
}

int ABS_FM::SizeInByte_count()
{
	return TreeSizeInByte(root);
}

int ABS_FM::SizeInByte_locate()
{
	return TreeSizeInByte(root) + SAL->GetMemorySize();
}

int ABS_FM::SizeInByte_extract()
{
	return TreeSizeInByte(root) + RankL->GetMemorySize();
}

int ABS_FM::TreeNodeCount(BitMap * r)
{
	if(r==NULL)
		return 0;
	return TreeNodeCount(r->Left()) + TreeNodeCount(r->Right()) + 1;
}

int ABS_FM::TreeSizeInByte(BitMap * r)
{
	int size = 0;
	if(r->Left())
		size += TreeSizeInByte(r->Left());
	if(r->Right())
		size+=TreeSizeInByte(r->Right());
	size = size + r->SizeInByte();
	return size;
}

/********************************weidb*************************************/
vector<vector<unsigned int> > ABS_FM::CountingAllOneMisLmer(int index_length)
{
	//建立第一层索引
	build_index_table1(index_length);

	//对于first_index的每一个节点，分别建立一个长为(12 - index_length)的四叉树
	//然后求出文本中每一个l-mer的1-mismatch
	//文本中共有pow(4, index_length)个节点
	vector<u32> secondindex_count;

	//alllmer_1mismatch_count初始化
	vector<u32> temp_count;
	for (int i = 0; i < filebinary_firstindex.size(); ++i)
	{
		for (int j = 0; j < filebinary_firstindex[i].size(); ++j)
		{
			temp_count.push_back(0);
		}
		alllmer_1mismatch_count.push_back(temp_count);
		temp_count.clear();
	}

	for (int i = 0; i < vfirstbinary_indexinfo.size(); ++i)
	{
		secondindex_count = build_index_table2(12 - index_length, vcs_firstindexrange[i].left, vcs_firstindexrange[i].right);
		for (int j = 0; j < vfirstbinary_indexinfo[i].size(); ++j)
		{
			u_1mismatchflag = (vfirstbinary_indexinfo[i][j] << 31) >> 31;
			u_1mismatchrow = (vfirstbinary_indexinfo[i][j] >> 9);
			u_1mismatchcolumn = (vfirstbinary_indexinfo[i][j] << 23) >>24;

			
			if (u_1mismatchflag == 0)
			{
				alllmer_1mismatch_count[u_1mismatchrow][u_1mismatchcolumn] += secondindex_count[filebinary_secondindex[u_1mismatchrow][u_1mismatchcolumn]];
				v_1mismatch = onemismatch(filebinary_secondindex[u_1mismatchrow][u_1mismatchcolumn], 12 - index_length);


				for (int l = 0; l < v_1mismatch.size(); ++l)
				{
					alllmer_1mismatch_count[u_1mismatchrow][u_1mismatchcolumn] += secondindex_count[v_1mismatch[l]];
				}
			}
			else
			{
				alllmer_1mismatch_count[u_1mismatchrow][u_1mismatchcolumn] += secondindex_count[filebinary_secondindex[u_1mismatchrow][u_1mismatchcolumn]];
			}
		}
	}
	return alllmer_1mismatch_count;
}

void ABS_FM::build_index_table1(int index_length)
{
	//initializing vcs_firstindexrange
	cs.str = " ";
	cs.left = 0;
	cs.right = 0;

	for (int i = 0; i < (int)pow(4, index_length); i++)
	{
		vcs_firstindexrange.push_back(cs);
	}
	
	stack<CountingString> stk;
	// push A, G, C, T to stk
	for (int i = 0; i < vDNA.size(); i++)
	{
		cs.str = vDNA[i];
		cs.left = 0;
		cs.right = 0;
		DrawBackSearch_index1(cs);
		stk.push(cs);
		if (cs.str.size() == index_length)
		{
			strtoint = StrtoInt(cs.str);
			vcs_firstindexrange[strtoint] = cs;

		}
	}

	while (stk.size() != 0)
	{
		CountingString temp = stk.top();
		stk.pop();
		if (temp.str.size() < index_length)
		{
			for (int i = 0; i < vDNA.size(); i++)
			{
				cs.str = vDNA[i] + temp.str;
				cs.left = temp.left;
				cs.right = temp.right;
				DrawBackSearch_index1(cs);
				if (cs.str.size() == index_length)
				{
					strtoint = StrtoInt(cs.str);
					vcs_firstindexrange[strtoint] = cs;
				}
				stk.push(cs);
			}


		}
	}

	//存储第一个查询表的具体表项
	ifstream inFile(fileName.c_str());
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

	//将DNA文本转换成二进制
	u32 str_2bit = 0;
	u32 str_result = 0;
	vector<unsigned int> vuint;
	//存在错误：最后一个文本的长度为200，而不是200+1
	//目前解决方法，在最后一行加上一个空格，使其长度为201
	for (int i = 0; i < filestring.size(); ++i)
	{
		for (int j = 0; j < filestring[i].size() - 1; ++j)		
		{
			switch (filestring[i][j])
			{
				case 'A': // 00
		       	 	str_2bit = 0;
					break;

				case 'C': // 01
		    	   	str_2bit = 1;
					break;

				case 'G': // 10
		    	    str_2bit = 2;
					break;

				case 'T': // 11
		    	    str_2bit = 3;
					break;

				default:
					cout<<"The input sequences contain other characters except for ACGT."<<endl;
					break;

			}
			if((j + 1) % 16 != 0)
			{
				str_result = (str_result << 2) + str_2bit;
				if(j + 2 == filestring[i].size())
				{
					vuint.push_back(str_result);
					str_result = 0;
				}
				
			}
			else
			{
				str_result = (str_result << 2) + str_2bit;
				vuint.push_back(str_result);
				str_result = 0;
			}

		}
		filebinary.push_back(vuint);
		vuint.clear();		
	}

	//得到文本中所有模式的1-mismatch
	u32 patternbinary;//模式长为12的二进制字符串
	u32 currentbinary;
	vector<unsigned int> temp_firstindex;//二进制文件,长度为length(8)
	vector<unsigned int> temp_secondindex;//二进制文件，长度为16-length(8)

	for (int i = 0; i < filebinary.size(); ++i)
	{
		patternbinary = filebinary[i][0] >> 8;

		//长度为12的模式分别赋值给firstindex和secondindex
		temp_firstindex.push_back((patternbinary << 16) >> 16);
		temp_secondindex.push_back(patternbinary >>16);

		for (int j = 0; j < filebinary[i].size(); ++j)
		{
			if (j == 0)
			{
				for (int l = 0; l < (32 - 24); l = l + 2)
				{
					currentbinary = (filebinary[i][j] << (24 + l)) >> 30;
					//update patternbinary,using overlap
					patternbinary = (patternbinary << 10) >> 8;
					patternbinary = patternbinary + currentbinary;

					//长度为12的模式分别赋值给firstindex和secondindex
					temp_firstindex.push_back((patternbinary << 16) >> 16);
					temp_secondindex.push_back(patternbinary >>16);
				}
			}
			else
			{
				if (j < filebinary[i].size() - 1)
				{
					for (int l = 0; l < 32; l = l + 2)
					{
						currentbinary = (filebinary[i][j] << l) >> 30;
						//update patternbinary,using overlap
						patternbinary = (patternbinary << 10) >> 8;
						patternbinary = patternbinary + currentbinary;

						//长度为12的模式分别赋值给firstindex和secondindex
						temp_firstindex.push_back((patternbinary << 16) >> 16);
						temp_secondindex.push_back(patternbinary >>16);
					}
				}

				else
				{
					for (int l = 0; l < ((filestring[i].size() - 1) % 16) * 2; l = l + 2)
					{
						currentbinary = (filebinary[i][j] << (l + 16)) >> 30;
						//update patternbinary,using overlap
						patternbinary = (patternbinary << 10) >> 8;
						patternbinary = patternbinary + currentbinary;

						//长度为12的模式分别赋值给firstindex和secondindex
						temp_firstindex.push_back((patternbinary << 16) >> 16);
						temp_secondindex.push_back(patternbinary >>16);
					}

				}
			}
		}
		filebinary_firstindex.push_back(temp_firstindex);
		filebinary_secondindex.push_back(temp_secondindex);
		//清空vector
		temp_firstindex.clear();
		temp_secondindex.clear();
	}

	//对于每一个文本中的每一个长为index_length的模式，找到其1-mismatch的模式串
	//共有t * n(3 * l + 1)个模式

	//vfirstbinary_indexinfo初始化
	vector<u32> temp;
	for (int i = 0; i < pow(4, index_length); ++i)
	{
		vfirstbinary_indexinfo.push_back(temp);
	}

	for (int i = 0; i < filebinary_firstindex.size(); ++i)
	{
		for (int j = 0; j < filebinary_firstindex[i].size(); ++j)
		{
			u_1mismatchflag = 0;
			//行占23位，列占8位，失配标记占1位
			u_1mismatchrow = (i << 9);
			u_1mismatchcolumn = (j << 1);
			u_1mismatchinfo = u_1mismatchrow + u_1mismatchcolumn + u_1mismatchflag;			
			vfirstbinary_indexinfo[filebinary_firstindex[i][j]].push_back(u_1mismatchinfo);

			//失配标记位为1
			u_1mismatchflag = 1;
			v_1mismatch = onemismatch(filebinary_firstindex[i][j], index_length);
			for (int l = 0; l < v_1mismatch.size(); ++l)
			{
				u_1mismatchinfo = u_1mismatchrow + u_1mismatchcolumn + u_1mismatchflag;
				vfirstbinary_indexinfo[v_1mismatch[l]].push_back(u_1mismatchinfo);

			}

		}

	}

}

vector<u32> ABS_FM::build_index_table2(int index_length, int left, int right)
{
	stack<CountingString> stk;
	// push A, G, C, T to stk

	cs.str = " ";
	cs.left = 0;
	cs.right = 0;

	vector<CountingString> vcs_2;//第二个查询表的频率范围
	for (int i = 0; i < (int)pow(4, index_length); i++)
	{
		vcs_2.push_back(cs);
	}
	
	for (int i = 0; i < vDNA.size(); i++)
	{
		cs.str = vDNA[i];
		cs.left = left;
		cs.right = right;
		DrawBackSearch_index2(cs);
		stk.push(cs);
		if (cs.str.size() == index_length)
		{
			strtoint = StrtoInt(cs.str);
			vcs_2[strtoint] = cs;
		}
	}

	while (stk.size() != 0)
	{
		CountingString temp = stk.top();
		stk.pop();
		if (temp.str.size() < index_length)
		{
			for (int i = 0; i < vDNA.size(); i++)
			{
				cs.str = vDNA[i] + temp.str;
				cs.left = temp.left;
				cs.right = temp.right;
				DrawBackSearch_index2(cs);
				if (cs.str.size() == index_length)
				{
					strtoint = StrtoInt(cs.str);
					vcs_2[strtoint] = cs;
				}
				stk.push(cs);
			}


		}
	}

	vector<u32> v_count;
	for (int i = 0;  i < vcs_2.size();  i++)
	{
		if (vcs_2[i].str.size() == index_length)
		{
			v_count.push_back(vcs_2[i].right - vcs_2[i].left + 1);
		}
		
	}
	return v_count;

}

vector<unsigned int> ABS_FM::onemismatch(unsigned int numberbinary, int length)
{
	u32 currentbinary;
	u32 tempbinary;
	u32 temponemis_1;
	u32 temponemis_2;
	u32 temponemis_3;
	vector<unsigned int> v_onemismatch;
//	cout<<numberbinary<<endl;
	for (int i = 0; i < length; i++)
	{
		if (i < length - 1)
		{
			currentbinary = (numberbinary << ((32 - 2 * length) + (2 * i))) >> 30;
//			cout<<numberbinary<<endl;
			//得到一个具体模式串的所有1-mismatch结果
			tempbinary = ((currentbinary + 1) << 30) >> 30;
//			v_onemismatch.push_back((numberbinary >> (16 - 2 * i)) << (16 - 2 * i) + (numberbinary << (18 + 2 * i)) >> (18 + 2 * i) + tempbinary << (14 -2 * i));
			temponemis_1 = (numberbinary >> (2 * length - 2 * i)) << (2 * length - 2 * i);
			temponemis_2 = (numberbinary << (32 -(2 * (length - 1)) + 2 * i)) >> (32 -(2 * (length - 1)) + 2 * i);
			temponemis_3 = tempbinary << ((2 * (length - 1)) -2 * i);
			v_onemismatch.push_back(temponemis_1 + temponemis_2 + temponemis_3);

			tempbinary = ((currentbinary + 2) << 30) >> 30;
//			v_onemismatch.push_back(((numberbinary >> (16 - 2 * i)) << (16 - 2 * i)) + ((numberbinary << (18 + 2 * i)) >> (18 + 2 * i)) + tempbinary << (14 -2 * i));
			temponemis_1 = (numberbinary >> (2 * length - 2 * i)) << (2 * length - 2 * i);
			temponemis_2 = (numberbinary << (32 -(2 * (length - 1)) + 2 * i)) >> (32 -(2 * (length - 1)) + 2 * i);
			temponemis_3 = tempbinary << ((2 * (length - 1)) -2 * i);
			v_onemismatch.push_back(temponemis_1 + temponemis_2 + temponemis_3);


			tempbinary = ((currentbinary + 3) << 30) >> 30;
//			v_onemismatch.push_back(((numberbinary >> (16 - 2 * i)) << (16 - 2 * i)) + ((numberbinary << (18 + 2 * i)) >> (18 + 2 * i)) + tempbinary << (14 -2 * i));
			temponemis_1 = (numberbinary >> (2 * length - 2 * i)) << (2 * length - 2 * i);
			temponemis_2 = (numberbinary << (32 -(2 * (length - 1)) + 2 * i)) >> (32 -(2 * (length - 1)) + 2 * i);
			temponemis_3 = tempbinary << ((2 * (length - 1)) -2 * i);
			v_onemismatch.push_back(temponemis_1 + temponemis_2 + temponemis_3);
		}
		else
		{
			currentbinary = (numberbinary << ((32 - 2 * length) + (2 * i))) >> 30;
//			cout<<numberbinary<<endl;
			//得到一个具体模式串的所有1-mismatch结果
			tempbinary = ((currentbinary + 1) << 30) >> 30;
//			v_onemismatch.push_back((numberbinary >> (16 - 2 * i)) << (16 - 2 * i) + (numberbinary << (18 + 2 * i)) >> (18 + 2 * i) + tempbinary << (14 -2 * i));
			temponemis_1 = (numberbinary >>2 ) << 2;
			v_onemismatch.push_back(temponemis_1 + tempbinary);

			tempbinary = ((currentbinary + 2) << 30) >> 30;
//			v_onemismatch.push_back(((numberbinary >> (16 - 2 * i)) << (16 - 2 * i)) + ((numberbinary << (18 + 2 * i)) >> (18 + 2 * i)) + tempbinary << (14 -2 * i));
			temponemis_1 = (numberbinary >>2 ) << 2;
			v_onemismatch.push_back(temponemis_1 + tempbinary);


			tempbinary = ((currentbinary + 3) << 30) >> 30;
//			v_onemismatch.push_back(((numberbinary >> (16 - 2 * i)) << (16 - 2 * i)) + ((numberbinary << (18 + 2 * i)) >> (18 + 2 * i)) + tempbinary << (14 -2 * i));
			temponemis_1 = (numberbinary >>2 ) << 2;
			v_onemismatch.push_back(temponemis_1 + tempbinary);
		}
		
	}
	return v_onemismatch;
}

//第二层索引后向搜索
void ABS_FM::DrawBackSearch_index1(CountingString & c)
{
	int len = c.str.size();
	int occ_left=0;
	int occ_right=0;
	
	int i = len -1;
	//当len = 1时
	if ((c.left <= c.right) and (i == 0))
	{
		unsigned char uc = c.str[i];//字符c为每次需要匹配的字符
		int coding = code[uc];
		if (coding == -1)
		{
			c.left = 1;
			c.right = 0;
			return;
		}
		c.left = C[coding];
		c.right = C[coding + 1] - 1;
	}

	
	//当len > 1时
	if ((c.left <= c.right) and (i >= 1))
	{
		unsigned char uc = c.str[0];
//		cout<< c.str <<" "<<uc<<" "<<c.left<<" "<<c.right<<endl;
		int coding = code[uc];
		if(coding == -1)
		{
			c.left = 1;
			c.right = 0;
			return;
		}
	/*	
	    Left = C[coding]+Occ(c,Left-1);
		Right =C[coding]+Occ(c,Right)-1;
		i=i-1;
	*/
		Occ(uc,c.left - 1, c.right, occ_left,occ_right);
		c.left = C[coding]+occ_left;
		c.right = C[coding]+occ_right-1;
	}
	if(c.right < c.left)
	{
		c.left = 1;
		c.right = 0;
		return ;
	}
	return;
} 

//第二层索引后向搜索
void ABS_FM::DrawBackSearch_index2(CountingString & c)
{
	int len = c.str.size();
	int occ_left=0;
	int occ_right=0;
	
	int i = len -1;
	//当len = 1时
	if ((c.left <= c.right) and (i == 0))
	{
		unsigned char uc = c.str[i];//字符c为每次需要匹配的字符
		int coding = code[uc];
		if (coding == -1)
		{
			c.left = 1;
			c.right = 0;
			return;
		}
		Occ(uc,c.left - 1, c.right, occ_left,occ_right);
		c.left = C[coding] + occ_left;
		c.right = C[coding] + occ_right-1;
	}

	//当len > 1时
	if ((c.left <= c.right) and (i >= 1))
	{
		unsigned char uc = c.str[0];
		int coding = code[uc];
		if(coding == -1)
		{
			c.left = 1;
			c.right = 0;
			return;
		}
	/*	
	    Left = C[coding]+Occ(c,Left-1);
		Right =C[coding]+Occ(c,Right)-1;
		i=i-1;
	*/
		Occ(uc,c.left - 1, c.right, occ_left,occ_right);
		c.left = C[coding]+occ_left;
		c.right = C[coding]+occ_right-1;
	}
	if(c.right < c.left)
	{
		c.left = 1;
		c.right = 0;
		return ;
	}
	return;
} 

int ABS_FM::StrtoInt(string str)
{
	u32 str_2bit = 0;
	u32 str_result = 0;
	for (int i = 0; i < str.size(); ++i)
	{
		switch (str[i])
		{
			case 'A': // 00
		        str_2bit = 0;
				break;

			case 'C': // 01
		        str_2bit = 1;
				break;

			case 'G': // 10
		        str_2bit = 2;
				break;

			case 'T': // 11
		        str_2bit = 3;
				break;

			default:
				cout<<"The input sequences contain other characters except for ACGT."<<endl;
				break;
		}

		str_result = (str_result << 2) + str_2bit;
	}
	return str_result;
}
/********************************weidb*************************************/

void ABS_FM::counting(const char * pattern)
{

}


/*
 * this functon is BackSearch
*/
void ABS_FM::DrawBackSearch(const char * pattern,int & Left,int &Right)
{
	int len = strlen(pattern);
	int occ_left=0;
	int occ_right=0;
	if(len <=0)
	{
		Left =1;
		Right = 0;
		return;
	}
	int i = len -1;
	unsigned char c = pattern[i];
	int coding = code[c];
	if (coding==-1)//coding =-1 is stand for the character is not exit in the zifubiao
	{
		Left = 1;
		Right = 0;
		return ;
	}
	Left = C[coding];
	Right = C[coding+1]-1;
	i=i-1;
	while ((Left <= Right) and (i>=0))
	{
		c = pattern[i];
		coding = code[c];
		if(coding == -1)
		{
			Left = 1;
			Right = 0;
			return;
		}
	/*	
	    Left = C[coding]+Occ(c,Left-1);
		Right =C[coding]+Occ(c,Right)-1;
		i=i-1;
	*/
		Occ(c,Left-1,Right,occ_left,occ_right);
		Left = C[coding]+occ_left;
		Right = C[coding]+occ_right-1;
		i=i-1;
	}
	if(Right < Left)
	{
		Left = 1;
		Right = 0;
		return ;
	}
	return;
}


int * ABS_FM::Locating(const char * pattern,int &num)
{
	int Left=1;
	int Right = 0;
	DrawBackSearch(pattern,Left,Right);
	if(Right < Left)
		return NULL;
	num = Right - Left + 1;
	int *pos =new int[num];
	for(int i=0;i<num;i++)
		pos[i]=Lookup(Left + i);
	return pos;
}


unsigned char* ABS_FM::Extracting(int pos,int len)
{
	if(pos + len > n-1 || pos <0)
	{
		cout<<pos<<"  "<<len<<endl;
		cout<<pos+len<<" "<<n-1<<" "<<pos<<endl;
		cout<<"ABS_FM::Extracting  error parmaters"<<endl;
		return NULL;
	}

	unsigned char * sequence=new unsigned char[len+1];
	sequence[len]='\0';
	int end = pos + len -1;
	int anchor = 0;
	int overloop = 0;
	int step = this->D*16;
	overloop = (n-2-end)%step;
	anchor = (n-2-end)/step;

	int i= RankL->GetValue(anchor);
	for(int j=0;j<overloop;j++)
		i = LF(i);

	for(int j=0;j<len;j++)
	{
		sequence[len-1-j]=L(i);
		i = LF(i);
	}
	return sequence;
}


int ABS_FM::Lookup(int i)
{
	int step = 0;
	int D = this->D;
	while(i%D!=0)
	{
		i=LF(i);
		step =step +1;
	}
	i=i/D;
	return (SAL->GetValue(i)+step)%n;
}

//返回L串中c字符在位置pos_left 和pos_right之前出现的次数，结果由rank_left 和rank_right带回.
void ABS_FM::Occ(unsigned char c,int pos_left,int pos_right,int &rank_left,int &rank_right)
{
	BitMap *r = root;
	int level=0;
	char code = '0';
	while(r->Left())
	{
		code = codeTable[c][level];
		
		if(code == '1')//编码是1,走右分支
		{
			if(pos_left>-1 && pos_right >-1) //left right 都有待查找
			{
			
				r->Rank(pos_left,pos_right,rank_left,rank_right);
				pos_left = rank_left -1;
				pos_right = rank_right -1;
			
			/*	pos_left = r->Rank(pos_left)-1;
				pos_right=r->Rank(pos_right)-1;
			*/
			}
			else if(pos_right > -1)//只查右分支
			{
				pos_right=r->Rank(pos_right)-1;
			}
			else//该跳出循环了,此时pos_left 和pos_right都是-1.
			{
				break;
			}
			//====kkzone-bai=======
			if(pos_right==pos_left)
			{
				rank_left = pos_left+1;
 				rank_right= pos_right+1;
				return ;
			}
			//====kkzone-bai=======
			r= r->Right();
		}
		else //编码是0,走左分支.
		{
			if(pos_left>-1 && pos_right >-1)
			{
				
				r->Rank(pos_left,pos_right,rank_left,rank_right);
				pos_left = (pos_left+1) - rank_left-1;
				pos_right= (pos_right+1)- rank_right-1;
			/*
				pos_left = (pos_left+1)-r->Rank(pos_left)-1;
				pos_right= (pos_right+1)-r->Rank(pos_right)-1;
			*/
			}
			else if(pos_right > -1)
			{
				pos_right = (pos_right+1)-r->Rank(pos_right)-1;
			}
			else
			{
				break;
			}
			r=r->Left();
			//====kkzone-bai=======
			if(pos_right==pos_left)
			{
				rank_left = pos_left+1;
				rank_right= pos_right+1;
				return ;
			}
			//====kkzone-bai=======
		}
		level++;
	}
	rank_left = pos_left+1;
	rank_right= pos_right+1;
	return ;
}
int ABS_FM::Occ(unsigned char c,int pos)
{
	BitMap * r = root;
	int level = 0;
	char code ='0';
	while(r->Left() && pos > -1)
	{
		code = codeTable[c][level];
		if(code == '1')
		{
			pos = r->Rank(pos) - 1;
			r = r->Right();
		}
		else
		{
			pos = (pos+1) - r->Rank(pos)-1;
			r = r->Left();
		}
		level = level +1;
	}
	return pos+1;
}


int ABS_FM::LF(int i)
{
/*
	unsigned char c = L(i);
//	cout<<"L(i)  "<<c<<endl;zerostable = zerostable;
	int coding = code[c];
	return C[coding]+Occ(c,i)-1;
*/
	
	int occ =0;
	unsigned char label =0;
	Occ(occ,label,i);
	int coding = code[label];
	return occ + C[coding] -1;

}


unsigned char ABS_FM::L(int i)
{
	BitMap * r = root;
	int bit =0;
	int rank = 0;
	
	while(r->Left())
	{
	//	bit = r->Rank(i) - r->Rank(i-1);
		rank = r->Rank(i,bit);
		if(bit ==1)
		{
			//i = r->Rank(i)-1;
			i = rank -1;
			r=r->Right();
		}
		else
		{
			//i = (i+1) - r->Rank(i)-1;
			i = (i+1) - rank -1;
			r=r->Left();
		}
	}
	return r->Label();

}

int ABS_FM::Occ(int & occ , unsigned char & label,int pos)
{
	BitMap * r = root;
	int bit =0;
	int rank =0;
	while(r->Left())
	{
		rank = r->Rank(pos,bit);
		if(bit==1)
		{
			pos = rank -1;
			r =r ->Right();
		}
		else
		{
			pos = (pos+1) - rank -1;
			r = r->Left();
		}
	}
	occ = pos +1;
	label = r->Label();
	return 0;
}


unsigned char * ABS_FM::Getfile(const char *filename)
{
	FILE * fp = fopen(filename,"r+");
	if(fp==NULL)
	{
		cout<<"Be sure the file is available "<<filename<<endl;
		exit(0);
	}
	fseek(fp,0,SEEK_END);
	this->n = ftell(fp)+1;
	unsigned char * T = new unsigned char[n];
	fseek(fp,0,SEEK_SET);

	int e=0;
	int num=0;
	while((e=fread(T+num,sizeof(uchar),n-1-num,fp))!=0)
		num = num +e;
	if(num!=n-1)
	{
		cout<<"Read source file failed"<<endl;
		exit(0);
	}
	T[n-1]=0;
	fclose(fp);

	memset(charFreq,0,256*sizeof(int));
	memset(charMap,0,256*sizeof(bool));
	for(int i=0;i<n;i++)
		charFreq[T[i]]++;
	this->alphabetsize = 0;
	for(i32 i=0;i<256;i++)
		if(charFreq[i])
		{
			this->alphabetsize++;
			this->charMap[i]=true;
		}
	this->code = new int[256];//kkzone-bai debug
	this->C = new int[alphabetsize+1];
	memset(C,0,(alphabetsize+1)*4);
	this->C[alphabetsize] = n;
	this->C[0] = 0;
	int k=1;
	i32 pre =0;
	for(int i=0;i<256;i++)
	{
		if(charFreq[i]!=0)
		{
			code[i]=k-1;
			C[k]=pre + charFreq[i];
			pre = C[k];
			k++;
		}
		else
			code[i]=-1;
	}
	return T;
}


int ABS_FM::BWT(unsigned char *T,int * SA,unsigned char * bwt,int len)
{
	int i=0;
	int index=0;
	for(i=0;i<len;i++)
	{
		index = (SA[i]-1+len)%len;
		bwt[i]=T[index];
	}
	return 0;
}


int ABS_FM::BuildTree(int speedlevel)
{
	int *SA = new int[n];
	divsufsort(T,SA,n);		
	//SA和Rank数组的采样
	int step1 =this->D;
	int step2 =this->D*16;
	SAL=new InArray(n/step1+1,blog(n));
	RankL=new InArray(n/step2+1,blog(n));

	int i=0;
	int j=0;
	for(i=0,j=0;i<n;i=i+step1,j++)
		SAL->SetValue(j,SA[i]);
	
	for(i=0;i<n;i++)
	{
		if(SA[i]==0)
			continue;
		if((n-2-(SA[i]-1))%step2 == 0)
		{
			RankL->SetValue((n-2-(SA[i]-1))/step2,i);
		}
	}

	bwt = new unsigned char[n];
	BWT(T,SA,bwt,n);

	double runs=0.0;
	for(int i=0;i<n-1;i++)
		if(bwt[i]!=bwt[i+1])
			runs++;
	runs=n/runs;
//	cout<<runs<<endl;
	int a=0;
	int b=0;
	if(speedlevel<0 || speedlevel >2)
	{
		cerr<<"speedlevel error"<<endl;
		exit(0);
	}
	switch(speedlevel)
	{
		case 0:a=2;b=10;break;
		case 1:a=4;b=20;break;
		case 2:a=10;b=50;break;
		default:a=4;b=20;break;
	}
	
	if(runs<a)
		block_size=block_size*1;
	else if(runs<b)
		block_size=block_size*2;
	else
		block_size=block_size*4;
	
//	cout<<"block_size: "<<block_size<<endl;
	TreeCode();
	root=CreateWaveletTree(bwt,n);
//	cout<<"CreatWaveletTree"<<endl;
	
	//Test_L();
	//Test_Occ();
	//Test_Shape(root);

	delete [] T;
	T=NULL;
	delete [] SA;
	SA=NULL;
	delete[] bwt;
	bwt=NULL;

	return 0;
}

void ABS_FM::Test_Shape(BitMap * r)
{
	if(r->Left() && r->Right())
	{
		Test_Shape(r->Left());
		Test_Shape(r->Right());
	}
	else if(r->Left() || r->Right())
	{
		cout<<"one child"<<endl;
	}
}
	
void ABS_FM::Test_Occ()
{
	int count[256]={0};
	unsigned long long int mis=0;
	for(int i=0;i<n;i++)
	{
		count[bwt[i]]++;
		if(Occ(bwt[i],i) != count[bwt[i]])
		{
			cout<<count[bwt[i]]<<" "<<Occ(bwt[i],i)<<" "<<(int)bwt[i]<<" "<<(char)bwt[i]<<" "<<i<<endl;
			cout<<codeTable[bwt[i]]<<endl;
			mis++;
		}
		if(i%100000==0)
			cout<<(i*1.0)/n<<endl;
	}
	cout<<"missing :"<<mis<<endl;

}

void ABS_FM::Test_L()
{
	int mis=0;
	for(int i=0;i<n;i++)
	{
		if(bwt[i]!=L(i))
		{
			cout<<bwt[i]<<" "<<L(i)<<" "<<i<<endl;
			mis++;;
		}
	}
	cout<<"mis: "<<mis<<endl;
}


BitMap * ABS_FM::CreateWaveletTree(unsigned char * bwt,int n)
{
	BitMap * root = NULL;

	root = FullFillWTNode(bwt,n,0);
	if(!root)
	{
		cout<<"FullfillWTNode failed"<<endl;
		exit(0);
	}
	return root;
}


BitMap * ABS_FM::FullFillWTNode(unsigned char * buff,int len,int level)
{
//	cout<<level<<endl;
	int CurrentLevel = level;
	unsigned int CurrentBitLen = len;
	unsigned char CurrentLabel = '\0';
	unsigned long long int *CurrentBitBuff = NULL;
	if ((int)strlen((const char*)codeTable[buff[0]])==level)
	{
		CurrentLabel = buff[0];
		CurrentBitBuff = NULL;
		//uchar * tables[5]={this->zerostable,this->R1,this->R2,this->R3,this->R4};
		uchar * tables[2] ={this->Z,this->R};
		BitMap * node = new BitMap(CurrentBitBuff,CurrentBitLen,CurrentLevel,block_size,CurrentLabel,tables);
		node->Left(NULL);
		node->Right(NULL);
		return node;
	}
	
	int u64Len=0;
	if(len%64==0)
		u64Len = len/64+1;
	else
		u64Len = len/64+2;
	CurrentBitBuff = new unsigned long long int[u64Len];
	memset(CurrentBitBuff,0,u64Len*8);
	unsigned char * lptr=NULL;
	unsigned char * rptr=NULL;
	int leftLen=0;
	int rightLen=0;

	lptr = new unsigned char[len];
	rptr = new unsigned char[len];
	memset(lptr,0,len);
	memset(rptr,0,len);

	//computer bitvect;

	int i=0;
	int bytePos=0;
	int bitOffset=0;
	u64 last = 0;
	for(i=0;i<len;i++)
	{
		if(codeTable[buff[i]][level]=='1')
		{
			//CurrentBitBuff[bytePos] |= (0x01<<(7-bitOffset));
			CurrentBitBuff[bytePos] |= (0x01ull<<(63-bitOffset));
			//construct right data buff
			rptr[rightLen++]=buff[i];
			last = 0;
		}
		else
		{
			lptr[leftLen++]=buff[i];
			last = 1;
		}
		bitOffset++;
		//if(bitOffset == 8)
		if(bitOffset == 64)
		{
			bytePos++;
			bitOffset = 0;
		}
	}
	CurrentBitBuff[bytePos] |= (last<<(63-bitOffset));

	//uchar * tables[5] = {this->zerostable,this->R1,this->R2,this->R3,this->R4};
	uchar * tables[2] = {this->Z,this->R};
	BitMap * node = new BitMap(CurrentBitBuff,CurrentBitLen,CurrentLevel,block_size,CurrentLabel,tables);

	if(leftLen !=0)
	{
		BitMap * left =FullFillWTNode(lptr,leftLen,level+1);
		node->Left(left);
		delete [] lptr;
		lptr=NULL;
	}
	if(rightLen!=0)
	{
		BitMap * right = FullFillWTNode(rptr,rightLen,level+1);
		node->Right(right);
		delete [] rptr;
		rptr=NULL;
	}
	return node;
}


int ABS_FM::DestroyWaveletTree()
{
	delete root ;
	root=NULL;
	return 0;
}


int ABS_FM::blog(int x)
{
	int ans=0;
	while(x>0)
	{
		ans++;
		x=(x>>1);
	}
	return ans;
}


void ABS_FM::Inittable()
{
	this -> Z = new uchar[1<<8];
	int tablewidth = 8;
	for(int i=0;i<tablewidth;i++)
		for(int j=(1<<i);j<(2<<i);j++)
			Z[j] = tablewidth-1-i;
	Z[0] = tablewidth;
	
	u64 tablesize = (1<<16);
	R  = new uchar[tablesize<<2];
	
	//查找表的初始化：在16bits的0,1串上模拟gamma解码的过程，得到
	//这些表
	u64 B[2]={0xffffffffffffffffull,0xffffffffffffffffull};
	int sum =0;//gamma编码的和,含义为原串被编码的bits数目。
	int step=0;//16bits可以完整解码的bits数,该值不会大于16.
	int rank = 0;//16bits表示的几个完整的runs,假设第一个runs是1-runs,这几个runs的rank值。
	int runs = 0 ;//runs 个数.
	
	int x = 0;//工作变量，保存本次的gamma解码值.
	int prestep = 0;//前一次正确解码的bits数(累加),<=16.
	for(u64 i=0;i<tablesize;i++)
	{
		B[0] = (i<<48);
		sum  =0 ;
		step = 0;
		prestep=0;
		rank = 0;
		runs = 0;
		while(1)
		{
			x = GammaDecode(B,step,this);//step会联动.
			if(step > 16)
				break;
			sum = sum + x;
			prestep = step;
			runs ++;
			if(runs%2==1)
				rank = rank + x;
		}
		R[i<<2] = runs;//r4
		R[(i<<2)+1] = prestep;//r2
		R[(i<<2)+2] = sum; //r1;
		R[(i<<2)+3] = rank;//r3

	}
}

//递归保存节点的编号信息
int ABS_FM::SaveNodePosition(BitMap * r,u32 position,savekit &s)
{
	if(!r)
		return 1;
	s.writei32(position);
	SaveNodePosition(r->Left(), 2 * position,s);
	SaveNodePosition(r->Right(),2 * position +1,s);
	return 0;
}

//递归保存节点的数据信息
int ABS_FM::SaveNodeData(BitMap *r,savekit &s)
{
	if(!r)
		return 1 ;
	r->Save(s);
	SaveNodeData(r->Left(),s);
	SaveNodeData(r->Right(),s);
	return 0;
}

int ABS_FM::SaveWTTree(savekit &s)
{
	//保存编号信息
	//int nodecount = 2*alphabetsize -1;
	//s.writei32(nodecount);
	SaveNodePosition(root,1,s);

	//保存节点数据信息
	SaveNodeData(root,s);
	return 0;
}

int ABS_FM::LoadWTTree(loadkit &s,uchar **tables)
{
	//读取数据，map的int域对应该节点的位置
	int nodecount = 2*alphabetsize -1;
//	cout<<alphabetsize<<endl;
//	s.loadi32(nodecount);
	int * p = new int[nodecount];
	s.loadi32array(p,nodecount);
	map<int,BitMap * > pmap;
	BitMap * r=NULL;
	for(int i=0;i<nodecount;i++)
	{
		if(tables)
			r = new BitMap(tables);
		else
			r = new BitMap();
		r->Load(s);
		pmap[p[i]] = r;
	}
	//挂链
	map<int ,BitMap *>::iterator iter;
	map<int ,BitMap *>::iterator f_iter;
	for(iter = pmap.begin();iter!=pmap.end();iter++)
	{
		f_iter = pmap.find(2*iter->first);
		if(f_iter != pmap.end())
			iter->second->Left(f_iter->second);
		else
			iter->second->Left(NULL);
		
		f_iter = pmap.find(2*iter->first +1);
		if(f_iter!=pmap.end())
			iter->second->Right(f_iter->second);
		else
			iter->second->Right(NULL);
	}
//	cout<<"767"<<endl;	
	f_iter = pmap.find(1);
	if(f_iter !=pmap.end())
		this->root = f_iter->second;
	else
	{
		cerr<<"Load WTTree error"<<endl;
		this->root = NULL;
		exit(0);
	}
//	cout<<"778"<<endl;
	return 0;
}

int ABS_FM::Load(loadkit &s)
{
	s.loadi32(this->n);
	s.loadi32(this->alphabetsize);
	s.loadi32(this->D);

	//for C
	this->C = new int[alphabetsize+1];
	s.loadi32array(this->C,alphabetsize+1);
	//for code
	this->code = new int[256];//kkzone-bai debug
	s.loadi32array(this->code,256);
	//for codeTable;
	memset(codeTable,0,sizeof(codeTable));
	for(int i=0;i<256;i++)
	{
		uchar len=0;
		s.loadu8(len);
		if(len!=0)
		{
			int bytes = len%8?len/8+1:len/8;
			uchar * bits = new uchar[bytes];
			s.loadu8array(bits,bytes);
			int in_index =0;
			int off_index =0;
			for(int j=0;j<len;j++)
			{
				if(bits[off_index] & (0x01<<(7-in_index)))
					codeTable[i][j] = '1';
				else
					codeTable[i][j] = '0';
				in_index++;
				if(in_index==8)
				{
					in_index =0;
					off_index ++;
				}
			}
		}
	}
	
	//for SAL
	this->SAL = new InArray();
	this->SAL->load(s);
	//for Rankl
	this->RankL = new InArray();
	this->RankL->load(s);

	Inittable();
	uchar * par[2]={Z,R};
	//cout<<"cs"<<endl;
	LoadWTTree(s,par);
//	cout<<"835"<<endl;
	T=NULL;
	bwt=NULL;
	return 0;
}
int ABS_FM::Save(savekit &s)
{
	s.writei32(n);
	s.writei32(alphabetsize);
	s.writei32(D);//SA的采样率
	
	//C表
	//s.writei32(alphabetsize+1);
	s.writei32array(C,alphabetsize+1);
	//code表
	//s.writei32(256);
	s.writei32array(code,256);//kkzone-bai debug
	
	//codeTable
	for(int i=0;i<256;i++)
	{
		uchar len = strlen(codeTable[i]);
		s.writeu8(len);
		if(0!=len)
		{
			int bytes = len%8?len/8+1:len/8;
			uchar *bits = new uchar[bytes];
			memset(bits,0,bytes);
			int off_index=0;
			int in_index =0;
			for(int j=0;j<len;j++)
			{
				if(codeTable[i][j]=='1')
					bits[off_index] = bits[off_index]|(0x01<<(7-in_index));
				in_index++;
				if(8==in_index)
				{
					in_index=0;
					off_index++;
				}
			}
			s.writeu8array(bits,bytes);
		}
	}
	
	//for SAL
	SAL->write(s);
	//for RankL
	RankL->write(s);
	//for WT tree
//	cout<<"SaveWTTree"<<endl;
	SaveWTTree(s);
	return 0;
}
