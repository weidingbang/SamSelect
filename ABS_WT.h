/*============================================
# Filename: ABS_WT.h
# Ver 1.0 2014-06-08
# Copyright (C) 2014 ChenLonggang (chenlonggang.love@163.com)
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
  Abstract class for FM-index.the datamember root denotes
  the root bitmap of the tree.
=============================================*/
#ifndef ABS_FM_H
#define ABS_FM_H
#define CHAR_SET_SIZE 256
#define CODE_MAX_LEN 256
#include<string.h>
#include"BitMap.h"
#include"InArray.h"
#include"loadkit.h"
#include"savekit.h"
#include"divsufsort.h"
#include <vector>
#include <set>

class ABS_FM
{
	/********************************weidb*************************************/
	public:
		//第一层索引表
		struct CountingString
		{
			string str;
			int left;
			int right;
		};
		CountingString cs;
		vector<CountingString> vcs_firstindexrange;//第一个查询表的频率范围
		
	
		vector<string> vDNA;//DNA碱基
		u32 strtoint;//字符串转化为二进制数
		string fileName;//文件文件名
		vector<string> filestring;//字符串文件
		vector<vector<unsigned int> > filebinary;//二进制文件
		vector<vector<unsigned int> > filebinary_firstindex;//二进制文件,长度为length(8)
		vector<vector<u32> > filebinary_secondindex;//二进制文件，长度为16-length(8)
		vector<vector<unsigned int> > lmerbinary;//所有l-mer的二进制
		vector<vector<u32> > vfirstbinary_indexinfo;//第一个查询表的具体表项
		vector<vector<u32> > alllmer_1mismatch_count;//文本中所有二进制串的1-mismatch count

		vector<unsigned int> v_1mismatch;
		u32 u_1mismatchinfo;
		u32 u_1mismatchrow;//存储行信息，占23bit
		u32 u_1mismatchcolumn;//存储列信息，占8bit
		u32 u_1mismatchflag;//存储标记位，是否存在失配，占1bit

		//第一层索引
		void build_index_table1(int index_length);
		void DrawBackSearch_index1(CountingString & c);

		//第二层索引
		vector<u32> build_index_table2(int index_length, int left, int right);
		void DrawBackSearch_index2(CountingString & c);

		vector<vector<unsigned int> > CountingAllOneMisLmer(int index_length);
		vector<unsigned int> onemismatch(unsigned int numberbinary, int length);
		int StrtoInt(string str);
	/********************************weidb*************************************/


		//D:SA数组采样步长,Rank采样步长D*16
		ABS_FM(const char * filename,int block_size=256,int D=32);
		ABS_FM(){};
		virtual ~ABS_FM();

		void counting(const char * partten);
		int * Locating(const char * pattern,int &num);
		unsigned char* Extracting(int pos,int len);
		int Load(loadkit & s);
		int Save(savekit & s);
		int BuildTree(int speedlevel=1);

		int GetAlphabetsize(){return alphabetsize;}
		int GetN(){return n;}
		int SizeInByte();
		int SizeInByte_count();
		int SizeInByte_locate();
		int SizeInByte_extract();
	protected:
		BitMap * root;
		uchar * Z;
	 	uchar * R;

	 	const char *pat;

		void Inittable();
		//int BuildTree();

		unsigned char * T;
		unsigned char * bwt;
		unsigned char * filename;
		int  n;//file size in byte
		int block_size;//the block_size of wavelettree's node.default 1024
		int D;//step of SAL
		InArray * SAL;//the sampling of SA array
		InArray * RankL;// The sampling of the Rank Array
		
		bool charMap[256];
		// if charMap[i] is true,char[i] is a member of alphabet
		
		int * C;// Cumulative Frequency sum
		int *code;

		int charFreq[CHAR_SET_SIZE];
		int alphabetsize;
		char codeTable[CHAR_SET_SIZE][CODE_MAX_LEN];

		int Occ(int &rank,unsigned char &label,int pos);
		int Occ(unsigned char c,int pos);
		void Occ(unsigned char c,int pos_left,int pos_right,int &rank_left,int &rank_right);
		int LF(int i);
		unsigned char L(int i);
		void DrawBackSearch(const char * pattern,int &Left,int & Right);
		int Lookup(int i);

		virtual int TreeCode(){return -1;};

		int BWT(unsigned char * T,int * SA,unsigned char * bwt,int n);
		BitMap * CreateWaveletTree(unsigned char * bwt,int n);
		BitMap * FullFillWTNode(unsigned char * bwt,int len,int level);
		int DestroyWaveletTree();
		int blog(int);
		unsigned char * Getfile(const char * filename);

		int SaveNodePosition(BitMap *,u32, savekit &);
		int SaveNodeData(BitMap *,savekit &s);
		int SaveWTTree(savekit & s);
		int LoadWTTree(loadkit & s,uchar ** tables=NULL);
		
		int TreeNodeCount(BitMap * root);

		int TreeSizeInByte(BitMap * r);
		friend int GammaDecode(u64 * buff,int & index,ABS_FM * t);
		friend int Zeros(u16 x,ABS_FM * t);
		void Test_L();
		void Test_Occ();
		void Test_Shape(BitMap *);


};
#endif

