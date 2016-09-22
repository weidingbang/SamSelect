/*============================================
# Filename: SamSelect.h
# Ver 1.0 2016-08-26
# Copyright (C) 2016 weidingbang (weidingbang1992@163.com)
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
 
=============================================*/
#ifndef SamSelect_H
#define SamSelect_H
#include<string.h>
#include"FM.h"
#include"AP.h"
#include"BitMap.h"
#include"InArray.h"
#include"loadkit.h"
#include"savekit.h"
#include"divsufsort.h"
#include <vector>
#include <set>

class SamSelect
{
	public:
		void SeqFiltering(MotifConfig mc);
		void GenSubString();
		SamSelect();

		//by GenSubString
		float OccRx();
		float OccMx();

		//by SeqFiltering
		int get_charlen(const char* a);
		void CorrentnessBeforeAP();
		void CorrentnessAfterAP();
		int HamDisTwoDNum(vector<subStrContext> a, vector<subStrContext> b);
		float HamDisOneDNum(string str, vector<subStrContext> vs);
		int HamDisNum(string str_a, string str_b);
		void Rank(int top);


		//by OccRx();
		int Combine(int x, int y);


		
		subStrContext ssc;
		vector<subStrContext> vssc;
		vector<vector<subStrContext> > temp_result;
		vector<AP_Result> ap_result;

	private:
		int index_table1_length;//第一层索引的长度，默认值等于8
		
		vector<string> filestring;
		vector<vector<unsigned int> > alllmer_1mismatch_count;
		int start_seq;
		int end_seq;

		const char *m_fileName;
		int m_l;
		int m_d;
		int m_q;
		int m_n;
		int m_t;
		float m_qpercent;
		int m_top;
};
#endif

