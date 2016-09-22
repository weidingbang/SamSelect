#ifndef _STDAFX_H_
#define _STDAFX_H_ 1

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct subStrContext
{
	int row;
	string subStr;
};

struct AP_Result
{
	vector<subStrContext> sub_result;
	int center;
};

struct MotifConfig {
	int l, d, q_seq_num;
	string inputFile;
};

#endif /* _STDAFX_H_*/


