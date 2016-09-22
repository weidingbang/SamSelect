#include<stdlib.h>
#include<string.h>
#include"FM.h"
#include"SamSelect.h"
#include <time.h>
#include<ctime>
#include <cmath>
#include <string>
#include <vector>
#include <getopt.h>
#include <stdlib.h>
#include<fstream>
#include<iostream>
using namespace std;

void printHelp() {
	cout
			<< "Arguments: inputFile -l <motifLength> -d <#Mutations> [-q <minMotifSequence>]  [-h]"
			<< endl;
	cout << "Where:" << endl;
	cout << "  -l: length of motif" << endl;
	cout << "  -d: max changes for planted instances of the motif" << endl;
	cout << "  -q: min sequence of strings that have motif "<< endl;
	cout << "  -h: print this help" << endl;
}

void parseArgs(int argc, char **argv, MotifConfig &mc) {
	mc.l = -1;
	mc.d = -1;
	mc.q_seq_num = -1;
	mc.inputFile = " ";
	for (int c; (c = getopt(argc, argv, "l:d:q:h")) != -1;) {
		switch (c) {
		case 'h':
			printHelp();
			exit(1);
			break;
		case 'l':
			mc.l = atoi(optarg);
			break;
		case 'd':
			mc.d = atoi(optarg);
			break;
		case 'q':
			mc.q_seq_num = atoi(optarg);
			break;
		case '?':
			cerr << "Unknown option `-" << optopt
					<< "' or missing required argument." << endl;
			exit(1);
			break;
		default:
			exit(1);
		}
	}

	if (optind < argc)
		mc.inputFile = string(argv[optind]);
	else {
		cerr << "No input file specified!" << endl;
		exit(1);
	}

	if (mc.l < 0 || mc.d < 0) {
		cerr << "Arguments -l and -d must be specified" << endl;
		exit(1);
	}
}

int main(int argc, char* argv[])
{
	MotifConfig mc;
	parseArgs(argc, argv, mc);

	SamSelect Sam;
	Sam.SeqFiltering(mc);


	return 0;
}
