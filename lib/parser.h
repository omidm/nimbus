#ifndef _PARSER
#define _PARSER

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <boost/thread.hpp>


using namespace std;

typedef set<string> CmSet;

void parseCommand(const string str, const CmSet& cms, string& cm, vector<int>& args);

int parseCommandFile(const char * fname, CmSet& cs);






#endif
