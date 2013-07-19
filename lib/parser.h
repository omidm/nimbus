#include "parser.h"



void parseCommand(const string str, const CmSet& cms, string& cm, vector<int>& args)
{
  cm.clear();
  args.clear();
  set<string>::iterator it;
  for(it = cms.begin(); it != cms.end(); ++it)
  {
    if (str.find(*it) == 0)
    {
      cm = *it;
      int arg;
      stringstream ss;
      ss << str.substr(it->length(), string::npos);
      while(true)
      {
        ss >> arg;
        if (ss.fail())
          break;
        args.push_back(arg);
      }
      break;
    }
  }
  if (cm == "")
    cout << "wrong command! try again." << endl;
};





