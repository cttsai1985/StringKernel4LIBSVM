#ifndef __GetLongOpt_h
#define __GetLongOpt_h

#include <string>
#include <vector>
using namespace std;


class GetLongOpt
{
    private:
        int argc; 
        char** argv;
        
    public:
        GetLongOpt(int input_argc, char** input_argv);
        //valid process
        int CheckOpt(string option);
		int ValidOptType(string option, int index);

        //get argument
		int GetIntOpt(string option, int default_value);
        double GetFloatOpt(string option, double default_value);
        string GetStringOpt(string option, string default_string);
        int GetCommandOpt(string option);

		int GetListFormOpt(string option, vector<string>* stringlist);
 		string GetSpecStringOpt(int index);
        string GetLastStringOpt();

        //utility for code capacity old version
        //bool CheckFileExist(string filename); //true/false do/not exist 
		//string GetMainFilename(string filename);
        
		       
        //~GetLongOpt;
};

/*///////////////////////////////////////////////////////////
//utility
//
//
///////////////////////////////////////////////////////////*/

bool CheckFileExist(string filename); //true/false do/not exist 
string GetMainFilename(string filename);


#endif //__GetLongOpt_h
