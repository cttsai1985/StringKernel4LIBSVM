#include "GetLongOpt.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include <cstring>
using namespace std;


GetLongOpt::GetLongOpt(int input_argc, char** input_argv)
{
    argc = input_argc;
    argv = input_argv;
}

int GetLongOpt::CheckOpt(string option)
{
    int i;
    string temp;
    option.insert(0, "-");
    for(i = 1; i < argc; i++)
    {
        temp = argv[i];
        if(temp.compare(option) == 0)
        {
            return i;
        }
    }
    return -1;
}

int GetLongOpt::ValidOptType(string option, int index)
{
	int i, length;

	if(index < 1)
	{
		return index;
	}

	length = static_cast<int>(strlen(argv[index+1]));

    i=0;
	if(((argv[index+1][0] == '-') || (argv[index+1][0] == '+')) || (i == 0))
	{
		i = 1;
	}
	else
	{
		i = 0;
	}

	for(; i < length; i++)
	{
		if(isdigit(argv[index+1][i]) == 0)
		{
			if(argv[index+1][i] != '.')
			{
				cout << "has invalid syntex in argument: " << option << endl;
				abort();
				return -1;
			}
		}
	}
	return index;
}

/*///////////////////////////////////////////////////////////
//
//
///////////////////////////////////////////////////////////*/
int GetLongOpt::GetIntOpt(string option, int default_value)
{
    int i = GetLongOpt::CheckOpt(option);
    i = ValidOptType(option, i);
	if(i > 0)
    {
        return atoi(argv[i+1]);
    }
    else
    {
        return default_value;
    }
}

double GetLongOpt::GetFloatOpt(string option, double default_value)
{
    int i = GetLongOpt::CheckOpt(option);
	i = ValidOptType(option, i);
    if(i > 0)
    {
        return atof(argv[i+1]);
    }
    else
    {
        return default_value;
    }
}

string GetLongOpt::GetStringOpt(string option, string default_string)
{
    int i = GetLongOpt::CheckOpt(option);
    if(i > 0)
    {
        return argv[i+1];
    }
    else
    {
        return default_string;
    }
}

string GetLongOpt::GetSpecStringOpt(int index)
{
	string temp;
	if((index >= argc)||(index < 0))
		temp = argv[argc-1];
	else
		temp = argv[index];
	return temp;
}

string GetLongOpt::GetLastStringOpt()
{
	string temp;
	temp = argv[argc-1];
	return temp;
}

int GetLongOpt::GetCommandOpt(string option)
{
	int i = GetLongOpt::CheckOpt(option);
	return i;
}

int GetLongOpt::GetListFormOpt(string option, vector<string>* stringlist)
{
	int count = 0;
	int i = GetLongOpt::CheckOpt(option);

	i++;
	while((i > 0) && (i < argc))
    {
        if(argv[i][0] != '-')
		{
			stringlist->push_back(argv[i]);
			count++;

		}
		else
		{
			return count;
		}
		i++;
    }
	return count;
}

/*///////////////////////////////////////////////////////////
//utility
///////////////////////////////////////////////////////////*/
/*
bool GetLongOpt::CheckFileExist(string filename)
{
	return CheckFileExist(filename);
}


string GetLongOpt::GetMainFilename(string filename)
{
	return GetMainFilename(filename);
}
*/
/*///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////*/
bool CheckFileExist(string filename)
{
	ifstream *ifs;

	ifs = new ifstream;
	ifs->open(filename.c_str());
	if(!*ifs)//not exist
	{
		delete ifs;
		return false;
	}
	else //exitst
	{
		delete ifs;
		return true;
	}
}


string GetMainFilename(string filename)
{
	string::size_type temp;
	temp = filename.find_last_of('.');

	if(temp != string::npos)
	{
		return filename.substr(0, temp);
	}
	else
	{
		return filename;
	}
}

