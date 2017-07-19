#include "GetLongOpt.h"
#include "StringKernels.h"
#include <cstdio>
#include <string>

using namespace std;

//int buffer_size=65536;

//setting
int kernel_type=0;
int kmer;

int main(int argc, char** argv)
{
	int error_var=0;
	GetLongOpt GLO(argc, argv);
	kernel_type=GLO.GetIntOpt("t", 1);
	kmer=GLO.GetIntOpt("k", 5);


	string train_filename=GLO.GetStringOpt("train", "");
	string test_filename=GLO.GetStringOpt("test", "");

	KernelConfig kernel_config(kernel_type, kmer);
	kernel_config.SetCoef();

	DataSet* target_dataset;
	DataSet train_dataset, test_dataset;
	train_dataset.initialize();
	test_dataset.initialize();

	error_var=LoadDataset(train_filename, &train_dataset);
	if(error_var==-1)
		return 0;
	if(test_filename.empty()==false)
	{
		error_var=LoadDataset(test_filename, &test_dataset);
		target_dataset=&test_dataset;
		if(error_var==-1)
			return 0;
	}
	else
		target_dataset=&train_dataset;

	//PreComputeWDKernel(&kernel_config, &train_dataset, target_dataset);
	//WeightingKernelOuptputAB(&kernel_config, target_dataset);

	switch(kernel_config.kernel_type)
	{
		case  DSK:
			ComputeDKernel(&kernel_config, &train_dataset, target_dataset);
			break;
		case WDSK:
			ComputeWDKernel(&kernel_config, &train_dataset, target_dataset);
			break;
	}
    return 0;
}
