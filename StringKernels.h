#ifndef __StringKernels_h
#define __StringKernels_h

#include <string>
#include <vector>

using namespace std;

enum {DSK, WDSK};

struct DataSet
{
	int size;
	int count_compute_self;
	int count_compute_kernel;//amount of training vectors
	vector<string> seq;
	vector<double> value;
	vector<double> normalized_value;
	vector<double> self;


	double** precomputed;

	void initialize();
	void normalize();

	void copy(struct DataSet* from, struct DataSet* to);
	void erase_one_char_in_string(int no);
	void erase_chars_in_string(const int* no_array);
	void printout();
	void clear();
};

struct KernelConfig
{
	int option_rbf;
	int kernel_type, kmer;
	double gamma;

	double* weighting_coef;

	KernelConfig(int k_t, int k);
    void SetKmer(int k){
        kmer=k;}
	void SetCoef();
	void SetGamma(double value){
        gamma=value;}

};

int LoadDataset(string filename, struct DataSet* T_DATASET);
int ComputeDKernel(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target);
int ComputeWDKernel(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target);

int PreComputeDKernel(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target);
int PreComputeWDKernel(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target);

int AllocateKernelOuptputArray(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target);
int WeightingKernelOuptputAB(struct KernelConfig* KC, struct DataSet* target);

int ComputeWDKernelNormalizeCoefficient(struct KernelConfig* KC, struct DataSet* target);

double compute_distance_string_kernel(string* a, string* b);
double compute_weighted_degree_string_kernel(int* k, double* weighting_coef, string* a, string* b);

double compute_radial_basis_function(double gamma, double value);

#endif //StringKernels_h
