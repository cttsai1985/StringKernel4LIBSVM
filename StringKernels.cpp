#include "StringKernels.h"
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

int buffer_size=65536;

void DataSet::initialize()
{
	count_compute_self=0;
	count_compute_kernel=0;
}

void DataSet::normalize()
{
	int i;
	double max, min;
	max=min=value[0];
	for(i=1; i< size; i++)
		if(max<value[i])
			max=value[i];
		else if(min>value[i])
			min=value[i];
	for(i=0; i<size; i++)
		normalized_value.push_back((value[i]-min)/(max-min));
}

void DataSet::erase_one_char_in_string(int no)
{
    for(int i=0; i<size; i++)
        seq[i].erase(no, 1);
}

void DataSet::erase_chars_in_string(const int* no_array)
{
	int i, j, seq_length;
	//if==1 mask-out(delete), ==0, keeped
    seq_length=(int)(seq[0].length());
        for(i=0; i<size; i++)
            for(j=(seq_length-1); j>=0; j--)
                if(no_array[j]==1)
                    seq[i].erase(j, 1);
}

void DataSet::printout()
{
	int i;
	for(i=0; i<size; i++)
		printf("%s\t%lf\n",seq[i].c_str(),value[i]);
}

void DataSet::copy(struct DataSet* from, struct DataSet* to)
{
	*to=*from;
}

void DataSet::clear()
{
	seq.clear();
	value.clear();
	normalized_value.clear();
	self.clear();

	for(int i=0; i< size; i++)
		delete[] precomputed[i];
	delete[] precomputed;

}

KernelConfig::KernelConfig(int k_t, int k)
{
	kernel_type=k_t;
	SetKmer(k);
	SetGamma(1.0);
	option_rbf=0;//off
}

void KernelConfig::SetCoef()
{
	int i;
	weighting_coef=new double[kmer+1];

	weighting_coef[0]=0;
	for(i=1; i <= kmer; i++)
		weighting_coef[i]=(2.0*(double)(kmer-i+1)/double((kmer*(kmer+1))));

/*
	for(i=1; i <= kmer_a; i++)
		printf("%d: %lf\n", i, weighting_coef_a[i]);
	for(i=1; i <= kmer_b; i++)
		printf("%d: %lf\n", i, weighting_coef_b[i]);
*/
}


int LoadDataset(string filename, struct DataSet* T_DATASET)
{
	char* buffer=new char[buffer_size];
	char* temp_data=new char[buffer_size];
	string temp_data_str;
	double temp_value;

	FILE* fp=fopen(filename.c_str(), "r");
	if(fp)
		while(feof(fp)==0)
		{
			fgets(buffer, buffer_size-1, fp);
			if(sscanf(buffer, "%s\t%lf", temp_data, &temp_value)==2)
			{
				temp_data_str=temp_data;

				T_DATASET->seq.push_back(temp_data_str);
				T_DATASET->value.push_back(temp_value);
				//clear buffer
				temp_data_str.clear();
				sprintf(buffer,"\0");
			}
		}
	else
		return -1;
/*
	for(int i=0;i<(T_DATASET->value.size());i++)
		printf("%s\t%s\t%lf\n", T_DATASET->seq_a[i].c_str(),  T_DATASET->seq_b[i].c_str(), T_DATASET->value[i]);
*/
	delete buffer;
	delete temp_data;

	T_DATASET->size=T_DATASET->value.size();

	return T_DATASET->value.size();
}


int PreComputeDKernel(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target)
{
	int i, j;
	AllocateKernelOuptputArray(KC, ref, target);

	for(i=0; i< target->size; i++)
		for(j=0; j< ref->size; j++)
			target->precomputed[i][j]=compute_distance_string_kernel(&(target->seq[i]), &(ref->seq[j]))/(double)(ref->seq[j].length());
	target->count_compute_kernel=ref->size;
	return target->size*ref->size;
}


int PreComputeWDKernel(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target)
{
	int i, j;
	ComputeWDKernelNormalizeCoefficient(KC, target);
	ComputeWDKernelNormalizeCoefficient(KC, ref);
	AllocateKernelOuptputArray(KC, ref, target);

	for(i=0; i< target->size; i++)
		for(j=0; j< ref->size; j++)
			if((target->self[i]!=0)&&(ref->self[j]!=0))
				target->precomputed[i][j]=compute_weighted_degree_string_kernel(&(KC->kmer), KC->weighting_coef, &(target->seq[i]), &(ref->seq[j]))/sqrt(target->self[i]*ref->self[j]);
			else
				target->precomputed[i][j]=0;
/*//DEBUG//
	for(i=0; i< target->size; i++)
		for(j=0; j< ref->size; j++)
			printf("(%d, %d):(%lf, %lf)\n", i, j, target->precomputed_a[i][j], target->precomputed_b[i][j]);
*/
	target->count_compute_kernel=ref->size;
	return target->size*ref->size;
}


int AllocateKernelOuptputArray(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target)
{
	int i;
	target->precomputed=new double* [target->size];
	for(i=0; i< target->size; i++)
		target->precomputed[i]=new double [ref->size];
	return 0;
}


int WeightingKernelOuptputAB(struct KernelConfig* KC, struct DataSet* target)
{
	int i,j;
	for(i=0; i< target->size; i++)
	{
		printf("%lf 0:%d", target->value[i], i+1);
		for(j=0; j< target->count_compute_kernel; j++)
			printf(" %d:%lf", j+1, target->precomputed[i][j]);
		printf("\n");
	}
	return 0;
}


int ComputeDKernel(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target)
{
	int i, j;
	double score, value;
//	double score_upper_bound=KC->weighting_a*(double)(ref->seq_a[0].length())
//		+KC->weighting_b*(double)(ref->seq_b[0].length());

	for(i=0; i< target->size; i++)
	{
		printf("%lf 0:%d", target->value[i], i+1);
		for(j=0; j< ref->size; j++)
		{
			value=compute_distance_string_kernel(&(target->seq[i]), &(ref->seq[j]))/(double)(ref->seq[j].length());
			score=value;
			if(KC->option_rbf==1)
				score=compute_radial_basis_function(KC->gamma, score);
			printf(" %d:%lf", j+1, score);
		}
		printf("\n");
	}

	return target->size*ref->size;
}


int ComputeWDKernel(struct KernelConfig* KC, struct DataSet* ref, struct DataSet* target)
{
	int i, j;
	double score, value;
	ComputeWDKernelNormalizeCoefficient(KC, target);
	ComputeWDKernelNormalizeCoefficient(KC, ref);

	for(i=0; i< target->size; i++)
	{
		printf("%lf 0:%d", target->value[i], i+1);
		for(j=0; j< ref->size; j++)
		{
			if((target->self[i]!=0)&&(ref->self[j]!=0))
				value=compute_weighted_degree_string_kernel(&(KC->kmer), KC->weighting_coef, &(target->seq[i]), &(ref->seq[j]))
							/sqrt(target->self[i]*ref->self[j]);
			else
				value=0;
			score=value;
			if(KC->option_rbf==1)
				score=compute_radial_basis_function(KC->gamma, score);
			printf(" %d:%lf", j+1, score);
		}
		printf("\n");
	}
	return target->size*ref->size;
}

int ComputeWDKernelNormalizeCoefficient(struct KernelConfig* KC, struct DataSet* target)
{
	int i;
	if(target->count_compute_self==0)
		for(i=0; i< target->size; i++)
			target->self.push_back(compute_weighted_degree_string_kernel(&(KC->kmer), KC->weighting_coef, &(target->seq[i]), &(target->seq[i])));

//	for(i=0; i< target->size; i++)
//		printf("%d:%lf\n", i, target->self_a[i]);

	(target->count_compute_self)++;
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////

double compute_distance_string_kernel(string* a, string* b)
{
	int i, length, score=0;
	if(a->length() <= b->length())
		length=(int)(a->length());
	else
		length=(int)(b->length());
	//printf("\n%d\t%s\t%s\n", length, a->c_str(), b->c_str());

	for(i=0;i<length;i++)
		if(a->c_str()[i]==b->c_str()[i])
			score++;
	return score;
}

double compute_weighted_degree_string_kernel(int* k, double* weighting_coef, string* a, string* b)
{
	int i, index, length, count;
	double score;

	if(a->length() <= b->length())
		length=(int)(a->length());
	else
		length=(int)(b->length());

	for(i=1, score=0; i <= *k; i++)
	{
		if(length < i)
			break;
		else
			for(index=0, count=0; (index+i)<=length; index++)
				if(a->substr(index, i) == b->substr(index, i))
					count++;
		score+=(double)(count)*(weighting_coef[i]);
		//printf("(%d: %d, %lf)\n", i, count, score);
	}
	return score;
}

double compute_radial_basis_function(double gamma, double value)
{
	return exp((-gamma)*value*value);
}

