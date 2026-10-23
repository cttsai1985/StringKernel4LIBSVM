#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>


#include "StringKernels.h"
#include "OrthogonalArray.h"
#include "GetLongOpt.h"
#include <cfloat>
#include <string>
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//EXTENSION SECTION/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//extern svm function
#include "svm.h"
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
extern void print_null(const char *s);
extern void read_problem(const char *filename);
extern double do_cross_validation(int type, int *fold_start, int *perm);
extern void do_cross_validation(int *fold_start, int *perm);

extern struct svm_parameter param;
extern struct svm_problem prob;		// set by read_problem
extern struct svm_node *x_space;
extern int nr_fold;
//for gold generation//start
int *fold_start;
int *perm;
void generate_fold(int l, int fold)
{
	int i;
	fold_start = Malloc(int,fold+1);
	//int l = prob.l;
	perm = Malloc(int,l);
	//int nr_class;

	for(i=0;i<l;i++) perm[i]=i;
	for(i=0;i<l;i++)
	{
		int j = i+rand()%(l-i);
		swap(perm[i],perm[j]);
	}
	for(i=0;i<=nr_fold;i++)
		fold_start[i]=i*l/fold;
	if(nr_fold==1)
		fold_start[0]=fold_start[1]=0;
}

void destroy_fold()
{
	free(fold_start);
	free(perm);
}
//for fold generation//end

//for analysis method
int type_analysis_method;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//setting
int oa_seed=-1;
int enlarge=0, shift_a=0, shift_b=0;
int fitness_type=0;
FILE* fp;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//for string kernel//start
DataSet* ptr_dataset;
int string_kernel_type, kmer_a, kmer_b;
double weighting_a, weighting_b;
void PositionalOFAT(int seq_type, DataSet backup_dataset, KernelConfig* ptr_kernel_config);
void PairedOFAT(DataSet backup_dataset, KernelConfig* ptr_kernel_config);
void PositionalMED(int seq_type, DataSet backup_dataset, KernelConfig* ptr_kernel_config);
void PairedMED(DataSet backup_dataset, KernelConfig* ptr_kernel_config);

void initialize_problem(struct DataSet* target)
{
	int i, j, k;
	prob.l=target->size;
	int elements=(target->count_compute_kernel+2)*(target->size);
	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node,elements);

	for(i=0, j=0; i<prob.l; i++)
	{
		prob.x[i]=&x_space[j];
		prob.y[i]=target->value[i];
		x_space[j].index=0;
		x_space[j++].value=(double)(i+1);
		for(k=0; k<target->count_compute_kernel; k++, j++)
		{
			x_space[j].index=k+1;
			x_space[j].value=0;
		}
		x_space[j++].index = -1;
	}
/*//DEBUG
	j=1;
	printf("%d: ", j);
	for(i=0; i<elements; i++)
		if(x_space[i].index == -1)
			printf(" %d\n%d: ", x_space[i].index, j++);
		else
			printf(" %d", x_space[i].index);
*///
}

void modify_problem(double weighting_a, double weighting_b, struct DataSet* target)
{
	int i, j;
	for(i=0, j=0; i<prob.l; i++)
	{
		for(j=0; j<=target->count_compute_kernel; j++)
			prob.x[i][j+1].value=weighting_a*target->precomputed_a[i][j]+weighting_b*target->precomputed_b[i][j];
	}
/*//DEBUG
	for(i=0; i<prob.l; i++)
	{
		for(j=0; j<(target->count_compute_kernel+2); j++)
			printf("%d:%lf ", prob.x[i][j].index, prob.x[i][j].value);
		printf("\n");
	}*/
}
//for string kernel//end
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int preload_settings(string filename)
{
	double temp;
	char c;
	FILE *fp_i;

	//setting
	string_kernel_type=2;
	weighting_a=0.5;
	weighting_b=0.5;
	kmer_a=3;
	kmer_b=3;

	param.svm_type=4;
	param.nu = 0.5;
	param.C = 1;
	param.p = 0.1;

	//loading file
	if(filename.empty()==true)
		return 0;
	else
		fp_i=fopen(filename.c_str(), "r");
	if(fp_i==NULL)
		return 0;
	while(!feof(fp_i))
	{
		c=fgetc(fp_i);
		//printf("%c\n", c);
		fscanf(fp_i, "Cost=%lf (%lf)", &param.C, &temp);
		if(fscanf(fp_i, "epsilon=%lf (%lf)", &param.p, &temp)==2)
			param.svm_type=3;
		if(fscanf(fp_i, "nu=%lf", &param.nu)==1)
			param.svm_type=4;
		fscanf(fp_i, "wa=%lf, wb=%lf", &weighting_a, &weighting_b);
		fscanf(fp_i, "ka=%d, kb=%d", &kmer_a, &kmer_b);
	}
	printf("loading setting:\n");
	if(param.svm_type==3)
		printf("Cost=%lf, epsilon=%lf, wa=%lf, wb=%lf, ka=%d, kb=%d\n", param.C, param.p, weighting_a, weighting_b, kmer_a, kmer_b);
	if(param.svm_type==4)
		printf("Cost=%lf, nu=%lf, wa=%lf, wb=%lf, ka=%d, kb=%d\n", param.C, param.nu, weighting_a, weighting_b, kmer_a, kmer_b);
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//MISC////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_seperate_bar(int type, int limit)
{
	int i;
	switch(type)
	{
		case 0:
			for(i=0; i<limit; i++)
				printf("%d \t", i+1);
			break;
		case 1:
			for(i=0; i<limit; i++)
				printf("-------- \t");
			break;
	}
	printf("\n");
}


void print_compare(double a, double b)
{
	if(a>b)
		printf(" +\t");
	else if(a<b)
		printf(" -\t");
	else
		printf(" =\t");
}

void print_svm_dataset()
{
	int i, j;
	for(i=0, printf("no %d, ", i+1); i< prob.l; i++, printf("\nno %d, ", i+1))
		for(j=0; prob.x[i][j].index!=-1;j++)
			printf("%d:%lf\t", prob.x[i][j].index,prob.x[i][j].value);
}

double fitness_adjust(double score)
{
	switch(fitness_type)
	    {
	        case 0://Mean Squared Error
                return -score;
            case 1://Mean Absolute Percentage Error
                return -score;
            case 2://Mean Absolute Error
                return -score;
			case 3://Absolute Total Percentage Error
				return -score;
			case 4://Weighted Mean Squared Error
				return -score;
            case 5://Weighted Mean Absolute Rrror
				return -score;
			case 6://Root Mean Squared Error
				return -score;
			case 7://prediction
				return score;
			case 8://prediction std
				return score;
            default:
                return 0;
	    }
}
//END MISC SECTION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PSO MAIN SECTION///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int neval;
void exit_with_help()
{
	printf(
	"Usage: svmsk_analyser [options] training_set_file\n"
	"options:\n"
	"-F fitness_type: select fitness type (default 0)\n"
	"	0 -- Mean squared error\n"
	"	1 -- Mean absolute percentage error\n"
	"	2 -- Mean absolute error\n"
	"	3 -- Absolute percentage total error\n"
	"	4 -- Weighted squared error\n"
	"	5 -- Weighted absolute error\n"
	"	6 -- Root mean squared error\n"
	"	7 -- Prediction\n"
	"-S setting filename: preload settings from desinated file\n"
	"-A Factor Analysis Methods: select one method to perform (default 0)\n"
	"	0 -- One-Factor-At-a-Time method (OFAT)\n"
	"	1 -- Main Effect Difference (MED)\n"
	"	2 -- Main Effect Difference with Interactions (MED)\n"
	"-L enlarge: enlarge OA table to prevent potential OA problems (default 0)\n"
	"-W*: Shift OA starting factor to prevent potential OA problems (default 0)\n"
	"  *=a: from string a (default=0)\n"
	"  *=b: from string b (default=0)\n"
	"-R OA seed: Random seed of initializing OA table (default -1=inactive, >0=active)\n"
    "-----=====  LIBSVM  =====-----\n"
	"-s svm_type: set type of SVM (default 4)\n"
	"	3 -- epsilon-SVR\n"
	"	4 -- nu-SVR\n"
	"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
	"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
	"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
	"-t kernel_type: set type of string kernel function (default 2)\n"
	"	0 -- Distance String Kernel Function(Not Offical)\n"
	"	1 -- Distance String Kernel Function with Weight-Adjusting (Not Offical)\n"
	"	2 -- Weighted Degree String Kernel Function\n"
	"	3 -- Weighted Degree String Kernel Function with Weight-Adjusting\n"
	"-k*: set the max length of k-mer\n"
	"  *=a : k-mer of string a (default=3)\n"
	"  *=b : k-mer of string b (default=3)\n"
	"-w*: set the weighting of designated string for kernel value\n"
	"  *=a : weighting of string a (default=0.5)\n"
	"  *=b : weighting of string b (default=0.5)\n"
	"-m cachesize: set cache memory size in MB (default 100)\n"
	"-e epsilon: set tolerance of termination criterion (default=0.001)\n"
	"-h shrinking: whether to use the shrinking heuristics, 0 or 1 (default=1)\n"
	"-b probability_estimates: whether to train a SVR model for probability estimates, 0 or 1 (default=0)\n"
    "-v n: n-fold cross validation mode (default=10)\n"
    "   1 -- self train and test mode\n"
	"  -1 -- automatic leave-one-out cross validation\n"
	);
	//system("pause");
	exit(1);
}


// Rewrite the fitness function as your wish
float fitness(/*float *x, int dim*/)
{
    //int i;
    float s = 0;
    modify_problem(weighting_a, weighting_b, ptr_dataset);
    s=do_cross_validation(fitness_type, fold_start, perm);
    neval++;
    return s;
}

int main(int argc, char* argv[])
{
    //--CTTSAI--
    //read parameter
    if(argc==1)
        exit_with_help();
    GetLongOpt GLO(argc, argv);

	//preloading settings
	string settings_filename=GLO.GetStringOpt("S","");
	preload_settings(settings_filename);

    //svm parameters
    param.svm_type = GLO.GetIntOpt("s", param.svm_type);
    param.kernel_type = PRECOMPUTED;
    param.degree = 2;
	param.gamma = 0;	// 1/k //can be tuned
	param.coef0 = 0;
	param.nu = GLO.GetFloatOpt("n", param.nu);
	param.cache_size = GLO.GetIntOpt("m", 100);
	param.C = GLO.GetFloatOpt("c", param.C);
	param.eps = GLO.GetFloatOpt("e", 1e-3);
	param.p = GLO.GetFloatOpt("p", param.p);
	param.shrinking = GLO.GetIntOpt("h", 1);
	param.probability = GLO.GetIntOpt("b", 0);
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;

	//set analysis method
	type_analysis_method=GLO.GetIntOpt("A", 0);
	enlarge=GLO.GetIntOpt("L", 0);
	shift_a=GLO.GetIntOpt("Wa", 0);
	shift_b=GLO.GetIntOpt("Wb", 0);
	oa_seed=GLO.GetIntOpt("R", -1);

	//set string kernel
	string_kernel_type=GLO.GetIntOpt("t", string_kernel_type);
	weighting_a=GLO.GetFloatOpt("wa", weighting_a);
	weighting_b=GLO.GetFloatOpt("wb", weighting_b);
	kmer_a=GLO.GetIntOpt("ka", kmer_a);
	kmer_b=GLO.GetIntOpt("kb", kmer_b);
	KernelConfig kernel_config(string_kernel_type, kmer_a, kmer_b, weighting_a, weighting_b);
	kernel_config.SetCoef();

    //filename
    string dataset_filename, log_filename;
    dataset_filename=argv[argc-1];
    log_filename=GetMainFilename(dataset_filename)+"_log.txt";
    //fp=fopen(log_filename.c_str(), "w");

    fitness_type=GLO.GetIntOpt("F", 0); //set fitness
	switch(string_kernel_type)
    {
		case 0:
			printf("using Distance String Kernel Function (Not Offical)\n");
			break;
		case 1:
			printf("using Distance String Kernel Function with Weight-Adjusting (Not Offical)\n");
			break;
		case 2:
			printf("using Weighted Degree String Kernel Function\n");
			break;
		case 3:
			printf("using Weighted Degree String Kernel Function with Weight-Adjusting\n");
			break;
    }

    time_t t_start, t_stop;
    time(&t_start);
    printf("start at %s", ctime(&t_start));
    int seed = time(NULL);    // using time clock as the seed of the random number generator
    srand(seed); //move to front

	//read dataset
    DataSet train_dataset;
    DataSet backup_dataset, using_dataset;
	train_dataset.initialize();
	int error_var=LoadDataset(dataset_filename, &train_dataset);
	if(error_var==-1)
		return 0;

	ptr_dataset=&train_dataset;
	PreComputeDKernel(&kernel_config, ptr_dataset, ptr_dataset);
	initialize_problem(ptr_dataset);

	//cv
    nr_fold=GLO.GetIntOpt("v", 10);
    if(nr_fold==-1)
        nr_fold=prob.l;
    generate_fold(prob.l,nr_fold);//generate fold for svm
    //check input error
    const char *error_msg;
	error_msg = svm_check_parameter(&prob,&param);
	if(error_msg)
	{
		fprintf(stderr,"Error: %s\n",error_msg);
		exit(1);
	}
	svm_set_print_string_function(&print_null);//silence svm output

	printf("using setting:\n");
	if(param.svm_type==3)
		printf("Cost=%lf, epsilon=%lf, wa=%lf, wb=%lf, ka=%d, kb=%d\n", param.C, param.p, weighting_a, weighting_b, kmer_a, kmer_b);
	if(param.svm_type==4)
		printf("Cost=%lf, nu=%lf, wa=%lf, wb=%lf, ka=%d, kb=%d\n", param.C, param.nu, weighting_a, weighting_b, kmer_a, kmer_b);

	switch(type_analysis_method)
	{
		case 0:
			printf("Performing OFAT method: \n");
			PositionalOFAT(1, train_dataset, &kernel_config);
			PositionalOFAT(2, train_dataset, &kernel_config);
			PairedOFAT(train_dataset, &kernel_config);
			break;
		case 1:
			printf("Performing MED method: \n");
			PositionalMED(1, train_dataset, &kernel_config);
			PositionalMED(2, train_dataset, &kernel_config);
			break;
		case 2:
			printf("Performing MED method for interaction: \n");
			PairedMED(train_dataset, &kernel_config);
			break;
	}

    time(&t_stop);
    int hour, min, sec;
    sec = t_stop - t_start;
    min = sec / 60;
    sec %= 60;
    hour = min / 60;
    min %= 60;

    printf("\n");
    printf("Analysis finished at %s", ctime(&t_stop));
    printf("execution time: %3dh %2dm %2ds.\n", hour, min, sec);

    // print final results
    printf("\n");
    printf("Total number of function evaluations: %d\n", neval);
    /*
    for(i=0;i<dim;i++)
        printf("%f ", ptr->p[i]);
    printf("\n");*/
    //fclose(fp); //end of log

	//free
	svm_destroy_param(&param);
	free(prob.y);
	free(prob.x);
	free(x_space);
	destroy_fold();

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void PositionalOFAT(int seq_type, DataSet backup_dataset, KernelConfig* ptr_kernel_config)
{
	//seq_a:seq_type=1, seq_a:seq_type=2,
	int i, str_length=0;
	double base_line, temp, fitness_max, fitness_min, fitness_avg;
	DataSet using_dataset;

	if(seq_type==1)
		str_length=(int)(backup_dataset.seq_a[0].length());
	else if(seq_type==2)
		str_length=(int)(backup_dataset.seq_b[0].length());

	ptr_dataset=&backup_dataset;
	if((string_kernel_type==0)||(string_kernel_type==1))
		PreComputeDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
	if((string_kernel_type==2)||(string_kernel_type==3))
		PreComputeWDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
	base_line=fitness();

	if(seq_type==1)
		printf("string a \tfitness \tdiff    \tr_diff(%%): \n");
	else if(seq_type==2)
		printf("string b \tfitness \tdiff    \tr_diff(%%): \n");
	printf("-----------\t---------\t---------\t----------\n");

	for(i=0, fitness_max=-DBL_MAX, fitness_min=DBL_MAX, fitness_avg=0; i<str_length; i++)
	{
		using_dataset=backup_dataset;
		ptr_dataset=&using_dataset;
		using_dataset.erase_one_char_in_string(seq_type, i);
		if((string_kernel_type==0)||(string_kernel_type==1))
			PreComputeDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
		if((string_kernel_type==2)||(string_kernel_type==3))
			PreComputeWDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
		temp=fitness();
		ptr_dataset->clear();
		if(temp>fitness_max)
			fitness_max=temp;
		else if (temp<fitness_min)
			fitness_min=temp;
		fitness_avg+=temp;
		printf("position %d:\t%lf\t%lf\t%lf\n",i+1, temp, temp-base_line, 100.0*(temp-base_line)/base_line);
	}
	printf("-----------\t---------\t---------\t----------\n");
	printf("base=%lf, max=%lf, min=%lf, avg=%lf\n", base_line, fitness_max, fitness_min, fitness_avg/(double)(str_length));
	printf("dif_max=%lf, dif_min=%lf\n\n", fitness_max-base_line, fitness_min-base_line);
}


void PairedOFAT(DataSet backup_dataset, KernelConfig* ptr_kernel_config)
{
	int i, j;
	double temp, base_line, fitness_max, fitness_min, fitness_avg;
	DataSet using_dataset;

	int str_length_a=(int)(backup_dataset.seq_a[0].length());
	int str_length_b=(int)(backup_dataset.seq_b[0].length());

	ptr_dataset=&backup_dataset;
	if((string_kernel_type==0)||(string_kernel_type==1))
		PreComputeDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
	if((string_kernel_type==2)||(string_kernel_type==3))
		PreComputeWDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
	base_line=fitness();

	printf("pairs \t");
	print_seperate_bar(0, str_length_b);
	print_seperate_bar(1, str_length_b);

	for(i=0, fitness_max=-DBL_MAX, fitness_min=DBL_MAX, fitness_avg=0; i<str_length_a; i++)
	{
		printf("%d: \t", i+1);
		for(j=0; j<str_length_b; j++)
		{
			using_dataset=backup_dataset;
			ptr_dataset=&using_dataset;
			using_dataset.erase_one_char_in_string(1,i);
			using_dataset.erase_one_char_in_string(2,j);
			if((string_kernel_type==0)||(string_kernel_type==1))
				PreComputeDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
			if((string_kernel_type==2)||(string_kernel_type==3))
				PreComputeWDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
			temp=fitness();
			ptr_dataset->clear();
			printf("%lf\t", temp);
			if(temp>fitness_max)
				fitness_max=temp;
			else if (temp<fitness_min)
				fitness_min=temp;
			fitness_avg+=temp;
		}
		printf("\n");
	}
	print_seperate_bar(1, str_length_b);
	printf("base=%lf, max=%lf, min=%lf, avg=%lf\n", base_line, fitness_max, fitness_min, fitness_avg/(double)(str_length_a*str_length_b));
}


void PositionalMED(int seq_type, DataSet backup_dataset, KernelConfig* ptr_kernel_config)
{
	//seq_a:seq_type=1, seq_a:seq_type=2,
	int i, j, str_length=0, shift=0;
	double score, fitness_max, fitness_min, fitness_avg;
	DataSet using_dataset;

	switch(seq_type)
	{
		case 1:
			str_length=(int)(backup_dataset.seq_a[0].length());
			shift=shift_a;
			break;
		case 2:
			str_length=(int)(backup_dataset.seq_b[0].length());
			shift=shift_b;
			break;
	}

	double* position_score_add=new double[str_length];
	double* position_score_remove=new double[str_length];
	OrthogonalArray OA(2, str_length+enlarge);
	if(OA.factors < str_length+shift)
		shift=OA.factors-str_length;
	printf("using OA: OA factors=%d, length=%d, enlarge=%d, shift=%d\n\n", OA.factors, str_length, enlarge, shift);

	if(oa_seed>=0)
		OA.ShuffleColumn(oa_seed);

	for(j=0;j<str_length;j++)
		position_score_add[j]=position_score_remove[j]=0;

	if(seq_type==1)
		printf("string a \t(+)score \t(-)score \tdiff \tabs_diff(%%) \t+/-: \n");
	else if(seq_type==2)
		printf("string b \t(+)score \t(-)score \tdiff \tabs_diff(%%) \t+/-: \n");
	printf("-------- \t-------- \t-------- \t------ \t--------- \t---- \n");

	for(i=0, fitness_max=-DBL_MAX, fitness_min=DBL_MAX, fitness_avg=0; i<OA.rows; i++)
	{
		using_dataset=backup_dataset;
		ptr_dataset=&using_dataset;
		using_dataset.erase_chars_in_string(seq_type, &OA[i][shift]);
		if((string_kernel_type==0)||(string_kernel_type==1))
			PreComputeDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
		if((string_kernel_type==2)||(string_kernel_type==3))
			PreComputeWDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
		score=fitness();
		fitness_avg+=score;
		if(score>fitness_max)
			fitness_max=score;
		else if(score<fitness_min)
			fitness_min=score;

		//into MED section
		for(j=0;j<str_length;j++)
			if(OA[i][shift+j]==0)
				position_score_add[j]+=score;
			else if(OA[i][shift+j]==1)
				position_score_remove[j]+=score;
	}

	for(j=0;j<str_length;j++)
	{
		printf("position %d: \t%lf\t%lf\t%.4lf\t%.6lf\t", j+1,
				position_score_add[j]*2.0/(double)(OA.rows),
				position_score_remove[j]*2.0/(double)(OA.rows),
				fitness_adjust((position_score_add[j]-position_score_remove[j]))*2.0/(double)(OA.rows),
				100.0*fabs((position_score_add[j]-position_score_remove[j])/position_score_add[j]));
		print_compare(fitness_adjust(position_score_add[j]), fitness_adjust(position_score_remove[j]));
		printf("\n");
	}
	printf("-------- \t-------- \t-------- \t------ \t--------- \t---- \n");
	printf("average=%lf, max=%lf, min=%lf, differ=(Pscore_add-Pscore_remove)\n\n",
			fitness_avg/(double)(OA.rows), fitness_max, fitness_min);
}

void PairedMED(DataSet backup_dataset, KernelConfig* ptr_kernel_config)
{
	int i, j, k, a, b;
	double score, fitness_avg, fitness_max, fitness_min;
	DataSet using_dataset;

	int str_length_a=(int)(backup_dataset.seq_a[0].length());
	int str_length_b=(int)(backup_dataset.seq_b[0].length());
	int str_length=str_length_a+str_length_b;

	double* position_score_add=new double[str_length];
	double* position_score_remove=new double[str_length];
	for(j=0;j<str_length;j++)
		position_score_add[j]=position_score_remove[j]=0;

	OrthogonalArray OA(2, str_length+enlarge);
	if(OA.factors < str_length+shift_a+shift_b)
		shift_b=(OA.factors-str_length);
	printf("using OA: factors=%d, enlarge=%d, length_a=%d, shift_a=%d, length_b=%d, shift_b=%d\n\n",
				OA.factors, enlarge, str_length_a, shift_a, str_length_b, shift_b);

	if(oa_seed>=0)
		OA.ShuffleColumn(oa_seed);

	double** pairs_score_add=new double* [str_length_a];
	double** pairs_score_remove=new double* [str_length_a];

	for(a=0;a<str_length_a;a++)
	{
		pairs_score_add[a]=new double [str_length_b];
		pairs_score_remove[a]=new double [str_length_b];
		for(b=0; b<str_length_b; b++)
			pairs_score_add[a][b]=pairs_score_remove[a][b]=0;
	}

	//OA.print();
	printf("string a \t(+)score \t(-)score \tdiff \tabs_diff(%%) \t+/-: \n");
	printf("-------- \t-------- \t-------- \t------ \t--------- \t---- \n");
	for(i=0, fitness_avg=0, fitness_max=-DBL_MAX, fitness_min=DBL_MAX; i<OA.rows; i++)
	{
		using_dataset=backup_dataset;
		ptr_dataset=&using_dataset;
		using_dataset.erase_chars_in_string(1, &OA[i][shift_a]);
		using_dataset.erase_chars_in_string(2, &OA[i][str_length_a+shift_a+shift_b]);
		if((string_kernel_type==0)||(string_kernel_type==1))
			PreComputeDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
		if((string_kernel_type==2)||(string_kernel_type==3))
			PreComputeWDKernel(ptr_kernel_config, ptr_dataset, ptr_dataset);
		score=fitness();
		fitness_avg+=score;
		if(score>fitness_max)
			fitness_max=score;
		else if(score<fitness_min)
			fitness_min=score;

		//into MED section
		for(j=0;j<str_length;j++)
			if(OA[i][shift_a+j]==0)
				position_score_add[j]+=score;
			else if(OA[i][shift_a+j]==1)
				position_score_remove[j]+=score;

		for(a=0; a<str_length_a; a++)
			for(b=0; b<str_length_b; b++)
				if(OA[i][shift_a+a]==OA[i][shift_a+shift_b+str_length_a+b])
					pairs_score_add[a][b]+=score;
				else
					pairs_score_remove[a][b]+=score;
	}

	for(j=0, k=0;j<str_length;j++, k++)
	{
		if(j==str_length_a)
		{
			printf("-------- \t-------- \t-------- \t------ \t--------- \t---- \n\n"
					"string b \t(+)score \t(-)score \tdiff \tabs_diff(%%) \t+/-: \n"
					"-------- \t-------- \t-------- \t------ \t--------- \t---- \n");
			k=0;
		}
		printf("position %d: \t%lf\t%lf\t%.4lf\t%.6lf\t",k+1,
				position_score_add[j]*2.0/(double)(OA.rows),
				position_score_remove[j]*2.0/(double)(OA.rows),
				fitness_adjust((position_score_add[j]-position_score_remove[j]))*2.0/(double)(OA.rows),
				100.0*fabs((position_score_add[j]-position_score_remove[j])/position_score_add[j]));
		print_compare(position_score_add[j], position_score_remove[j]);
			printf("\n");
	}
	printf("-------- \t-------- \t-------- \t------ \t--------- \t---- \n");
	printf("average=%lf, max=%lf, min=%lf, differ=(Pscore_add-Pscore_remove)\n\n",
		fitness_avg/(double)(OA.rows), fitness_max, fitness_min);

	////
	printf("pairs \t");
	print_seperate_bar(0, str_length_b);
	print_seperate_bar(1, str_length_b);
	for(a=0; a<str_length_a; a++)
	{
		printf("%d: \t", a+1);
		for(b=0; b<str_length_b; b++)
			printf("%lf\t",
				fitness_adjust((pairs_score_add[a][b]-pairs_score_remove[a][b]))*2.0/(double)(OA.rows));
		printf("\n");
	}
	print_seperate_bar(1, str_length_b);

	////
	printf("pairs \t");
	print_seperate_bar(0, str_length_b);
	print_seperate_bar(1, str_length_b);
	for(a=0; a<str_length_a; a++)
	{
		printf("%d: \t", a+1);
		for(b=0; b<str_length_b; b++)
			print_compare(fitness_adjust(pairs_score_add[a][b]), fitness_adjust(pairs_score_remove[a][b]));
		printf("\n");
	}
	printf("\n");
	print_seperate_bar(1, str_length_b);
}

