
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* --- new data type --- */
typedef float *Vector;		/* Vector[1,...,size] */
typedef float **Matrix;		/* Matrix[1,...,nrow][1,...,ncol] */
typedef float ***Tensor;	/* Tensor[1,...,dim_x][1,...,dim_y][1,...,dim_z] */
typedef int *iArray; 
typedef int **iMatrix;
#define MaxProb 1.0e60
#define MinProb -1.0e60
#define SQUARE(A)   ((A)*(A))

static const double LOGZERO  = -1.0e60 ;            /* ~log(0) */
static const double LSMALL	 = -0.5e55 ;            /* log values < LSMALL are set to LOGZERO */
static const double LOG2PI   = 1.837877066 ;        /* log(2 * PI) */
static const double SQRT2PI  = 2.506628275 ;        /* sqrt (2 * PI) */
static const double PI       = 3.14159265358979 ;   /* PI */
static const double TPI      = 6.28318530717959 ;   /* PI*2 */
#define D 84240
#define MATH_INFINITE	1e30
#define MATH_PI			3.14159265358979			// PI
#define MATH_TPI		6.28318530717959			// 2*PI
#define MATH_LOG2PI		1.837877066					// log(2*PI)
#define MATH_SQRT2PI	2.506628275 				// sqrt(2*PI)
/* --- global variables --- */ 
/* Feature-related */
int FeaDim;
int MaxObsNum;
int MaxFrameNum;
int MaxPhoneNum;		/* Maximum phone number within a single utterance */

/* Model-related */
int MaxStateNum;		/* including non-emitting state */
int MaxMixNum;
int ModelNum;			/* number of total HMMs */
int VecSize;
char **phone_table;		/* model names */
int TotalGaussianNo;
int SP_NO;
/* IO path */
char AdaptFileList[200];
char CofPath[200];
char RecPath[200];
char FeaExtension[20];
char LabExtension[20];
char FeaType[20];
char OldHmmPath[200];
char NewHmmPath[200];
char ResultPath[200];
char TreeFileName[200];
char ModelList[200];
float TH;
float TAU;

struct Model
{
	char phone_name[20];
	int phone_id;

	int state_no;		/* emitting state */
	iArray mix_no;		/* mix_no[j] */

	Matrix TransP;		/* TransP[i][j]: state_no+2 x state_no+2 */ /* NOTE!!! */
	Matrix MixWeight;       /* MixWeight[j][m]*/
	Matrix GConst;		/* GConst[j][m] */
	Tensor Mean;		/* Mean[j][m][k] */
	Tensor Var;		/* Var[j][m][k] */
	
	Matrix GammaAcc;	/* GammaAcc[j][m] */
	Tensor GammaOtAcc;	/* GammaOtAcc[j][m][k] */
	Tensor GammaOt2Acc;	/* GammaOt2Acc[j][m][k] */
};
struct Model *hmm;

struct Observation
{
	int obs_id;
	int CLEARED;

	/* --- for transcription --- */
	int phone_no;		/* number of phoneme within the transcription */
	iArray phone;
	iArray state_no_this_phone;
	
	/* --- for model parameters --- */
	int state_no;

	iArray phone_belong;	/* state profile: which hmm_ID this state belongs to */
	iArray state_belong;	/* state profile: which state within the hmm this state belongs to */
	iArray mix_no;		/* state profile: number of mixture within this state */

	/* Tensor Mean;	*/	/* Mean[j][m][k]. Actually MatrixOfVector */
	/* Tensor Var; */	/* Var[j][m][k]. Actually MatrixOfVector */
	/* Matrix MixWeight; */	/* MixWeight[j][m]. Actually VectorOfVector */
	/* Matrix GConst; */	/* GConst[j][m]. Actually VectorOfVector */
	Matrix TransP;	/* TransP[i][j]. Actually VectorOfVector */

	/* --- for feature --- */
	int T;
	Matrix feature;		/* feature[t][k] */
	iArray s;		/* for signaling the result of Viterbi decoding */	
	/* Tensor gamma; */	/* gamma[t][j][m] */
	
};
struct Observation sentence;
/* NOTE: activate one sentence each time */

struct HMM_struct 
{			/* for reading-in HTK model    */
/*
   short phone_ID;
   short state_no;
   short mix_no[StateNum];
   float mixweight[StateNum][MaxMixNum];
   float mean[StateNum][MaxMixNum][VectorSize];
   float var[StateNum][MaxMixNum][VectorSize];
   float transP[StateNum][StateNum];
   float Gconst[StateNum][MaxMixNum];
*/
/* NOTE: including non-emitting states */

	short phone_ID;
	short state_no;         /* including non-emitting states */
	short *mix_no;		/* mix_no[MaxStateNum] */
	iMatrix mix_flag;	/* mix_flag[MaxStateNum][MaxMixNum] */
	Matrix mixweight;	/* mixweight[MaxStateNum][MaxMixNum] */
	Tensor mean;		/* mean[MaxStateNum][MaxMixNum][VecSize] */
	Tensor var;		/* var[MaxStateNum][MaxMixNum][VecSize] */
	Matrix transP;		/* transP[MaxStateNum][MaxStateNum] */
	Matrix Gconst;		/* Gconst[MaxStateNum][MaxMixNum] */
};

struct featype {      /* for HTK Feature Header */
   long number;
   long shift;
   short size;
   short kind;
};

void ClearMAPAcc(void);
void AccMAPAcc(void);
void AllocateMemSentence(void);
void FreeMemSentence(void);
int GetAdaptData(char *file_name, int obs_id);
void GetSenModel(int obs_id);

int LoadModel(void);
void WriteModel(void);
void AllocateMemSentence(void);
void FreeMemSentence(void);
int GetAdaptData(char *file_name, int obs_id);
void GetSenModel(int obs_id);
int ReadHMM(char *model_name, struct HMM_struct *model, int RCDmodel_index, int VecSize);
int LoadModel(void);
void WriteModel(void);
float TakeLog(float value);

float CalStateLikeli(int LogKey, int M, Vector data, Matrix Mean, Matrix Var, Vector C, Vector GConst);
float CalMixLikeli(int LogKey, Vector data, Vector Mean, Vector Var, float GConst);

Vector CreateVector(int dim);
void FreeVector(Vector v);
Matrix CreateMatrix(int n_row, int n_col);
void FreeMatrix(Matrix m, int n_row);
Tensor CreateTensor(int dim_x, int dim_y, int dim_z);
void FreeTensor(Tensor t, int dim_x, int dim_y);
void ZeroVector(Vector v, int dim);
void ZeroMatrix(Matrix m, int dim_x, int dim_y);
void ZeroTensor(Tensor t, int dim_x, int dim_y, int dim_z);
iArray CreateiArray(int dim);
void FreeiArray(iArray b);
iMatrix CreateiMatrix(int nrow, int ncol);
void FreeiMatrix(iMatrix c, int nrow);
void ZeroiArray(iArray b, int dim);
void ZeroiMatrix(iMatrix c, int n_row, int n_col);
Matrix CreateVectorOfVector(int dim);
Tensor CreateMatrixOfVector(int dim_x, int dim_y);
void FreeVectorOfVector(Matrix m);
void FreeMatrixOfVector(Tensor t, int dim_x);
void Abort(char *err_msg1, char *err_msg2);
void MAP_Initialize(char *filename);
void InitTempModel(struct HMM_struct *model);
void InitHmm(int RCDmodel_index, int state_no, short *mix_no);
float ViterbiDecode(void);
