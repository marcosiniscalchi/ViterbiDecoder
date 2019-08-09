/*-------------------------------------------------------*/
/*                   Initialize_MAP.c                    */
/*           Module: Memory Allocation Operations        */
/*                                                       */
/* (1) Read in global constants from SA.ini              */
/* (2) Memory allocation for basic data structure        */
/*                                                       */
/*                             Last Update: 2000/05/26   */
/*                             Kenneth, Kuan-Ting Chen   */
/* Bugs Report: please contact kenneth@iis.sinica.edu.tw */
/*-------------------------------------------------------*/
#include "BasicDef_MAP.h"
/*-------------------------- Create Vectors ------------------------------*/
Vector CreateVector(int dim)
{
   Vector v;
   v = (float *)calloc(dim, sizeof(float));
   return (v);
}

/*--------------------------- Free Vectors -------------------------------*/
void FreeVector(Vector v)
{
   free(v);
}

/*------------------------- Create Matrices ------------------------------*/
Matrix CreateMatrix(int n_row, int n_col)
{
   Matrix m;
   int i;
   m = (float **)calloc(n_row, sizeof(float *));
   for (i = 0; i < n_row; i++)
      m[i] = CreateVector(n_col);
   return (m);
}

/*-------------------------- Free Matrices -------------------------------*/
void FreeMatrix(Matrix m, int n_row)
{
   int i;
   for (i = 0; i < n_row; i++)
      FreeVector(m[i]);
   free(m);
}

/*-------------------------- Create Tensors ------------------------------*/
Tensor CreateTensor(int dim_x, int dim_y, int dim_z)
{
   Tensor t;
   int i;
   t = (float ***)calloc(dim_x, sizeof(float **));
   for (i = 0; i < dim_x; i++) {
      t[i] = CreateMatrix(dim_y, dim_z);
   }
   return (t);
}

/*--------------------------- Free Tensors -------------------------------*/
void FreeTensor(Tensor t, int dim_x, int dim_y)
{
   int i;
   for (i = 0; i < dim_x; i++) {
      FreeMatrix(t[i], dim_y);
   }
   free(t);
}

/* --------------------- Zero Vector/Matrix/Tensor -----------------------*/
void ZeroVector(Vector v, int dim) {
   int i;
   for (i=0; i<dim; i++)
      v[i] = 0.0;
}

void ZeroMatrix(Matrix m, int dim_x, int dim_y) {
   int i, j;
   for (i=0; i<dim_x; i++)
      for (j=0; j<dim_y; j++)
	 m[i][j] = 0.0;
}

void ZeroTensor(Tensor t, int dim_x, int dim_y, int dim_z) {
   int i, j, k;
   for (i=0; i<dim_x; i++)
      for (j=0; j<dim_y; j++)
	 for (k=0; k<dim_z; k++)
	    t[i][j][k] = 0.0;
}

/* ------------------ Create/Free iArray/iMatrix -------------------------*/
iArray CreateiArray(int dim)
{
	iArray b;
	b = (int *)calloc(dim, sizeof(int));
	return (b);
}

void FreeiArray(iArray b)
{
	free(b);
}

iMatrix CreateiMatrix(int nrow, int ncol)
{
	int j;
	iMatrix c;

	c = (int **)calloc(nrow, sizeof(int *));
	for (j=0; j<nrow; j++)
		c[j] = CreateiArray(ncol);
	return (c);
}

void FreeiMatrix(iMatrix c, int nrow)
{
	int j;
	
	for (j=0; j<nrow; j++)
		FreeiArray(c[j]);
	free(c);
}

void ZeroiArray(iArray b, int dim)
{
	int i;
	for (i=0; i<dim; i++)
		b[i] = 0;
}

void ZeroiMatrix(iMatrix c, int n_row, int n_col)
{
	int i, j;
	for (i=0; i<n_row; i++)
		for (j=0; j<n_col; j++)
			c[i][j] = 0;
}

/* --------------- kenneth_2000/05/26 --------------- */
/* For sentence.Mean ... etc. */
Matrix CreateVectorOfVector(int dim)
{
	Matrix m;
	
	m=(float **)calloc(dim, sizeof(float *));
	return(m);	
}

Tensor CreateMatrixOfVector(int dim_x, int dim_y)
{
	Tensor t;				
	int i;

	t=(float ***)calloc(dim_x, sizeof(float **));
	for (i=0; i<dim_x; i++)
		t[i]=CreateVectorOfVector(dim_y);
	
	return(t);
}
	
void FreeVectorOfVector(Matrix m)
{
	free(m);
}

void FreeMatrixOfVector(Tensor t, int dim_x)
{
	int i;
	for (i=0; i<dim_x; i++)
		FreeVectorOfVector(t[i]);
	free(t);
}
	
/* -------------------------  Abort Message ------------------------------*/
void Abort(char *err_msg1, char *err_msg2)
{
	printf("[%s] %s\n", err_msg1, err_msg2);
	exit(-1);
}

/* --------------------------- Read SA.ini -------------------------------*/
void MAP_Initialize(char *filename)
{
	FILE *f_ini, *f_mlist;
	char tmpstr[300];
	char *p1;
	int flag[30];
	int i, j, k;
	int temp;
	int count=0;
	int count1;
	int sentence_no;
    int MAX=300;
	for (i=0; i<17; i++)
		flag[i] = 0;
	
	if ((f_ini = fopen(filename, "r"))==NULL)
		Abort("ERROR Opening File", filename);
	fscanf(f_ini, "%s", tmpstr);
	while (!feof(f_ini))
	{
		if (strcmp(tmpstr, "[Feature]")==0)
			;
		else if (strcmp(tmpstr, "[Model]")==0)
			;
		else if (strcmp(tmpstr, "[MAP]")==0)
			;
		else if (strcmp(tmpstr, "[Tree]")==0)
			;
		else if (strcmp(tmpstr, "FeaDim")==0)
		{
			fscanf(f_ini, "%d", &FeaDim);
			flag[0] = 1;
		}
		else if (strcmp(tmpstr, "MaxObsNum")==0)
		{
			fscanf(f_ini, "%d", &MaxObsNum);
			flag[1] = 1;
		}
		else if (strcmp(tmpstr, "MaxFrameNum")==0)
		{
			fscanf(f_ini, "%d", &MaxFrameNum);
			flag[2] = 1;
		}
		else if (strcmp(tmpstr, "MaxPhoneNum")==0)
		{
			fscanf(f_ini, "%d", &MaxPhoneNum);
			flag[3] = 1;
		}
		else if (strcmp(tmpstr, "MaxStateNum")==0)
		{
			fscanf(f_ini, "%d", &MaxStateNum);
			flag[4] = 1;
		}		
		else if (strcmp(tmpstr, "MaxMixNum")==0)
		{
			fscanf(f_ini, "%d", &MaxMixNum);
			flag[5] = 1;
		}
		else if (strcmp(tmpstr, "ModelNum")==0)
		{
			fscanf(f_ini, "%d", &ModelNum);
			flag[6] = 1;
		}
		else if (strcmp(tmpstr, "VecSize")==0)
		{
			fscanf(f_ini, "%d", &VecSize);
			flag[7] = 1;
			if (VecSize!=FeaDim)
				Abort("ERROR initializing", "VecSize!=FeaDim");
		}
		
		else if (strcmp(tmpstr, "CofPath")==0)
		{
			fscanf(f_ini, "%s", CofPath);
			flag[8] = 1;
		}
		else if (strcmp(tmpstr, "RecPath")==0)
		{
			fscanf(f_ini, "%s", RecPath);
			flag[9] = 1;
		}
		else if (strcmp(tmpstr, "FeaExtension")==0)
		{
			fscanf(f_ini, "%s",FeaExtension);
			flag[10] = 1;
		}
		else if (strcmp(tmpstr, "LabExtension")==0)
		{
			fscanf(f_ini, "%s", LabExtension);
			flag[11] = 1;
		}
		else if (strcmp(tmpstr, "TAU")==0)
		{
			fscanf(f_ini, "%f", &TH);
			flag[12] = 1;
		}
		else if (strcmp(tmpstr, "Model_list")==0)
		{
			fscanf(f_ini, "%s", ModelList);
			flag[13] = 1;
		}
		else if (strcmp(tmpstr, "AdaptFileList")==0)
		{
			fscanf(f_ini, "%s", AdaptFileList);
			flag[14] = 1;
		}
		else if (strcmp(tmpstr, "SP_NO")==0)
		{
			fscanf(f_ini, "%d", &SP_NO);
			flag[15] = 1;
		}
		else if (strcmp(tmpstr, "[EndSA]")==0)
		{
			for (i=0; i<16; i++)
			{
				if (flag[i]!=1)
					Abort("ERROR", "Environmental constants incompletely specified");
			}
			break;
		}
		else
			Abort("ERROR undefined", tmpstr);

		fscanf(f_ini, "%s", tmpstr);
	}
	fclose(f_ini);
	phone_table = (char **)calloc(ModelNum, sizeof(char *));
	for (j=0; j<ModelNum; j++)
		phone_table[j] = (char *)calloc(20, sizeof(char));

	if ((f_mlist = fopen(ModelList, "r"))==NULL)
		Abort("ERROR Opening File", ModelList);	

	k=0;
	while(fscanf(f_mlist, "%s\n", tmpstr)!=EOF)
	{

			strcpy(phone_table[k], tmpstr);
			k++;
	}
	fclose(f_mlist);
	if (k!=ModelNum)
		Abort("ERROR", "ModelNum and item number in ModelList don't match");

	/* --- Allocate momory for model parameters --- */
	hmm = (struct Model *)calloc(ModelNum, sizeof(struct Model));
		
}

void InitTempModel(struct HMM_struct *model)
{
	model->mix_no = (short *)calloc(MaxStateNum, sizeof(short));
	model->mix_flag = CreateiMatrix(MaxStateNum, MaxMixNum);
	model->mixweight = CreateMatrix(MaxStateNum, MaxMixNum);
	model->mean = CreateTensor(MaxStateNum, MaxMixNum, VecSize);
	model->var = CreateTensor(MaxStateNum, MaxMixNum, VecSize);
	model->transP = CreateMatrix(MaxStateNum, MaxStateNum);
	model->Gconst = CreateMatrix(MaxStateNum, MaxMixNum);
}

void InitHmm(int RCDmodel_index, int state_no, short *mix_no)
{				/* including non-emitting states */
	int j, max_M;
	int mix_no_this_state;
	
	max_M = 0;
	for (j=1; j<=state_no-2; j++)	/* e.g. state_no=5, j=(0),1,2,3,(4) */
	{
		mix_no_this_state = mix_no[j];
		if (mix_no_this_state > max_M)
			max_M = mix_no_this_state;
	}
	if (max_M==0)
		Abort("ERROR", "max_M == 0 in InitHmm()\n");

	hmm[RCDmodel_index].phone_id = RCDmodel_index;
	hmm[RCDmodel_index].state_no = state_no - 2;					/* emitting state */
	hmm[RCDmodel_index].mix_no = CreateiArray(hmm[RCDmodel_index].state_no);	/* mix_no[j] */

	hmm[RCDmodel_index].TransP = CreateMatrix(hmm[RCDmodel_index].state_no+2, hmm[RCDmodel_index].state_no+2);	/* TransP[i][j]: state_no+2 x state_no+2 */
	hmm[RCDmodel_index].MixWeight = CreateMatrix(hmm[RCDmodel_index].state_no, max_M);  /* MixWeight[j][m]*/
	hmm[RCDmodel_index].GConst = CreateMatrix(hmm[RCDmodel_index].state_no, max_M);	    /* GConst[j][m] */
	hmm[RCDmodel_index].Mean = CreateTensor(hmm[RCDmodel_index].state_no, max_M, VecSize);		/* Mean[j][m][k] */
	hmm[RCDmodel_index].Var = CreateTensor(hmm[RCDmodel_index].state_no, max_M, VecSize);		/* Var[j][m][k] */
	
	hmm[RCDmodel_index].GammaAcc = CreateMatrix(hmm[RCDmodel_index].state_no, max_M);	/* GammaAcc[j][m] */
	hmm[RCDmodel_index].GammaOtAcc = CreateTensor(hmm[RCDmodel_index].state_no, max_M, VecSize);	/* GammaOtAcc[j][m][k] */
	hmm[RCDmodel_index].GammaOt2Acc = CreateTensor(hmm[RCDmodel_index].state_no, max_M, VecSize);	/* GammaOt2Acc[j][m][k] */
}
