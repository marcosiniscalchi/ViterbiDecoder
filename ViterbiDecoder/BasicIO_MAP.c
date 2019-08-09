/*-------------------------------------------------------*/
/*                    BasicIO_MAP.c                      */
/*              Module: Basic Definitions                */
/*                                                       */
/* (1) Read Model Parameters                             */
/* (2) Write New Model Parameters                        */
/* (3) Read Adaptation Data Feature                      */
/*                                                       */
/*                             Last Update: 2000/05/27   */
/*                             Kenneth, Kuan-Ting Chen   */
/* Bugs Report: please contact kenneth@iis.sinica.edu.tw */
/*-------------------------------------------------------*/
#include "BasicDef_MAP.h"
/*------------- Allocate Memory for each single sentence -----------------*/
/*Swap for int or float*/
float float_swap(float value){
    union v {
        float       f;
        unsigned int    i;
    };
     
    union v val;
     
    val.f = value;
    val.i = htonl(val.i);
                   
    return val.f;
}

void Swap32 (void* a )
{

  long *Long = (long *)a;
  *Long = ((*Long&0x000000ffL)<<24 )| \
          ((*Long&0x0000ff00L)<<8 )| \
          ((*Long&0x00ff0000L)>>8 )| \
          ((*Long&0xff000000L)>>24 );
}


/*Swap for short*/
void Swap16 ( short *Short )
{
  *Short  = ((*Short & 0x00ff) << 8) | \
      ((*Short & 0xff00) >> 8);
}


void AllocateMemSentence(void)
{
	int i;
	int T = sentence.T;
	int phone_no = sentence.phone_no;
	int new_N = sentence.state_no;	/* total state number within that sentence */
	int max_M = MaxMixNum;		/* maximum mixture number within a state */
	int total_state = phone_no * 2 + new_N;

   	if (sentence.CLEARED != 1)
		Abort("ERROR", "Sentence not cleared!");
   	else
		sentence.CLEARED = 0;

	sentence.feature = CreateMatrix(T, FeaDim);
	sentence.s = CreateiArray(T);
	/* sentence.gamma = CreateTensor(T, new_N, max_M); 
	ZeroTensor(sentence.gamma, T, new_N, max_M); */

	sentence.phone = CreateiArray(phone_no);
	sentence.state_no_this_phone = CreateiArray(phone_no);

	sentence.phone_belong = CreateiArray(new_N);
	sentence.state_belong = CreateiArray(new_N);
	sentence.mix_no = CreateiArray(new_N);

	/* sentence.Mean = CreateMatrixOfVector(new_N, max_M);
	sentence.Var = CreateMatrixOfVector(new_N, max_M);
	sentence.MixWeight = CreateVectorOfVector(new_N);
	sentence.GConst = CreateVectorOfVector(new_N); */

	sentence.TransP = CreateVectorOfVector(total_state);
	/* ZeroMatrix(sentence.TransP, total_state, MaxStateNum); */
}


/*---------- free memory of observations: for ALL observ_seq -------------*/
void FreeMemSentence(void)
{
	int obs_id = sentence.obs_id;
	int i;
	int T = sentence.T;
	int new_N = sentence.state_no;
	int max_M = MaxMixNum;
	int total_state = (sentence.phone_no)*2+new_N;

	if (sentence.CLEARED == 1)
		Abort("ERROR", "sentence already cleared");

	FreeMatrix(sentence.feature, T);
	FreeiArray(sentence.s);
	/* FreeTensor(sentence.gamma, T, new_N); */
	
	FreeiArray(sentence.phone);
	FreeiArray(sentence.state_no_this_phone);

	FreeiArray(sentence.phone_belong);
	FreeiArray(sentence.state_belong);
	FreeiArray(sentence.mix_no);

	/* FreeMatrixOfVector(sentence.Mean, new_N);
	FreeMatrixOfVector(sentence.Var, new_N);
	FreeVectorOfVector(sentence.MixWeight);
	FreeVectorOfVector(sentence.GConst); */
	FreeVectorOfVector(sentence.TransP);

	sentence.CLEARED = 1;

}

/*------------------------ Get Adaptation Data ---------------------------*/
int GetAdaptData(char *file_name, int obs_id)
{                    /* e.g. 100001_01.cof */
	struct featype header;        /* for read_in features */
	FILE *f_cof, *f_lab;          /* for open aflist and .cof, .lab file */
	char phone_name[10];          /* specifying phone name */
	char cof_name[150];           /* for specifying .cof file */
	char lab_name[150];           /* for specifying .lab file */
	char model_name[10];
	char *tmpstr1, *tmpstr2,temp_str[50];
	int x, i, j;
	int k, now_phone;
	int phone_no;                 /* number of phones within an utterance */
	int frame_num;                /* number of total frames */
	long data_length;              /* total data length */
	int begin, end, tmp1, tmp2;   /* begin and end time index */
	float tmp;
	float temp_data[100000];             /* for temporary readed feature */
	iArray tmp_sub_id;

  


	/* --- initialize --- */
	tmp_sub_id = CreateiArray(MaxPhoneNum);

	/* --- file name processing --- */
	tmpstr1 = strtok(file_name, ".");
	printf("\n===== file: %s =====\n", tmpstr1);
	sprintf(cof_name,"%s/%s.%s", CofPath, tmpstr1, FeaExtension);
	sprintf(lab_name,"%s/%s.%s", RecPath, tmpstr1, LabExtension);
	printf("cof_name = %s, lab_name = %s\n",cof_name, lab_name);/*marco's comment*/
	
	/* open *.cof file */
	fprintf(stdout, "[%03d] %s: ", obs_id, cof_name);
	f_cof = fopen(cof_name,"rb");   /* get oooooo_xx.cof */
	if (f_cof == NULL) 
		Abort("ERROR opening file", cof_name);

	/* read feature header */
	fseek(f_cof, 0, 2);
	data_length = ftell(f_cof) / ( sizeof(float) * FeaDim );
	rewind(f_cof);

	fread(&header, sizeof(struct featype), 1, f_cof);
	begin = 0;
	fseek(f_cof, begin/50000*sizeof(float)*FeaDim,1);

  Swap32(&header.number); /*  htk format, reverse byte order  */
  Swap32(&header.shift);  /*  htk format, reverse byte order  */
  Swap16(&header.size); /*  htk format, reverse byte order  */
  Swap16(&header.kind); /*  htk format, reverse byte order  */


	/* Error checking */
	
	printf("header.number = %ld\n", header.number);
	printf("data_length = %ld\n", data_length);
	if (header.number!=data_length) 
	{
		fprintf(stderr, "ERROR! data_length disagree!\n");
		fprintf(stderr, "data_length = %ld\n", data_length);
		fprintf(stderr, "header.number = %ld\n", header.number);
		exit (-1);
	}
	fprintf(stdout, "total %ld frames\n", data_length);

	/* open .rec file */
	fprintf(stdout, "Marco: [%03d] %s: ", obs_id, lab_name);
	f_lab = fopen(lab_name, "r");
	if (f_lab == NULL) 
	{
		Abort("ERROR opening file", lab_name);	
	}

	/* get transcription */
	j = 0;   /* j-th sub-model */

	fscanf(f_lab, "%s\n", temp_str);
	fscanf(f_lab, "%s\n", temp_str);

	while (fscanf(f_lab, "%s\n", temp_str)!=EOF) 
	{

	    if (strcmp(temp_str,".")==0)
		{
			break;
		}
		else 
		{
			fscanf(f_lab, "%s\n", temp_str);
	        fscanf(f_lab, "%s\n", model_name);
		for (i=0; i<ModelNum; i++) 
		{
			if (strcmp(hmm[i].phone_name, model_name) == 0) 
			{  /* got hmm[i] */
				tmp_sub_id[j] = i;
				/*printf("\nmodel name, hmm.phone_name [%d %d], model_name=%s\n", j, i,model_name); marco's comment*/
				break;
			}
		}
		j++;
		}
	}  /* end while, indicating EOF of current .lab */

	/* get feature */
	sentence.obs_id = obs_id;
	sentence.T = data_length;
	sentence.phone_no = j;
	sentence.state_no = 0;
	for (k=0; k<sentence.phone_no; k++)
	{
		now_phone = tmp_sub_id[k];
		sentence.state_no += hmm[now_phone].state_no;
	}

	fprintf(stdout, "total %d phones, %d states\n", sentence.phone_no, sentence.state_no);
	
	AllocateMemSentence();

	i = fread(temp_data, sizeof(float), data_length*FeaDim, f_cof);
  for(i = 0; i < data_length*FeaDim; i++)
     temp_data[i] = float_swap(temp_data[i]);

 for(i = 0; i<1; i++)
  {
    for(j=0; j<FeaDim; j++){


      printf("f[%d][%d]=%f  ",i,j,temp_data[j]);}
                        printf("\n");
  }

	for (i=0; i<sentence.T; i++)
		for (j=0; j<FeaDim; j++)
			sentence.feature[i][j] = temp_data[i*FeaDim+j];

	for (i=0; i < sentence.phone_no; i++) 
	{
		now_phone = tmp_sub_id[i];
		sentence.phone[i] = now_phone;
		sentence.state_no_this_phone[i] = hmm[now_phone].state_no;
	}

	fclose(f_cof);
	fclose(f_lab);
	return obs_id;
}

/* --------------------- Get Cascaded HMM parameters -------------------- */
void GetSenModel(int obs_id)
{
	int i, j;
   	int current_phone;
   	int n, m, d, n2, n3;
	int state_index;
	int state_index2;
	int N, M;
	int total_state;
	FILE *f_debug;

	if ((f_debug = fopen("debug_BasicIO.txt", "w"))==NULL)
		Abort("ERROR", "can not open \"debug_BasicIO.txt\"");	

   	if (obs_id != sentence.obs_id) 
	{
      		fprintf(stderr, "ERROR in GetSenModel()!\n");
      		exit(-1);
   	}

	state_index = 0;
	state_index2 = 0;
	for (i=0; i< sentence.phone_no; i++) 
	{  	/* for each sub model */
      
		current_phone= sentence.phone[i];   /* current phone index */
		N = sentence.state_no_this_phone[i];
		if (N!=hmm[current_phone].state_no)
			Abort("ERROR", "N!=hmm[current_phone].state_no");		

		for (n=0; n<N; n++) 
		{
			sentence.phone_belong[state_index] = current_phone;
			sentence.state_belong[state_index] = n;
			sentence.mix_no[state_index] = hmm[current_phone].mix_no[n]; 
			M = sentence.mix_no[state_index];				

			fprintf(f_debug, "[%3d] current phone:[%3d]%s \tstate:%d M:%d\n", state_index, sentence.phone_belong[state_index], hmm[current_phone].phone_name, sentence.state_belong[state_index], sentence.mix_no[state_index]);
					
			state_index++;

      		} /* end loop for state index */
		
		for (n2=0; n2<N+2; n2++)
		{
			for (n3=0; n3<N+2; n3++)			
			{
				sentence.TransP[state_index2] = hmm[current_phone].TransP[n2];
			}
			state_index2++;
		}

   	} /* end loop for sub-sentence model index */

	if (state_index!=sentence.state_no)
		Abort("ERROR", "state_index!=sentence.state_no");
	total_state = sentence.state_no + 2 * sentence.phone_no;
	if (state_index2!=total_state)
		Abort("ERROR", "state_index2!=total_state");

	fclose(f_debug);
}   

/*--------------------- Read in Models from HMM Files --------------------*/
int ReadHMM(char *model_name, struct HMM_struct *model, int RCDmodel_index, int VecSize)
{
	char tmpstr[30], HMM_name[100];
    char debug_name[30];
	FILE *f_HMM;
	FILE *f_debug;
	int state, now_state, mix, now_mix, dim, i, j;
	float gg;
	float temp;
    int flag;
	/* open file for HMM definition */
	strcpy(HMM_name, OldHmmPath);
	strcat(HMM_name, "/");
	strcat(HMM_name, model_name);     /* specifying the path to open HMM file */
	sprintf(debug_name,"debug");
    
	f_debug = fopen(debug_name,"w");
	if (f_debug == NULL) 
	{  /* cannot open HMM file */
		Abort("ERROR opening Debug", debug_name);
	}

	f_HMM = fopen(HMM_name,"r");
	if (f_HMM == NULL) 
	{  /* cannot open HMM file */
 		Abort("ERROR opening HMM", HMM_name);
	}

	/* scan f_HMM for HMM detail */
	now_mix = 1;
	for (i=0; i<MaxStateNum; i++) 
	{       /* useful for the case mix_no = 1 */
		model->mix_no[i] = now_mix;
	}

	model->phone_ID = RCDmodel_index;    /* get phone_ID */

	fscanf(f_HMM, "%s", tmpstr);         /* scan for 1st string */

	while (!feof(f_HMM)) 
	{               
		if ((strcmp(tmpstr,"<BEGINHMM>") == 0)||(strcmp(tmpstr,"<BeginHMM>") == 0)) 
		{
	 		;
      	}
      		else if ((strcmp(tmpstr,"<NUMSTATES>") == 0)||(strcmp(tmpstr,"<NUMSTATES>") == 0)) 
		{
	 		fscanf(f_HMM, "%d ", &state);     /* get number of states */
	 		model->state_no = state;          /* and assign it to state_no */
      		}
  
      		else if (strcmp(tmpstr,"<STREAMINFO>") == 0)
		{
	 		fscanf(f_HMM, "%s ", tmpstr);   /* skip twice the items following */
	 		fscanf(f_HMM, "%s ", tmpstr);   /* <StreamInfo> 1 30 assumed */
      		}
 
      	   else if ((strcmp(tmpstr,"<STATE>") == 0)||(strcmp(tmpstr,"<State>") == 0) )
		{
	 		fscanf(f_HMM, "%d ", &now_state);   /* get current_state */
			mix=0;
		//	if (now_state==4)
		//		printf("error here!!");
      		}
		   else if (strcmp(tmpstr,"<VECSIZE>") == 0) 
		{
	 		fscanf(f_HMM, "%s ", tmpstr);   /* get current_state */
		//	if (now_state==4)
		//		printf("error here!!");
      		}
      		else if ((strcmp(tmpstr,"<NUMMIXES>") == 0)||(strcmp(tmpstr,"<NumMixes>") == 0))
		{
	 		fscanf(f_HMM, "%d ", &now_mix);   /* get current_mix */
			model->mix_no[now_state-1] = now_mix;
      		}
      	else if ((strcmp(tmpstr,"<MIXTURE>") == 0)||(strcmp(tmpstr,"<Mixture>") == 0) )
		{
	 		mix++;
			fscanf(f_HMM, "%d ", &now_mix);   /* get current_mix */
	 		fscanf(f_HMM, "%f ", &gg);        /* get mixweight */
	 		model->mixweight[now_state-1][mix-1] = gg;
			model->mix_flag[now_state-1][mix-1] = 1;
      		}				    /*********/
		else if ((strcmp(tmpstr,"<MEAN>") == 0) ||(strcmp(tmpstr,"<Mean>") == 0))
		{
	 		fscanf(f_HMM, "%d ", &dim);         /* get dim for Mean */
	 		if (dim != VecSize) 
			{               /* check if dim == VectorSize */
	    			printf("ERROR! Vector Sizes Disagree!\n");
	    			return (-1);
	 		}
	 		for (i=0; i<VecSize; i++) 
			{  		/* get mean value for current Gaussian */
	    			fscanf(f_HMM, "%f ", &gg);
	    			model->mean[now_state-1][mix-1][i] = gg;
					fprintf(f_debug,"Mean %d %f\n",i,model->mean[now_state-1][mix-1][i]);
	 		}
      		}
      		else if ((strcmp(tmpstr,"<VARIANCE>") == 0)||(strcmp(tmpstr,"<Variance>") == 0)) 
	    	{
			flag=0;
	 		fscanf(f_HMM, "%d ", &dim);         /* get dim for Var */
	 		if (dim != VecSize) 
			{            /* check if dim == VectorSize */
	    			printf("ERROR! Vector Sizes Disagree!\n");
	    			return (-1);
	 		}
	 		for (i=0; i<VecSize; i++) 
			{  		/* get var value for current Gaussian */
	    			fscanf(f_HMM, "%f ", &gg);
					model->var[now_state-1][mix-1][i] = 1.0 / gg;
					fprintf(f_debug,"Variance %d %f\n",i,model->var[now_state-1][mix-1][i]);
				//	if (gg<1e-007)
				//	{
				//	   model->var[now_state-1][mix-1][i] = model->var[now_state-1][mix-1][i-1];
				//	}
	    			
	 		}     
			for (i=0; i<VecSize; i++) 
			{ 
				temp=1/(model->var[now_state-1][mix-1][i]);
				model->Gconst[now_state-1][mix-1] = model->Gconst[now_state-1][mix-1]+0.5*log(temp);
			}
			   model->Gconst[now_state-1][mix-1] = model->Gconst[now_state-1][mix-1]+(VecSize/2.0)*log(2.0*MATH_PI);
			   model->Gconst[now_state-1][mix-1] = 2*model->Gconst[now_state-1][mix-1];
      		}
      		else if (strcmp(tmpstr,"<GCONST>") == 0) 
		    {
	 		fscanf(f_HMM, "%f ", &gg);        /* get GConst and assign it to */
	 		model->Gconst[now_state-1][mix-1] = gg;  /* note:now_state begins at 2! */
      		}
      		else if ((strcmp(tmpstr,"<TRANSP>") == 0)||(strcmp(tmpstr,"<TransP>") == 0)) 
		{
	 		fscanf(f_HMM, "%s ", tmpstr);     /* skip the following next 1 figure */
	 		for (i=0; i<state; i++) 
			{
	    			for (j=0; j<state; j++) 
				{   /* get TransP values for the model */
	       				fscanf(f_HMM, "%e ", &gg);
	       				if (gg < 0) 
					{
		  				gg = 1.0e-30;
	       				}
	       				model->transP[i][j] = gg;
	    			}     		            /*^^^^^^^*/
	 		}
      		}
      		else if ((strcmp(tmpstr,"<ENDHMM>") == 0)||(strcmp(tmpstr,"<EndHMM>") == 0)) 
			{
	 			break;
      		}
      		else 
		    {
	 			/*printf("%s Not Defined!\n",tmpstr); marco's comment*/
      		}
      		fscanf(f_HMM, "%s", tmpstr);
   	} /* END OF while */
	fclose(f_debug);
   	fclose(f_HMM);	
   return (1);
} /* END OF Read_HMM */
/*----------------------- Get Model Parameters ---------------------------*/
int LoadModel(void)
{
	/* using temp_model to call read_HMM, storing the model parameter in hmm[] */
   	struct HMM_struct temp_model;    /* declare temp_model for read_HMM() */
   	int state_index, state_index_2, mix_index, vec_index, RCDmodel_index;
   	int n, n2; /*N_index*/
   	int m; /*M_index*/
	int N;
//	int M;
   	int READ_HMM_SUCCESS;
/*
	printf("Initializing temp_model...\n");
*/
	TotalGaussianNo = 0;
	InitTempModel(&temp_model);
   	printf("Loading HMM: ");
   	for (RCDmodel_index=0; RCDmodel_index < ModelNum; RCDmodel_index++) 
	{
	/* Initialization (clear old value) of temp_model, which is of HMM_struct */


   	    	temp_model.phone_ID = 0;
      		temp_model.state_no = 0;
      		for (state_index=0; state_index < MaxStateNum; state_index++) 
			{
 	 		temp_model.mix_no[state_index] = 0;
	 		for (state_index_2=0; state_index_2 < MaxStateNum; state_index_2++) 
			{
	    			temp_model.transP[state_index][state_index_2] = 0.0;
	 		}
	 		for (mix_index=0; mix_index < MaxMixNum; mix_index++) 
			{
	    			temp_model.mixweight[state_index][mix_index] = 0.0;
	    			temp_model.Gconst[state_index][mix_index] = 0.0;
				    temp_model.mix_flag[state_index][mix_index] = 0;
	    			for (vec_index=0; vec_index < VecSize; vec_index++) 
				{
	       				temp_model.mean[state_index][mix_index][vec_index] = 0.0;
	       				temp_model.var[state_index][mix_index][vec_index] = 0.0;
	    		} /* end loop for vec_index */
			} /* end loop for mix_index */
			} /* end loop for state_index */

      /* Call read_HMM(), storing parameters in temp_model */
   		READ_HMM_SUCCESS = ReadHMM(phone_table[RCDmodel_index], &temp_model, RCDmodel_index,vec_index);
    //ReadHMM(phone_table[RCDmodel_index], &temp_model, RCDmodel_index);
      		if (!READ_HMM_SUCCESS) 
		{
	 		fprintf(stderr, "ERROR in LoadHMM()!\n");
	 		exit (-1);
      		} 
		if (RCDmodel_index==0)
			printf("%3d/%3d", RCDmodel_index+1, ModelNum);
      		else
	 		printf("\b\b\b\b\b\b\b%3d/%3d", RCDmodel_index+1, ModelNum);


      /* Allocate memory for hmm[RCDmodel_index] */
     		strcpy(hmm[RCDmodel_index].phone_name, phone_table[RCDmodel_index]);
     		hmm[RCDmodel_index].state_no = temp_model.state_no - 2; /* emitting states */
		InitHmm(RCDmodel_index, temp_model.state_no, temp_model.mix_no); 

      /* Store parameters in temp_model into hmm[RCDmodel_index] */
		N = hmm[RCDmodel_index].state_no;
		if (N!=(temp_model.state_no-2))
			Abort("ERROR", "N != temp_model.state_no-2");
 
		for (n=0; n<N+2; n++)	/* n actually means state (n+2) */
		{
			for (n2=0; n2<N+2; n2++)
			{
				hmm[RCDmodel_index].TransP[n][n2] = temp_model.transP[n][n2];			   			
			}
		}			
		for (n=0; n<N;n++)
		{
			hmm[RCDmodel_index].mix_no[n] = temp_model.mix_no[n+1];
			mix_index = 0;
			for (m=0; m<temp_model.mix_no[n+1]; m++) 
			{
	    		hmm[RCDmodel_index].MixWeight[n][mix_index] = temp_model.mixweight[n+1][m];
				hmm[RCDmodel_index].GConst[n][mix_index] = temp_model.Gconst[n+1][m];
			
	    			for (vec_index=0; vec_index < VecSize; vec_index++) 
				{
	       				hmm[RCDmodel_index].Mean[n][mix_index][vec_index] = temp_model.mean[n+1][m][vec_index];
	       				hmm[RCDmodel_index].Var[n][mix_index][vec_index] = temp_model.var[n+1][m][vec_index];
				} /* end loop for vec_index */
				mix_index++;

	 		} /* end loop for m */
      		} /* end loop for n */
   	} /* end loop for RCDmodel_index */
  
	fprintf(stdout, "\nTotal %d Gaussian Loaded\n", TotalGaussianNo);
	return (RCDmodel_index);
} /* END OF LoadModel() */

/*------------------------- Write Adapted Models -------------------------*/
void WriteModel(void)
{
   	char new_file_name[150];
   	FILE *f_new;
   	int model_index;
   	int n, m, i, n2;
	int N, M;
	int StateNum;

   	fprintf(stdout, "\nWriting New HMM Files...\n");

   	for (model_index=0; model_index < ModelNum; model_index++) 
	{
		sprintf(new_file_name, "%s/%s", NewHmmPath, hmm[model_index].phone_name);
		f_new = fopen(new_file_name, "w");
		if (f_new == NULL) 
		{
	 		fprintf(stderr, "ERROR! cannot create HMM file %s\n", new_file_name);
	 		exit (-1);
      		}


		StateNum = hmm[model_index].state_no + 2;
		N = StateNum - 2;
      		fprintf(f_new, "~o\n");
      		fprintf(f_new, "<STREAMINFO> 1 39\n");
			fprintf(f_new, "<VECSIZE> 39<NULLD><MFCC_E_D_A><DIAGC>\n");
			fprintf(f_new, "~h \"%s\"\n",hmm[model_index].phone_name);
            fprintf(f_new, "<BEGINHMM>\n");
			fprintf(f_new, "<NUMSTATES> %d\n",StateNum);
      		for (n=0; n<N; n++)  /* state by state */ 
		{   
			fprintf(f_new, "<STATE> %d\n", n+2);
			M = hmm[model_index].mix_no[n];
	 		fprintf(f_new, "<NUMMIXES> %d\n", M);
	 		for (m=0; m<M; m++)  /* mixture by mixture */
			{  
	    			fprintf(f_new, "<MIXTURE> %d %e\n", m+1, hmm[model_index].MixWeight[n][m]);

				fprintf(f_new, "<MEAN> %d\n", VecSize);
	    			for (i=0; i<VecSize; i++)
	       				fprintf(f_new, " %e", hmm[model_index].Mean[n][m][i]);
	    			fprintf(f_new, "\n");

	    			fprintf(f_new, "<VARIANCE> %d\n", VecSize);
	    			for (i=0; i<VecSize; i++)
	       				fprintf(f_new, " %e", 1.0 / hmm[model_index].Var[n][m][i]);
	    			fprintf(f_new, "\n");

	    			fprintf(f_new, "<GCONST> %e\n", hmm[model_index].GConst[n][m]);
	 		} /* end loop for single mixture */
      		} /* end loop for single state */

      		fprintf(f_new, "<TRANSP> %d\n", StateNum);
		for (n=0; n<StateNum; n++)
			for (n2=0; n2<StateNum; n2++)
			{
				fprintf(f_new, "%e", hmm[model_index].TransP[n][n2]);
				if (n2==StateNum-1)
					fprintf(f_new, "\n");
				else
					fprintf(f_new, " ");
			}

      		fprintf(f_new, "<ENDHMM> \n");

      		fclose(f_new);
   	}
   	fprintf(stdout, "New HMM definitions written in %s...\n", NewHmmPath);
}
