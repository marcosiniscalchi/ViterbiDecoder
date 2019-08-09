/*-------------------------------------------------------*/
/*                 ViterbiDecoding.c                     */
/*         Module: Viterbi Decoding Algorithm            */
/*                                                       */
/* (1) Do Viterbi Decoding to find optimal               */
/*     state sequence                                    */
/*                                                       */
/*                             Last Update: 2000/05/27   */
/*                             Kenneth, Kuan-Ting Chen   */
/* Bugs Report: please contact kenneth@iis.sinica.edu.tw */
/*-------------------------------------------------------*/
#include "BasicDef_MAP.h"
/*---------------- Do Viterbi Docoding to find state alignment -----------*/
float ViterbiDecode(void)  /* return LogP */
{
	FILE *f_debug;
   	int i, j, k, index, m, n ,t;
	iArray segment;
	int current_phone;
	int current_state;
   	float max_prob;
   	float alternative_prob;
   	float log_likeli;
   	float LogP;
   	int frame_no;
   	int new_N = sentence.state_no;
	int total_state = sentence.phone_no * 2 + new_N;
	int phone_index, state_index;
	int ENTRY=-100;
	int EXIT=100;
	int NUL=-500;
	int entry_state;
	iMatrix state_table;	/* (e.g.) state_table[0][0] = -100  */
				/*        state_table[1][0] = 0     */
				/*        state_table[2][0]= 1      */
				/*        state_table[3][0] = 100   */
				/*        state_table[4][0] = -100  */
				/*        state_table[5][0] = 0     */
/* state_table[state_index][0]: specifying the state property */						
/* state_table[state_index][1]: specifying the ID of the phone belonged to */
/* state_table[state_index][2]: specifying the sequential NO. within the transcription of the phone */ 
/* state_table[state_index][3]: specifying the corresponding ENTRY state index */
/* state_table[state_index][4]: specifying the sequential NO. of the emitting state */

	float trans_prob;
   	int max_index;    	/* previous state */
	Matrix likeli_table;	/* likeli_table[t][j] */
	Matrix delta;		/* delta[t][j]: for forward probability */
	iMatrix psy;		/* psy[t][j] = 	argmax (delta[t-1][j]*a[i][j]) */
	       			/*	       1<=i<=N    		       */

	f_debug=fopen("debug_vitedecode.txt","w");

	frame_no = sentence.T;
	segment = CreateiArray(frame_no);
	likeli_table = CreateMatrix(frame_no, new_N);
	state_table = CreateiMatrix(total_state, 5);
	delta = CreateMatrix(frame_no, total_state);
	psy = CreateiMatrix(frame_no, total_state);

	/* --- Preparation for the required information --- */
	state_index = 0;
	phone_index = 0;
	for (j=0; j<new_N; j++)
	{

		if (sentence.state_belong[j]==0)
		{
			state_table[state_index][0] = ENTRY;
			state_table[state_index][1] = sentence.phone_belong[j];
			state_table[state_index][2] = phone_index;
			state_table[state_index][3] = NUL;
			state_table[state_index][4] = NUL;
			entry_state = state_index;
			

			state_index++;
		
			state_table[state_index][0] = sentence.state_belong[j];
			state_table[state_index][1] = sentence.phone_belong[j];
			state_table[state_index][2] = phone_index;
			state_table[state_index][3] = entry_state;
			state_table[state_index][4] = j;
			

			
			/* NOTE: when the model contains only 1 state */
			if (sentence.state_belong[j] == (hmm[sentence.phone_belong[j]].state_no-1))
			{
			
				state_index++;
				state_table[state_index][0] = EXIT;
				state_table[state_index][1] = sentence.phone_belong[j];
				state_table[state_index][2] = phone_index;
				state_table[state_index][3] = entry_state;			
				state_table[state_index][4] = NUL;
				

				phone_index++;

			}
		}
		else if (sentence.state_belong[j] == (hmm[sentence.phone_belong[j]].state_no-1))
		{
			state_table[state_index][0] = sentence.state_belong[j];
			state_table[state_index][1] = sentence.phone_belong[j];
			state_table[state_index][2] = phone_index;
			state_table[state_index][3] = entry_state;
			state_table[state_index][4] = j;
			

			state_index++;
		
			state_table[state_index][0] = EXIT;
			state_table[state_index][1] = sentence.phone_belong[j];
			state_table[state_index][2] = phone_index;
			state_table[state_index][3] = entry_state;			
			state_table[state_index][4] = NUL;

			phone_index++;
			
			
		}
		else
		{
		
			
			state_table[state_index][0] = sentence.state_belong[j];
			state_table[state_index][1] = sentence.phone_belong[j];
			state_table[state_index][2] = phone_index;
			state_table[state_index][3] = entry_state;
			state_table[state_index][4] = j;
						
		}
		
		state_index++;	
	} 
  fprintf(f_debug, "state_index: %d, total_states: %d\n",state_index,total_state);
/* Marco's comment
	for (j=0; j<state_index; j++)
	{
		fprintf(f_debug, "[%03d] %4d %4d(%s)\t%4d %4d %4d\n", j, state_table[j][0], state_table[j][1], phone_table[state_table[j][1]], state_table[j][2], state_table[j][3], state_table[j][4]);
	}
*/

	/* --- Pre-processing --- */
	if (state_index!=total_state)
		Abort("ERROR", "total_state invalid");
	if (phone_index!=sentence.phone_no)
		Abort("ERROR", "total_phone_no invalid");

	for (t=0; t<frame_no; t++)
	{
		for (j=0; j<new_N; j++)
		{
			current_phone = sentence.phone_belong[j];
			current_state = sentence.state_belong[j];
		
			likeli_table[t][j] = CalStateLikeli(1, sentence.mix_no[j], sentence.feature[t], hmm[current_phone].Mean[current_state], hmm[current_phone].Var[current_state], hmm[current_phone].MixWeight[current_state], hmm[current_phone].GConst[current_state]);
		}
		
	}

	/* --- Initialization --- */
	
	t=0;
	for (state_index=0; state_index<total_state; state_index++)
	{
		if (state_index==0)
		{
			delta[t][state_index] = 0.0;
			psy[t][state_index] = -1;
			
		}
		else
		{
			if (state_table[state_index][0]==ENTRY)
			{
				
				index = state_index - 1;	/* previous EXIT state */
				entry_state = state_table[index][3];	/* corresponding ENTRY state */
				
				phone_index = state_table[index][1];	/* corresponding phone ID */	
				n = hmm[phone_index].state_no;  	/* e.g. 3 */
				trans_prob = TakeLog(sentence.TransP[entry_state][n+2-1]);
	
				delta[t][state_index] = delta[t][entry_state] + trans_prob;
				psy[t][state_index] = entry_state;
				
			}
			else if (state_table[state_index][0]==EXIT)
			{
				entry_state = state_table[state_index][3];
				phone_index = state_table[state_index][1];
				n = hmm[phone_index].state_no;				

				k=entry_state+1;
				max_index = k;
				trans_prob = TakeLog(sentence.TransP[k][n+2-1]);
				max_prob = delta[t][k] + trans_prob;

				for (k=entry_state+2; k<state_index; k++)
				{	/* find optimal previous state */
					trans_prob = TakeLog(sentence.TransP[k][n+2-1]);
					alternative_prob = delta[t][k] + trans_prob;
					if (alternative_prob > max_prob)
					{
						max_index = k;
						max_prob = alternative_prob;
					}
				}
				delta[t][state_index] = max_prob;				
				psy[t][state_index] = max_index;
				
			}	
			else	/* emitting states */
			{ 
				n = state_table[state_index][0];
				j = state_table[state_index][4];
				entry_state = state_table[state_index][3];
				
				
				trans_prob = TakeLog(sentence.TransP[entry_state][n+1]);
				delta[t][state_index] = delta[t][entry_state] + trans_prob + likeli_table[t][j];
				psy[t][state_index] = entry_state;
				
			} 
		}
		
	}
	
	

	/* --- Recursion --- */
for (t=1; t<frame_no; t++)
{

	for (state_index=0; state_index<total_state; state_index++)
	{
		if (state_index==0)
		{
			delta[t][state_index] = MinProb;
			psy[t][state_index] = -1;
			
		}
		else
		{
			if (state_table[state_index][0]==ENTRY)
			{
				index = state_index - 1;	/* previous EXIT state */
				entry_state = state_table[index][3];	/* corresponding ENTRY state */
				phone_index = state_table[index][1];	/* corresponding phone ID */	
				n = hmm[phone_index].state_no;

				trans_prob = TakeLog(sentence.TransP[entry_state][n+2-1]);
				max_index = entry_state;	/* previous ENTRY state */
				max_prob = delta[t][entry_state] + trans_prob;				
				
				alternative_prob = delta[t-1][index] + 0.0;
				if (alternative_prob>max_prob)
				{
					max_index = index;	/* prebious EXIT state */
					max_prob = alternative_prob;
				}
				delta[t][state_index] = max_prob;
				psy[t][state_index] = max_index;
				
			}
			else if (state_table[state_index][0]==EXIT)
			{
			
				entry_state = state_table[state_index][3];
				phone_index = state_table[state_index][1];
				n = hmm[phone_index].state_no;				
				
				k=entry_state+1;
				max_index = k;
				trans_prob = TakeLog(sentence.TransP[k][n+2-1]);
				max_prob = delta[t][k] + trans_prob;
								
				for (k=entry_state+2; k<state_index; k++)
				{	/* find optimal previous state */
					
					trans_prob = TakeLog(sentence.TransP[k][n+2-1]);
					alternative_prob = delta[t][k] + trans_prob;
					if (alternative_prob > max_prob)
					{
						max_index = k;
						max_prob = alternative_prob;
					}
				}
				
				delta[t][state_index] = max_prob;				
				psy[t][state_index] = max_index;
				
			}	
			else	/* emitting states */
			{ 
				n = state_table[state_index][0];
				j = state_table[state_index][4];
				entry_state = state_table[state_index][3];
				
				trans_prob = TakeLog(sentence.TransP[entry_state][n+1]);
				max_index = entry_state;				
				max_prob = delta[t][entry_state] + trans_prob;
				for (k=entry_state+1; state_table[k][0]!=EXIT; k++)
				{	            /*^^^^^^^^^^^^^^^^^^^^^^^*/  	
					trans_prob = TakeLog(sentence.TransP[k][n+1]);
					alternative_prob = delta[t-1][k] + trans_prob;
					if (alternative_prob>max_prob)
					{
						max_index = k;
						max_prob = alternative_prob;
					}
				}
				delta[t][state_index] = max_prob + likeli_table[t][j];
				psy[t][state_index] = max_index;
				} 
		}
	}	
		
}


	/* --- Backtracking --- */
	
	t=frame_no-1;
	/* find optimal terminal state (which EXIT state) */
	max_index = total_state - 1;
	max_prob = delta[t][max_index];

	for (state_index = 0; state_index<(total_state-1); state_index++)
	{
		if (state_table[state_index][0]!=EXIT)
			continue;		
		alternative_prob = delta[t][state_index];
		if (alternative_prob>max_prob)
		{
			max_index = state_index;
			max_prob = alternative_prob;
		}
	}
	state_index = max_index;
	LogP = max_prob;
			
	/* back-trace */
	index = psy[t][state_index];            /* must be emitting state */	
	segment[t] = index;			/* segment[t] records the physical state ID */
	sentence.s[t] = state_table[index][4];  /* sentence.s[t] records the emitting state ID */
	if (sentence.s[t]==NUL)
		Abort("ERROR", "NUL appears in sentence.s[t]");	

	for (t=frame_no-2; t>=0; t--)
	{
		k=0;
		state_index = segment[t+1];	 /* next state: MUST BE EMITTING STATE */
		index = psy[t+1][state_index];   /* this state */
		
		/* pass the ENTRY states until this state is not ENTRY */
		while(state_table[index][0]==ENTRY)  
		{
			k = 1;
			j = psy[t+1][index];
			index = j;
		}
		/* if this state is EXIT, find preceeding state */
		if (state_table[index][0]==EXIT)
		{
			k = 0;
			j = psy[t][index];
			index = j;
		}
		else 
		{
			if (k==1)
				Abort("ERROR", "invalid state transition: ENTRY NOT from EXIT");
		}
		if (state_table[index][0]==ENTRY || state_table[index][0]==EXIT)
			Abort("ERROR", "invalid state transition: can't find emitting states");
		segment[t] = index;
		sentence.s[t] = state_table[index][4];
	}

	for (t=1; t<frame_no; t++)
	{
		if ((sentence.s[t]!=sentence.s[t-1]) && (sentence.s[t]!=(sentence.s[t-1]+1)))
			printf("SKIP!\n");
		fprintf(f_debug, "%d", sentence.s[t]);
		fprintf(f_debug, (t%20==0)? "\n":" ");
	}

	FreeiArray(segment);
	FreeiMatrix(psy, frame_no);
	FreeiMatrix(state_table, total_state);
	FreeMatrix(likeli_table, frame_no);
	FreeMatrix(delta, frame_no);


	fclose(f_debug);	
	
	return(LogP);

	

}

