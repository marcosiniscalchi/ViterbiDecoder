
#include "BasicDef_MAP.h"

int main(int argc, char **argv)
{
	FILE *f_list;
	char file_name[20];
	int q, Q;
	float LogP;

	int temp[1000];

/*	if (argc!=4)
		Abort("USAGE", "EXEC InitFileName ModelListName OldHmmPath NewHmmPath RecPath FileList");
*/
	MAP_Initialize(argv[1]);
	strcpy(OldHmmPath, argv[2]);
	strcpy(NewHmmPath, argv[3]);

        printf("-------------- Initialization -------------\n");
	printf("FeaDim:      %5d\n", FeaDim);
	printf("MaxObsNum:   %5d\n", MaxObsNum);
	printf("MaxFrameNum: %5d\n", MaxFrameNum);
	printf("MaxPhoneNum: %5d\n", MaxPhoneNum);
	printf("MaxStateNum: %5d\n", MaxStateNum);
	printf("MaxMixNum:   %5d\n", MaxMixNum);
	printf("ModelNum:    %5d\n", ModelNum);
	printf("VecSize:     %5d\n", VecSize);
	printf("AdaptFileList: %s\n", AdaptFileList);
	printf("CofPath:       %s\n", CofPath);
	printf("RecPath:       %s\n", RecPath);
	printf("OldHmmPath:    %s\n", OldHmmPath);
	printf("NewHmmPath:    %s\n", NewHmmPath);
	printf("TAU: %4.4f\n", TH);
        printf("-------------------------------------------\n");

        LoadModel();

	/*ClearMAPAcc(); marco's comment*/ 

	f_list = fopen(AdaptFileList, "r");
	if (f_list == NULL)
		Abort("ERROR opening file", AdaptFileList);

	q = 0;	
	sentence.CLEARED = 1;
	while(fscanf(f_list, "%s", file_name)!=EOF)
	{
		printf("file_name = %s",file_name);/*marco's comment*/
		GetAdaptData(file_name, q);
		GetSenModel(q);

		LogP = ViterbiDecode();
		printf("===== Average LogP: %.6f\n", LogP/sentence.T);
		FreeMemSentence();
		q++;
	}
	printf("\n");
}
