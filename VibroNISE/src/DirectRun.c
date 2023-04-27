#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "ToolsMod.h"
#include "ParmsHandleMod.h"

#include <string.h>


void GenerateFileName(char *FileName, const ParmStruct *Parms);
void CreateInputFile(ParmStruct *Parms);


int main(int ArgsTotal, char *Args[])
{
	ParmStruct *Parms;
	int FilePos = 0, ScriptToggle = 0, ScriptTotal;
	char FileName[100], Command[200] = "", CommandLine[100];
	FILE *FileID1, *FileID2, *FileID3, *FileID4;
	
	if (ArgsTotal > 3) sscanf(Args[3], "%d", &ScriptTotal);
	else ScriptTotal = 1;
	
	if (ArgsTotal > 2) sprintf(Command, "%s", Args[2]);
	else ScriptTotal = false;
	
	InitializeLog(false);
	
	PrintTitle(false);
	PrintMessage("                // Direct Run //", "", "", false, false);
	PrintMessage("", "", "", false, false);
	
	Parms = (ParmStruct *) calloc(1, sizeof(ParmStruct));
	InitializeParms(Parms);
	
	if (ScriptTotal)		FileID1 = fopen("Script1.sh", "w"), fprintf(FileID1, "#!/bin/sh\nsource ~/.bashrc\n");
	if (ScriptTotal > 1)	FileID2 = fopen("Script2.sh", "w"), fprintf(FileID2, "#!/bin/sh\nsource ~/.bashrc\n");
	if (ScriptTotal > 2)	FileID3 = fopen("Script3.sh", "w"), fprintf(FileID3, "#!/bin/sh\nsource ~/.bashrc\n");
	if (ScriptTotal > 3)	FileID4 = fopen("Script4.sh", "w"), fprintf(FileID4, "#!/bin/sh\nsource ~/.bashrc\n");
		
	while (FilePos != -1) {
		
		ReadInputFile(Args, Parms, &FilePos);
		GenerateFileName(FileName, Parms);
		
		if (ScriptTotal) {
			if (ScriptToggle == 0 || ScriptToggle == 2 * ScriptTotal - 1)		fprintf(FileID1, "cd %s\n%s\ncd ..\n", FileName, Command);
			else if (ScriptToggle == 1 || ScriptToggle == 2 * ScriptTotal - 2)	fprintf(FileID2, "cd %s\n%s\ncd ..\n", FileName, Command);
			else if (ScriptToggle == 2 || ScriptToggle == 2 * ScriptTotal - 3)	fprintf(FileID3, "cd %s\n%s\ncd ..\n", FileName, Command);
			else																fprintf(FileID4, "cd %s\n%s\ncd ..\n", FileName, Command);
			ScriptToggle = (ScriptToggle + 1) % (2 * ScriptTotal);
		}
		
		sprintf(CommandLine, "mkdir %s", FileName);
		system(CommandLine);
		
		CreateInputFile(Parms);
		
		sprintf(CommandLine, "mv _Input.dat %s/Input.dat", FileName);
		system(CommandLine);
	}
	
	if (ScriptTotal)		system("chmod +x Script1.sh"), fclose(FileID1);
	if (ScriptTotal > 1)	system("chmod +x Script2.sh"), fclose(FileID2);
	if (ScriptTotal > 2)	system("chmod +x Script3.sh"), fclose(FileID3);
	if (ScriptTotal > 3)	system("chmod +x Script4.sh"), fclose(FileID4);
		
	free(Parms);
	
	return 0;
}


void GenerateFileName(char *FileName, const ParmStruct *Parms)
{
	sprintf(FileName, "");
	sprintf(FileName, "%sU_%d", FileName, Parms->UnitTotal);
	//sprintf(FileName, "%sD_%d_T_%d", FileName, (int) Parms->BathStaticDev, (int) Parms->Temperature);
}


void CreateInputFile(ParmStruct *Parms)
{
	FILE *FileID;
	
	FileID = fopen("_Input.dat", "w");
	
	if (Parms->Mode == ModeAnnihilation)	fprintf(FileID, "Mode\nAnnihilation\n");
	
	if (Parms->Periodic)					fprintf(FileID, "Periodic\n");
	
	if (Parms->SampleTotal != 1)			fprintf(FileID, "SampleTotal\n%d\n", Parms->SampleTotal);
	if (Parms->UnitTotal != 1)				fprintf(FileID, "UnitTotal\n%d\n", Parms->UnitTotal);
	if (Parms->XCouplingNN != 0)			fprintf(FileID, "XCouplingNN\n%f\n", Parms->XCouplingNN);
	if (Parms->XPointDipoleRadius != 1)		fprintf(FileID, "XPointDipoleRadius\n%d\n", Parms->XPointDipoleRadius);
	if (Parms->BathStaticDev != 0)			fprintf(FileID, "BathStaticDev\n%f\n", Parms->BathStaticDev);
	if (Parms->Temperature != 0)			fprintf(FileID, "Temperature\n%f\n", Parms->Temperature);
	
	fclose(FileID);
}
