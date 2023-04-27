#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "ToolsMod.h"
#include "ParmsHandleMod.h"


void SaveResponse2D(const char *FileName, const ComplexType *Response, const ParmStruct *Parms, const int MemSkip);


int main(int ArgsTotal, char *Args[])
{
	typedef struct { ComplexType GB, SE, EA; } ResponseStruct;
	
	ParmStruct *Parms;
	
	ResponseStruct *ResponseNR, *ResponseRR, *ResponseTotal;
	
	ComplexType *Pulse;
	
	double Dummy, ReadReal, ReadImag;
	
	int PriTimeCount, SecTimeCount, WaitTimeCount, PriPulseCount, SecPulseCount, PulseTimeMin, ResponseTimeMin, PriTimeSkip;
	int PriPulseOffset, SecPulseOffset, TerPulseOffset, PulsiveEntry, ImpulsiveEntry, PriPulseEntry, SecPulseEntry, TerPulseEntry;
	
	FILE *FileID;
	
	char FileName[100];
	
	
	PrintTitle(false);
	PrintMessage("                // Pulse Shaper //", "", "", false, false);
	PrintMessage("", "", "", false, false);
	
	Parms = (ParmStruct *) calloc(1, sizeof(ParmStruct));
	InitializeParms(Parms);
	ReadInputFile(Args, Parms, NULL);
	DerivedParms(Parms);
	
	
	PulseTimeMin = -Parms->PulseTimeTotal + 1;
	ResponseTimeMin = -Parms->ResponseTimeTotal + 1;
	PriTimeSkip = (2 * Parms->PulseTimeTotal - 1) * Parms->ResponseTimeTotal;
	
	Pulse = malloc((2 * Parms->PulseTimeTotal - 1) * sizeof(ComplexType));
	ResponseNR = calloc(IntPow2(Parms->ResponseTimeTotal), sizeof(ResponseStruct));
	ResponseRR = calloc(IntPow2(Parms->ResponseTimeTotal), sizeof(ResponseStruct));
	ResponseTotal = calloc((2 * Parms->ResponseTimeTotal - 1) * PriTimeSkip, sizeof(ResponseStruct));
	
	FileID = fopen(Parms->PulseFile, "r");
	for (PriTimeCount = 0; PriTimeCount < Parms->PulseTimeTotal; PriTimeCount++) {
		fscanf(FileID, "%lf %lf %lf\n", &Dummy, &ReadReal, &ReadImag);
		Pulse[IntToArrayEntry(PriTimeCount, PulseTimeMin)] = (ReadReal + I * ReadImag);
		Pulse[IntToArrayEntry(-PriTimeCount, PulseTimeMin)] = (ReadReal - I * ReadImag);
	}
	fclose(FileID);
	
	for (WaitTimeCount = PulseTimeMin; WaitTimeCount < Parms->PulseTimeTotal; WaitTimeCount++) {
		
		sprintf(FileName, "%s%d%s", "../../WT_", Parms->WaitTimeTotal + WaitTimeCount, "/ResponseGBNR.dat");
		FileID = fopen(FileName, "r");
		for (PriTimeCount = 0; PriTimeCount < Parms->ResponseTimeTotal; PriTimeCount++) {
			for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
				fscanf(FileID, "%lf %lf %lf %lf\n", &Dummy, &Dummy, &ReadImag, &ReadReal);
				ImpulsiveEntry = IntToArrayEntry(PriTimeCount, ResponseTimeMin) * PriTimeSkip + IntToArrayEntry(WaitTimeCount, PulseTimeMin) * Parms->ResponseTimeTotal + SecTimeCount;
				ResponseTotal[ImpulsiveEntry].GB = (ReadReal + I * ReadImag);
			}
		}
		fclose(FileID);
		
		sprintf(FileName, "%s%d%s", "../../WT_", Parms->WaitTimeTotal + WaitTimeCount, "/ResponseGBRR.dat");
		FileID = fopen(FileName, "r");
		for (PriTimeCount = 0; PriTimeCount < Parms->ResponseTimeTotal; PriTimeCount++) {
			for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
				fscanf(FileID, "%lf %lf %lf %lf\n", &Dummy, &Dummy, &ReadImag, &ReadReal);
				ImpulsiveEntry = IntToArrayEntry(-PriTimeCount, ResponseTimeMin) * PriTimeSkip + IntToArrayEntry(WaitTimeCount, PulseTimeMin) * Parms->ResponseTimeTotal + SecTimeCount;
				ResponseTotal[ImpulsiveEntry].GB = (ReadReal + I * ReadImag);
			}
		}
		fclose(FileID);
		
		sprintf(FileName, "%s%d%s", "../../WT_", Parms->WaitTimeTotal + WaitTimeCount, "/ResponseSENR.dat");
		FileID = fopen(FileName, "r");
		for (PriTimeCount = 0; PriTimeCount < Parms->ResponseTimeTotal; PriTimeCount++) {
			for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
				fscanf(FileID, "%lf %lf %lf %lf\n", &Dummy, &Dummy, &ReadImag, &ReadReal);
				ImpulsiveEntry = IntToArrayEntry(PriTimeCount, ResponseTimeMin) * PriTimeSkip + IntToArrayEntry(WaitTimeCount, PulseTimeMin) * Parms->ResponseTimeTotal + SecTimeCount;
				ResponseTotal[ImpulsiveEntry].SE = (ReadReal + I * ReadImag);
			}
		}
		fclose(FileID);
		
		sprintf(FileName, "%s%d%s", "../../WT_", Parms->WaitTimeTotal + WaitTimeCount, "/ResponseSERR.dat");
		FileID = fopen(FileName, "r");
		for (PriTimeCount = 0; PriTimeCount < Parms->ResponseTimeTotal; PriTimeCount++) {
			for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
				fscanf(FileID, "%lf %lf %lf %lf\n", &Dummy, &Dummy, &ReadImag, &ReadReal);
				ImpulsiveEntry = IntToArrayEntry(-PriTimeCount, ResponseTimeMin) * PriTimeSkip + IntToArrayEntry(WaitTimeCount, PulseTimeMin) * Parms->ResponseTimeTotal + SecTimeCount;
				ResponseTotal[ImpulsiveEntry].SE = (ReadReal + I * ReadImag);
			}
		}
		fclose(FileID);
		
		sprintf(FileName, "%s%d%s", "../../WT_", Parms->WaitTimeTotal + WaitTimeCount, "/ResponseEANR.dat");
		FileID = fopen(FileName, "r");
		for (PriTimeCount = 0; PriTimeCount < Parms->ResponseTimeTotal; PriTimeCount++) {
			for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
				fscanf(FileID, "%lf %lf %lf %lf\n", &Dummy, &Dummy, &ReadImag, &ReadReal);
				ImpulsiveEntry = IntToArrayEntry(PriTimeCount, ResponseTimeMin) * PriTimeSkip + IntToArrayEntry(WaitTimeCount, PulseTimeMin) * Parms->ResponseTimeTotal + SecTimeCount;
				ResponseTotal[ImpulsiveEntry].EA = (ReadReal + I * ReadImag);
			}
		}
		fclose(FileID);
		
		sprintf(FileName, "%s%d%s", "../../WT_", Parms->WaitTimeTotal + WaitTimeCount, "/ResponseEARR.dat");
		FileID = fopen(FileName, "r");
		for (PriTimeCount = 0; PriTimeCount < Parms->ResponseTimeTotal; PriTimeCount++) {
			for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
				fscanf(FileID, "%lf %lf %lf %lf\n", &Dummy, &Dummy, &ReadImag, &ReadReal);
				ImpulsiveEntry = IntToArrayEntry(-PriTimeCount, ResponseTimeMin) * PriTimeSkip + IntToArrayEntry(WaitTimeCount, PulseTimeMin) * Parms->ResponseTimeTotal + SecTimeCount;
				ResponseTotal[ImpulsiveEntry].EA = (ReadReal + I * ReadImag);
			}
		}
		fclose(FileID);
	}
	
	PrintMessage("Finished loading!", "", "", true, false);
	
	for (PriPulseCount = 0; PriPulseCount < Parms->ResponseTimeTotal; PriPulseCount++)
		for (SecPulseCount = 0; SecPulseCount < Parms->ResponseTimeTotal; SecPulseCount++) {
			PulsiveEntry = PriPulseCount * Parms->ResponseTimeTotal + SecPulseCount;
			
			for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
				
				TerPulseOffset = SecPulseCount - SecTimeCount;
				if (abs(TerPulseOffset) > Parms->PulseTimeTotal - 1) continue;
				
				for (WaitTimeCount = PulseTimeMin; WaitTimeCount < Parms->PulseTimeTotal; WaitTimeCount++) {
					
					SecPulseOffset = TerPulseOffset - WaitTimeCount; // (SecPulseCount + Parms->WaitTimeTotal) - (SecTimeCount + Parms-WaitTimeTotal + WaitTimeCount);
					if (abs(SecPulseOffset) + abs(TerPulseOffset) > Parms->PulseTimeTotal - 1) continue;
					
					for (PriTimeCount = ResponseTimeMin; PriTimeCount < Parms->ResponseTimeTotal; PriTimeCount++) {
						
						PriPulseOffset = SecPulseOffset + PriPulseCount - PriTimeCount;
						// (SecPulseCount + Parms->WaitTimeTotal + PriPulseCount) - (SecTimeCount + Parms->WaitTimeTotal + WaitTimeCount + PriTimeCount);
						if (abs(PriPulseOffset) + abs(SecPulseOffset) + abs(TerPulseOffset) > Parms->PulseTimeTotal - 1) continue;
						
						PriPulseEntry = IntToArrayEntry(PriPulseOffset, PulseTimeMin);
						SecPulseEntry = IntToArrayEntry(SecPulseOffset, PulseTimeMin);
						TerPulseEntry = IntToArrayEntry(TerPulseOffset, PulseTimeMin);
						
						ImpulsiveEntry = IntToArrayEntry(-PriTimeCount, ResponseTimeMin) * PriTimeSkip + IntToArrayEntry(WaitTimeCount, PulseTimeMin) * Parms->ResponseTimeTotal + SecTimeCount;
						
						ResponseRR[PulsiveEntry].GB += ResponseTotal[ImpulsiveEntry].GB * conj(Pulse[PriPulseEntry]) * Pulse[SecPulseEntry] * Pulse[TerPulseEntry];
						ResponseRR[PulsiveEntry].SE += ResponseTotal[ImpulsiveEntry].SE * conj(Pulse[PriPulseEntry]) * Pulse[SecPulseEntry] * Pulse[TerPulseEntry];
						ResponseRR[PulsiveEntry].EA += ResponseTotal[ImpulsiveEntry].EA * conj(Pulse[PriPulseEntry]) * Pulse[SecPulseEntry] * Pulse[TerPulseEntry];
						
						ImpulsiveEntry = IntToArrayEntry(PriTimeCount, ResponseTimeMin) * PriTimeSkip + IntToArrayEntry(WaitTimeCount, PulseTimeMin) * Parms->ResponseTimeTotal + SecTimeCount;
						
						ResponseNR[PulsiveEntry].GB += ResponseTotal[ImpulsiveEntry].GB * Pulse[PriPulseEntry] * conj(Pulse[SecPulseEntry]) * Pulse[TerPulseEntry];
						ResponseNR[PulsiveEntry].SE += ResponseTotal[ImpulsiveEntry].SE * Pulse[PriPulseEntry] * conj(Pulse[SecPulseEntry]) * Pulse[TerPulseEntry];
						ResponseNR[PulsiveEntry].EA += ResponseTotal[ImpulsiveEntry].EA * Pulse[PriPulseEntry] * conj(Pulse[SecPulseEntry]) * Pulse[TerPulseEntry];
					}
					
				}
			}
		}

	SaveResponse2D("ResponseGBRR", &(*ResponseRR).GB, Parms, sizeof(ResponseStruct) / sizeof(ResponseRR[0].GB));
	SaveResponse2D("ResponseGBNR", &(*ResponseNR).GB, Parms, sizeof(ResponseStruct) / sizeof(ResponseRR[0].GB));
	SaveResponse2D("ResponseSERR", &(*ResponseRR).SE, Parms, sizeof(ResponseStruct) / sizeof(ResponseRR[0].GB));
	SaveResponse2D("ResponseSENR", &(*ResponseNR).SE, Parms, sizeof(ResponseStruct) / sizeof(ResponseRR[0].GB));
	SaveResponse2D("ResponseEARR", &(*ResponseRR).EA, Parms, sizeof(ResponseStruct) / sizeof(ResponseRR[0].GB));
	SaveResponse2D("ResponseEANR", &(*ResponseNR).EA, Parms, sizeof(ResponseStruct) / sizeof(ResponseRR[0].GB));
	
	PrintMessage("Finished!", "", "", true, false);
	
	free(Pulse);
	free(ResponseNR); free(ResponseRR); free(ResponseTotal);
	free(Parms);
	
	return 0;
}


void SaveResponse2D(const char *FileName, const ComplexType *Response, const ParmStruct *Parms, const int MemSkip)
{
	int PriCount, SecCount;
	char ModFileName[20];
	FILE *FileID;
	
	CreateFileName(ModFileName, FileName, Parms->Embarrassing); FileID = fopen(ModFileName, "w");
	for (PriCount = 0; PriCount < Parms->ResponseTimeTotal; PriCount++)
		for (SecCount = 0; SecCount < Parms->ResponseTimeTotal; SecCount++)
			fprintf(FileID, "%lf %lf %lf %lf\n", PriCount * Parms->TimeStep, SecCount * Parms->TimeStep, cimag(Response[(PriCount * Parms->ResponseTimeTotal + SecCount) * MemSkip]), creal(Response[(PriCount * Parms->ResponseTimeTotal + SecCount) * MemSkip]));
	fclose(FileID);
}