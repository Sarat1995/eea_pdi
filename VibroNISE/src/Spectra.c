#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "ToolsMod.h"
#include "ParmsHandleMod.h"

#include <math.h>
#include <fftw3.h>


static void Generate1DSpec(const ParmStruct *Parms);
static void Generate2DSpec(ComplexType *TotalSpectrum, ComplexType *PriSubSpectrum, ComplexType *SecSubSpectrum, const char *IdChar, const ParmStruct *Parms, const int Rephasing);
static void FourierTransform1D(ComplexType *Output, const ComplexType *Input, const int DimensionIn, const int DimensionOut);
static void FourierTransform2D(ComplexType *Output, const ComplexType *Input, const int DimensionIn, const int DimensionOut);
static void SaveSpectrum(const ComplexType *Spectrum, const char *FileName, const ParmStruct *Parms, const int Rephasing);


int main(int ArgsTotal, char *Args[])
{
	ParmStruct *Parms;
	ComplexType *Spectrum2D, *SpectrumRR, *SpectrumNR, *SpectrumGB, *SpectrumSE, *SpectrumEA;
	
	PrintTitle(false);
	PrintMessage("              // Spectra Generator //", "", "", false, false);
	PrintMessage("", "", "", false, false);
	
	Parms = calloc(1, sizeof(ParmStruct));
	InitializeParms(Parms);
	ReadInputFile(Args, Parms, NULL);
	DerivedParms(Parms);
	
	Parms->SpectrumTotal = 2 * Parms->ResponseTimeTotal;
	
	if (Parms->Mode == ModeResponse1D) Generate1DSpec(Parms);
	else if (Parms->Mode == ModeResponse2D) {
		
		Spectrum2D = (ComplexType *) calloc(IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumRR = (ComplexType *) calloc(IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumNR = (ComplexType *) calloc(IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumGB = (ComplexType *) calloc(IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumSE = (ComplexType *) calloc(IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumEA = (ComplexType *) calloc(IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		
		Generate2DSpec(Spectrum2D, SpectrumRR, SpectrumGB, "GBRR", Parms, true);
		Generate2DSpec(Spectrum2D, SpectrumNR, SpectrumGB, "GBNR", Parms, false);
		Generate2DSpec(Spectrum2D, SpectrumRR, SpectrumSE, "SERR", Parms, true);
		Generate2DSpec(Spectrum2D, SpectrumNR, SpectrumSE, "SENR", Parms, false);
		Generate2DSpec(Spectrum2D, SpectrumRR, SpectrumEA, "EARR", Parms, true);
		Generate2DSpec(Spectrum2D, SpectrumNR, SpectrumEA, "EANR", Parms, false);
		
		SaveSpectrum(Spectrum2D, "Spectrum2D.dat", Parms, false);
		SaveSpectrum(SpectrumRR, "SpectrumRR.dat", Parms, true);
		SaveSpectrum(SpectrumNR, "SpectrumNR.dat", Parms, false);
		SaveSpectrum(SpectrumGB, "SpectrumGB.dat", Parms, false);
		SaveSpectrum(SpectrumSE, "SpectrumSE.dat", Parms, false);
		SaveSpectrum(SpectrumEA, "SpectrumEA.dat", Parms, false);
		
		free(Spectrum2D); free(SpectrumRR); free(SpectrumNR); free(SpectrumGB); free(SpectrumSE); free(SpectrumEA);
		
	}
	
	PrintMessage("Finished!", "", "", true, false);
	
	free(Parms);
	
	return 0;
}


static void Generate1DSpec(const ParmStruct *Parms)
{
	double Dummy, ReadReal, ReadImag;
	RealType Energy;
	int Count;
	FILE *FileID;
	
	ComplexType *Response = malloc(Parms->ResponseTimeTotal * sizeof(ComplexType)), *Spectrum = malloc(Parms->SpectrumTotal * sizeof(ComplexType));
	
	FileID = fopen("Response1D.dat", "r");
	for (Count = 0; Count < Parms->ResponseTimeTotal; Count++) {
		fscanf(FileID, "%lf %lf %lf\n", &Dummy, &ReadReal, &ReadImag);
		Response[Count] = ReadReal + I * ReadImag;
		if (Parms->XLifeTime) Response[Count] *= exp(-Parms->TimeStep * Count / (2 * Parms->XLifeTime));
	}
	fclose(FileID);
	
	FourierTransform1D(Spectrum, Response, Parms->ResponseTimeTotal, Parms->SpectrumTotal);
	
	NormalizeCArray(Spectrum, Parms->SpectrumTotal);
	
	FileID = fopen("Spectrum1D.dat", "w");
	for (Count = 0; Count < Parms->SpectrumTotal; Count++) {
		Energy = (Count - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->XEnergyMean;
		if (Energy > Parms->EnergyMin && Energy < Parms->EnergyMax) fprintf(FileID, "%f %f\n", Energy, creal(Spectrum[Count]));
	}
	fclose(FileID);
	
	free(Response); free(Spectrum);
}


static void Generate2DSpec(ComplexType *TotalSpectrum, ComplexType *PriSubSpectrum, ComplexType *SecSubSpectrum, const char *IdChar, const ParmStruct *Parms, const int Rephasing)
{
	ComplexType *Response = malloc(IntPow2(Parms->ResponseTimeTotal) * sizeof(ComplexType)), *Spectrum = malloc(IntPow2(Parms->SpectrumTotal) * sizeof(ComplexType));
	ComplexType SpectralValue;
	double Dummy, ReadReal, ReadImag;
	RealType PriEnergy, SecEnergy;
	int PriCount, SecCount;
	char FileName[100];
	FILE *FileID;
	
	sprintf(FileName, "%s%s%s", "Response", IdChar, ".dat");
	
	FileID = fopen(FileName, "r");
	for (PriCount = 0; PriCount < Parms->ResponseTimeTotal; PriCount++)
		for (SecCount = 0; SecCount < Parms->ResponseTimeTotal; SecCount++) {
			fscanf(FileID, "%lf %lf %lf %lf\n", &Dummy, &Dummy, &ReadImag, &ReadReal);
			Response[PriCount * Parms->ResponseTimeTotal + SecCount] = (-ReadReal - I * ReadImag);
			if (Parms->XLifeTime) Response[PriCount * Parms->ResponseTimeTotal + SecCount] *= exp(-Parms->TimeStep * (PriCount + SecCount) / (2 * Parms->XLifeTime));
		}
	fclose(FileID);
	
	FourierTransform2D(Spectrum, Response, Parms->ResponseTimeTotal, Parms->SpectrumTotal);
	
	sprintf(FileName, "%s%s%s", "Spectrum", IdChar, ".dat");
	FileID = fopen(FileName, "w");
	for (PriCount = 0; PriCount < Parms->SpectrumTotal; PriCount++) {
		
		PriEnergy = (PriCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->XEnergyMean;
		if (PriEnergy < Parms->EnergyMin || PriEnergy > Parms->EnergyMax) continue;
		
		for (SecCount = 0; SecCount < Parms->SpectrumTotal; SecCount++) {
			SecEnergy = (SecCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->XEnergyMean;
			if (SecEnergy < Parms->EnergyMin || SecEnergy > Parms->EnergyMax) continue;
			
			SpectralValue = Spectrum[(Rephasing * Parms->SpectrumTotal + Minus1Pow(Rephasing) * PriCount) * Parms->SpectrumTotal + SecCount];
			TotalSpectrum[PriCount * Parms->SpectrumTotal + SecCount] += SpectralValue;
			PriSubSpectrum[PriCount * Parms->SpectrumTotal + SecCount] += SpectralValue;
			SecSubSpectrum[PriCount * Parms->SpectrumTotal + SecCount] += SpectralValue;
			fprintf(FileID, "%f %f %f %f\n", Minus1Pow(Rephasing) * PriEnergy, SecEnergy, -cimag(SpectralValue), -creal(SpectralValue));
		}
	}
	fclose(FileID);
	
	free(Response); free(Spectrum);
}


static void FourierTransform1D(ComplexType *Output, const ComplexType *Input, const int DimensionIn, const int DimensionOut)
{
	fftw_complex *In, *Out;
	fftw_plan Plan;
	int Count;
	
	In = fftw_malloc(sizeof(fftw_complex) * DimensionOut);
	Out = fftw_malloc(sizeof(fftw_complex) * DimensionOut);
	Plan = fftw_plan_dft_1d(DimensionOut, In, Out, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	for (Count = 0; Count < DimensionOut; Count++) In[Count] = 0;
	for (Count = 0; Count < DimensionIn; Count++) {
		In[Count] = Input[Count];
		if (Count > 0) In[DimensionOut - Count] = conj(In[Count]);
	}
	
	fftw_execute(Plan);
	
	for (Count = 0; Count < DimensionOut; Count++)
		Output[Count] = Out[(Count < (int) ceil(.5 * DimensionOut)) * (Count + (int) floor(.5 * DimensionOut)) + (Count > (int) floor(.5 * DimensionOut)) * (Count - (int) ceil(.5 * DimensionOut))];
	
	fftw_destroy_plan(Plan);
	fftw_free(In); fftw_free(Out);
}


static void FourierTransform2D(ComplexType *Output, const ComplexType *Input, const int DimensionIn, const int DimensionOut)
{
	fftw_complex *In, *Out;
	fftw_plan Plan;
	int PriCount, SecCount;
	
	In = fftw_malloc(sizeof(fftw_complex) * IntPow2(DimensionOut));
	Out = fftw_malloc(sizeof(fftw_complex) * IntPow2(DimensionOut));
	Plan = fftw_plan_dft_2d(DimensionOut, DimensionOut, In, Out, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	for (PriCount = 0; PriCount < IntPow2(DimensionOut); PriCount++) In[PriCount] = 0;
	for (PriCount = 0; PriCount < DimensionIn; PriCount++)
		for (SecCount = 0; SecCount < DimensionIn; SecCount++) In[PriCount * DimensionOut + SecCount] = Input[PriCount * DimensionIn + SecCount];
	
	fftw_execute(Plan);
	
	for (PriCount = 0; PriCount < DimensionOut; PriCount++)
		for (SecCount = 0; SecCount < DimensionOut; SecCount++)
			Output[PriCount * DimensionOut + SecCount] = Out[((PriCount < (int) ceil(.5 * DimensionOut)) * (PriCount + (int) floor(.5 * DimensionOut)) + (PriCount > (int) floor(.5 * DimensionOut)) * (PriCount - (int) ceil(.5 * DimensionOut))) * DimensionOut + (SecCount < (int) ceil(.5 * DimensionOut)) * (SecCount + (int) floor(.5 * DimensionOut)) + (SecCount > (int) floor(.5 * DimensionOut)) * (SecCount - (int) ceil(.5 * DimensionOut))];
	
	fftw_destroy_plan(Plan);
	fftw_free(In); fftw_free(Out);
}


static void SaveSpectrum(const ComplexType *Spectrum, const char *FileName, const ParmStruct *Parms, const int Rephasing)
{
	FILE *FileID;
	RealType PriEnergy, SecEnergy;
	int PriCount, SecCount;
	
	FileID = fopen(FileName, "w");
	for (PriCount = 0; PriCount < Parms->SpectrumTotal; PriCount++) {
		PriEnergy = (PriCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->XEnergyMean;
		if (PriEnergy < Parms->EnergyMin || PriEnergy > Parms->EnergyMax) continue;
		for (SecCount = 0; SecCount < Parms->SpectrumTotal; SecCount++) {
			SecEnergy = (SecCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->XEnergyMean;				
			if (SecEnergy < Parms->EnergyMin || SecEnergy > Parms->EnergyMax) continue;
			fprintf(FileID, "%f %f %f %f\n", Minus1Pow(Rephasing) * PriEnergy, SecEnergy, -cimag(Spectrum[PriCount * Parms->SpectrumTotal + SecCount]), -creal(Spectrum[PriCount * Parms->SpectrumTotal + SecCount]));
		}
	}
	fclose(FileID);
}