#include "GlobalsMod.h"
#include "OverlapsMod.h"

#include <math.h>


static void SaveOverlaps(const char *FileName, const RealType *Overlaps, const int Dimension);

static RealType Overlap(const int Shifted, const int Unshifted, const RealType Coupling, const RealType CouplingSquared);
static RealType Laguerre(const int K, const RealType X);
static RealType Function1(const int N, const int M);
static RealType Function2(const int N, const int M, const RealType X);
static long int Factorial(const int LowBound, const int UpBound);


/////////////
// Globals //
/////////////


void BuildOverlaps(const char *FileName, RealType **Overlaps, const RealType Coupling, const int Dimension, const int Embarrassing)
{
	RealType CouplingSquared = Coupling * Coupling;
	int Shifted, Unshifted;
	
	*Overlaps = malloc(Dimension * Dimension * sizeof(RealType));
	
	for (Shifted = 0; Shifted < Dimension; Shifted++)
		for (Unshifted = 0; Unshifted < Dimension; Unshifted++)
			(*Overlaps)[Shifted * Dimension + Unshifted] = Overlap(Shifted, Unshifted, Coupling, CouplingSquared);
	
	if (!Embarrassing) SaveOverlaps(FileName, *Overlaps, Dimension);
}


////////////
// Locals //
////////////


static RealType Overlap(const int Shifted, const int Unshifted, const RealType Coupling, const RealType CouplingSquared)
{
	RealType Output = exp(-CouplingSquared / 2);
	
	if (Shifted == Unshifted)
		Output *= Laguerre(Unshifted, CouplingSquared);
	else if (Shifted > Unshifted)
		Output *= pow(Coupling, Shifted - Unshifted) * Function1(Unshifted, Shifted) * Function2(Unshifted, Shifted - Unshifted, CouplingSquared);
	else
		Output *= pow(-Coupling, Unshifted - Shifted) * Function1(Shifted, Unshifted) * Function2(Shifted, Unshifted - Shifted, CouplingSquared);
	
	return Output;
}


static RealType Laguerre(const int K, const RealType X)
{
	RealType Output, Dummy1, Dummy2;
	int Count;
	
	if (K == 0)			Output = 1;
	else if (K == 1)	Output = 1 - X;
	else {
		Dummy2 = 1;
		Dummy1 = 1 - X;
		for (Count = 2; Count <= K; Count++) {
			Output = ((2. * Count - 1 - X) * Dummy1 - (Count - 1) * Dummy2) / Count;
			Dummy2 = Dummy1;
			Dummy1 = Output;
		}		
	}
	
	return Output;
}


static RealType Function1(const int N, const int M)
{
	RealType Output;
	
	if (N >= (M - N))	Output = sqrt(Factorial(N + 1, M)) / Factorial(1, M - N);
	else if (N > 0)		Output = sqrt((RealType) Factorial(M - N + 1, M) / Factorial(1, N) / Factorial(1, M - N));
	else				Output = 1 / sqrt(Factorial(1, M));
	
	return Output;
}


static RealType Function2(const int N, const int M, const RealType X)
{
	RealType Output = 1;
	int Count;
	
	for (Count = N; Count >= 1; Count--) Output = 1 - (RealType) (N - Count + 1) / (Count * (M + Count)) * X * Output;
	
	return Output;
}


static long int Factorial(const int LowBound, const int UpBound)
{
	long int Output = LowBound;
	int Count;
	
	for (Count = LowBound + 1; Count <= UpBound; Count++) Output *= Count;
	
	return Output;
}


static void SaveOverlaps(const char *FileName, const RealType *Overlaps, const int Dimension)
{
	char ModFileName[20];
	int RowCount, ColumnCount;
	FILE *FileID;
	
	sprintf(ModFileName, "%s%s", FileName, ".dat");
	
	FileID = fopen(ModFileName, "w");
	for (RowCount = 0; RowCount < Dimension; RowCount++) {
		for (ColumnCount = 0; ColumnCount < Dimension; ColumnCount++) fprintf(FileID, "%f, ", Overlaps[RowCount * Dimension + ColumnCount]);
		fprintf(FileID, "\n");
	}
	fclose(FileID);
}