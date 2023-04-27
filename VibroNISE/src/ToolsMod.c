#include "GlobalsMod.h"
#include "ToolsMod.h"

#include <time.h>
#include <math.h>


#ifdef DoublePrec
extern void dsyev_(char *Jobz, char *UpLo, int *N, double *A, int *LDA, double *W, double *Work, int *LWork, int *Info);
extern void dsyevd_(char *Jobz, char *UpLo, int *N, double *A, int *LDA, double *W, double *Work, int *LWork, int *IWork, int *LIWork, int *Info);
extern void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *Info);
extern void dgetri_(int *N, double *A, int *LDA, int *IPIV, double *Work, int *LWork, int *Info);
#else
extern void ssyev_(char *Jobz, char *UpLo, int *N, float *A, int *LDA, float *W, float *Work, int *LWork, int *Info);
extern void ssyevd_(char *Jobz, char *UpLo, int *N, float *A, int *LDA, float *W, float *Work, int *LWork, int *IWork, int *LIWork, int *Info);
extern void sgetrf_(int *M, int *N, float *A, int *LDA, int *IPIV, int *Info);
extern void sgetri_(int *N, float *A, int *LDA, int *IPIV, float *Work, int *LWork, int *Info);
#endif

RealType *Work;
int *IWork;
int LIWork, LWork;



void InitializeLog(const int Embarrassing)
{
	FILE *FileID;
	char FileName[20];
	
	CreateFileName(FileName, "Log", Embarrassing); FileID = fopen(FileName, "w");
	fclose(FileID);
}


void PrintTitle(const int Embarrassing)
{
	
#ifdef DoublePrec
	char PrecisionMode[] = "double";
#else
	char PrecisionMode[] = "single";
#endif
	
	PrintMessage("", "", "", false, Embarrassing);
	PrintMessage("     // // // ////  ////   ///  /  // //  //// ////           ", "", "", false, Embarrassing);
	PrintMessage("    // // // // // // // // // // // // //    //              ", "", "", false, Embarrassing);
	PrintMessage("   // // // ////  ////  // // ///// //  ///  ///              ", "", "", false, Embarrassing);
	PrintMessage("   ///  // // // // // // // // // //    // //                ", "", "", false, Embarrassing);
	PrintMessage("   /   // ////  //  // ///  //  / // ////  ////  Version 3.02 ", "", "", false, Embarrassing);
	PrintMessage("", "", "", false, Embarrassing);
	PrintMessage("   Running in ", PrecisionMode, " precision mode.", false, Embarrassing);
	PrintMessage("", "", "", false, Embarrassing);
}


void PrintProgress(int *Percent, int *Benchmark, const int Total, const int Embarrassing)
{
	char DisplayNumber[3];
	int NewBenchmark;
	
	IntToString(DisplayNumber, *Percent);
	PrintMessage("Progress ", DisplayNumber, "%.", true, Embarrassing);
	
	do {
		(*Percent)++;
		NewBenchmark = Total * *Percent / 100;
	} while (NewBenchmark == *Benchmark);
	*Benchmark = NewBenchmark;
}


void CreateFileName(char *FileName, const char *Name, const int Embarrassing)
{
	sprintf(FileName, "%s", Name);
	if (Embarrassing) sprintf(FileName, "%s_%d", FileName, Embarrassing);
	sprintf(FileName, "%s%s", FileName, ".dat");
}


void PrintMessage(const char *Left, const char *Middle, const char *Right, const int DisplayTime, const int Embarrassing)
{
	time_t TimeRawFormat;
	struct tm *GetTime;
	char Time[25];
	time (&TimeRawFormat);
	char FileName[20];
	FILE *FileID;
	
	CreateFileName(FileName, "Log", Embarrassing); FileID = fopen(FileName, "a");
	
	if (DisplayTime) {
		GetTime = localtime(&TimeRawFormat);
		strftime(Time, 25, "%Y:%m:%d %H:%M:%S", GetTime);
		printf("%s - ", Time);
		fprintf(FileID, "%s - ", Time);
	}
	
	printf("%s%s%s\n", Left, Middle, Right);
	fprintf(FileID, "%s%s%s\n", Left, Middle, Right);
	fclose(FileID);
}


long int SetTimer(void)
{
	time_t TimeRawFormat;
	struct tm *GetTime;
	char ContrChar[25];
	time (&TimeRawFormat);
	int ContrInt;
	long int Timer;
	
	GetTime = localtime(&TimeRawFormat);
	
	strftime(ContrChar, 25, "%S", GetTime);
	sscanf(ContrChar, "%d", &ContrInt);
	Timer = ContrInt;
	
	strftime(ContrChar, 25, "%M", GetTime);
	sscanf(ContrChar, "%d", &ContrInt);
	Timer += ContrInt * 60;
	
	strftime(ContrChar, 25, "%H", GetTime);
	sscanf(ContrChar, "%d", &ContrInt);
	Timer += ContrInt * 3600;
	
	strftime(ContrChar, 25, "%d", GetTime);
	sscanf(ContrChar, "%d", &ContrInt);
	Timer += ContrInt * 3600 * 24;
	
	return Timer;
}


int Minus1Pow(const int Input)
{
	return Input % 2 == 0 ? 1 : -1;
}


int IntPow2(const int Input)
{
	return Input * Input;
}


int IntPowN(const int Input, const int N)
{
	return pow(Input, N);
}


RealType RealPow2(const RealType Input)
{
	return Input * Input;
}


RealType RealPowN(const RealType Input, const int N)
{
	return pow(Input, N);
}


RealType RealSqrt(const RealType Input)
{
	return sqrt(Input);
}


RealType RealExp(const RealType Input)
{
	return exp(Input);
}


int IntMin(const int PriInput, const int SecInput)
{
	return (PriInput < SecInput) ? PriInput : SecInput;
}


int IntMax(const int PriInput, const int SecInput)
{
	return (PriInput > SecInput) ? PriInput : SecInput;
}


int Separation(const int PriSite, const int SecSite, const int Length, const int Periodic)
{
	int Distance = abs(PriSite - SecSite);
	
	return Periodic ? IntMin(Length - Distance, Distance) : Distance;
}


int SignedSeparation(const int PriSite, const int SecSite, const int Length, const int Periodic)
{
	int Distance = abs(PriSite - SecSite);
	
	if( Distance == 0 )
		return 0;
	
    if( Periodic )
    {
        int length_minus_distance = Length - Distance;
        
        if( length_minus_distance < Distance )
        {
            if( SecSite > PriSite )
                return -length_minus_distance;
            else
                return length_minus_distance;
        }
        else
        {
            if( SecSite > PriSite )
                return Distance;
            else
                return -Distance;
        }
    }
    else
    {
        if( SecSite > PriSite )
            return Distance;
        else
            return -Distance;
    }
}


int IntToArrayEntry(const int Int, const int MinValue)
{
	return Int - MinValue;
}


RealType Lorentzian(const RealType Value, const RealType Width)
{
	return 1 / (RealPow2(Value) + RealPow2(Width));
}


RealType Gaussian(const RealType Value, const RealType StandardDeviation)
{
	return exp(-.5 * RealPow2(Value) / RealPow2(StandardDeviation));
}


int SumIArray(const int *Array, const int Dimension)
{
	int Count, Value = 0;
	const int *ArrayPos = Array;
	
	for (Count = 0; Count < Dimension; Count++) Value += (*ArrayPos++);
	
	return Value;
}


void CumulateRArray(RealType *Array, const int Dimension)
{
	int Count;
	
	for (Count = 1; Count < Dimension; Count++) Array[Count] += Array[Count - 1];
}


ComplexType Trace(const ComplexType *Matrix, const int Dimension)
{
	int Count;
	ComplexType Value = 0;
	
	for (Count = 0; Count < Dimension; Count++) Value += Matrix[Count * (Dimension + 1)];
	
	return Value;
}


void FillRArray(RealType *Array, const RealType Value, const int Dimension)
{
	RealType *ArrayElement = &Array[0];
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) = Value;
}


void FillCArray(ComplexType *Array, const ComplexType Value, const int Dimension)
{
	ComplexType *ArrayElement = &Array[0];
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) = Value;
}


void NormalizeDArray(double *Array, const int Dimension)
{
	RealType Maximum = 0;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++)
		if (fabs(Array[Count]) > Maximum) Maximum = fabs(Array[Count]);
	
	if (Maximum > 0) for (Count = 0; Count < Dimension; Count++) Array[Count] /= Maximum;
}


void NormalizeCArray(ComplexType *Array, const int Dimension)
{
	RealType Maximum = 0;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++)
		if (cabs(Array[Count]) > Maximum) Maximum = cabs(Array[Count]);
	
	if (Maximum > 0) for (Count = 0; Count < Dimension; Count++) Array[Count] /= Maximum;
}


void CreateMatrixMask(int *Mask, const RealType *Matrix, const int Dimension)
{
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) Mask[Count] = -1 * (!Matrix[Count]);
}


void RealPlusRArray(RealType *Array, const RealType Addition, int Dimension)
{
	RealType *ArrayPos = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayPos++) += Addition;
}


void FactorCArrays(ComplexType *Array, const ComplexType *FactorArray, const int Dimension)
{
	const ComplexType *FactorArrayElement = FactorArray;
	ComplexType *ArrayElement = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) *= (*FactorArrayElement++);
}


void RealTimesRArray(RealType *Array, const RealType Factor, int Dimension)
{
	RealType *ArrayPos = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayPos++) *= Factor;
}


void RealTimesDArray(double *Array, const RealType Factor, int Dimension)
{
	double *ArrayPos = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayPos++) *= Factor;
}


void RealTimesCArray(ComplexType *Array, const RealType Factor, int Dimension)
{
	ComplexType *ArrayPos = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayPos++) *= Factor;
}


void RealTimesDCArray(double complex *Array, const RealType Factor, int Dimension)
{
	double complex *ArrayPos = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayPos++) *= Factor;
}


RealType RInnerProduct(const RealType *PriInput, const RealType *SecInput, const int Dimension)
{
	const RealType *PriElement = PriInput, *SecElement = SecInput;
	RealType Value = 0;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) Value += (*PriElement++) * (*SecElement++);
	return Value;
}


ComplexType CInnerProduct(const ComplexType *Bra, const ComplexType *Ket, const int Dimension)
{
	const ComplexType *BraEntry = &Bra[0], *KetEntry = &Ket[0];
	ComplexType Value = 0;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) Value += conj(*BraEntry++) * (*KetEntry++);
	return Value;
}


ComplexType RCInnerProduct(const RealType *Bra, const ComplexType *Ket, const int Dimension)
{
	const RealType *BraEntry = &Bra[0];
	const ComplexType *KetEntry = &Ket[0];
	ComplexType Value = 0;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) Value += (*BraEntry++) * (*KetEntry++);
	return Value;
}


void RMatrixTimesCVector(ComplexType *Product, const RealType *Matrix, const ComplexType *Vector, const int Dimension)
{
	const RealType *MatrixElement = Matrix;
	const ComplexType *VectorElement;
	ComplexType *ProductElement = Product;
	ComplexType Contribution;
	int PriCount, SecCount;
	
	for (PriCount = 0; PriCount < Dimension; PriCount++) {
		Contribution = 0;
		VectorElement = Vector;
		for(SecCount = 0; SecCount < Dimension; SecCount++) Contribution += (*MatrixElement++) * (*VectorElement++);
		(*ProductElement++) = Contribution;
	}
}


void RMatrixTimesCVectorMasked(ComplexType *Product, const RealType *Matrix, const ComplexType *Vector, const int *Mask, const int Dimension)
{
	int PriCount, SmartCount, SecCount;
	ComplexType Contribution;
	
	for (PriCount = 0; PriCount < Dimension; PriCount++) {
		Contribution = 0;
		SmartCount = PriCount * Dimension;
		for (SecCount = 0; SecCount < Dimension; SecCount++)
			if (Mask[SmartCount + SecCount] != -1) Contribution += Matrix[SmartCount + SecCount] * Vector[SecCount];
		Product[PriCount] = Contribution;
	}
}


void CMatrixTimesCVector(ComplexType *Product, const ComplexType *Matrix, const ComplexType *Vector, const int Dimension)
{
	const ComplexType *MatrixElement = Matrix, *VectorElement;
	ComplexType *ProductElement = Product;
	ComplexType Contribution;
	int PriCount, SecCount;
	
	for (PriCount = 0; PriCount < Dimension; PriCount++) {
		Contribution = 0;
		// Can be modified to "= Vector;"
		VectorElement = &Vector[0];
		for(SecCount = 0; SecCount < Dimension; SecCount++) Contribution += (*MatrixElement++) * (*VectorElement++);
		(*ProductElement++) = Contribution;
	}
}


void Conjugate(ComplexType *Conjugate, const ComplexType *Array, const int Dimension)
{
	ComplexType *ConjugateElement = Conjugate;
	int PriCount, SecCount;
	
	for (PriCount = 0; PriCount < Dimension; PriCount++)
		for (SecCount = 0; SecCount < Dimension; SecCount++)
			(*ConjugateElement++) = conj(Array[SecCount * Dimension + PriCount]);
}


void RealToComplex(ComplexType *Output, const RealType *Input, const int Dimension)
{
	ComplexType *OutputElement = Output;
	const RealType *InputElement = Input;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*OutputElement++) = (*InputElement++);
}


void IntToString(char *Output, const int Input)
{
	sprintf(Output, "%d", Input);
}


void RealToString(char *Output, const RealType Input)
{
	sprintf(Output, "%7.4E", Input);
}


void CreateGrid(RealType *Grid, const RealType Min, const RealType Max, const int Steps)
{
	int Count;
	
	for (Count = 0; Count < Steps; Count++) Grid[Count] = Min + (Max - Min) * (RealType) Count / (Steps - 1);
}


void PrintIArray(const int *Array, const int RowTotal, const int ColumnTotal)
{
	int RowCount, ColumnCount;
	
	for (RowCount = 0; RowCount < RowTotal; RowCount++) {
		for (ColumnCount = 0; ColumnCount < ColumnTotal; ColumnCount++) printf("%d, ", Array[RowCount * ColumnTotal + ColumnCount]);
		printf("\n");
	}
}


void PrintRArray(const RealType *Array, const int RowTotal, const int ColumnTotal)
{
	int RowCount, ColumnCount;
	
	for (RowCount = 0; RowCount < RowTotal; RowCount++) {
		for (ColumnCount = 0; ColumnCount < ColumnTotal; ColumnCount++) printf("%f, ", Array[RowCount * ColumnTotal + ColumnCount]);
		printf("\n");
	}
}


void PrintCArray(const ComplexType *Array, const int RowTotal, const int ColumnTotal)
{
	int RowCount, ColumnCount;
	
	for (RowCount = 0; RowCount < RowTotal; RowCount++) {
		for (ColumnCount = 0; ColumnCount < ColumnTotal; ColumnCount++) printf("%f + %fi, ", creal(Array[RowCount * ColumnTotal + ColumnCount]), cimag(Array[RowCount * ColumnTotal + ColumnCount]));
		printf("\n");
	}
}


void SaveRArray(const char *FileName, const RealType *Array, const int RowTotal, const int ColumnTotal, const RealType *RowIdVector, const int Embarrassing)
{
	int RowCount, ColumnCount;
	char ModFileName[20];
	FILE *FileID;
	
	CreateFileName(ModFileName, FileName, Embarrassing); FileID = fopen(ModFileName, "w");
	for (RowCount = 0; RowCount < RowTotal; RowCount++) {
		if (RowIdVector != NULL) fprintf(FileID, "%f, ", RowIdVector[RowCount]);
		for (ColumnCount = 0; ColumnCount < ColumnTotal; ColumnCount++) fprintf(FileID, "%f, ", Array[RowCount * ColumnTotal + ColumnCount]);
		fprintf(FileID, "\n");
	}
	fclose(FileID);
}


void SaveCArray(const char *FileName, const ComplexType *Array, const int RowTotal, const int ColumnTotal, const RealType *RowIdVector, const int Embarrassing)
{
	int RowCount, ColumnCount;
	char ModFileName[20];
	FILE *FileID;
	
	CreateFileName(ModFileName, FileName, Embarrassing); FileID = fopen(ModFileName, "w");
	for (RowCount = 0; RowCount < RowTotal; RowCount++) {
		if (RowIdVector != NULL) fprintf(FileID, "%f, ", RowIdVector[RowCount]);
		for (ColumnCount = 0; ColumnCount < ColumnTotal; ColumnCount++) fprintf(FileID, "%f, %f, ", creal(Array[RowCount * ColumnTotal + ColumnCount]), cimag(Array[RowCount * ColumnTotal + ColumnCount]));
		fprintf(FileID, "\n");
	}
	fclose(FileID);
}


void SaveDVector(const char *FileName, const double *Vector, const RealType *IdVector, const int Dimension, const int Embarrassing)
{
	int Count;
	char ModFileName[30];
	FILE *FileID;
	
	CreateFileName(ModFileName, FileName, Embarrassing); FileID = fopen(ModFileName, "w");
	for (Count = 0; Count < Dimension; Count++) fprintf(FileID, "%lf %lf\n", IdVector[Count], Vector[Count]);
	fclose(FileID);
}


void SaveCVector(const char *FileName, const ComplexType *Vector, const RealType *IdVector, const int Dimension, const int Embarrassing)
{
	int Count;
	char ModFileName[30];
	FILE *FileID;
	
	CreateFileName(ModFileName, FileName, Embarrassing); FileID = fopen(ModFileName, "w");
	for (Count = 0; Count < Dimension; Count++) fprintf(FileID, "%lf %lf %lf\n", IdVector[Count], creal(Vector[Count]), cimag(Vector[Count]));
	fclose(FileID);
}


void SaveDCVector(const char *FileName, const double complex *Vector, const RealType *IdVector, const int Dimension, const int Embarrassing)
{
	int Count;
	char ModFileName[30];
	FILE *FileID;
	
	CreateFileName(ModFileName, FileName, Embarrassing); FileID = fopen(ModFileName, "w");
	for (Count = 0; Count < Dimension; Count++) fprintf(FileID, "%lf %lf %lf\n", IdVector[Count], creal(Vector[Count]), cimag(Vector[Count]));
	fclose(FileID);
}


void InitializeDiag(RealType *Matrix, RealType *EigenValues, int Dimension, const int DAndC)
{
	int Info;
	FILE *FileID;
	
	LWork = -1;
	Work = calloc(1, sizeof(RealType));
	if (DAndC) {
		LIWork = -1;
		IWork = calloc(1, sizeof(int));
#ifdef DoublePrec
		dsyevd_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, IWork, &LIWork, &Info);			
#else
		ssyevd_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, IWork, &LIWork, &Info);
#endif
		LIWork = IWork[0];
		free(IWork);
		IWork = calloc(LIWork, sizeof(int));
	}
	else {
#ifdef DoublePrec
		dsyev_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, &Info);
#else
		ssyev_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, &Info);
#endif
	}
	LWork = Work[0];
	free(Work);
	if (LWork < 1 + 6 * Dimension + 2 * IntPow2(Dimension)) LWork = 1 + 6 * Dimension + 2 * IntPow2(Dimension);
	Work = calloc(LWork, sizeof(RealType));
}


void DiagonalizeMatrix(RealType *Matrix, RealType *EigenValues, int Dimension, const int DAndC)
{
	int Info;
	
	if (Dimension == 0) return;
	
	if (DAndC) {
#ifdef DoublePrec
		dsyevd_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, IWork, &LIWork, &Info);
#else
		ssyevd_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, IWork, &LIWork, &Info);
#endif
	}
	else {
#ifdef DoublePrec
		dsyev_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, &Info);
#else
		ssyev_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, &Info);
#endif
	}
}


void CloseDiag(const int DAndC)
{
	free(Work);
	if (DAndC) free(IWork);
}


void InvertMatrix(RealType *Matrix, int Dimension)
{
	int Info;
	int *PivotArray = calloc(Dimension, sizeof(int));
	RealType *Work = calloc(IntPow2(Dimension), sizeof(RealType));
	
#ifdef DoublePrec
	dgetrf_(&Dimension, &Dimension, Matrix, &Dimension, PivotArray, &Info);
	dgetri_(&Dimension, Matrix, &Dimension, PivotArray, Work, &Dimension, &Info);
#else
	sgetrf_(&Dimension, &Dimension, Matrix, &Dimension, PivotArray, &Info);
	sgetri_(&Dimension, Matrix, &Dimension, PivotArray, Work, &Dimension, &Info);
#endif
	
	free(Work); free(PivotArray);
}
