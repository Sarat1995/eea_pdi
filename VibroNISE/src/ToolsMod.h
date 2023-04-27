#ifndef ToolsMod
#define ToolsMod

void InitializeLog(const int Embarrassing);
void PrintTitle(const int Embarrassing);
void PrintProgress(int *Percent, int *Benchmark, const int Total, const int Embarrassing);
void CreateFileName(char *FileName, const char *Name, const int Embarrassing);
void PrintMessage(const char *Left, const char *Middle, const char *Right, const int DisplayTime, const int Embarrassing);
void CumulateRArray(RealType *Array, const int Dimension);
void FillRArray(RealType *Array, const RealType Value, const int Dimension);
void FillCArray(ComplexType *Array, const ComplexType Value, const int Dimension);
void NormalizeDArray(double *Array, const int Dimension);
void NormalizeCArray(ComplexType *Array, const int Dimension);
void CreateMatrixMask(int *Mask, const RealType *Matrix, const int Dimension);
void RealPlusRArray(RealType *Array, const RealType Addition, int Dimension);
void FactorCArrays(ComplexType *Array, const ComplexType *FactorArray, const int Dimension);
void RealTimesRArray(RealType *Array, const RealType Factor, int Dimension);
void RealTimesDArray(double *Array, const RealType Factor, int Dimension);
void RealTimesCArray(ComplexType *Array, const RealType Factor, int Dimension);
void RealTimesDCArray(double complex *Array, const RealType Factor, int Dimension);
void RMatrixTimesCVector(ComplexType *Product, const RealType *Matrix, const ComplexType *Vector, const int Dimension);
void RMatrixTimesCVectorMasked(ComplexType *Product, const RealType *Matrix, const ComplexType *Vector, const int *Mask, const int Dimension);
void CMatrixTimesCVector(ComplexType *Product, const ComplexType *Matrix, const ComplexType *Vector, const int Dimension);
void Conjugate(ComplexType *Conjugate, const ComplexType *Array, const int Dimension);
void RealToComplex(ComplexType *Output, const RealType *Input, const int Dimension);
void IntToString(char *Output, const int Input);
void RealToString(char *Output, const RealType Input);
void CreateGrid(RealType *Grid, const RealType Min, const RealType Max, const int Steps);
void PrintIArray(const int *Array, const int RowTotal, const int ColumnTotal);
void PrintRArray(const RealType *Array, const int RowTotal, const int ColumnTotal);
void PrintCArray(const ComplexType *Array, const int RowTotal, const int ColumnTotal);
void SaveRArray(const char *FileName, const RealType *Array, const int RowTotal, const int ColumnTotal, const RealType *RowIdVector, const int Embarrassing);
void SaveCArray(const char *FileName, const ComplexType *Array, const int RowTotal, const int ColumnTotal, const RealType *RowIdVector, const int Embarrassing);
void SaveDVector(const char *FileName, const double *Vector, const RealType *IdVector, const int Dimension, const int Embarrassing);
void SaveCVector(const char *FileName, const ComplexType *Vector, const RealType *IdVector, const int Dimension, const int Embarrassing);
void SaveDCVector(const char *FileName, const double complex *Vector, const RealType *IdVector, const int Dimension, const int Embarrassing);
void InitializeDiag(RealType *Matrix, RealType *EigenValues, int Dimension, const int DAndC);
void DiagonalizeMatrix(RealType *Matrix, RealType *EigenValues, int Dimension, const int DAndC);
void CloseDiag(const int DAndC);
void InvertMatrix(RealType *Matrix, int Dimension);

long int SetTimer(void);
int Minus1Pow(const int Input);
int IntPow2(const int Input);
int IntPowN(const int Input, const int N);
RealType RealPow2(const RealType Input);
RealType RealPowN(const RealType Input, const int N);
RealType RealSqrt(const RealType Input);
RealType RealExp(const RealType Input);
int IntMin(const int PriInput, const int SecInput);
int IntMax(const int PriInput, const int SecInput);
int Separation(const int PriSite, const int SecSite, const int Length, const int Periodic);
int SignedSeparation(const int PriSite, const int SecSite, const int Length, const int Periodic);
int IntToArrayEntry(const int Int, const int MinValue);
RealType Lorentzian(const RealType Value, const RealType Width);
RealType Gaussian(const RealType Value, const RealType StandardDeviation);
int SumIArray(const int *Array, const int Dimension);
ComplexType Trace(const ComplexType *Matrix, const int Dimension);
RealType RInnerProduct(const RealType *PriInput, const RealType *SecInput, const int Dimension);
ComplexType CInnerProduct(const ComplexType *Bra, const ComplexType *Ket, const int Dimension);
ComplexType RCInnerProduct(const RealType *Bra, const ComplexType *Ket, const int Dimension);

#endif
