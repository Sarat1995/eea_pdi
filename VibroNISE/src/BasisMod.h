#ifndef BasisMod
#define BasisMod

#define XTotal 2
#define ETotal 2
#define HTotal ETotal
#define GTotal 3


typedef struct {
	
	int U;
	int VX, VE, VH, VG;
	
} ParticleStruct;


typedef struct {
	
	ParticleStruct X[XTotal], E[ETotal], H[HTotal], G[GTotal];
	
} BasisStruct;


typedef struct {
	
	RealType *XG, *XE, *XH, *EG, *HG;
	
} VOverlapStruct;


typedef struct {
	
	int Ket, Bra;
	
} MatrixBasisStruct;


void BuildBases(BasisStruct **GBasis, BasisStruct **SBasis, BasisStruct **DBasis, int *GBasisBound, int *SBasisBound, int *DBasisBound, const int UnitTotal, const int EHRadius, const int HDomain, const int VTotal, const int VRadius, const RealType VSlope, const int VExtra, const int Periodic, const int Embarrassing);
void BuildMatrixBasis(const char *FileName, MatrixBasisStruct **Basis, int *BasisTotal, const BasisStruct *KetBasis, const BasisStruct *BraBasis, const int KetBasisTotal, const int BraBasisTotal, const int Embarrassing);
void BuildHierarchy(int **Hierarchy, int **HierarchyCouplings, int *HierarchyTotal, const int UnitTotal, const int Depth, const int Embarrassing);
void SaveOperator(const char *FileName, const RealType *Operator, const BasisStruct *PriBasis, const BasisStruct *SecBasis, const int PriDimension, const int SecDimension, const int Embarrassing);
void PrintBasisState(const BasisStruct BS);

int CreateX(BasisStruct *BS, const int UnitId);
int CreateE(BasisStruct *BS, const int UnitId);
int CreateH(BasisStruct *BS, const int UnitId);
int AnnihilateX(BasisStruct *BS, const int UnitId);
int AnnihilateE(BasisStruct *BS, const int UnitId);
int AnnihilateH(BasisStruct *BS, const int UnitId);
int IdentifyBasisState(RealType *EffectiveVOverlap, const BasisStruct BS, const BasisStruct RefBS, const VOverlapStruct VOverlaps, const int VTotal);
int IdentifyXEHBasisState(const BasisStruct BS, const BasisStruct RefBS);
int SaveBasisState(FILE *FileID, const BasisStruct BS);

#endif