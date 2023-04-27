#include "GlobalsMod.h"
#include "ToolsMod.h"
#include "BasisMod.h"

#include <string.h>


static void BuildGBasis(BasisStruct *GBasis, int *GBasisBound, const int UnitTotal, const int VTotal, const int VRadius, const RealType VSlope, const int VExtra, const int Periodic);
static void BuildSBasis(BasisStruct *SBasis, int *SBasisBound, const int UnitTotal, const int EHRadius, const int HDomain, const int VTotal, const int VRadius, const RealType VSlope, const int VExtra, const int Periodic);
static void BuildDBasis(BasisStruct *DBasis, int *DBasisBound, const int UnitTotal, const int VTotal, const int VRadius, const RealType VSlope, const int VExtra, const int Periodic);
static void ClearBasisState(BasisStruct *BS);
static void ClearParticle(ParticleStruct *Particle);
static void AssignBasisState(BasisStruct *Basis, const BasisStruct BS, int *Id);
static void SortBasisState(BasisStruct *BS);
static void SortParticles(ParticleStruct *Particles, const int Total);
//static void BuildHierarchyLocal(int *Hierarchy, int *HierarchyTotal, const int XEHBasisTotal, const int Depth);
static void BuildHierarchyLocal(int *Hierarchy, int *HierarchyTotal, const int UnitTotal, const int Depth);
static void SaveGBasis(const BasisStruct *GBasis, const int *GBasisBound);
static void SaveSBasis(const BasisStruct *SBasis, const int *SBasisBound);
static void SaveDBasis(const BasisStruct *DBasis, const int *DBasisBound);
static void SaveMatrixBasis(const char *FileName, const MatrixBasisStruct *MatrixBasis, const BasisStruct *KetBasis, const BasisStruct *BraBasis, const int Dimension);
static void BasisStateToString(char *Output, const BasisStruct BS);
static void ParticleToString(char *Output, int *CharCount, int *Printed, const char *Character, const int U, const int V);

static int ValidBasisState(const BasisStruct BS, const int UnitTotal, const int VTotal, const int VRadius, const RealType VSlope, const int Periodic);
static int ValidVibron(const BasisStruct BS, const ParticleStruct Vibron, const int UnitTotal, const int VTotal, const int VRadius, const RealType VSlope, const int Periodic);
static int IncreaseHierarchyIndex(int *Index, const int UnitId, const int UnitTotal, const int Depth);
static int IdentifyHierarchyIndex(const int *Index, const int *RefIndex, int Units);


/////////////
// Globals //
/////////////


void BuildBases(BasisStruct **GBasis, BasisStruct **SBasis, BasisStruct **DBasis, int *GBasisBound, int *SBasisBound, int *DBasisBound, const int UnitTotal, const int EHRadius, const int HDomain, const int VTotal, const int VRadius, const RealType VSlope, const int VExtra, const int Periodic, const int Embarrassing)
{
	PrintMessage("Constructing the bases.", "", "", true, Embarrassing);
	
	if (GBasis != NULL) {
		
		BuildGBasis(NULL, GBasisBound, UnitTotal, VTotal, VRadius, VSlope, VExtra, Periodic);
		*GBasis = malloc(*GBasisBound * sizeof(BasisStruct));
		BuildGBasis(*GBasis, NULL, UnitTotal, VTotal, VRadius, VSlope, VExtra, Periodic);
		
		if (!Embarrassing) SaveGBasis(*GBasis, GBasisBound);
	}
	
	if (SBasis != NULL) {
		
		BuildSBasis(NULL, SBasisBound, UnitTotal, EHRadius, HDomain, VTotal, VRadius, VSlope, VExtra, Periodic);
		*SBasis = malloc(*SBasisBound * sizeof(BasisStruct));
		BuildSBasis(*SBasis, NULL, UnitTotal, EHRadius, HDomain, VTotal, VRadius, VSlope, VExtra, Periodic);
		
		if (!Embarrassing) SaveSBasis(*SBasis, SBasisBound);
	}
	
	if (DBasis != NULL) {
		
		BuildDBasis(NULL, DBasisBound, UnitTotal, VTotal, VRadius, VSlope, VExtra, Periodic);
		*DBasis = malloc(*DBasisBound * sizeof(BasisStruct));
		BuildDBasis(*DBasis, NULL, UnitTotal, VTotal, VRadius, VSlope, VExtra, Periodic);
		
		if (!Embarrassing) SaveDBasis(*DBasis, DBasisBound);
	}
}


void BuildMatrixBasis(const char *FileName, MatrixBasisStruct **Basis, int *BasisTotal, const BasisStruct *KetBasis, const BasisStruct *BraBasis, const int KetBasisTotal, const int BraBasisTotal, const int Embarrassing)
{
	int KetCount, BraCount;
	
	PrintMessage("Constructing ", FileName, ".", true, Embarrassing);
	
	*BasisTotal = KetBasisTotal * BraBasisTotal;
	
	*Basis = malloc(*BasisTotal * sizeof(MatrixBasisStruct));
	
	for (KetCount = 0; KetCount < KetBasisTotal; KetCount++)
		for (BraCount = 0; BraCount < BraBasisTotal; BraCount++) {
			(*Basis)[KetCount * BraBasisTotal + BraCount].Ket = KetCount;
			(*Basis)[KetCount * BraBasisTotal + BraCount].Bra = BraCount;
		}
	
	if (!Embarrassing) SaveMatrixBasis(FileName, *Basis, KetBasis, BraBasis, *BasisTotal);
}


int CreateX(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < XTotal; Count++) if ((*BS).X[Count].U == UnitId) return false;
	for (Count = 0; Count < ETotal; Count++) if ((*BS).E[Count].U == UnitId || (*BS).H[Count].U == UnitId) return false;
	
	for (Count = 0; Count < GTotal; Count++)
		if ((*BS).G[Count].U == UnitId) {
			memcpy(&(*BS).X[XTotal - 1], &(*BS).G[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).G[Count]);
			SortBasisState(BS);
			return true;
		}
	
	(*BS).X[XTotal - 1].U = UnitId;
	(*BS).X[XTotal - 1].VG = 0;
	SortBasisState(BS);
	return true;
}


int CreateE(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < XTotal; Count++) if ((*BS).X[Count].U == UnitId) return false;
	for (Count = 0; Count < ETotal; Count++) if ((*BS).E[Count].U == UnitId) return false;
	
	for (Count = 0; Count < HTotal; Count++)
		if ((*BS).H[Count].U == UnitId) {
			memcpy(&(*BS).X[XTotal - 1], &(*BS).H[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).H[Count]);
			SortBasisState(BS);
			return true;
		}
	
	for (Count = 0; Count < GTotal; Count++)
		if ((*BS).G[Count].U == UnitId) {
			memcpy(&(*BS).E[ETotal - 1], &(*BS).G[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).G[Count]);
			SortBasisState(BS);
			return true;
		}
	
	(*BS).E[ETotal - 1].U = UnitId;
	(*BS).E[ETotal - 1].VG = 0;
	SortBasisState(BS);
	return true;
}


int CreateH(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < XTotal; Count++) if ((*BS).X[Count].U == UnitId) return false;
	for (Count = 0; Count < HTotal; Count++) if ((*BS).H[Count].U == UnitId) return false;
	
	for (Count = 0; Count < ETotal; Count++)
		if ((*BS).E[Count].U == UnitId) {
			memcpy(&(*BS).X[XTotal - 1], &(*BS).E[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).E[Count]);
			SortBasisState(BS);
			return true;
		}
	
	for (Count = 0; Count < GTotal; Count++)
		if ((*BS).G[Count].U == UnitId) {
			memcpy(&(*BS).H[HTotal - 1], &(*BS).G[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).G[Count]);
			SortBasisState(BS);
			return true;
		}
	
	(*BS).H[HTotal - 1].U = UnitId;
	(*BS).H[HTotal - 1].VG = 0;
	SortBasisState(BS);
	return true;
}


int AnnihilateX(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < XTotal; Count++)
		if ((*BS).X[Count].U == UnitId) {
			memcpy(&(*BS).G[GTotal - 1], &(*BS).X[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).X[Count]);
			SortBasisState(BS);
			return true;
		}
	
	return false;
}


int AnnihilateE(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < XTotal; Count++)
		if ((*BS).X[Count].U == UnitId) {
			memcpy(&(*BS).H[HTotal - 1], &(*BS).X[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).X[Count]);
			SortBasisState(BS);
			return true;
		}
	
	for (Count = 0; Count < ETotal; Count++)
		if ((*BS).E[Count].U == UnitId) {
			memcpy(&(*BS).G[GTotal - 1], &(*BS).E[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).E[Count]);
			SortBasisState(BS);
			return true;
		}
	
	return false;
}


int AnnihilateH(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < XTotal; Count++)
		if ((*BS).X[Count].U == UnitId) {
			memcpy(&(*BS).E[ETotal - 1], &(*BS).X[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).X[Count]);
			SortBasisState(BS);
			return true;
		}
	
	for (Count = 0; Count < HTotal; Count++)
		if ((*BS).H[Count].U == UnitId) {
			memcpy(&(*BS).G[GTotal - 1], &(*BS).H[Count], sizeof(ParticleStruct));
			ClearParticle(&(*BS).H[Count]);
			SortBasisState(BS);
			return true;
		}
	
	return false;
}


int IdentifyBasisState(RealType *EffectiveVOverlap, const BasisStruct BS, const BasisStruct RefBS, const VOverlapStruct VOverlaps, const int VTotal)
{
	int Count, Offset = 0;
	
	*EffectiveVOverlap = 1;
	
	for (Count = 0; Count < XTotal; Count++) {
		if (BS.X[Count].U != RefBS.X[Count].U) return false;
		else if (BS.X[Count].VX == -1) {
			if		(BS.X[Count].VG != -1)		(*EffectiveVOverlap) *= VOverlaps.XG[RefBS.X[Count].VX * VTotal + BS.X[Count].VG];
			else if	(BS.X[Count].VE != -1)		(*EffectiveVOverlap) *= VOverlaps.XE[RefBS.X[Count].VX * VTotal + BS.X[Count].VE];
			else if	(BS.X[Count].VH != -1)		(*EffectiveVOverlap) *= VOverlaps.XH[RefBS.X[Count].VX * VTotal + BS.X[Count].VH];
		}
		else if (BS.X[Count].VX != RefBS.X[Count].VX) return false;
	}
	
	for (Count = 0; Count < ETotal; Count++) {
		if (BS.E[Count].U != RefBS.E[Count].U) return false;
		else if (BS.E[Count].VE == -1) {
			if		(BS.E[Count].VX != -1)		(*EffectiveVOverlap) *= VOverlaps.XE[BS.E[Count].VX * VTotal + RefBS.E[Count].VE];
			else if	(BS.E[Count].VG != -1)		(*EffectiveVOverlap) *= VOverlaps.EG[RefBS.E[Count].VE * VTotal + BS.E[Count].VG];
		}
		else if (BS.E[Count].VE != RefBS.E[Count].VE) return false;
	}
	
	for (Count = 0; Count < HTotal; Count++) {
		if (BS.H[Count].U != RefBS.H[Count].U) return false;
		else if (BS.H[Count].VH == -1) {
			if		(BS.H[Count].VX != -1)		(*EffectiveVOverlap) *= VOverlaps.XH[BS.H[Count].VX * VTotal + RefBS.H[Count].VH];
			else if	(BS.H[Count].VG != -1)		(*EffectiveVOverlap) *= VOverlaps.HG[RefBS.H[Count].VH * VTotal + BS.H[Count].VG];
		}
		else if (BS.H[Count].VH != RefBS.H[Count].VH) return false;
	}
	
	for (Count = 0; Count < GTotal; Count++) {
		if (BS.G[Count].U != RefBS.G[Count - Offset].U) {
			if (BS.G[Count].U != -1 && BS.G[Count].VG == -1) {
				if		(BS.G[Count].VX != -1)	(*EffectiveVOverlap) *= VOverlaps.XG[BS.G[Count].VX * VTotal];
				else if	(BS.G[Count].VE != -1)	(*EffectiveVOverlap) *= VOverlaps.EG[BS.G[Count].VE * VTotal];
				else							(*EffectiveVOverlap) *= VOverlaps.HG[BS.G[Count].VH * VTotal];
				Offset++;
			}
			else return false;
		}
		else if (BS.G[Count].VG == -1) {
			if		(BS.G[Count].VX != -1)		(*EffectiveVOverlap) *= VOverlaps.XG[BS.G[Count].VX * VTotal + RefBS.G[Count - Offset].VG];
			else if	(BS.G[Count].VE != -1)		(*EffectiveVOverlap) *= VOverlaps.EG[BS.G[Count].VE * VTotal + RefBS.G[Count - Offset].VG];
			else if (BS.G[Count].VH != -1)		(*EffectiveVOverlap) *= VOverlaps.HG[BS.G[Count].VH * VTotal + RefBS.G[Count - Offset].VG];
		}
		else if (BS.G[Count].VG != RefBS.G[Count - Offset].VG) return false;
	}
	
	return true;
}


int IdentifyXEHBasisState(const BasisStruct BS, const BasisStruct RefBS)
{
	int Count;
	const int *BSEntry = &BS.X[0].U, *RefBSEntry = &RefBS.X[0].U;
	
	for (Count = 0; Count < XTotal + ETotal + HTotal; Count++) if ((*BSEntry + sizeof(ParticleStruct)) != (*RefBSEntry + sizeof(ParticleStruct))) return false;
	return true;
}


/*void BuildHierarchy(int **Hierarchy, int **HierarchyCouplings, int *HierarchyTotal, const int XEHBasisTotal, const int Depth, const int Embarrassing)
{
	int *Index = calloc(XEHBasisTotal, sizeof(int));
	int PriCount, SecCount, XEHBasisCount;
	
	PrintMessage("Constructing the hierarchy.", "", "", true, Embarrassing);
	
	BuildHierarchyLocal(NULL, HierarchyTotal, XEHBasisTotal, Depth);
	*Hierarchy = malloc(*HierarchyTotal * XEHBasisTotal * sizeof(int));
	BuildHierarchyLocal(*Hierarchy, NULL, XEHBasisTotal, Depth);
	
	*HierarchyCouplings = malloc(*HierarchyTotal * 2 * XEHBasisTotal * sizeof(int));
	memset(*HierarchyCouplings, -1, *HierarchyTotal * 2 * XEHBasisTotal * sizeof(int));
	
	for (PriCount = 0; PriCount < *HierarchyTotal; PriCount++)
		for (XEHBasisCount = 0; XEHBasisCount < XEHBasisTotal; XEHBasisCount++) {
			memcpy(Index, &(*Hierarchy)[PriCount * XEHBasisTotal], XEHBasisTotal * sizeof(int));
			
			if (Index[XEHBasisCount] != 0) {
				Index[XEHBasisCount]--;
				for (SecCount = PriCount - 1; SecCount >= 0; SecCount--)
					if (IdentifyHierarchyIndex(Index, &(*Hierarchy)[SecCount * XEHBasisTotal], XEHBasisTotal)) {
						(*HierarchyCouplings)[PriCount * 2 * XEHBasisTotal + XEHBasisCount] = SecCount;
						break;
					}
				Index[XEHBasisCount]++;
			}
			
			if (SumIArray(Index, XEHBasisTotal) != Depth) {
				Index[XEHBasisCount]++;
				for (SecCount = PriCount + 1; SecCount < *HierarchyTotal; SecCount++)
					if (IdentifyHierarchyIndex(Index, &(*Hierarchy)[SecCount * XEHBasisTotal], XEHBasisTotal)) {
						(*HierarchyCouplings)[(PriCount * 2 + 1) * XEHBasisTotal + XEHBasisCount] = SecCount;
						break;
					}
			}
		}
	
	free(Index);
}*/


void BuildHierarchy(int **Hierarchy, int **HierarchyCouplings, int *HierarchyTotal, const int UnitTotal, const int Depth, const int Embarrassing)
{
	int *Index = calloc(UnitTotal, sizeof(int));
	int PriCount, SecCount, UnitCount;
	
	PrintMessage("Constructing the hierarchy.", "", "", true, Embarrassing);
	
	BuildHierarchyLocal(NULL, HierarchyTotal, UnitTotal, Depth);
	*Hierarchy = malloc(*HierarchyTotal * UnitTotal * sizeof(int));
	BuildHierarchyLocal(*Hierarchy, NULL, UnitTotal, Depth);
	
	*HierarchyCouplings = malloc(*HierarchyTotal * 2 * UnitTotal * sizeof(int));
	memset(*HierarchyCouplings, -1, *HierarchyTotal * 2 * UnitTotal * sizeof(int));
	
	for (PriCount = 0; PriCount < *HierarchyTotal; PriCount++)
		for (UnitCount = 0; UnitCount < UnitTotal; UnitCount++) {
			memcpy(Index, &(*Hierarchy)[PriCount * UnitTotal], UnitTotal * sizeof(int));
			
			if (Index[UnitCount] != 0) {
				Index[UnitCount]--;
				for (SecCount = PriCount - 1; SecCount >= 0; SecCount--)
					if (IdentifyHierarchyIndex(Index, &(*Hierarchy)[SecCount * UnitTotal], UnitTotal)) {
						(*HierarchyCouplings)[PriCount * 2 * UnitTotal + UnitCount] = SecCount;
						break;
					}
				Index[UnitCount]++;
			}
			
			if (SumIArray(Index, UnitTotal) != Depth) {
				Index[UnitCount]++;
				for (SecCount = PriCount + 1; SecCount < *HierarchyTotal; SecCount++)
					if (IdentifyHierarchyIndex(Index, &(*Hierarchy)[SecCount * UnitTotal], UnitTotal)) {
						(*HierarchyCouplings)[(PriCount * 2 + 1) * UnitTotal + UnitCount] = SecCount;
						break;
					}
			}
		}
	
	free(Index);
}


void SaveOperator(const char *FileName, const RealType *Operator, const BasisStruct *PriBasis, const BasisStruct *SecBasis, const int PriDimension, const int SecDimension, const int Embarrassing)
{
	char ModFileName[20];
	int PriCount, SecCount, TagCount;
	FILE *FileID;
	
	if (Embarrassing) return;
	
	sprintf(ModFileName, "%s%s", FileName, ".dat");
	
	FileID = fopen(ModFileName, "w");
	
	if (IntMax(PriDimension, SecDimension) > 70) fprintf(FileID, "Operator too large for storage!");
	else if (IntMin(PriDimension, SecDimension) == 0) fprintf(FileID, "Operator empty!");
	else {
		
		fprintf(FileID, "                ");
		for (PriCount = 0; PriCount < PriDimension; PriCount++) {
			if (PriBasis == NULL) fprintf(FileID, "%2d               ", PriCount);
			else {
				fprintf(FileID, "|");
				TagCount = SaveBasisState(FileID, PriBasis[PriCount]);
				fprintf(FileID, ">");
				for (; TagCount < 15; TagCount++) fprintf(FileID, " ");
			}
		}
		fprintf(FileID, "\n");
		for (SecCount = 0; SecCount < SecDimension; SecCount++) {
			if (SecBasis == NULL) fprintf(FileID, "%2d               ", SecCount);
			else {
				fprintf(FileID, "|");
				TagCount = SaveBasisState(FileID, SecBasis[SecCount]);
				fprintf(FileID, ">");
				for (; TagCount < 15; TagCount++) fprintf(FileID, " ");
			}
			for (PriCount = 0; PriCount < PriDimension; PriCount++) fprintf(FileID, "%10.3f       ", Operator[PriCount * SecDimension + SecCount]);
			fprintf(FileID, "\n");
		}
		
	}
	
	fclose(FileID);
}


void PrintBasisState(const BasisStruct BS)
{
	char Tag[100];
	
	BasisStateToString(Tag, BS);
	printf("%s\n", Tag);
}


int SaveBasisState(FILE *FileID, const BasisStruct BS)
{
	char Tag[100];
	
	BasisStateToString(Tag, BS);
	fprintf(FileID, "%s", Tag);
	
	return strlen(Tag);
}


////////////
// Locals //
////////////


static void BuildGBasis(BasisStruct *GBasis, int *GBasisBound, const int UnitTotal, const int VTotal, const int VRadius, const RealType VSlope, const int VExtra, const int Periodic)
{
	BasisStruct BS;
	int Id = 0;
	
	ClearBasisState(&BS);
	
	AssignBasisState(GBasis, BS, &Id);
	
	for (BS.G[0].U = 0; BS.G[0].U < UnitTotal; BS.G[0].U++)
		for (BS.G[0].VG = 1; BS.G[0].VG < VTotal; BS.G[0].VG++)
			AssignBasisState(GBasis, BS, &Id);
	
	if (GBasisBound != NULL) GBasisBound[1] = Id;
	
	for (BS.G[0].U = 1; BS.G[0].U < UnitTotal; BS.G[0].U++)
		for (BS.G[1].U = 0; BS.G[1].U < BS.G[0].U; BS.G[1].U++)
			for (BS.G[0].VG = 1; BS.G[0].VG < VTotal - 1; BS.G[0].VG++)
				for (BS.G[1].VG = 1; BS.G[1].VG < VTotal - BS.G[0].VG; BS.G[1].VG++)
					if (ValidBasisState(BS, UnitTotal, VTotal, VRadius, VSlope, Periodic)) AssignBasisState(GBasis, BS, &Id);
	
	if (GBasisBound != NULL) GBasisBound[2] = Id;
	
	if (VExtra)
		for (BS.G[0].U = 2; BS.G[0].U < UnitTotal; BS.G[0].U++)
			for (BS.G[1].U = 1; BS.G[1].U < BS.G[0].U; BS.G[1].U++)
				for (BS.G[2].U = 0; BS.G[2].U < BS.G[1].U; BS.G[2].U++)
					for (BS.G[0].VG = 1; BS.G[0].VG < VTotal - 2; BS.G[0].VG++)
						for (BS.G[1].VG = 1; BS.G[1].VG < VTotal - BS.G[0].VG - 1; BS.G[1].VG++)
							for (BS.G[2].VG = 1; BS.G[2].VG < VTotal - BS.G[0].VG - BS.G[1].VG; BS.G[2].VG++)
								if (ValidBasisState(BS, UnitTotal, VTotal, VRadius, VSlope, Periodic)) AssignBasisState(GBasis, BS, &Id);
	
	if (GBasisBound != NULL) GBasisBound[3] = Id;
	
	if (GBasisBound != NULL) GBasisBound[0] = Id;
}


static void BuildSBasis(BasisStruct *SBasis, int *SBasisBound, const int UnitTotal, const int EHRadius, const int HDomain, const int VTotal, const int VRadius, const RealType VSlope, const int VExtra, const int Periodic)
{
	BasisStruct BS;
	int Id = 0;
	
	ClearBasisState(&BS);
	
	for (BS.X[0].U = 0; BS.X[0].U < (HDomain ? HDomain : UnitTotal); BS.X[0].U++)
		for (BS.X[0].VX = 0; BS.X[0].VX < VTotal; BS.X[0].VX++)
			AssignBasisState(SBasis, BS, &Id);
	
	if (SBasisBound != NULL) SBasisBound[1] = Id;
	
	for (BS.X[0].U = 0; BS.X[0].U < (HDomain ? HDomain : UnitTotal); BS.X[0].U++)
		for (BS.G[0].U = 0; BS.G[0].U < UnitTotal; BS.G[0].U++) {
			if (BS.G[0].U == BS.X[0].U) continue;
			for (BS.X[0].VX = 0; BS.X[0].VX < VTotal - 1; BS.X[0].VX++)
				for (BS.G[0].VG = 1; BS.G[0].VG < VTotal - BS.X[0].VX; BS.G[0].VG++)
					if (ValidBasisState(BS, UnitTotal, VTotal, VRadius, VSlope, Periodic)) AssignBasisState(SBasis, BS, &Id);
		}
	
	if (SBasisBound != NULL) SBasisBound[2] = Id;
	
	if (VExtra)
		for (BS.X[0].U = 0; BS.X[0].U < (HDomain ? HDomain : UnitTotal); BS.X[0].U++)
			for (BS.G[0].U = 1; BS.G[0].U < UnitTotal; BS.G[0].U++) {
				if (BS.G[0].U == BS.X[0].U) continue;
				for (BS.G[1].U = 0; BS.G[1].U < BS.G[0].U; BS.G[1].U++) {
					if (BS.G[1].U == BS.X[0].U) continue;
					for (BS.X[0].VX = 0; BS.X[0].VX < VTotal - 2; BS.X[0].VX++)
						for (BS.G[0].VG = 1; BS.G[0].VG < VTotal - BS.X[0].VX - 1; BS.G[0].VG++)
							for (BS.G[1].VG = 1; BS.G[1].VG < VTotal - BS.X[0].VX - BS.G[0].VG; BS.G[1].VG++)
								if (ValidBasisState(BS, UnitTotal, VTotal, VRadius, VSlope, Periodic)) AssignBasisState(SBasis, BS, &Id);
				}
			}
	
	if (SBasisBound != NULL) SBasisBound[3] = Id;
	
	ClearBasisState(&BS);
	
	for (BS.E[0].U = HDomain; BS.E[0].U < UnitTotal; BS.E[0].U++)
		for (BS.H[0].U = 0; BS.H[0].U < (HDomain ? HDomain : UnitTotal); BS.H[0].U++) {
			if (BS.H[0].U == BS.E[0].U || Separation(BS.E[0].U, BS.H[0].U, UnitTotal, Periodic) > EHRadius) continue;
			for (BS.E[0].VE = 0; BS.E[0].VE < VTotal; BS.E[0].VE++)
				for (BS.H[0].VH = 0; BS.H[0].VH < VTotal - BS.E[0].VE; BS.H[0].VH++)
					if (ValidBasisState(BS, UnitTotal, VTotal, VRadius, VSlope, Periodic)) AssignBasisState(SBasis, BS, &Id);
		}
	
	if (SBasisBound != NULL) SBasisBound[4] = Id;
	
	if (SBasisBound != NULL) SBasisBound[0] = Id;
}


static void BuildDBasis(BasisStruct *DBasis, int *DBasisBound, const int UnitTotal, const int VTotal, const int VRadius, const RealType VSlope, const int VExtra, const int Periodic)
{
	BasisStruct BS;
	int Id = 0;
	
	ClearBasisState(&BS);
	
	for (BS.X[0].U = 1; BS.X[0].U < UnitTotal; BS.X[0].U++)
		for (BS.X[1].U = 0; BS.X[1].U < BS.X[0].U; BS.X[1].U++)
			for (BS.X[0].VX = 0; BS.X[0].VX < VTotal; BS.X[0].VX++)
				for (BS.X[1].VX = 0; BS.X[1].VX < VTotal - BS.X[0].VX; BS.X[1].VX++)
					AssignBasisState(DBasis, BS, &Id);
	
	if (DBasisBound != NULL) DBasisBound[1] = Id;
	
	for (BS.X[0].U = 1; BS.X[0].U < UnitTotal; BS.X[0].U++)
		for (BS.X[1].U = 0; BS.X[1].U < BS.X[0].U; BS.X[1].U++)
			for (BS.G[0].U = 0; BS.G[0].U < UnitTotal; BS.G[0].U++) {
				if (BS.G[0].U == BS.X[0].U || BS.G[0].U == BS.X[1].U) continue;
				for (BS.X[0].VX = 0; BS.X[0].VX < VTotal - 1; BS.X[0].VX++)
					for (BS.X[1].VX = 0; BS.X[1].VX < VTotal - BS.X[0].VX - 1; BS.X[1].VX++)
						for (BS.G[0].VG = 1; BS.G[0].VG < VTotal - BS.X[0].VX - BS.X[1].VX; BS.G[0].VG++)
							if (ValidBasisState(BS, UnitTotal, VTotal, VRadius, VSlope, Periodic)) AssignBasisState(DBasis, BS, &Id);
			}
	
	if (DBasisBound != NULL) DBasisBound[2] = Id;
	
	if (VExtra)
		for (BS.X[0].U = 1; BS.X[0].U < UnitTotal; BS.X[0].U++)
			for (BS.X[1].U = 0; BS.X[1].U < BS.X[0].U; BS.X[1].U++)
				for (BS.G[0].U = 1; BS.G[0].U < UnitTotal; BS.G[0].U++) {
					if (BS.G[0].U == BS.X[0].U || BS.G[0].U == BS.X[1].U) continue;
					for (BS.G[1].U = 0; BS.G[1].U < BS.G[0].U; BS.G[1].U++) {
						if (BS.G[1].U == BS.X[0].U || BS.G[1].U == BS.X[1].U) continue;
						for (BS.X[0].VX = 0; BS.X[0].VX < VTotal - 2; BS.X[0].VX++)
							for (BS.X[1].VX = 0; BS.X[1].VX < VTotal - BS.X[0].VX - 2; BS.X[1].VX++)
								for (BS.G[0].VG = 1; BS.G[0].VG < VTotal - BS.X[0].VX - BS.X[1].VX; BS.G[0].VG++)
									for (BS.G[1].VG = 1; BS.G[1].VG < VTotal - BS.X[0].VX - BS.X[1].VX - BS.G[0].VG; BS.G[1].VG++)
										if (ValidBasisState(BS, UnitTotal, VTotal, VRadius, VSlope, Periodic)) AssignBasisState(DBasis, BS, &Id);
					}
				}
	
	if (DBasisBound != NULL) DBasisBound[3] = Id;
	
	// Room for charge transfer states!
	
	if (DBasisBound != NULL) DBasisBound[0] = Id;
}


static void ClearBasisState(BasisStruct *BS)
{
	memset(BS, -1, sizeof(BasisStruct));
}


static void ClearParticle(ParticleStruct *Particle)
{
	memset(Particle, -1, sizeof(ParticleStruct));
}


static void AssignBasisState(BasisStruct *Basis, const BasisStruct BS, int *Id)
{
	if (Basis != NULL) memcpy(&Basis[*Id], &BS, sizeof(BasisStruct));
	(*Id)++;
}


static void SortBasisState(BasisStruct *BS)
{
	SortParticles(&(*BS).X[0], XTotal);
	SortParticles(&(*BS).E[0], ETotal);
	SortParticles(&(*BS).H[0], HTotal);
	SortParticles(&(*BS).G[0], GTotal);
}


static void SortParticles(ParticleStruct *Particles, const int Total)
{
	ParticleStruct *WorkParticle = malloc(sizeof(ParticleStruct));
	int PriCount, SecCount;
	
	for (PriCount = 1; PriCount < Total; PriCount++) {
		memcpy(WorkParticle, &Particles[PriCount], sizeof(ParticleStruct));
		
		for (SecCount = PriCount - 1; SecCount >= 0 && Particles[SecCount].U < (*WorkParticle).U; SecCount--)
			memcpy(&Particles[SecCount + 1], &Particles[SecCount], sizeof(ParticleStruct));
		
		memcpy(&Particles[SecCount + 1], WorkParticle, sizeof(ParticleStruct));
	}
	
	free(WorkParticle);
}


static int ValidBasisState(const BasisStruct BS, const int UnitTotal, const int VTotal, const int VRadius, const RealType VSlope, const int Periodic)
{
	int Count;
	
	for (Count = 0; Count < GTotal && BS.G[Count].U != -1; Count++)
		if (!ValidVibron(BS, BS.G[Count], UnitTotal, VTotal, VRadius, VSlope, Periodic)) return false;
	
	return true;
}


static int ValidVibron(const BasisStruct BS, const ParticleStruct Vibron, const int UnitTotal, const int VTotal, const int VRadius, const RealType VSlope, const int Periodic)
{
	int Distance = UnitTotal, Count;
	
	for (Count = 0; Count < XTotal; Count++)
		if (BS.X[Count].U != -1) Distance = IntMin(Distance, Separation(BS.X[Count].U, Vibron.U, UnitTotal, Periodic));
	
	for (Count = 0; Count < ETotal; Count++)
		if (BS.E[Count].U != -1) Distance = IntMin(Distance, Separation(BS.E[Count].U, Vibron.U, UnitTotal, Periodic));
	
	for (Count = 0; Count < HTotal; Count++)
		if (BS.H[Count].U != -1) Distance = IntMin(Distance, Separation(BS.H[Count].U, Vibron.U, UnitTotal, Periodic));
	
	if (Distance < UnitTotal) {
		
		if (Distance > VRadius) return false;
		if (VSlope != 0 && !(Vibron.VG < VTotal - Distance / VSlope)) return false;
		
	}
	else {
		
		for (Count = 0; Count < GTotal; Count++)
			if (BS.G[Count].U != -1 && BS.G[Count].U != Vibron.U) Distance = IntMin(Distance, Separation(BS.G[Count].U, Vibron.U, UnitTotal, Periodic));
		
		if (Distance > VRadius) return false;
		if (VSlope != 0 && !(VTotal - Distance / VSlope) > 0) return false;
		
	}
	
	return true;
}


/*static void BuildHierarchyLocal(int *Hierarchy, int *HierarchyTotal, const int XEHBasisTotal, const int Depth)
{
	int *Index = calloc(XEHBasisTotal, sizeof(int));
	int Count = 0;
	
	do {
		if (Hierarchy != NULL) memcpy(&Hierarchy[Count * XEHBasisTotal], Index, XEHBasisTotal * sizeof(int));
		Count++;
	} while (IncreaseHierarchyIndex(Index, 0, XEHBasisTotal, Depth));
	
	if (HierarchyTotal != NULL) *HierarchyTotal = Count;
	
	free(Index);
}*/


static void BuildHierarchyLocal(int *Hierarchy, int *HierarchyTotal, const int UnitTotal, const int Depth)
{
	int *Index = calloc(UnitTotal, sizeof(int));
	int Count = 0;
	
	do {
		if (Hierarchy != NULL) memcpy(&Hierarchy[Count * UnitTotal], Index, UnitTotal * sizeof(int));
		Count++;
	} while (IncreaseHierarchyIndex(Index, 0, UnitTotal, Depth));
	
	if (HierarchyTotal != NULL) *HierarchyTotal = Count;
	
	free(Index);
}


static int IncreaseHierarchyIndex(int *Index, const int UnitId, const int UnitTotal, const int Depth)
{
	if (SumIArray(Index, UnitTotal) == Depth) {
		if (UnitId + 1 == UnitTotal) return false;
		Index[UnitId] = 0;
		return IncreaseHierarchyIndex(Index, UnitId + 1, UnitTotal, Depth);
	}
	Index[UnitId]++;
	return true;
}


static int IdentifyHierarchyIndex(const int *Index, const int *RefIndex, int Units)
{
	for (Units--; Units >= 0; Units--) if (Index[Units] != RefIndex[Units]) return false;
	
	return true;
}


static void SaveGBasis(const BasisStruct *GBasis, const int *GBasisBound)
{
	int Count;
	FILE *FileID;
	
	FileID = fopen("GBasis.dat", "w");
	
	fprintf(FileID, "Total: %d\n", *GBasisBound);
	fprintf(FileID, "\n   Vacuum state\n");
	for (Count = 0; Count < *GBasisBound; Count++) {
		if (Count == 1)						fprintf(FileID, "\n   Single vibration: %d\n", GBasisBound[1] - 1);
		else if (Count == GBasisBound[1])	fprintf(FileID, "\n   Two vibrations: %d\n", GBasisBound[2] - GBasisBound[1]);
		else if (Count == GBasisBound[2])	fprintf(FileID, "\n   Three vibrations: %d\n", GBasisBound[3] - GBasisBound[2]);
		fprintf(FileID, "%8d: |", Count);
		SaveBasisState(FileID, GBasis[Count]);
		fprintf(FileID, ">\n");
	}
	
	fclose(FileID);
}


static void SaveSBasis(const BasisStruct *SBasis, const int *SBasisBound)
{
	int Count;
	FILE *FileID;
	
	FileID = fopen("SBasis.dat", "w");
	
	fprintf(FileID, "Total: %d\n", *SBasisBound);
	fprintf(FileID, "\n   Purely vibronic: %d\n", SBasisBound[1]);
	for (Count = 0; Count < *SBasisBound; Count++) {
		if (Count == SBasisBound[3])		fprintf(FileID, "\n   Charge-transfer, purely vibronic: %d\n", SBasisBound[4] - SBasisBound[3]);
		else if (Count == SBasisBound[2])	fprintf(FileID, "\n   Vibronic + two vibrations: %d\n", SBasisBound[3] - SBasisBound[2]);
		else if (Count == SBasisBound[1])	fprintf(FileID, "\n   Vibronic + single vibration: %d\n", SBasisBound[2] - SBasisBound[1]);
		fprintf(FileID, "%8d: |", Count);
		SaveBasisState(FileID, SBasis[Count]);
		fprintf(FileID, ">\n");
	}
	
	fclose(FileID);
}


static void SaveDBasis(const BasisStruct *DBasis, const int *DBasisBound)
{
	int Count;
	FILE *FileID;
	
	FileID = fopen("DBasis.dat", "w");
	
	fprintf(FileID, "Total: %d\n", *DBasisBound);
	fprintf(FileID, "\n   Purely vibronic: %d\n", DBasisBound[1]);
	for (Count = 0; Count < *DBasisBound; Count++) {
		if (Count == DBasisBound[1])		fprintf(FileID, "\n   Vibronic + single vibration: %d\n", DBasisBound[2] - DBasisBound[1]);
		else if (Count == DBasisBound[2])	fprintf(FileID, "\n   Vibronic + two vibrations: %d\n", DBasisBound[3] - DBasisBound[2]);
		fprintf(FileID, "%8d: |", Count);
		SaveBasisState(FileID, DBasis[Count]);
		fprintf(FileID, ">\n");
	}
	
	fclose(FileID);
}


static void SaveMatrixBasis(const char *FileName, const MatrixBasisStruct *MatrixBasis, const BasisStruct *KetBasis, const BasisStruct *BraBasis, const int Dimension)
{
	char ModFileName[20];
	int Count;
	FILE *FileID;
	
	sprintf(ModFileName, "%s%s", FileName, ".dat");
	
	FileID = fopen(ModFileName, "w");
	
	for (Count = 0; Count < Dimension; Count++) {
		fprintf(FileID, "%8d: |", Count);
		SaveBasisState(FileID, KetBasis[MatrixBasis[Count].Ket]);
		fprintf(FileID, "><");
		SaveBasisState(FileID, BraBasis[MatrixBasis[Count].Bra]);
		fprintf(FileID, "|\n");
	}
	
	fclose(FileID);
}


static void BasisStateToString(char *Output, const BasisStruct BS)
{
	int Count, CharCount = 2, Printed = false;
	
	sprintf(Output, "");
	
	for (Count = 0; Count < XTotal; Count++)
		if (BS.X[Count].U != -1) ParticleToString(Output, &CharCount, &Printed, "x", BS.X[Count].U, BS.X[Count].VX);
	for (Count = 0; Count < ETotal; Count++)
		if (BS.E[Count].U != -1) ParticleToString(Output, &CharCount, &Printed, "e", BS.E[Count].U, BS.E[Count].VE);
	for (Count = 0; Count < HTotal; Count++)
		if (BS.H[Count].U != -1) ParticleToString(Output, &CharCount, &Printed, "h", BS.H[Count].U, BS.H[Count].VH);
	for (Count = 0; Count < GTotal; Count++)
		if (BS.G[Count].U != -1) ParticleToString(Output, &CharCount, &Printed, "g", BS.G[Count].U, BS.G[Count].VG);
	
	if (!Printed) sprintf(Output, "%sVAC", Output), CharCount += 3;
	
	sprintf(Output, "%s", Output);
}


static void ParticleToString(char *Output, int *CharCount, int *Printed, const char *Character, const int U, const int V)
{
	if (*Printed) sprintf(Output, "%s;", Output), (*CharCount)++;
	sprintf(Output, "%s%s%d,%d", Output, Character, U, V), (*CharCount) += 4;
	if (U > 9) (*CharCount)++;
	*Printed = true;
}