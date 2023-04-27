#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "RandomMod.h"
#include "ToolsMod.h"
#include "BasisMod.h"
#include "SubsMod.h"

#include <string.h>
#include <math.h>


static void BuildPropagator(ComplexType *Propagator, const RealType *EigenVectors, const RealType *EigenValues, const ParmStruct *Parms, const int Dimension);
static void BuildHamiltonianDiag(RealType *Diag, const RealType *BasicHamiltonian, const BathStruct *Bath, const BasisStruct *Basis, const int Dimension);
static void ApplyTrotter(ComplexType *Coef, const RealType *HamiltonianDiag, const int Dimension, const ParmStruct *Parms);
static void CoefDerivative(ComplexType *CoefDerivative, const ComplexType *Coef, const RealType *Liouville, const int *LiouvilleMask, const int *Hieperator, const int *AntiHieperator, const int BasisTotal, const int *Hierarchy, const int *HierarchyCouplings, const ParmStruct *Parms);
static void BuildQuantumForce(RealType **QuantumForce, const RealType *FeedbackVector, const BasisStruct *Basis, const ParmStruct *Parms);

static RealType LineStrength(const int GStateId, const RealType *GSDipoles, const int *GSDipolesMask, const RealType *PermDipoles, const RealType *SEigenVector, const PolStruct *Pols, const ParmStruct *Parms);
static RealType TransitionDipole(const int UnitId, const int GStateId, const RealType *GSDipoles, const int *GSDipolesMask, const RealType *SEigenVector, const ParmStruct *Parms);
static RealType BrownianOscillator(const RealType BathDynamicValue, const RealType Deviation, const ParmStruct *Parms);
static RealType DiagonalEnergy(const BasisStruct BasisState, const RealType *XEnergies, const ParmStruct *Parms);
static RealType TotalBath(const BasisStruct BasisState, const BathStruct *Bath);


RealType RandomValue;


/////////////
// Globals //
/////////////


void BuildXEnergies(RealType **XEnergies, const ParmStruct *Parms)
{
	RealType *XCouplings;
	int PriCount, SecCount, MolCount;
	double ReadDouble;
	FILE *FileID;
	
	*XEnergies = calloc(IntPow2(Parms->UnitTotal), sizeof(RealType));
	
	if (strncmp(Parms->XEnergiesFile, "|", 1) != 0) {

        printf("Reading XEnergies from file.");
		
		FileID = fopen(Parms->XEnergiesFile, "r");
		for (PriCount = 0; PriCount < Parms->UnitTotal; PriCount++) {
			for (SecCount = 0; SecCount < Parms->UnitTotal; SecCount++) {
				fscanf(FileID, "%lf ", &ReadDouble);
				(*XEnergies)[PriCount * Parms->UnitTotal + SecCount] = ReadDouble;
			}
			(*XEnergies)[PriCount * (Parms->UnitTotal  + 1)] -= Parms->XEnergyMean;
		}
		fclose(FileID);
		
	}
	else {
		
		XCouplings = calloc(Parms->UnitTotal * Parms->MolPerCell, sizeof(RealType));
		
		if (strncmp(Parms->XCouplingsFile, "|", 1) != 0) {
			
			FileID = fopen(Parms->XCouplingsFile, "r");
			for (PriCount = 1; PriCount < Parms->UnitTotal; PriCount++) 
            {
                for( MolCount = 0; MolCount < Parms->MolPerCell; MolCount++ )
                {
				    if (fscanf(FileID, "%lf ", &ReadDouble) == EOF) break;
				    XCouplings[PriCount + MolCount * Parms->UnitTotal] = ReadDouble;
                }
			}
			fclose(FileID);
			
		}
		else for (PriCount = 1; PriCount < IntMin(Parms->UnitTotal, Parms->XPointDipoleRadius + 1); PriCount++) XCouplings[PriCount] = Parms->XCouplingNN / IntPowN(PriCount, 3);
		
		for (PriCount = 0; PriCount < Parms->UnitTotal; PriCount++)
			for (SecCount = 0; SecCount < Parms->UnitTotal; SecCount++)
			{
				int separation = SignedSeparation(PriCount, SecCount, Parms->UnitTotal, Parms->Periodic);
				
				if( ( separation > 0 && (PriCount % Parms->MolPerCell == 0) ) || ( separation < 0 && (PriCount % Parms->MolPerCell == 1) ) || Parms->MolPerCell == 1 )
					(*XEnergies)[PriCount * Parms->UnitTotal + SecCount] = XCouplings[abs(separation)];
				else
					(*XEnergies)[PriCount * Parms->UnitTotal + SecCount] = XCouplings[abs(separation) + Parms->UnitTotal];
			}
				
		free(XCouplings);
		
		if (strncmp(Parms->XDiagsFile, "|", 1) != 0) {
			
			FileID = fopen(Parms->XDiagsFile, "r");
			for (PriCount = 0; PriCount < Parms->UnitTotal; PriCount++) {
				fscanf(FileID, "%lf ", &ReadDouble);
				(*XEnergies)[PriCount * (Parms->UnitTotal  + 1)] = ReadDouble;
			}
			fclose(FileID);
			
		}
		else
			for (PriCount = 0; PriCount < Parms->UnitTotal; PriCount++)
				(*XEnergies)[PriCount * (Parms->UnitTotal  + 1)] = Parms->XEnergyOffset * (PriCount - (RealType) (Parms->UnitTotal - 1) / 2);
		
	}
	
	if (!Parms->Embarrassing) SaveRArray("XEnergies", *XEnergies, Parms->UnitTotal, Parms->UnitTotal, NULL, false);
}


void BuildPermDipoles(RealType **PermDipoles, PolStruct *Pols, const ParmStruct *Parms)
{
	int UnitCount, PolCount;
	double ReadDouble;
	FILE *FileID;
	
	*PermDipoles = calloc(3 * Parms->UnitTotal, sizeof(RealType));
	
	Pols->Active[0] = true;
	
	if (Parms->HierarchyDepth) for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++) (*PermDipoles)[UnitCount] = 1;
	
	else if (strncmp(Parms->PermDipolesFile, "|", 1) == 0) {
		for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++)
			(*PermDipoles)[UnitCount] = cos(Parms->TwistAngle * UnitCount), (*PermDipoles)[Parms->UnitTotal + UnitCount] = sin(Parms->TwistAngle * UnitCount);
		if (fabs(Parms->TwistAngle != 0)) Pols->Active[1] = true;
	}
	
	else {
		FileID = fopen(Parms->PermDipolesFile, "r");
		for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++)
			for (PolCount = 0; PolCount < 3; PolCount++) {
				fscanf(FileID, "%lf ", &ReadDouble);
				(*PermDipoles)[PolCount * Parms->UnitTotal + UnitCount] = ReadDouble;
				if (ReadDouble) Pols->Active[PolCount] = true;
			}
		fclose(FileID);
	}
}


void BuildSOS(int **SOS, const char *SOSChar, const int Dimension)
{
	char Dummy[1000];
	int Count = 0;
	
	*SOS = malloc(Dimension * sizeof(int));
	
	sprintf(Dummy, "%s,|", SOSChar);
	
	if (strncmp(Dummy, ",|", 2) == 0) for (Count = 0; Count < Dimension; Count++) (*SOS)[Count] = Count;
	else {
		while (strncmp(Dummy, "|", 1) != 0) sscanf(Dummy, "%d,%s", &(*SOS)[Count++], Dummy);
		(*SOS)[Count] = -1;
	}
}


void BuildEnergyGrid(RealType **EnergyGrid, const ParmStruct *Parms)
{
	int Count;
	
	*EnergyGrid = malloc(Parms->SpectrumTotal * sizeof(RealType));
	
	for (Count = 0; Count < Parms->SpectrumTotal; Count++) (*EnergyGrid)[Count] = Parms->EnergyMin + (Parms->EnergyMax - Parms->EnergyMin) * (Count + 1) / (Parms->SpectrumTotal + 1) - Parms->XEnergyMean;
}


void BuildBasicHamiltonian(const char *FileName, RealType **BasicHamiltonian, int *VelocityMap, const RealType *XEnergies, const VOverlapStruct VOverlaps, const BasisStruct *Basis, const int *BasisBound, const ParmStruct *Parms)
{
	BasisStruct PriBasisState, SecBasisState;
	RealType EffectiveVOverlap;
	int PriBasisCount, SecBasisCount, PriUnitCount, SecUnitCount;
	
	PrintMessage("Constructing ", FileName, ".", true, Parms->Embarrassing);
	
	*BasicHamiltonian = calloc(IntPow2(*BasisBound), sizeof(RealType));
	
	for (PriBasisCount = 0; PriBasisCount < *BasisBound; PriBasisCount++) (*BasicHamiltonian)[PriBasisCount * (*BasisBound + 1)] = DiagonalEnergy(Basis[PriBasisCount], XEnergies, Parms);
	
	for (PriBasisCount = 0; PriBasisCount < BasisBound[3]; PriBasisCount++)
		for (PriUnitCount = 0; PriUnitCount < Parms->UnitTotal; PriUnitCount++) {
			
			memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
			if (!AnnihilateX(&PriBasisState, PriUnitCount)) continue;
			
			for (SecUnitCount = 0; SecUnitCount < Parms->UnitTotal; SecUnitCount++) {
				if (PriUnitCount == SecUnitCount) continue;
				
				memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
				if (!CreateX(&SecBasisState, SecUnitCount)) continue;
				
				for (SecBasisCount = 0; SecBasisCount < PriBasisCount; SecBasisCount++)
					if (IdentifyBasisState(&EffectiveVOverlap, SecBasisState, Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
						(*BasicHamiltonian)[PriBasisCount * *BasisBound + SecBasisCount] += XEnergies[PriUnitCount * Parms->UnitTotal + SecUnitCount] * EffectiveVOverlap;
						(*BasicHamiltonian)[SecBasisCount * *BasisBound + PriBasisCount] += XEnergies[PriUnitCount * Parms->UnitTotal + SecUnitCount] * EffectiveVOverlap;
						if (VelocityMap != NULL) VelocityMap[PriBasisCount * *BasisBound + SecBasisCount] = 2 * (PriUnitCount - SecUnitCount);
					}
			}
		}

	for (PriBasisCount = BasisBound[3]; PriBasisCount < *BasisBound; PriBasisCount++)
		for (PriUnitCount = 0; PriUnitCount < Parms->UnitTotal; PriUnitCount++) {
			
			memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
			if (AnnihilateE(&PriBasisState, PriUnitCount)) 
            {
				for (SecUnitCount = 0; SecUnitCount < Parms->UnitTotal; SecUnitCount++) 
                {
                    int separation = SignedSeparation(PriUnitCount, SecUnitCount, Parms->UnitTotal, Parms->Periodic);
					if (abs(separation) != 1) continue;
					
					memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
					if (!CreateE(&SecBasisState, SecUnitCount)) continue;
					
					for (SecBasisCount = 0; SecBasisCount < PriBasisCount; SecBasisCount++)
						if (IdentifyBasisState(&EffectiveVOverlap, SecBasisState, Basis[SecBasisCount], VOverlaps, Parms->VTotal)) 
                        {
                            if( Parms->MolPerCell > 1 )
                            {
				                if( ( separation > 0 && (PriUnitCount % Parms->MolPerCell == 0) ) || ( separation < 0 && (PriUnitCount % Parms->MolPerCell == 1) ) )
                                {
							        (*BasicHamiltonian)[PriBasisCount * *BasisBound + SecBasisCount] += Parms->ECouplingNNRight * EffectiveVOverlap;
							        (*BasicHamiltonian)[SecBasisCount * *BasisBound + PriBasisCount] += Parms->ECouplingNNRight * EffectiveVOverlap;
                                }
                                else
                                {
							        (*BasicHamiltonian)[PriBasisCount * *BasisBound + SecBasisCount] += Parms->ECouplingNNLeft * EffectiveVOverlap;
							        (*BasicHamiltonian)[SecBasisCount * *BasisBound + PriBasisCount] += Parms->ECouplingNNLeft * EffectiveVOverlap;
                                }
                            }
                            else
                            {
							    (*BasicHamiltonian)[PriBasisCount * *BasisBound + SecBasisCount] += Parms->ECouplingNN * EffectiveVOverlap;
							    (*BasicHamiltonian)[SecBasisCount * *BasisBound + PriBasisCount] += Parms->ECouplingNN * EffectiveVOverlap;
                            }
							if (VelocityMap != NULL) VelocityMap[PriBasisCount * *BasisBound + SecBasisCount] = PriUnitCount - SecUnitCount;
						}
				}
			}
			
			memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
			if (AnnihilateH(&PriBasisState, PriUnitCount)) {
				for (SecUnitCount = 0; SecUnitCount < Parms->UnitTotal; SecUnitCount++) {
					if (Separation(PriUnitCount, SecUnitCount, Parms->UnitTotal, Parms->Periodic) != 1) continue;
					
                    int separation = SignedSeparation(PriUnitCount, SecUnitCount, Parms->UnitTotal, Parms->Periodic);
					if (abs(separation) != 1) continue;
					
					memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
					if (!CreateH(&SecBasisState, SecUnitCount)) continue;
					
					for (SecBasisCount = 0; SecBasisCount < PriBasisCount; SecBasisCount++)
						if (IdentifyBasisState(&EffectiveVOverlap, SecBasisState, Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
                            if( Parms->MolPerCell > 1 )
                            {
				                if( ( separation > 0 && (PriUnitCount % Parms->MolPerCell == 0) ) || ( separation < 0 && (PriUnitCount % Parms->MolPerCell == 1) ) )
                                {
							        (*BasicHamiltonian)[PriBasisCount * *BasisBound + SecBasisCount] += Parms->HCouplingNNRight * EffectiveVOverlap;
							        (*BasicHamiltonian)[SecBasisCount * *BasisBound + PriBasisCount] += Parms->HCouplingNNRight * EffectiveVOverlap;
                                }
                                else
                                {
							        (*BasicHamiltonian)[PriBasisCount * *BasisBound + SecBasisCount] += Parms->HCouplingNNLeft * EffectiveVOverlap;
							        (*BasicHamiltonian)[SecBasisCount * *BasisBound + PriBasisCount] += Parms->HCouplingNNLeft * EffectiveVOverlap;
                                }
                            }
                            else
                            {
							    (*BasicHamiltonian)[PriBasisCount * *BasisBound + SecBasisCount] += Parms->HCouplingNN * EffectiveVOverlap;
							    (*BasicHamiltonian)[SecBasisCount * *BasisBound + PriBasisCount] += Parms->HCouplingNN * EffectiveVOverlap;
                            }
							if (VelocityMap != NULL) VelocityMap[PriBasisCount * *BasisBound + SecBasisCount] = PriUnitCount - SecUnitCount;
						}
				}
			}
		}
	
	SaveOperator(FileName, *BasicHamiltonian, Basis, Basis, *BasisBound, *BasisBound, Parms->Embarrassing);
}


void BuildDipoles(const char *FileName, RealType **Dipoles, int **DipolesMask, const BasisStruct *LoBasis, const BasisStruct *HiBasis, const int LoBasisTotal, const int HiBasisTotal, const VOverlapStruct VOverlaps, const ParmStruct *Parms)
{
	BasisStruct HiBasisState, LoBasisState;
	RealType EffectiveVOverlap;
	int LoBasisCount, HiBasisCount, UnitCount;
	
	// Start of playground!
	
	/*memcpy(&HiBasisState, &HiBasis[0], sizeof(BasisStruct));
	PrintBasisState(HiBasisState);
	
	AnnihilateX(&HiBasisState, 0);
	PrintBasisState(HiBasisState);
	
	PrintBasisState(LoBasis[0]);
	
	if (IdentifyBasisState(&EffectiveVOverlap, HiBasisState, LoBasis[0], VOverlaps, Parms->VTotal)) printf("%f\n", EffectiveVOverlap);
	
	printf("%f\n", VOverlaps.XG[0]);
	
	exit(0);*/
	
	// End of playground!
	
	PrintMessage("Constructing ", FileName, ".", true, Parms->Embarrassing);
	
	*Dipoles = calloc(HiBasisTotal * LoBasisTotal, sizeof(RealType));
	*DipolesMask = malloc(HiBasisTotal * LoBasisTotal * sizeof(int)),
	memset(*DipolesMask, -1, HiBasisTotal * LoBasisTotal * sizeof(int));
	
	for (HiBasisCount = 0; HiBasisCount < HiBasisTotal; HiBasisCount++)
		for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++) {
			memcpy(&HiBasisState, &HiBasis[HiBasisCount], sizeof(BasisStruct));
			if (!AnnihilateX(&HiBasisState, UnitCount)) continue;
			
			memcpy(&LoBasisState, &HiBasisState, sizeof(BasisStruct));
			for (LoBasisCount = 0; LoBasisCount < LoBasisTotal; LoBasisCount++)
				if (IdentifyBasisState(&EffectiveVOverlap, LoBasisState, LoBasis[LoBasisCount], VOverlaps, Parms->VTotal)) {
					// Are the two lines below consistent? [+= vs. =]
					(*Dipoles)[LoBasisCount * HiBasisTotal + HiBasisCount] += EffectiveVOverlap;
					(*DipolesMask)[LoBasisCount * HiBasisTotal + HiBasisCount] = UnitCount;
				}
		}
		
	SaveOperator(FileName, *Dipoles, LoBasis, HiBasis, LoBasisTotal, HiBasisTotal, Parms->Embarrassing);
}


void BuildMatrixOperator(const char *FileName, RealType **MatrixOperator, int **MatrixOperatorMask, const RealType *KetOperator, const RealType *BraOperator, const MatrixBasisStruct *PriMatrixBasis, const MatrixBasisStruct *SecMatrixBasis, const int PriKetBasisTotal, const int PriBraBasisTotal, const int SecKetBasisTotal, const int SecBraBasisTotal, const int Embarrassing)
{
	int PriCount, SecCount, PriDimension = PriKetBasisTotal * PriBraBasisTotal, SecDimension = SecKetBasisTotal * SecBraBasisTotal;
	
	PrintMessage("Constructing ", FileName, ".", true, Embarrassing);
	
	*MatrixOperator = calloc(PriDimension * SecDimension, sizeof(RealType));
	
	for (PriCount = 0; PriCount < PriDimension; PriCount++)
		for (SecCount = 0; SecCount < SecDimension; SecCount++) {
			if (PriMatrixBasis[PriCount].Bra == SecMatrixBasis[SecCount].Bra && KetOperator != NULL)
				(*MatrixOperator)[PriCount * SecDimension + SecCount] += KetOperator[PriMatrixBasis[PriCount].Ket * SecKetBasisTotal + SecMatrixBasis[SecCount].Ket];
			if (PriMatrixBasis[PriCount].Ket == SecMatrixBasis[SecCount].Ket && BraOperator != NULL)
				(*MatrixOperator)[PriCount * SecDimension + SecCount] -= BraOperator[PriMatrixBasis[PriCount].Bra * SecBraBasisTotal + SecMatrixBasis[SecCount].Bra];
		}
	
	if (MatrixOperatorMask != NULL)
		*MatrixOperatorMask = malloc(PriDimension * SecDimension * sizeof(int)),
		CreateMatrixMask(*MatrixOperatorMask, *MatrixOperator, PriDimension * SecDimension);
	
	SaveOperator(FileName, *MatrixOperator, NULL, NULL, PriDimension, SecDimension, Embarrassing);
}


void BuildGPropagator(ComplexType *GPropagator, const BasisStruct *GBasis, const ParmStruct *Parms)
{
	RealType PreCalc = InvCmToFreq * Parms->TimeStep;
	int Count;
	
	for (Count = 0; Count < *Parms->GBasisBound; Count++) GPropagator[Count] = cexp(-I * DiagonalEnergy(GBasis[Count], NULL, Parms) * PreCalc);
}


void PrepareTrotter(RealType *Matrix, ComplexType *Propagator, const int Dimension, const ParmStruct *Parms)
{
	int Count;
	
	RealType *MatrixCopy = malloc(IntPow2(Dimension) * sizeof(RealType)), *EigenValues = malloc(Dimension * sizeof(RealType));
	
	memcpy(MatrixCopy, Matrix, IntPow2(Dimension) * sizeof(RealType));
	
	for (Count = 0; Count < Dimension; Count++) Matrix[Count * (Dimension + 1)] = 0;
	
	DiagonalizeMatrix(Matrix, EigenValues, Dimension, Parms->DAndC);
	BuildPropagator(Propagator, Matrix, EigenValues, Parms, Dimension);
	
	memcpy(Matrix, MatrixCopy, IntPow2(Dimension) * sizeof(RealType));
	
	free(MatrixCopy); free(EigenValues);
}


void SetPolarizations(PolStruct *Pols)
{
	Pols->Weight = .2;
	if (Pols->Id < 0 && Pols->Active[0])	{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 0; Pols->Pulses[2] = 0; Pols->Pulses[3] = 0; Pols->Id = 0;  return; }
	if (Pols->Id < 1 && Pols->Active[1])	{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 1; Pols->Pulses[2] = 1; Pols->Pulses[3] = 1; Pols->Id = 1;  return; }
	if (Pols->Id < 2 && Pols->Active[2])	{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 2; Pols->Pulses[2] = 2; Pols->Pulses[3] = 2; Pols->Id = 2;  return; }
	
	Pols->Weight = .0666667;
	if (Pols->Active[0] && Pols->Active[1]) {
		if (Pols->Id < 3)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 0; Pols->Pulses[2] = 1; Pols->Pulses[3] = 1; Pols->Id = 3;  return; }
		if (Pols->Id < 4)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 1; Pols->Pulses[2] = 0; Pols->Pulses[3] = 0; Pols->Id = 4;  return; }
		if (Pols->Id < 5)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 1; Pols->Pulses[2] = 0; Pols->Pulses[3] = 1; Pols->Id = 5;  return; }
		if (Pols->Id < 6)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 0; Pols->Pulses[2] = 1; Pols->Pulses[3] = 0; Pols->Id = 6;  return; }
		if (Pols->Id < 7)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 1; Pols->Pulses[2] = 1; Pols->Pulses[3] = 0; Pols->Id = 7;  return; }
		if (Pols->Id < 8)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 0; Pols->Pulses[2] = 0; Pols->Pulses[3] = 1; Pols->Id = 8;  return; }
	}
	if (Pols->Active[0] && Pols->Active[2]) {
		if (Pols->Id < 9)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 0; Pols->Pulses[2] = 2; Pols->Pulses[3] = 2; Pols->Id = 9;  return; }
		if (Pols->Id < 10)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 2; Pols->Pulses[2] = 0; Pols->Pulses[3] = 0; Pols->Id = 10; return; }
		if (Pols->Id < 11)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 2; Pols->Pulses[2] = 0; Pols->Pulses[3] = 2; Pols->Id = 11; return; }
		if (Pols->Id < 12)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 0; Pols->Pulses[2] = 2; Pols->Pulses[3] = 0; Pols->Id = 12; return; }
		if (Pols->Id < 13)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 2; Pols->Pulses[2] = 2; Pols->Pulses[3] = 0; Pols->Id = 13; return; }
		if (Pols->Id < 14)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 0; Pols->Pulses[2] = 0; Pols->Pulses[3] = 2; Pols->Id = 14; return; }
	}
	if (Pols->Active[1] && Pols->Active[2]) {
		if (Pols->Id < 15)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 1; Pols->Pulses[2] = 2; Pols->Pulses[3] = 2; Pols->Id = 15; return; }
		if (Pols->Id < 16)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 2; Pols->Pulses[2] = 1; Pols->Pulses[3] = 1; Pols->Id = 16; return; }
		if (Pols->Id < 17)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 2; Pols->Pulses[2] = 1; Pols->Pulses[3] = 2; Pols->Id = 17; return; }
		if (Pols->Id < 18)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 1; Pols->Pulses[2] = 2; Pols->Pulses[3] = 1; Pols->Id = 18; return; }
		if (Pols->Id < 19)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 2; Pols->Pulses[2] = 2; Pols->Pulses[3] = 1; Pols->Id = 19; return; }
		if (Pols->Id < 20)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 1; Pols->Pulses[2] = 1; Pols->Pulses[3] = 2; Pols->Id = 20; return; }
	}
	Pols->Id = 21;
}


void InitBath(BathStruct *Bath, const ParmStruct *Parms)
{
	int Count;
	
	if (Parms->BathMass > 0) {
		for (Count = 0; Count < Parms->UnitTotal; Count++) {
			Bath[Count].BathPos = RandomGaussian(0, sqrt(Parms->Temperature / Parms->BathSpring));
			Bath[Count].BathVel = RandomGaussian(0, sqrt(Parms->Temperature / Parms->BathMass));
			Bath[Count].BathVelBra = Bath[Count].BathVel;
			Bath[Count].XBathDynamic = Parms->BathCoupling * Bath[Count].BathPos;
		}
	}
	else if (Parms->BathDynamicDev) {
		if (Parms->EHRadius)
			for (Count = 0; Count < Parms->UnitTotal; Count++) {
				Bath[Count].EBathDynamic = RandomGaussian(0, SqrtOfHalf * Parms->BathDynamicDev);
				Bath[Count].HBathDynamic = RandomGaussian(0, SqrtOfHalf * Parms->BathDynamicDev);
				Bath[Count].XBathDynamic = Bath[Count].EBathDynamic + Bath[Count].HBathDynamic;
			}
		else for (Count = 0; Count < Parms->UnitTotal; Count++) Bath[Count].XBathDynamic = RandomGaussian(0, Parms->BathDynamicDev);
	}
	
	if (Parms->BathStaticDev) {
		if (Parms->EHRadius)
			for (Count = 0; Count < Parms->UnitTotal; Count++) {
				Bath[Count].EBathStatic = RandomGaussian(0, SqrtOfHalf * Parms->BathStaticDev);
				Bath[Count].HBathStatic = RandomGaussian(0, SqrtOfHalf * Parms->BathStaticDev);
				Bath[Count].XBathStatic = Bath[Count].EBathStatic + Bath[Count].HBathStatic;
			}
		else for (Count = 0; Count < Parms->UnitTotal; Count++) Bath[Count].XBathStatic = RandomGaussian(0, Parms->BathStaticDev);
	}
}


void UpdateBath(BathStruct *Bath, const RealType *FeedbackVector, const RealType *FeedbackVectorBra, const BasisStruct *Basis, const ParmStruct *Parms, const int EnableVelBra)
{
	RealType *BathForce, *QuantumForce, *QuantumForceBra;
	int Count;
	
	if (!Parms->BathDynamics) return;
	
	if (Parms->BathMass > 0) {
		
		BathForce = malloc(Parms->UnitTotal * sizeof(RealType));
		
		for (Count = 0; Count < Parms->UnitTotal; Count++)
			BathForce[Count] = Parms->BathSpring * Bath[Count].BathPos + RandomGaussian(0, Parms->BathTempDev);
		
		BuildQuantumForce(&QuantumForce, FeedbackVector, Basis, Parms); BuildQuantumForce(&QuantumForceBra, FeedbackVectorBra, Basis, Parms);
		
		if (EnableVelBra)
			for (Count = 0; Count < Parms->UnitTotal; Count++) {
				Bath[Count].BathVel -= ((BathForce[Count] + QuantumForce[Count]) / Parms->BathMass + Parms->BathDamping * Bath[Count].BathVel) * Parms->TimeStep;
				Bath[Count].BathVelBra -= ((BathForce[Count] + QuantumForceBra[Count]) / Parms->BathMass + Parms->BathDamping * Bath[Count].BathVelBra) * Parms->TimeStep;
				Bath[Count].BathPos += .5 * (Bath[Count].BathVel + Bath[Count].BathVelBra) * Parms->TimeStep;
			}
		else
			for (Count = 0; Count < Parms->UnitTotal; Count++) {
				Bath[Count].BathVel -= ((BathForce[Count] + QuantumForce[Count]) / Parms->BathMass + Parms->BathDamping * Bath[Count].BathVel) * Parms->TimeStep;
				Bath[Count].BathPos += Bath[Count].BathVel * Parms->TimeStep;
			}
		
		/*if (FeedbackVectorBra == NULL) for (Count = 0; Count < Parms->UnitTotal; Count++) Bath[Count].BathPos += Bath[Count].BathVel * Parms->TimeStep;
		else
			for (Count = 0; Count < Parms->UnitTotal; Count++) {
				Bath[Count].BathVelBra -= ((Parms->BathSpring * Bath[Count].BathPos + QuantumForceBra[Count] + ThermalForce[Count]) / Parms->BathMass + Parms->BathDamping * Bath[Count].BathVelBra) * Parms->TimeStep;
				Bath[Count].BathPos += .5 * (Bath[Count].BathVel + Bath[Count].BathVelBra) * Parms->TimeStep;
			}*/
		
		free(BathForce); free(QuantumForce); free(QuantumForceBra);
		
		for (Count = 0; Count < Parms->UnitTotal; Count++) Bath[Count].XBathDynamic = Parms->BathCoupling * Bath[Count].BathPos;
		
	}
	else if (Parms->EHRadius)
		for (Count = 0; Count < Parms->UnitTotal; Count++) {
			Bath[Count].EBathDynamic = BrownianOscillator(Bath[Count].EBathDynamic, SqrtOfHalf * Parms->BathDynamicDev, Parms);
			Bath[Count].HBathDynamic = BrownianOscillator(Bath[Count].HBathDynamic, SqrtOfHalf * Parms->BathDynamicDev, Parms);
			Bath[Count].XBathDynamic = Bath[Count].EBathDynamic + Bath[Count].HBathDynamic;
		}
	else for (Count = 0; Count < Parms->UnitTotal; Count++) Bath[Count].XBathDynamic = BrownianOscillator(Bath[Count].XBathDynamic, Parms->BathDynamicDev, Parms);
}


void Initialize(RealType *Matrix, RealType *EigenValues, ComplexType *Propagator, const BathStruct *Bath, const RealType *BasicHamiltonian, const BasisStruct *Basis, const ParmStruct *Parms, const int Dimension, const int TrotterNow)
{
	if (TrotterNow) return;
	
	BuildHamiltonian(Matrix, BasicHamiltonian, Bath, Basis, Dimension);
	DiagonalizeMatrix(Matrix, EigenValues, Dimension, Parms->DAndC);
	BuildPropagator(Propagator, Matrix, EigenValues, Parms, Dimension);
}


void Update(RealType *Matrix, RealType *EigenValues, ComplexType *Propagator, const BathStruct *Bath, const RealType *BasicHamiltonian, const BasisStruct *Basis, const ParmStruct *Parms, const int Dimension, const int TrotterNow)
{	
	if (!Parms->BathDynamics) return;
	
	if (TrotterNow) BuildHamiltonianDiag(Matrix, BasicHamiltonian, Bath, Basis, Dimension);
	else {		
		BuildHamiltonian(Matrix, BasicHamiltonian, Bath, Basis, Dimension);
		DiagonalizeMatrix(Matrix, EigenValues, Dimension, Parms->DAndC);
		BuildPropagator(Propagator, Matrix, EigenValues, Parms, Dimension);
	}
}


void Propagate(ComplexType *Coef, ComplexType *CoefOld, const RealType *Matrix, const ComplexType *Propagator, const int Dimension, const ParmStruct *Parms, const int TrotterNow)
{
	if (TrotterNow) ApplyTrotter(Coef, Matrix, Dimension, Parms);
	
	memcpy(CoefOld, Coef, Dimension * sizeof(ComplexType));
	CMatrixTimesCVector(Coef, Propagator, CoefOld, Dimension);
	
	if (TrotterNow == 1) ApplyTrotter(Coef, Matrix, Dimension, Parms);
}


void BuildHamiltonian(RealType *Hamiltonian, const RealType *BasicHamiltonian, const BathStruct *Bath, const BasisStruct *Basis, const int Dimension)
{
	int Count;
	
	memcpy(Hamiltonian, BasicHamiltonian, IntPow2(Dimension) * sizeof(RealType));
	for (Count = 0; Count < Dimension; Count++) Hamiltonian[Count * (Dimension + 1)] += TotalBath(Basis[Count], Bath);
}


void UpdateFeedbackState(int *FeedbackState, BathStruct *Bath, RealType *OldEigenVectors, const ComplexType *Coef, const RealType *EigenVectors, const RealType *EigenValues, const BasisStruct *Basis, const ParmStruct *Parms, const int FeedbackMode)
{
	ComplexType *CoefAdiabat;
	RealType *NACoupling;
	
	const RealType *FeedbackVector = &EigenVectors[*FeedbackState * *Parms->SBasisBound], *TargetVector;
	RealType InnerProduct, Contribution, Discriminant, VelScaling, ProbMagn = 0, QuantityA = 0, QuantityB = 0;
	int TargetState, UnitCount, BasisCount, Bra = (FeedbackMode == FMBra || FeedbackMode == FMBraOnly);
	
	if (*OldEigenVectors == -2) {
		if (FeedbackMode != FMKet) memcpy(OldEigenVectors, EigenVectors, IntPow2(*Parms->SBasisBound) * sizeof(RealType));
		return;
	}
	
	CoefAdiabat = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	NACoupling = calloc(Parms->UnitTotal, sizeof(RealType));
	
	if (FeedbackMode != FMBra) RandomValue = RandomUniform();
	
	RMatrixTimesCVector(CoefAdiabat, EigenVectors, Coef, *Parms->SBasisBound);
	
	for (TargetState = 0; TargetState < *Parms->SBasisBound; TargetState++) {
		if (TargetState == *FeedbackState) continue;
		
		InnerProduct = RInnerProduct(&EigenVectors[TargetState * *Parms->SBasisBound], &OldEigenVectors[*FeedbackState * *Parms->SBasisBound], *Parms->SBasisBound);
		if (InnerProduct > .9) continue;
		
		if ((Contribution = 2 * creal(CoefAdiabat[TargetState] / CoefAdiabat[*FeedbackState]) * InnerProduct) > 0)
			if ((ProbMagn += Contribution) > RandomValue) {
				
				TargetVector = &EigenVectors[TargetState * *Parms->SBasisBound];
				for (BasisCount = 0; BasisCount < Parms->SBasisBound[3]; BasisCount++) NACoupling[Basis[BasisCount].X[0].U] += (*FeedbackVector++) * (*TargetVector++);
				
				for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++)
					QuantityA += RealPow2(NACoupling[UnitCount]),
					QuantityB += (Bra ? Bath[UnitCount].BathVelBra : Bath[UnitCount].BathVel) * NACoupling[UnitCount];
				
				QuantityA /= 2 * Parms->BathMass;
				Discriminant = RealPow2(QuantityB) + 4 * QuantityA * (EigenValues[*FeedbackState] - EigenValues[TargetState]) * InvCmToKJPMol;
				
				VelScaling = QuantityB;
				if (Discriminant < 0) VelScaling /= QuantityA;
				else {
					if (QuantityB < 0)	VelScaling += sqrt(Discriminant);
					else				VelScaling += -sqrt(Discriminant);
					VelScaling /= 2 * QuantityA;
					*FeedbackState = TargetState;
				}
				
				if (Bra) for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++) Bath[UnitCount].BathVelBra -= VelScaling * NACoupling[UnitCount] / Parms->BathMass;
				else for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++) Bath[UnitCount].BathVel -= VelScaling * NACoupling[UnitCount] / Parms->BathMass;
				
				break;
			}
	}
	
	if (FeedbackMode != FMKet) memcpy(OldEigenVectors, EigenVectors, IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	
	free(CoefAdiabat); free(NACoupling);
}


void UpdateSign(int *Sign, const RealType *EigenVectors, const RealType *OldEigenVectors, const int Feedback, const int Dimension)
{
	if (*OldEigenVectors == -2) return;
	if (RInnerProduct(&EigenVectors[Feedback * Dimension], &OldEigenVectors[Feedback * Dimension], Dimension) < .9) *Sign = -*Sign;
}


void DirectAbsorption(double *Absorption, const RealType *EnergyGrid, const int *SOS, const RealType *PermDipoles, const RealType *SEigenVectors, const RealType *SEigenValues, const RealType *GSDipoles, const int *GSDipolesMask, const PolStruct *Pols, const ParmStruct *Parms)
{
	RealType *LineStrengths = malloc(*Parms->SBasisBound * sizeof(RealType));
	int SpectrumCount, SCount;
	
	for (SCount = 0; SCount < *Parms->SBasisBound && SOS[SCount] != -1; SCount++)
		LineStrengths[SCount] = LineStrength(0, GSDipoles, GSDipolesMask, PermDipoles, &SEigenVectors[SOS[SCount] * *Parms->SBasisBound], Pols, Parms) * (SEigenValues[SOS[SCount]] + Parms->MonomerShift - Parms->CrystalShift);
	
	for (SpectrumCount = 0; SpectrumCount < Parms->SpectrumTotal; SpectrumCount++)
		for (SCount = 0; SCount < *Parms->SBasisBound && SOS[SCount] != -1; SCount++)
            if (Parms->BroadFunc == Gauss)
			    Absorption[SpectrumCount] += LineStrengths[SCount] * Gaussian(EnergyGrid[SpectrumCount] - SEigenValues[SOS[SCount]], Parms->HBroad);
            else if (Parms->BroadFunc == Lorentz)
			    Absorption[SpectrumCount] += LineStrengths[SCount] * Lorentzian(EnergyGrid[SpectrumCount] - SEigenValues[SOS[SCount]], Parms->HBroad);
	
	free(LineStrengths);
}


void DirectFluorescence(double *Fluorescence, const RealType *EnergyGrid, const int *SOS, const BasisStruct *GBasis, const RealType *PermDipoles, const RealType *SEigenVectors, const RealType *SEigenValues, const RealType *GSDipoles, const int *GSDipolesMask, const PolStruct *Pols, const ParmStruct *Parms)
{
	RealType *LineStrengths = malloc(Parms->VTotal * sizeof(RealType));
	RealType PartitionFunction = 0, BoltzmannFactor = 1;
	int SpectrumCount, GCount, SCount, VCount;
	
	if (Parms->Temperature != 0)
		for (SCount = 0; SCount < *Parms->SBasisBound && BoltzmannFactor > .000001; SCount++)
			PartitionFunction += BoltzmannFactor = exp(-(SEigenValues[SCount] - SEigenValues[0]) / (Parms->Temperature * KJPMolToInvCm));
	
	for (SCount = 0; SCount < *Parms->SBasisBound && SOS[SCount] != -1; SCount++) {
		
		FillRArray(LineStrengths, 0, Parms->VTotal);
		
		BoltzmannFactor = (Parms->Temperature == 0) ? (SOS[SCount] == 0) : exp(-(SEigenValues[SOS[SCount]] - SEigenValues[0]) / (Parms->Temperature * KJPMolToInvCm)) / PartitionFunction;
		
		if (BoltzmannFactor < .000001) break;
		
		for (GCount = 0; GCount < *Parms->GBasisBound; GCount++) {
			VCount = IntMax(GBasis[GCount].G[0].VG, 0) + IntMax(GBasis[GCount].G[1].VG, 0) + IntMax(GBasis[GCount].G[2].VG, 0);
			LineStrengths[VCount] += LineStrength(GCount, GSDipoles, GSDipolesMask, PermDipoles, &SEigenVectors[SOS[SCount] * *Parms->SBasisBound], Pols, Parms);
		}
		
		for (VCount = 0; VCount < Parms->VTotal; VCount++)
			for (SpectrumCount = 0; SpectrumCount < Parms->SpectrumTotal; SpectrumCount++)
                if (Parms->BroadFunc == Gauss)
				    Fluorescence[SpectrumCount] += RealPowN(EnergyGrid[SpectrumCount] + Parms->MonomerShift, 3) * BoltzmannFactor * LineStrengths[VCount] * Gaussian(EnergyGrid[SpectrumCount] - SEigenValues[SOS[SCount]] + VCount * Parms->VEnergy, Parms->HBroad);
                else if (Parms->BroadFunc == Lorentz)
				    Fluorescence[SpectrumCount] += BoltzmannFactor * LineStrengths[VCount] * Lorentzian(EnergyGrid[SpectrumCount] - SEigenValues[SOS[SCount]] + VCount * Parms->VEnergy, Parms->HBroad);
	}
	
	free(LineStrengths);
}


void InitGCoef(ComplexType *GCoef, const int Dimension)
{
	FillCArray(GCoef, 0, Dimension);
	GCoef[0] = 1;
}


void Excite(ComplexType *HiCoef, const ComplexType *LoCoef, const RealType *PermDipoles, const RealType *LoHiDipoles, const int *LoHiDipolesMask, const int HiTotal, const int LoTotal, const int HierarchyTotal)
{
	ComplexType WorkValue;
	int HierarchyCount, HiCount, LoCount, HiSmart, LoSmart;
	
	for (HierarchyCount = 0; HierarchyCount < HierarchyTotal; HierarchyCount++) {
		HiSmart = HierarchyCount * HiTotal;
		LoSmart = HierarchyCount * LoTotal;
		
		for (HiCount = 0; HiCount < HiTotal; HiCount++) {
			WorkValue = 0;
			for (LoCount = 0; LoCount < LoTotal; LoCount++)
				if (LoHiDipolesMask[LoCount * HiTotal + HiCount] != -1) WorkValue += (LoHiDipoles[LoCount * HiTotal + HiCount] * PermDipoles[LoHiDipolesMask[LoCount * HiTotal + HiCount]]) * LoCoef[LoSmart + LoCount];
			HiCoef[HiSmart + HiCount] = WorkValue;
		}
	}
}


void DeExcite(ComplexType *LoCoef, const ComplexType *HiCoef, const RealType *PermDipoles, const RealType *LoHiDipoles, const int *LoHiDipolesMask, const int LoTotal, const int HiTotal, const int HierarchyTotal)
{
	ComplexType WorkValue;
	int HierarchyCount, LoCount, HiCount, LoSmart, HiSmart;
	
	for (HierarchyCount = 0; HierarchyCount < HierarchyTotal; HierarchyCount++) {
		LoSmart = HierarchyCount * LoTotal;
		HiSmart = HierarchyCount * HiTotal;
		
		for (LoCount = 0; LoCount < LoTotal; LoCount++) {
			WorkValue = 0;
			for (HiCount = 0; HiCount < HiTotal; HiCount++)
				if (LoHiDipolesMask[LoCount * HiTotal + HiCount] != -1) WorkValue += (LoHiDipoles[LoCount * HiTotal + HiCount] * PermDipoles[LoHiDipolesMask[LoCount * HiTotal + HiCount]]) * HiCoef[HiSmart + HiCount];
			LoCoef[LoSmart + LoCount] = WorkValue;
		}
	}
}


ComplexType ExciteAdiabat(ComplexType *SCoef, const ComplexType *GCoef, const RealType *GSDipoles, const int *GSDipolesMask, const RealType *EigenVector, const ParmStruct *Parms)
{
	int GCount, SCount, SmartCount;
	ComplexType Factor = 0;
	
	for (SCount = 0; SCount < Parms->SBasisBound[3]; SCount++) {
		SmartCount = SCount * *Parms->GBasisBound;
		for (GCount = 0; GCount < *Parms->GBasisBound; GCount++)
			if (GSDipolesMask[SmartCount + GCount] != -1) Factor += GSDipoles[SmartCount + GCount] * EigenVector[SCount] * GCoef[GCount];
	}
	
	if (SCoef != NULL) RealToComplex(SCoef, EigenVector, *Parms->SBasisBound);
	
	return Factor;
}


void DeExciteAdiabat(ComplexType *GCoef, const RealType *EigenVector, const RealType *GSDipoles, const int *GSDipolesMask, const ParmStruct *Parms)
{
	int GCount, SCount;
	
	for (GCount = 0; GCount < *Parms->GBasisBound; GCount++) {
		GCoef[GCount] = 0;
		for (SCount = 0; SCount < Parms->SBasisBound[3]; SCount++)
			if (GSDipolesMask[SCount * *Parms->GBasisBound + GCount] != -1) GCoef[GCount] += GSDipoles[SCount * *Parms->GBasisBound + GCount] * EigenVector[SCount];
	}
}


void ExciteMatrixAdiabat(ComplexType *SSCoef, const RealType *EigenVectors, const int PriAdiabat, const int SecAdiabat, const int Dimension)
{
	int PriCount, SecCount;
	
	for (PriCount = 0; PriCount < Dimension; PriCount++)
		for (SecCount = 0; SecCount < Dimension; SecCount++)
			SSCoef[PriCount * Dimension + SecCount] = EigenVectors[PriAdiabat * Dimension + PriCount] * EigenVectors[SecAdiabat * Dimension + SecCount];
}


ComplexType DeExciteMatrixAdiabat(const ComplexType *SSCoef, const RealType *EigenVectors, const int PriAdiabat, const int SecAdiabat, const int Dimension)
{
	ComplexType Value = 0;
	int PriCount, SecCount;
	
	for (PriCount = 0; PriCount < Dimension; PriCount++)
		for (SecCount = 0; SecCount < Dimension; SecCount++)
			Value += SSCoef[PriCount * Dimension + SecCount] * EigenVectors[PriAdiabat * Dimension + PriCount] * EigenVectors[SecAdiabat * Dimension + SecCount];
	
	return Value;
}


void CalculateXPopulations(RealType *Populations, const ComplexType *SCoef, const BasisStruct *SBasis, const int *SBasisBound)
{
	const ComplexType *SCoefPos = SCoef;
	int Count;
	
	for (Count = 0; Count < SBasisBound[3]; Count++) Populations[SBasis[Count].X[0].U] += RealPow2(cabs(*SCoefPos++));
}


void SSCoefToXPopulations(RealType *Populations, const ComplexType *SSCoef, const BasisStruct *SBasis, const int *SBasisBound)
{
	int Count;
	
	for (Count = 0; Count < SBasisBound[3]; Count++) Populations[SBasis[Count].X[0].U] = creal(SSCoef[Count * (*SBasisBound + 1)]);
}


void CalculateAllPopulations(RealType *Populations, const ComplexType *SCoef, const BasisStruct *SBasis, const int *SBasisBound, const int UnitTotal)
{
	const ComplexType *SCoefPos = SCoef;
	int Count;
	
	for (Count = 0; Count < SBasisBound[3]; Count++) Populations[SBasis[Count].X[0].U * (UnitTotal + 1)] += RealPow2(cabs(*SCoefPos++));
	for (Count = SBasisBound[3]; Count < *SBasisBound; Count++) Populations[SBasis[Count].E[0].U * UnitTotal + SBasis[Count].H[0].U] += RealPow2(cabs(*SCoefPos++));
}


void SSCoefToAllPopulations(RealType *Populations, const ComplexType *SSCoef, const BasisStruct *SBasis, const int *SBasisBound, const int UnitTotal)
{
	int Count;
	
	for (Count = 0; Count < SBasisBound[3]; Count++) Populations[SBasis[Count].X[0].U * (UnitTotal + 1)] = creal(SSCoef[Count * (*SBasisBound + 1)]);
	for (Count = SBasisBound[3]; Count < *SBasisBound; Count++) Populations[SBasis[Count].E[0].U * UnitTotal + SBasis[Count].H[0].U] = creal(SSCoef[Count * (*SBasisBound + 1)]);
}


void Trap(RealType *PopulationTrapValue, ComplexType *SCoef, const BasisStruct *SBasis, const int *SBasisBound, const RealType TrapPart, const int TrapUnit)
{
	RealType SqrtOneMinusTrapPart = sqrt(1 - TrapPart);
	int Count;
	
	for (Count = 0; Count < SBasisBound[3]; Count++)
		if (SBasis[Count].X[0].U == TrapUnit) {
			(*PopulationTrapValue) += TrapPart * RealPow2(cabs(SCoef[Count]));
			SCoef[Count] *= SqrtOneMinusTrapPart;
		}
}


void CalculateVelocity(RealType *Velocity, const ComplexType *SCoef, const RealType *BasicHamiltonian, const int *VelocityMap, const int BasisTotal)
{
	RealType VelocityValue = 0;
	int PriCount, SecCount;
	
	for (PriCount = 0; PriCount < BasisTotal; PriCount++)
		for (SecCount = 0; SecCount < PriCount; SecCount++)
			VelocityValue += cimag(conj(SCoef[PriCount]) * SCoef[SecCount]) * BasicHamiltonian[PriCount * BasisTotal + SecCount] * VelocityMap[PriCount * BasisTotal + SecCount];
	*Velocity += VelocityValue;
}


void Propagate_HEOM(ComplexType *Coef, ComplexType *CoefWork, const RealType *Liouville, const int *LiouvilleMask, const int *Hieperator, const int *AntiHieperator, const int BasisTotal, const int *Hierarchy, const int *HierarchyCouplings, const ParmStruct *Parms)
{
	ComplexType *CoefPos, *PriCoefWorkPos, *SecCoefWorkPos, *TerCoefWorkPos, *QuaCoefWorkPos;
	int Iteration, Count, Dimension = BasisTotal * Parms->HierarchyTotal;
	
	CoefDerivative(CoefWork, Coef, Liouville, LiouvilleMask, Hieperator, AntiHieperator, BasisTotal, Hierarchy, HierarchyCouplings, Parms); // k1
	
	for (Iteration = 0; Iteration < 3; Iteration++) {
		
		CoefPos = Coef; PriCoefWorkPos = &CoefWork[(2 * Iteration + 1) * Dimension]; SecCoefWorkPos = &CoefWork[2 * Iteration * Dimension];
		for (Count = 0; Count < Dimension; Count++) (*PriCoefWorkPos++) = (*CoefPos++) + (.5 + .5 * (Iteration == 2)) * (*SecCoefWorkPos++); // y + c_I * k_I, where I = Iteration + 1
		
		CoefDerivative(&CoefWork[(2 * Iteration + 2) * Dimension], &CoefWork[(2 * Iteration + 1) * Dimension], Liouville, LiouvilleMask, Hieperator, AntiHieperator, BasisTotal, Hierarchy, HierarchyCouplings, Parms); // k_(I+1)
	}
	
	CoefPos = Coef; PriCoefWorkPos = CoefWork; SecCoefWorkPos = &CoefWork[2 * Dimension]; TerCoefWorkPos = &CoefWork[4 * Dimension]; QuaCoefWorkPos = &CoefWork[6 * Dimension];
	for (Count = 0; Count < Dimension; Count++) (*CoefPos++) += 0.16666666667 * ((*PriCoefWorkPos++) + 2 * (*SecCoefWorkPos++) + 2 * (*TerCoefWorkPos++) + (*QuaCoefWorkPos++));
}


void BuildHieperators(int **Hieperator, int **AntiHieperator, const MatrixBasisStruct *MatrixBasis, const int BasisTotal, const BasisStruct *KetBasis, const BasisStruct *BraBasis, const ParmStruct *Parms)
{
	BasisStruct BasisState;
	int UnitCount, BasisCount;
	
	PrintMessage("Constructing hierarchy operators.", "", ".", true, Parms->Embarrassing);
	
	*Hieperator = malloc(Parms->UnitTotal * BasisTotal * sizeof(int)); *AntiHieperator = malloc(Parms->UnitTotal * BasisTotal * sizeof(int));
	
	for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++)
		for (BasisCount = 0; BasisCount < BasisTotal; BasisCount++) {
			
			memcpy(&BasisState, &KetBasis[MatrixBasis[BasisCount].Ket], sizeof(BasisStruct));
			(*Hieperator)[UnitCount * BasisTotal + BasisCount] = AnnihilateE(&BasisState, UnitCount);
			(*Hieperator)[UnitCount * BasisTotal + BasisCount] += AnnihilateH(&BasisState, UnitCount);
			(*AntiHieperator)[UnitCount * BasisTotal + BasisCount] = (*Hieperator)[UnitCount * BasisTotal + BasisCount];
			
			memcpy(&BasisState, &BraBasis[MatrixBasis[BasisCount].Bra], sizeof(BasisStruct));
			(*Hieperator)[UnitCount * BasisTotal + BasisCount] -= AnnihilateE(&BasisState, UnitCount);
			(*Hieperator)[UnitCount * BasisTotal + BasisCount] -= AnnihilateH(&BasisState, UnitCount);
			memcpy(&BasisState, &BraBasis[MatrixBasis[BasisCount].Bra], sizeof(BasisStruct));
			(*AntiHieperator)[UnitCount * BasisTotal + BasisCount] += AnnihilateE(&BasisState, UnitCount);
			(*AntiHieperator)[UnitCount * BasisTotal + BasisCount] += AnnihilateH(&BasisState, UnitCount);
		}
}


void CalculateBoltzmannFactors(RealType **BoltzmannFactor, const RealType *Hamiltonian, const RealType Temperature, const int Dimension)
{
	RealType PartitionFunction = 0;
	int Count;
	
	*BoltzmannFactor = malloc(Dimension * sizeof(RealType));
	
	for (Count = 0; Count < Dimension; Count++) {
		(*BoltzmannFactor)[Count] = exp(-Hamiltonian[Count * (Dimension + 1)] / (Temperature + .001));
		PartitionFunction += (*BoltzmannFactor)[Count];
		printf("%f\n", Hamiltonian[Count * (Dimension + 1)]);
	}
	
	for (Count = 0; Count < Dimension; Count++) (*BoltzmannFactor)[Count] /= PartitionFunction;
}


/*void CalculateDOccupation(RealType *DOccupation, const RealType *DEigenState, const ParmStruct *DBasis, const ParmStruct *Parms)
{
	int BasisCount;
	
	int PriCount, SecCount;
	
	for (BasisCount = 0; BasisCount < Parms->DBasisTotal; BasisCount++)
		DOccupation[DBasis[BasisCount].X[0].U * Parms->UnitTotal, DBasis[BasisCount].X[1].U] += RealPow2(DEigenState[BasisCount]);
}*/


void SaveResponse2D(const char *FileName, const ComplexType *Response, const ParmStruct *Parms, const int MemSkip, const int Minus)
{
	int PriCount, SecCount;
	char ModFileName[20];
	FILE *FileID;
	
	CreateFileName(ModFileName, FileName, Parms->Embarrassing); FileID = fopen(ModFileName, "w");
	for (PriCount = 0; PriCount < Parms->ResponseTimeTotal; PriCount++)
		for (SecCount = 0; SecCount < Parms->ResponseTimeTotal; SecCount++)
			fprintf(FileID, "%lf %lf %lf %lf\n", PriCount * Parms->TimeStep, SecCount * Parms->TimeStep, Minus1Pow(1 + Minus) * cimag(Response[(PriCount * Parms->ResponseTimeTotal + SecCount) * MemSkip]), Minus1Pow(1 + Minus) * creal(Response[(PriCount * Parms->ResponseTimeTotal + SecCount) * MemSkip]));
	fclose(FileID);
}


////////////
// Locals //
////////////


static void BuildPropagator(ComplexType *Propagator, const RealType *EigenVectors, const RealType *EigenValues, const ParmStruct *Parms, const int Dimension)
{
	ComplexType Exponent;
	int StateCount, PriCount, SecCount, SmartCount;
	
	FillCArray(Propagator, 0, IntPow2(Dimension));
	
	for (StateCount = 0; StateCount < Dimension; StateCount++) {
		SmartCount = StateCount * Dimension;
		Exponent = cexp(-I * EigenValues[StateCount] * InvCmToFreq * Parms->TimeStep);
		for (PriCount = 0; PriCount < Dimension; PriCount++)
			for (SecCount = 0; SecCount < Dimension; SecCount++)
				Propagator[SecCount * Dimension + PriCount] += EigenVectors[SmartCount + PriCount] * Exponent * EigenVectors[SmartCount + SecCount];
	}
}


static RealType LineStrength(const int GStateId, const RealType *GSDipoles, const int *GSDipolesMask, const RealType *PermDipoles, const RealType *SEigenVector, const PolStruct *Pols, const ParmStruct *Parms)
{
	RealType Contribution, Output = 0;
	int PolCount, UnitCount;
	
	for (PolCount = 0; PolCount < 3; PolCount++) {
		if (!Pols->Active[PolCount]) continue;
		Contribution = 0;
		
		for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++)
        {
            Contribution += PermDipoles[PolCount * Parms->UnitTotal + UnitCount] * TransitionDipole(UnitCount, GStateId, GSDipoles, GSDipolesMask, SEigenVector, Parms);
            //printf("%f\n", PermDipoles[PolCount * Parms->UnitTotal + UnitCount]);
        }		

		Output += RealPow2(Contribution);
	}
	
	return Output;
}


static RealType TransitionDipole(const int UnitId, const int GStateId, const RealType *GSDipoles, const int *GSDipolesMask, const RealType *SEigenVector, const ParmStruct *Parms)
{
	RealType Output = 0;
	int SCount, SCountLBound, SCountUBound;
	
	if (GStateId == 0)								SCountUBound = Parms->SBasisBound[1], SCountLBound = 0; // Single particles
	else if (GStateId < Parms->GBasisBound[1])		SCountUBound = Parms->SBasisBound[2], SCountLBound = 0; // Single/two particles
	else											SCountUBound = Parms->SBasisBound[3], SCountLBound = Parms->SBasisBound[1]; // Three particles
	
	for (SCount = SCountLBound; SCount < SCountUBound; SCount++)
		if (GSDipolesMask[GStateId * *Parms->SBasisBound + SCount] == UnitId)
			Output += GSDipoles[GStateId * *Parms->SBasisBound + SCount] * SEigenVector[SCount];
	
	return Output;
}


static RealType BrownianOscillator(const RealType BathDynamicValue, const RealType Deviation, const ParmStruct *Parms)
{
	return Parms->BathPreCalc[0] * BathDynamicValue + Parms->BathPreCalc[1] * RandomGaussian(0, Deviation);
}


static void BuildHamiltonianDiag(RealType *Diag, const RealType *BasicHamiltonian, const BathStruct *Bath, const BasisStruct *Basis, const int Dimension)
{
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) Diag[Count] = BasicHamiltonian[Count * (Dimension + 1)] + TotalBath(Basis[Count], Bath);
}


static RealType DiagonalEnergy(const BasisStruct BasisState, const RealType *XEnergies, const ParmStruct *Parms)
{
	RealType ReturnValue;
	int Count, VSum = 0;
	
	for (Count = 0; Count < XTotal; Count++) if (BasisState.X[Count].U != -1) VSum += BasisState.X[Count].VX;
	for (Count = 0; Count < ETotal; Count++) if (BasisState.E[Count].U != -1) VSum += BasisState.E[Count].VE;
	for (Count = 0; Count < HTotal; Count++) if (BasisState.H[Count].U != -1) VSum += BasisState.H[Count].VH;
	for (Count = 0; Count < GTotal; Count++) if (BasisState.G[Count].U != -1) VSum += BasisState.G[Count].VG;
	
	ReturnValue = Parms->VEnergy * VSum;
	
	if (XEnergies == NULL) return ReturnValue;
	
	for (Count = 0; Count < XTotal; Count++)
		if (BasisState.X[Count].U != -1) ReturnValue += XEnergies[BasisState.X[Count].U * (Parms->UnitTotal + 1)];
	
	if (BasisState.E[0].U != -1) ReturnValue += Parms->EHEnergyInf - Parms->EHEnergyNN / Separation(BasisState.E[0].U, BasisState.H[0].U, Parms->UnitTotal, Parms->Periodic);
	
	return ReturnValue;
}


static RealType TotalBath(const BasisStruct BasisState, const BathStruct *Bath)
{
	RealType ReturnValue = 0; 
	
	if (BasisState.X[0].U != -1) ReturnValue += Bath[BasisState.X[0].U].XBathStatic + Bath[BasisState.X[0].U].XBathDynamic;
	if (BasisState.X[1].U != -1) ReturnValue += Bath[BasisState.X[1].U].XBathStatic + Bath[BasisState.X[1].U].XBathDynamic;
	if (BasisState.E[0].U != -1) ReturnValue += Bath[BasisState.E[0].U].EBathStatic + Bath[BasisState.E[0].U].EBathDynamic;
	if (BasisState.H[0].U != -1) ReturnValue += Bath[BasisState.H[0].U].HBathStatic + Bath[BasisState.H[0].U].HBathDynamic;
	
	return ReturnValue;
}


static void BuildQuantumForce(RealType **QuantumForce, const RealType *FeedbackVector, const BasisStruct *Basis, const ParmStruct *Parms)
{
	const RealType *FeedbackVectorElement = FeedbackVector;
	int Count;
	
	*QuantumForce = calloc(Parms->UnitTotal, sizeof(RealType));
	
	if (FeedbackVector == NULL) return;
	
	for (Count = 0; Count < Parms->SBasisBound[3]; Count++) (*QuantumForce)[Basis[Count].X[0].U] += RealPow2(*FeedbackVectorElement++);
	RealTimesRArray(*QuantumForce, Parms->BathCoupling * InvCmToKJPMol, Parms->UnitTotal);
}


static void ApplyTrotter(ComplexType *Coef, const RealType *HamiltonianDiag, const int Dimension, const ParmStruct *Parms)
{
	const RealType *HamiltonianDiagElement = &HamiltonianDiag[0];
	ComplexType *CoefElement = &Coef[0];
	RealType PreCalc = InvCmToFreq * Parms->TimeStep;
	int Count;
	
	if (Parms->Trotter == 1) PreCalc *= .5;
	
	for (Count = 0; Count < Dimension; Count++) (*CoefElement++) *= cexp(-I * PreCalc * (*HamiltonianDiagElement++));
}


static void CoefDerivative(ComplexType *CoefDerivative, const ComplexType *Coef, const RealType *Liouville, const int *LiouvilleMask, const int *Hieperator, const int *AntiHieperator, const int BasisTotal, const int *Hierarchy, const int *HierarchyCouplings, const ParmStruct *Parms)
{
	const ComplexType *CoefPos;
	ComplexType *CoefDerivativePos;
	ComplexType PreFactor;
	int HierarchyCount, UnitCount, BasisCount, PriStartId, SecStartId, SmartCalc;
	
	for (HierarchyCount = 0; HierarchyCount < Parms->HierarchyTotal; HierarchyCount++) {
		PriStartId = HierarchyCount * BasisTotal;
		
		if (Liouville != NULL) RMatrixTimesCVectorMasked(&CoefDerivative[PriStartId], Liouville, &Coef[PriStartId], LiouvilleMask, BasisTotal);
		else FillCArray(&CoefDerivative[PriStartId], 0, BasisTotal);
		
		CoefPos = &Coef[PriStartId]; CoefDerivativePos = &CoefDerivative[PriStartId];
		for (BasisCount = 0; BasisCount < BasisTotal; BasisCount++) {
			(*CoefDerivativePos) *= -I;
			(*CoefDerivativePos++) -= Parms->HierarchyPreCalc[0] * SumIArray(&Hierarchy[HierarchyCount * Parms->UnitTotal], Parms->UnitTotal) * (*CoefPos++);
		}
		
		// Inter-hierarchy mixing!
		
		if (Hieperator == NULL) continue;
		
		for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++) {
			
			SmartCalc = UnitCount * BasisTotal;
			
			SecStartId = HierarchyCouplings[HierarchyCount * 2 * Parms->UnitTotal + UnitCount];
			if (SecStartId != -1) {
				SecStartId *= BasisTotal;
				
				for (BasisCount = 0; BasisCount < BasisTotal; BasisCount++) {
					PreFactor = 0;
					if (Hieperator[SmartCalc + BasisCount]) PreFactor += I * Parms->HierarchyPreCalc[1] * Hieperator[SmartCalc + BasisCount];
					if (AntiHieperator[SmartCalc + BasisCount]) PreFactor += Parms->HierarchyPreCalc[2] * AntiHieperator[SmartCalc + BasisCount];
					if (PreFactor) CoefDerivative[PriStartId + BasisCount] += Hierarchy[HierarchyCount * Parms->UnitTotal + UnitCount] * PreFactor * Coef[SecStartId + BasisCount];
				}
			}
			
			SecStartId = HierarchyCouplings[(HierarchyCount * 2 + 1) * Parms->UnitTotal + UnitCount];
			if (SecStartId != -1) {
				SecStartId *= BasisTotal;
				
				for (BasisCount = 0; BasisCount < BasisTotal; BasisCount++)
					if (Hieperator[SmartCalc + BasisCount]) CoefDerivative[PriStartId + BasisCount] += I * .5 * Hieperator[SmartCalc + BasisCount] * Coef[SecStartId + BasisCount] * InvCmToFreq;
			}
		}
	}
}
