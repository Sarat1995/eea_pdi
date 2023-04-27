typedef struct {
	
	int Pulses[4];
	int Active[3];
	int Id;
	RealType Weight;
	
} PolStruct;


typedef struct {
	
	RealType XBathStatic, EBathStatic, HBathStatic, XBathDynamic, EBathDynamic, HBathDynamic, BathPos, BathVel, BathVelBra;
	
} BathStruct;


enum FeedbackModeList { FMSerial, FMKet, FMBra, FMKetOnly, FMBraOnly };


void BuildXEnergies(RealType **XEnergies, const ParmStruct *Parms);
void BuildPermDipoles(RealType **PermDipoles, PolStruct *Pols, const ParmStruct *Parms);
void BuildSOS(int **SOS, const char *SOSChar, const int Dimension);
void BuildEnergyGrid(RealType **EnergyGrid, const ParmStruct *Parms);
void BuildBasicHamiltonian(const char *FileName, RealType **BasicHamiltonian, int *VelocityMap, const RealType *XEnergies, const VOverlapStruct VOverlaps, const BasisStruct *Basis, const int *BasisBound, const ParmStruct *Parms);
void BuildDipoles(const char *FileName, RealType **Dipoles, int **DipolesMap, const BasisStruct *LoBasis, const BasisStruct *HiBasis, const int LoBasisTotal, const int HiBasisTotal, const VOverlapStruct VOverlaps, const ParmStruct *Parms);
void BuildMatrixOperator(const char *FileName, RealType **MatrixOperator, int **MatrixOperatorMask, const RealType *KetOperator, const RealType *BraOperator, const MatrixBasisStruct *PriMatrixBasis, const MatrixBasisStruct *SecMatrixBasis, const int PriKetBasisTotal, const int PriBraBasisTotal, const int SecKetBasisTotal, const int SecBraBasisTotal, const int Embarrassing);
void BuildGPropagator(ComplexType *GPropagator, const BasisStruct *GBasis, const ParmStruct *Parms);
void PrepareTrotter(RealType *Matrix, ComplexType *Propagator, const int Dimension, const ParmStruct *Parms);
void SetPolarizations(PolStruct *Pols);
void InitBath(BathStruct *Bath, const ParmStruct *Parms);
void UpdateBath(BathStruct *Bath, const RealType *FeedbackVector, const RealType *FeedbackVectorBra, const BasisStruct *Basis, const ParmStruct *Parms, const int EnableVelBra);
void Initialize(RealType *Matrix, RealType *EigenValues, ComplexType *Propagator, const BathStruct *Bath, const RealType *BasicHamiltonian, const BasisStruct *Basis, const ParmStruct *Parms, const int Dimension, const int TrotterNow);
void Update(RealType *Matrix, RealType *EigenValues, ComplexType *Propagator, const BathStruct *Bath, const RealType *BasicHamiltonian, const BasisStruct *Basis, const ParmStruct *Parms, const int Dimension, const int TrotterNow);
void Propagate(ComplexType *Coef, ComplexType *CoefOld, const RealType *Matrix, const ComplexType *Propagator, const int Dimension, const ParmStruct *Parms, const int TrotterNow);
void BuildHamiltonian(RealType *Hamiltonian, const RealType *BasicHamiltonian, const BathStruct *Bath, const BasisStruct *Basis, const int Dimension);
void UpdateFeedbackState(int *FeedbackState, BathStruct *Bath, RealType *OldEigenVectors, const ComplexType *Coef, const RealType *EigenVectors, const RealType *EigenValues, const BasisStruct *Basis, const ParmStruct *Parms, const int FeedbackMode);
void UpdateSign(int *Sign, const RealType *EigenVectors, const RealType *OldEigenVectors, const int Feedback, const int Dimension);
void DirectAbsorption(double *Absorption, const RealType *EnergyGrid, const int *SOS, const RealType *PermDipoles, const RealType *SEigenVectors, const RealType *SEigenValues, const RealType *GSDipoles, const int *GSDipolesMask, const PolStruct *Pols, const ParmStruct *Parms);
void DirectFluorescence(double *Fluorescence, const RealType *EnergyGrid, const int *SOS, const BasisStruct *GBasis, const RealType *PermDipoles, const RealType *SEigenVectors, const RealType *SEigenValues, const RealType *GSDipoles, const int *GSDipolesMask, const PolStruct *Pols, const ParmStruct *Parms);
void InitGCoef(ComplexType *GCoef, const int Dimension);
void Excite(ComplexType *HiCoef, const ComplexType *LoCoef, const RealType *PermDipoles, const RealType *LoHiDipoles, const int *LoHiDipolesMask, const int HiTotal, const int LoTotal, const int HierarchyTotal);
void DeExcite(ComplexType *LoCoef, const ComplexType *HiCoef, const RealType *PermDipoles, const RealType *LoHiDipoles, const int *LoHiDipolesMask, const int LoTotal, const int HiTotal, const int HierarchyTotal);
void DeExciteAdiabat(ComplexType *GCoef, const RealType *EigenVector, const RealType *GSDipoles, const int *GSDipolesMask, const ParmStruct *Parms);
void ExciteMatrixAdiabat(ComplexType *SSCoef, const RealType *EigenVectors, const int PriAdiabat, const int SecAdiabat, const int Dimension);
void CalculateXPopulations(RealType *Populations, const ComplexType *SCoef, const BasisStruct *SBasis, const int *SBasisBound);
void SSCoefToXPopulations(RealType *Populations, const ComplexType *SSCoef, const BasisStruct *SBasis, const int *SBasisBound);
void CalculateAllPopulations(RealType *Populations, const ComplexType *SCoef, const BasisStruct *SBasis, const int *SBasisBound, const int UnitTotal);
void SSCoefToAllPopulations(RealType *Populations, const ComplexType *SSCoef, const BasisStruct *SBasis, const int *SBasisBound, const int UnitTotal);
void Trap(RealType *PopulationTrapValue, ComplexType *SCoef, const BasisStruct *SBasis, const int *SBasisBound, const RealType TrapPart, const int TrapUnit);
void CalculateVelocity(RealType *Velocity, const ComplexType *SCoef, const RealType *BasicHamiltonian, const int *VelocityMap, const int BasisTotal);
void Propagate_HEOM(ComplexType *Coef, ComplexType *CoefWork, const RealType *Liouville, const int *LiouvilleMask, const int *Hieperator, const int *AntiHieperator, const int BasisTotal, const int *Hierarchy, const int *HierarchyCouplings, const ParmStruct *Parms);
//void BuildHieperators(int **Hieperator, int **AntiHieperator, const MatrixBasisStruct *MatrixBasis, const int BasisTotal, const BasisStruct *KetBasis, const BasisStruct *BraBasis, const ParmStruct *Parms);
void BuildHieperators(int **Hieperator, int **AntiHieperator, const MatrixBasisStruct *MatrixBasis, const int BasisTotal, const BasisStruct *KetBasis, const BasisStruct *BraBasis, const ParmStruct *Parms);
void CalculateBoltzmannFactors(RealType **BoltzmannFactor, const RealType *Hamiltonian, const RealType Temperature, const int Dimension);
void SaveResponse2D(const char *FileName, const ComplexType *Response, const ParmStruct *Parms, const int MemSkip, const int Minus);


ComplexType ExciteAdiabat(ComplexType *SCoef, const ComplexType *GCoef, const RealType *GSDipoles, const int *GSDipolesMask, const RealType *EigenVector, const ParmStruct *Parms);
ComplexType DeExciteMatrixAdiabat(const ComplexType *SSCoef, const RealType *EigenVectors, const int PriAdiabat, const int SecAdiabat, const int Dimension);