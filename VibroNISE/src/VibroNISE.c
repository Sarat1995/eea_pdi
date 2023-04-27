#include "VibroNISE.h"


int main(int ArgsTotal, char *Args[])
{
	ParmStruct *Parms = calloc(1, sizeof(ParmStruct));
	PolStruct *Pols = calloc(1, sizeof(PolStruct));
	VOverlapStruct VOverlaps;
	RealType *XEnergies, *PermDipoles;
	long int InitialTimer, FinalTimer;
	char Timing[10], DisplayNumber[10];
	
	
	if (ArgsTotal < 3) Parms->Embarrassing = 0;
	else sscanf(Args[2], "%d", &Parms->Embarrassing);
	
	InitializeLog(Parms->Embarrassing);
	PrintTitle(Parms->Embarrassing);
	
	if (Parms->Embarrassing)
		IntToString(DisplayNumber, Parms->Embarrassing),
		PrintMessage("Running embarrassingly parallel session ", DisplayNumber, ".", true, Parms->Embarrassing);
	
	
	InitializeParms(Parms);
	ReadInputFile(Args, Parms, NULL);
	DerivedParms(Parms);

	BuildXEnergies(&XEnergies, Parms);
	
	BuildPermDipoles(&PermDipoles, Pols, Parms);
	
	// Change back to -sqrt(...) !!!
	BuildOverlaps("VOverlapsXG", &VOverlaps.XG, RealSqrt(Parms->XHuangRhys), Parms->VTotal, Parms->Embarrassing);
	BuildOverlaps("VOverlapsEG", &VOverlaps.EG, RealSqrt(Parms->EHuangRhys), Parms->VTotal, Parms->Embarrassing);
	BuildOverlaps("VOverlapsHG", &VOverlaps.HG, RealSqrt(Parms->HHuangRhys), Parms->VTotal, Parms->Embarrassing);
	BuildOverlaps("VOverlapsXE", &VOverlaps.XE, RealSqrt(Parms->XHuangRhys) - RealSqrt(Parms->EHuangRhys), Parms->VTotal, Parms->Embarrassing);
	BuildOverlaps("VOverlapsXH", &VOverlaps.XH, RealSqrt(Parms->XHuangRhys) - RealSqrt(Parms->HHuangRhys), Parms->VTotal, Parms->Embarrassing);
	
	
	RandomInitialise(Parms->PriRandomSeed, Parms->SecRandomSeed);
	
	InitialTimer = SetTimer();
	PrintMessage("Starting calculations!", "", "", true, Parms->Embarrassing);
	
	
	if (Parms->Mode == ModeDebug)				Debug(Parms, XEnergies, PermDipoles, Pols, VOverlaps);
	else if (Parms->Mode == ModeDirect1D)		Direct1D(Parms, XEnergies, PermDipoles, Pols, VOverlaps);
	else if (Parms->Mode == ModeResponse1D) {
		if (Parms->HierarchyDepth)				Response1D_HEOM(Parms, XEnergies, PermDipoles, VOverlaps);
		else if (Parms->QFeedback)				Response1D_SH(Parms, XEnergies, PermDipoles, VOverlaps);
		else									Response1D(Parms, XEnergies, PermDipoles, Pols, VOverlaps);
	}
	else if (Parms->Mode == ModeResponse2D) {
		if (Parms->QFeedback)					Response2D_SH(Parms, XEnergies, VOverlaps);
		else									Response2D_HEOM(Parms, XEnergies, PermDipoles, VOverlaps);
	}
	else if (Parms->Mode == ModePopTransfer) {
		if (Parms->HierarchyDepth)				PopTransfer_HEOM(Parms, XEnergies, VOverlaps);
		else if (Parms->QFeedback)				PopTransfer_SH(Parms, XEnergies, VOverlaps);
		else									PopTransfer(Parms, XEnergies, VOverlaps);
	}
	else if (Parms->Mode == ModeTempCoh) {
		if (Parms->HierarchyDepth) {
			if (Parms->SiteRep)					TempCoh_HEOM_Site(Parms, XEnergies, VOverlaps);
			else								TempCoh_HEOM(Parms, XEnergies, VOverlaps);
		}
		else if (Parms->QFeedback) {
			if (Parms->SiteRep)					TempCoh_SH_Site(Parms, XEnergies, VOverlaps);
			else								TempCoh_SH(Parms, XEnergies, VOverlaps);
		}
		else if (Parms->SiteRep)				TempCoh_Site(Parms, XEnergies, VOverlaps);
		else									TempCoh(Parms, XEnergies, VOverlaps);
		//else									TempCoh_(Parms, XEnergies, PermDipoles, VOverlaps);
	}
	else if (Parms->Mode == ModeAnnihilation)	Annihilation(Parms, XEnergies, VOverlaps);
	
	
	FinalTimer = SetTimer();
	
	IntToString(Timing, FinalTimer - InitialTimer);
	PrintMessage("Finished in ", Timing, " seconds!", true, Parms->Embarrassing);
	
	
	free(XEnergies); free(PermDipoles);
	free(VOverlaps.XG); free(VOverlaps.EG); free(VOverlaps.HG); free(VOverlaps.XE); free(VOverlaps.XH);
	free(Parms); free(Pols);
	
	return 0;
}


static void Debug(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const PolStruct *Pols, const VOverlapStruct VOverlaps)
{
	int TimeCount, FeedbackKet, FeedbackBra, ExciteState;
	
	BasisStruct *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *OldEigenVectors, *SEigenValues, *PopulationsKet, *PopulationsBra;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	ComplexType *SPropagator, *SCoefKet, *SCoefBra, *SCoefOld;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(NULL, &SBasis, NULL, NULL, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType)); OldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	SCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	PopulationsKet = calloc(Parms->WaitTimeTotal * *Parms->SBasisBound, sizeof(RealType));
	PopulationsBra = calloc(Parms->WaitTimeTotal * *Parms->SBasisBound, sizeof(RealType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	
	ExciteState = *Parms->SBasisBound - 1;
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(Bath, Parms);
		Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
		FeedbackKet = ExciteState; FeedbackBra = 0;
		RealToComplex(SCoefKet, &SMatrix[ExciteState * *Parms->SBasisBound], *Parms->SBasisBound);
		RealToComplex(SCoefBra, &SMatrix[ExciteState * *Parms->SBasisBound], *Parms->SBasisBound);
		
		OldEigenVectors[0] = -2;
		
		for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
			
			PopulationsKet[TimeCount * *Parms->SBasisBound + FeedbackKet]++;
			PopulationsBra[TimeCount * *Parms->SBasisBound + FeedbackBra]++;
			
			UpdateBath(Bath, &SMatrix[FeedbackKet * *Parms->SBasisBound], &SMatrix[FeedbackBra * *Parms->SBasisBound], SBasis, Parms, true);
			Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
			Propagate(SCoefKet, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
			Propagate(SCoefBra, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
			UpdateFeedbackState(&FeedbackKet, Bath, OldEigenVectors, SCoefKet, SMatrix, SEigenValues, SBasis, Parms, FMKet);
			UpdateFeedbackState(&FeedbackBra, Bath, OldEigenVectors, SCoefBra, SMatrix, SEigenValues, SBasis, Parms, FMBra);
		}
	}
	
	PrintCArray(SCoefKet, 1, *Parms->SBasisBound);
	PrintCArray(SCoefBra, 1, *Parms->SBasisBound);
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesRArray(PopulationsKet, 1.0 / Parms->SampleTotal, *Parms->SBasisBound * Parms->WaitTimeTotal);
	RealTimesRArray(PopulationsBra, 1.0 / Parms->SampleTotal, *Parms->SBasisBound * Parms->WaitTimeTotal);
	SaveRArray("PopulationsKet", PopulationsKet, Parms->WaitTimeTotal, *Parms->SBasisBound, TimeGrid, Parms->Embarrassing);
	SaveRArray("PopulationsBra", PopulationsBra, Parms->WaitTimeTotal, *Parms->SBasisBound, TimeGrid, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(Bath);
	free(SBasis);
	free(BasicSHamiltonian); free(SMatrix); free(OldEigenVectors); free(SEigenValues); free(PopulationsKet); free(PopulationsBra); free(TimeGrid);
	free(SPropagator); free(SCoefKet); free(SCoefBra); free(SCoefOld);
}


static void Direct1D(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const PolStruct *Pols, const VOverlapStruct VOverlaps)
{
	BasisStruct *GBasis, *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *SEigenValues, *GSDipoles, *EnergyGrid;
	
	double *Absorption = calloc(Parms->SpectrumTotal, sizeof(double)), *Fluorescence = calloc(Parms->SpectrumTotal, sizeof(double));
	
	int *GSDipolesMask, *SOS;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(&GBasis, &SBasis, NULL, Parms->GBasisBound, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildDipoles("GSDipoles", &GSDipoles, &GSDipolesMask, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildSOS(&SOS, Parms->SOSChar, *Parms->SBasisBound);
	BuildEnergyGrid(&EnergyGrid, Parms);
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(Bath, Parms);
		BuildHamiltonian(SMatrix, BasicSHamiltonian, Bath, SBasis, *Parms->SBasisBound);

	    SaveRArray("SHamiltonianArray", SMatrix, *Parms->SBasisBound, *Parms->SBasisBound, NULL, Parms->Embarrassing);

		DiagonalizeMatrix(SMatrix, SEigenValues, *Parms->SBasisBound, Parms->DAndC);

	    SaveRArray("SEigenVectors", SMatrix, *Parms->SBasisBound, *Parms->SBasisBound, NULL, Parms->Embarrassing);
	    SaveRArray("SEigenValues", SEigenValues, *Parms->SBasisBound, 1, NULL, Parms->Embarrassing);

		DirectAbsorption(Absorption, EnergyGrid, SOS, PermDipoles, SMatrix, SEigenValues, GSDipoles, GSDipolesMask, Pols, Parms);
		DirectFluorescence(Fluorescence, EnergyGrid, SOS, GBasis, PermDipoles, SMatrix, SEigenValues, GSDipoles, GSDipolesMask, Pols, Parms);
	}
	
	RealPlusRArray(EnergyGrid, Parms->XEnergyMean, Parms->SpectrumTotal);
	NormalizeDArray(Absorption, Parms->SpectrumTotal);
	SaveDVector("DirectAbs", Absorption, EnergyGrid, Parms->SpectrumTotal, Parms->Embarrassing);
	NormalizeDArray(Fluorescence, Parms->SpectrumTotal);
	SaveDVector("DirectFlu", Fluorescence, EnergyGrid, Parms->SpectrumTotal, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(GBasis); free(SBasis);
	free(Bath);
	free(BasicSHamiltonian); free(SMatrix); free(SEigenValues); free(GSDipoles); free(EnergyGrid);
	free(Absorption); free(Fluorescence);
	free(GSDipolesMask); free(SOS);
}


static void Response1D(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, PolStruct *Pols, const VOverlapStruct VOverlaps)
{
	int TimeCount;
	
	BasisStruct *GBasis, *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *SEigenValues, *GSDipoles;
	
	RealType *TimeGrid = malloc(Parms->ResponseTimeTotal * sizeof(RealType));
	
	ComplexType *SPropagator, *GCoef, *SCoef, *SCoefOld;
	
	double complex *Response = calloc(Parms->ResponseTimeTotal, sizeof(double complex));
	
	int *GSDipolesMask;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(&GBasis, &SBasis, NULL, Parms->GBasisBound, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildDipoles("GSDipoles", &GSDipoles, &GSDipolesMask, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	GCoef = malloc(*Parms->GBasisBound * sizeof(ComplexType));
	SCoef = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	if (Parms->Trotter) PrepareTrotter(BasicSHamiltonian, SPropagator, *Parms->SBasisBound, Parms);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		Pols->Id = -1;
		while (true) {
			SetPolarizations(Pols);
			
			if (Pols->Id > 2) break;
			InitBath(Bath, Parms);
			Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, Parms->Trotter);
			
			InitGCoef(GCoef, *Parms->GBasisBound);
			Excite(SCoef, GCoef, &PermDipoles[Pols->Pulses[0] * Parms->UnitTotal], GSDipoles, GSDipolesMask, *Parms->SBasisBound, *Parms->GBasisBound, 1);
			
			for (TimeCount = 0; TimeCount < Parms->ResponseTimeTotal; TimeCount++) {
				DeExcite(GCoef, SCoef, &PermDipoles[Pols->Pulses[1] * Parms->UnitTotal], GSDipoles, GSDipolesMask, *Parms->GBasisBound, *Parms->SBasisBound, 1);
				Response[TimeCount] += GCoef[0];
				
				UpdateBath(Bath, NULL, NULL, NULL, Parms, false);
				Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, Parms->Trotter);
				Propagate(SCoef, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, Parms->Trotter);
			}
		}
	}
	
	CreateGrid(TimeGrid, 0, (Parms->ResponseTimeTotal - 1) * Parms->TimeStep, Parms->ResponseTimeTotal);
	RealTimesDCArray(Response, 1.0 / Parms->SampleTotal, Parms->ResponseTimeTotal);
	SaveDCVector("Response1D", Response, TimeGrid, Parms->ResponseTimeTotal, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(GBasis); free(SBasis);
	free(Bath);
	free(BasicSHamiltonian); free(SMatrix); free(SEigenValues); free(GSDipoles);
	free(SPropagator); free(GCoef); free(SCoef); free(SCoefOld);
	free(Response);
	free(GSDipolesMask);
}


static void Response1D_SH(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const VOverlapStruct VOverlaps)
{
	int TimeCount, SCount, Feedback, Sign;
	
	BasisStruct *GBasis, *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *OldEigenVectors, *SEigenValues, *GSDipoles;
	
	RealType *TimeGrid = malloc(Parms->ResponseTimeTotal * sizeof(RealType));
	
	RealType PreCalc = InvCmToFreq * Parms->TimeStep;
	
	ComplexType *SPropagator, *GCoef, *SCoef, *SCoefOld;
	
	ComplexType Factor;
	
	double complex *Response = calloc(Parms->ResponseTimeTotal, sizeof(double complex));
	
	int *GSDipolesMask;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(&GBasis, &SBasis, NULL, Parms->GBasisBound, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildDipoles("GSDipoles", &GSDipoles, &GSDipolesMask, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType)); OldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	GCoef = malloc(*Parms->GBasisBound * sizeof(ComplexType));
	SCoef = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	if (Parms->Trotter) PrepareTrotter(BasicSHamiltonian, SPropagator, *Parms->SBasisBound, Parms);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		for (SCount = 0; SCount < *Parms->SBasisBound; SCount++) {
			
			InitBath(Bath, Parms);
			Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, Parms->Trotter);
			
			InitGCoef(GCoef, *Parms->GBasisBound);
			Factor = ExciteAdiabat(SCoef, GCoef, GSDipoles, GSDipolesMask, &SMatrix[SCount * *Parms->SBasisBound], Parms);
			Feedback = SCount;
			
			Sign = 1; *OldEigenVectors = -2;
			
			for (TimeCount = 0; TimeCount < Parms->ResponseTimeTotal; TimeCount++) {
				DeExciteAdiabat(GCoef, &SMatrix[Feedback * *Parms->SBasisBound], GSDipoles, GSDipolesMask, Parms);
				Response[TimeCount] += Sign * Factor * GCoef[0];
				
				UpdateBath(Bath, &SMatrix[Feedback * *Parms->SBasisBound], NULL, SBasis, Parms, true);
				Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
				UpdateSign(&Sign, SMatrix, OldEigenVectors, Feedback, *Parms->SBasisBound);
				Propagate(SCoef, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
				Factor *= cexp(-I * SEigenValues[Feedback] * PreCalc);
				UpdateFeedbackState(&Feedback, Bath, OldEigenVectors, SCoef, SMatrix, SEigenValues, SBasis, Parms, FMKetOnly);
			}
		}
	}
	
	CreateGrid(TimeGrid, 0, (Parms->ResponseTimeTotal - 1) * Parms->TimeStep, Parms->ResponseTimeTotal);
	RealTimesDCArray(Response, 1.0 / Parms->SampleTotal, Parms->ResponseTimeTotal);
	SaveDCVector("Response1D", Response, TimeGrid, Parms->ResponseTimeTotal, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(GBasis); free(SBasis);
	free(Bath);
	free(BasicSHamiltonian); free(SMatrix); free(OldEigenVectors); free(SEigenValues); free(GSDipoles);
	free(SPropagator); free(GCoef); free(SCoef); free(SCoefOld);
	free(Response);
	free(GSDipolesMask);
}


static void Response1D_HEOM(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const VOverlapStruct VOverlaps)
{
	int TimeCount, GGBasisTotal, SGBasisTotal;
	
	BasisStruct *GBasis, *SBasis;
	
	MatrixBasisStruct *GGBasis, *SGBasis;
	
	RealType *BasicSHamiltonian, *SGLiouville, *GSDipoles, *GGSGDipoles;
	
	RealType *TimeGrid = malloc(Parms->ResponseTimeTotal * sizeof(RealType));
	
	ComplexType *GGCoef, *SGCoef, *SGCoefWork;
	
	ComplexType *Response = malloc(Parms->ResponseTimeTotal * sizeof(ComplexType));
	
	int *Hierarchy, *HierarchyCouplings, *SGLiouvilleMask, *GSDipolesMask, *GGSGDipolesMask, *SGHieperator, *SGAntiHieperator;
	
	
	BuildHierarchy(&Hierarchy, &HierarchyCouplings, &Parms->HierarchyTotal, Parms->UnitTotal, Parms->HierarchyDepth, Parms->Embarrassing);
	
	BuildBases(&GBasis, &SBasis, NULL, Parms->GBasisBound, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	BuildMatrixBasis("GGBasis", &GGBasis, &GGBasisTotal, GBasis, GBasis, *Parms->GBasisBound, *Parms->GBasisBound, Parms->Embarrassing);
	BuildMatrixBasis("SGBasis", &SGBasis, &SGBasisTotal, SBasis, GBasis, *Parms->SBasisBound, *Parms->GBasisBound, Parms->Embarrassing);
	
	BuildHieperators(&SGHieperator, &SGAntiHieperator, SGBasis, SGBasisTotal, SBasis, GBasis, Parms);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildMatrixOperator("SGLiouville", &SGLiouville, &SGLiouvilleMask, BasicSHamiltonian, NULL, SGBasis, SGBasis, *Parms->SBasisBound, *Parms->GBasisBound, *Parms->SBasisBound, *Parms->GBasisBound, Parms->Embarrassing);
	RealTimesRArray(SGLiouville, InvCmToFreq * Parms->TimeStep, IntPow2(SGBasisTotal));
	
	BuildDipoles("GSDipoles", &GSDipoles, &GSDipolesMask, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildMatrixOperator("GGSGDipoles", &GGSGDipoles, &GGSGDipolesMask, GSDipoles, NULL, GGBasis, SGBasis, *Parms->GBasisBound, *Parms->GBasisBound, *Parms->SBasisBound, *Parms->GBasisBound, Parms->Embarrassing);
	
	
	GGCoef = calloc(GGBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType));
	SGCoef = calloc(SGBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType)); SGCoefWork = malloc(7 * SGBasisTotal * Parms->HierarchyTotal * sizeof(ComplexType));
	// Change calloc into malloc for SGCoef?
	
	GGCoef[0] = 1;
	Excite(SGCoef, GGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, SGBasisTotal, GGBasisTotal, 1);
	
	for (TimeCount = 0; TimeCount < Parms->ResponseTimeTotal; TimeCount++) {
		DeExcite(GGCoef, SGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, GGBasisTotal, SGBasisTotal, 1);
		Response[TimeCount] = GGCoef[0];
		Propagate_HEOM(SGCoef, SGCoefWork, SGLiouville, SGLiouvilleMask, SGHieperator, SGAntiHieperator, SGBasisTotal, Hierarchy, HierarchyCouplings, Parms);
	}
	
	CreateGrid(TimeGrid, 0, (Parms->ResponseTimeTotal - 1) * Parms->TimeStep, Parms->ResponseTimeTotal);
	SaveCVector("Response1D", Response, TimeGrid, Parms->ResponseTimeTotal, Parms->Embarrassing);
	
	
	free(GBasis); free(SBasis); free(GGBasis); free(SGBasis);
	free(BasicSHamiltonian); free(SGLiouville); free(GSDipoles); free(GGSGDipoles); free(TimeGrid);
	free(GGCoef); free(SGCoef); free(SGCoefWork); free(Response);
	free(Hierarchy); free(HierarchyCouplings); free(SGLiouvilleMask); free(GGSGDipolesMask); free(SGHieperator); free(SGAntiHieperator); free(GSDipolesMask);
}


static void Response2D_SH(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	typedef struct { ComplexType GBRR, GBNR, SERR, SENR, EARR, EANR; } ResponseStruct;
	
	int PriTimeCount, SecTimeCount, WaitTimeCount, PriSCount, SecSCount, Feedback, FeedbackKet, FeedbackBra, PriFeedback, Sign, PriSign;
	
	ResponseStruct *Response = malloc(IntPow2(Parms->ResponseTimeTotal) * sizeof(ResponseStruct));
	
	BasisStruct *GBasis, *SBasis, *DBasis;
	
	BathStruct *Bath = malloc(Parms->UnitTotal * sizeof(BathStruct)), *PriBath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	BathStruct *WaitBath = malloc(Parms->UnitTotal * sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *BasicDHamiltonian, *SMatrix, *DMatrix, *OldEigenVectors, *PriOldEigenVectors, *WaitOldEigenVectors, *SEigenValues, *DEigenValues;
	RealType *GSDipoles, *SDDipoles;
	
	RealType PreCalc = InvCmToFreq * Parms->TimeStep;
	
	ComplexType *GPropagator, *SPropagator, *DPropagator, *GCoef, *WaitSecGCoef, *SCoefKet, *SCoefBra, *SCoef, *PriSCoef, *SCoefOld;
	
	ComplexType Factor, PriFactor, WaitFactor;
	
	int *GSDipolesMask, *SDDipolesMask;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(&GBasis, &SBasis, &DBasis, Parms->GBasisBound, Parms->SBasisBound, Parms->DBasisBound, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildBasicHamiltonian("DHamiltonian", &BasicDHamiltonian, NULL, XEnergies, VOverlaps, DBasis, Parms->DBasisBound, Parms);
	BuildDipoles("GSDipoles", &GSDipoles, &GSDipolesMask, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildDipoles("SDDipoles", &SDDipoles, &SDDipolesMask, SBasis, DBasis, *Parms->SBasisBound, *Parms->DBasisBound, VOverlaps, Parms);
	
	GPropagator = (ComplexType *) malloc(*Parms->GBasisBound * sizeof(ComplexType));
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType)); SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	OldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType)); PriOldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	WaitOldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	
	DMatrix = malloc(IntPow2(*Parms->DBasisBound) * sizeof(RealType));
	DEigenValues = malloc(*Parms->DBasisBound * sizeof(RealType));
	DPropagator = malloc(IntPow2(*Parms->DBasisBound) * sizeof(ComplexType));
	
	GCoef = malloc(*Parms->GBasisBound * sizeof(ComplexType)); WaitSecGCoef = malloc(*Parms->GBasisBound * sizeof(ComplexType));
	SCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoef = malloc(*Parms->SBasisBound * sizeof(ComplexType)); PriSCoef = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	
	
	if (*Parms->DBasisBound > *Parms->SBasisBound)	InitializeDiag(BasicDHamiltonian, NULL, *Parms->DBasisBound, Parms->DAndC);
	else											InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	
	BuildGPropagator(GPropagator, GBasis, Parms);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		for (PriSCount = 0; PriSCount < *Parms->SBasisBound; PriSCount++)
			for (SecSCount = 0; SecSCount < *Parms->SBasisBound; SecSCount++) {
				
				InitBath(PriBath, Parms);
				Initialize(SMatrix, SEigenValues, SPropagator, PriBath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
				InitGCoef(GCoef, *Parms->GBasisBound);
				PriFactor = ExciteAdiabat(PriSCoef, GCoef, GSDipoles, GSDipolesMask, &SMatrix[PriSCount * *Parms->SBasisBound], Parms);
				
				PriFeedback = PriSCount; PriSign = 1; *PriOldEigenVectors = -2;
				
				for (PriTimeCount = 0; PriTimeCount < Parms->ResponseTimeTotal; PriTimeCount++) {
					
					DeExciteAdiabat(WaitSecGCoef, &SMatrix[PriFeedback * *Parms->SBasisBound], GSDipoles, GSDipolesMask, Parms);
					memcpy(WaitBath, PriBath, Parms->UnitTotal * sizeof(BathStruct));
					
					for (WaitTimeCount = 0; WaitTimeCount < Parms->WaitTimeTotal; WaitTimeCount++)
						UpdateBath(WaitBath, NULL, NULL, NULL, Parms, true),
						FactorCArrays(WaitSecGCoef, GPropagator, *Parms->GBasisBound);
					
					memcpy(Bath, WaitBath, Parms->UnitTotal * sizeof(BathStruct));
					Factor = PriFactor * ExciteAdiabat(SCoef, WaitSecGCoef, GSDipoles, GSDipolesMask, &SMatrix[SecSCount * *Parms->SBasisBound], Parms);
					
					Feedback = SecSCount; Sign = 1; *OldEigenVectors = -2;
					
					for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
						DeExciteAdiabat(GCoef, &SMatrix[Feedback * *Parms->SBasisBound], GSDipoles, GSDipolesMask, Parms);
						Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].GBNR += (PriSign * Sign) * Factor * GCoef[0];
						
						UpdateBath(Bath, &SMatrix[Feedback * *Parms->SBasisBound], NULL, SBasis, Parms, true);
						Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
						UpdateSign(&Sign, SMatrix, OldEigenVectors, Feedback, *Parms->SBasisBound);
						//if (*OldEigenVectors != -2  && RInnerProduct(&SMatrix[PriFeedback * *Parms->SBasisBound], &OldEigenVectors[PriFeedback * *Parms->SBasisBound], *Parms->SBasisBound) < .9) PriSign = -PriSign;
						Propagate(SCoef, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
						Factor *= cexp(-I * SEigenValues[Feedback] * PreCalc);
						UpdateFeedbackState(&Feedback, Bath, OldEigenVectors, SCoef, SMatrix, SEigenValues, SBasis, Parms, FMKetOnly);
					}
					
					memcpy(Bath, WaitBath, Parms->UnitTotal * sizeof(BathStruct));
					InitGCoef(GCoef, *Parms->GBasisBound);
					Factor = PriFactor * ExciteAdiabat(SCoef, GCoef, GSDipoles, GSDipolesMask, &SMatrix[SecSCount * *Parms->SBasisBound], Parms);
					
					Feedback = SecSCount; Sign = 1; *OldEigenVectors = -2;
					
					for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
						DeExciteAdiabat(GCoef, &SMatrix[Feedback * *Parms->SBasisBound], GSDipoles, GSDipolesMask, Parms);
						Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].GBRR += (PriSign * Sign) * Factor * CInnerProduct(WaitSecGCoef, GCoef, *Parms->GBasisBound);
						
						UpdateBath(Bath, NULL, &SMatrix[Feedback * *Parms->SBasisBound], SBasis, Parms, true);
						Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
						UpdateSign(&Sign, SMatrix, OldEigenVectors, Feedback, *Parms->SBasisBound);
						Propagate(SCoef, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
						Factor *= cexp(I * SEigenValues[Feedback] * PreCalc);
						UpdateFeedbackState(&Feedback, Bath, OldEigenVectors, SCoef, SMatrix, SEigenValues, SBasis, Parms, FMBraOnly);
						FactorCArrays(WaitSecGCoef, GPropagator, *Parms->GBasisBound);
					}
					
					memcpy(WaitBath, PriBath, Parms->UnitTotal * sizeof(BathStruct));
					memcpy(SCoefKet, PriSCoef, *Parms->SBasisBound);
					InitGCoef(GCoef, *Parms->GBasisBound);
					WaitFactor = PriFactor * ExciteAdiabat(SCoefBra, GCoef, GSDipoles, GSDipolesMask, &SMatrix[SecSCount * *Parms->SBasisBound], Parms);
					
					FeedbackKet = PriFeedback; FeedbackBra = SecSCount;
					memcpy(WaitOldEigenVectors, PriOldEigenVectors, IntPow2(*Parms->SBasisBound) * sizeof(RealType));
					
					for (WaitTimeCount = 0; WaitTimeCount < Parms->WaitTimeTotal; WaitTimeCount++) {
						UpdateBath(WaitBath, &SMatrix[FeedbackKet * *Parms->SBasisBound], &SMatrix[FeedbackBra * *Parms->SBasisBound], SBasis, Parms, true);
						Update(SMatrix, SEigenValues, SPropagator, WaitBath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
						Propagate(SCoefKet, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
						Propagate(SCoefBra, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
						WaitFactor *= cexp(-I * (SEigenValues[FeedbackKet] - SEigenValues[FeedbackBra]) * PreCalc);
						UpdateFeedbackState(&FeedbackKet, WaitBath, WaitOldEigenVectors, SCoefKet, SMatrix, SEigenValues, SBasis, Parms, FMKet);
						UpdateFeedbackState(&FeedbackBra, WaitBath, WaitOldEigenVectors, SCoefBra, SMatrix, SEigenValues, SBasis, Parms, FMBra);
					}
					
					memcpy(Bath, WaitBath, Parms->UnitTotal * sizeof(BathStruct));
					DeExciteAdiabat(WaitSecGCoef, &SMatrix[FeedbackBra * *Parms->SBasisBound], GSDipoles, GSDipolesMask, Parms);
					Factor = WaitFactor;
					
					Feedback = FeedbackKet; Sign = 1;
					memcpy(OldEigenVectors, WaitOldEigenVectors, IntPow2(*Parms->SBasisBound) * sizeof(RealType));
					
					for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
						DeExciteAdiabat(GCoef, &SMatrix[Feedback * *Parms->SBasisBound], GSDipoles, GSDipolesMask, Parms);
						Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].SENR += (PriSign * Sign) * Factor * CInnerProduct(GCoef, WaitSecGCoef, *Parms->GBasisBound);
						
						UpdateBath(Bath, &SMatrix[Feedback * *Parms->SBasisBound], NULL, SBasis, Parms, true);
						Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
						UpdateSign(&Sign, SMatrix, OldEigenVectors, Feedback, *Parms->SBasisBound);
						Propagate(SCoefKet, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
						Factor *= cexp(-I * SEigenValues[Feedback] * PreCalc);
						UpdateFeedbackState(&Feedback, Bath, OldEigenVectors, SCoefKet, SMatrix, SEigenValues, SBasis, Parms, FMKetOnly);
						FactorCArrays(WaitSecGCoef, GPropagator, *Parms->GBasisBound);
					}
					
					memcpy(Bath, WaitBath, Parms->UnitTotal * sizeof(BathStruct));
					DeExciteAdiabat(WaitSecGCoef, &SMatrix[FeedbackKet * *Parms->SBasisBound], GSDipoles, GSDipolesMask, Parms);
					Factor = WaitFactor;
					
					Feedback = FeedbackBra; Sign = 1;
					memcpy(OldEigenVectors, WaitOldEigenVectors, IntPow2(*Parms->SBasisBound) * sizeof(RealType));
					
					for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
						DeExciteAdiabat(GCoef, &SMatrix[Feedback * *Parms->SBasisBound], GSDipoles, GSDipolesMask, Parms);
						Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].SERR += (PriSign * Sign) * Factor * CInnerProduct(WaitSecGCoef, GCoef, *Parms->GBasisBound);
						
						UpdateBath(Bath, NULL, &SMatrix[Feedback * *Parms->SBasisBound], SBasis, Parms, true);
						Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
						UpdateSign(&Sign, SMatrix, OldEigenVectors, Feedback, *Parms->SBasisBound);
						Propagate(SCoefBra, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
						Factor *= cexp(I * SEigenValues[Feedback] * PreCalc);
						UpdateFeedbackState(&Feedback, Bath, OldEigenVectors, SCoefBra, SMatrix, SEigenValues, SBasis, Parms, FMBraOnly);
						FactorCArrays(WaitSecGCoef, GPropagator, *Parms->GBasisBound);
					}
					
					UpdateBath(PriBath, &SMatrix[PriFeedback * *Parms->SBasisBound], NULL, SBasis, Parms, true);
					Update(SMatrix, SEigenValues, SPropagator, PriBath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
					//if (PriTimeCount != 0 && RInnerProduct(&SMatrix[PriFeedback * *Parms->SBasisBound], &OldEigenVectors[PriFeedback * *Parms->SBasisBound], *Parms->SBasisBound) < .9) PriSign = -PriSign;
					UpdateSign(&PriSign, SMatrix, PriOldEigenVectors, PriFeedback, *Parms->SBasisBound);
					Propagate(PriSCoef, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
					PriFactor *= cexp(I * SEigenValues[PriFeedback] * PreCalc);
					UpdateFeedbackState(&PriFeedback, PriBath, PriOldEigenVectors, PriSCoef, SMatrix, SEigenValues, SBasis, Parms, FMKetOnly);
				}
			}
	}
	
	SaveResponse2D("ResponseGBRR", &(*Response).GBRR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseGBNR", &(*Response).GBNR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseSERR", &(*Response).SERR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseSENR", &(*Response).SENR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseEARR", &(*Response).EARR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), true);
	SaveResponse2D("ResponseEANR", &(*Response).EANR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), true);
	
	CloseDiag(Parms->DAndC);
	
	
	free(Response); free(Bath); free(PriBath); free(WaitBath);
	free(GBasis); free(SBasis); free(DBasis);
	free(BasicSHamiltonian); free(BasicDHamiltonian); free(SMatrix); free(DMatrix); free(OldEigenVectors); free(PriOldEigenVectors); free(WaitOldEigenVectors); free(SEigenValues); free(DEigenValues);
	free(GSDipoles); free(SDDipoles);
	free(GPropagator); free(SPropagator); free(DPropagator); free(GCoef); free(WaitSecGCoef); free(SCoefKet); free(SCoefBra); free(SCoef); free(PriSCoef); free(SCoefOld);
	free(GSDipolesMask); free(SDDipolesMask);
}


static void Response2D_HEOM(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const VOverlapStruct VOverlaps)
{
	typedef struct { ComplexType GBRR, GBNR, SERR, SENR, EARR, EANR; } ResponseStruct;
	
	int PriTimeCount, SecTimeCount, WaitTimeCount, GGBasisTotal, SGBasisTotal, SSBasisTotal, DSBasisTotal;
	
	// Temporarily!!!
	int HierarchyCount;
	
	ResponseStruct *Response = malloc(IntPow2(Parms->ResponseTimeTotal) * sizeof(ResponseStruct));
	
	BasisStruct *GBasis, *SBasis, *DBasis;
	
	MatrixBasisStruct *GGBasis, *SGBasis, *SSBasis, *DSBasis;
	
	RealType *BasicSHamiltonian, *BasicDHamiltonian, *SGLiouville, *SSLiouville, *DSLiouville;
	RealType *GSDipoles, *SDDipoles, *GGSGDipoles, *SGSSDipoles, *SSDSDipoles;
	
	ComplexType *GGCoef, *PriSGCoef, *SecSGCoef, *SSCoef, *DSCoef, *ConjGGCoef, *ConjSSCoef, *CoefWork;
	
	int *Hierarchy, *HierarchyCouplings, *SGLiouvilleMask, *SSLiouvilleMask, *DSLiouvilleMask;
	int *GSDipolesMask, *SDDipolesMask, *GGSGDipolesMask, *SGSSDipolesMask, *SSDSDipolesMask;
	int *SGHieperator, *SGAntiHieperator, *SSHieperator, *SSAntiHieperator, *DSHieperator, *DSAntiHieperator;
	
	
	BuildHierarchy(&Hierarchy, &HierarchyCouplings, &Parms->HierarchyTotal, Parms->UnitTotal, Parms->HierarchyDepth, Parms->Embarrassing);
	
	BuildBases(&GBasis, &SBasis, &DBasis, Parms->GBasisBound, Parms->SBasisBound, Parms->DBasisBound, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildMatrixBasis("GGBasis", &GGBasis, &GGBasisTotal, GBasis, GBasis, *Parms->GBasisBound, *Parms->GBasisBound, Parms->Embarrassing);
	BuildMatrixBasis("SGBasis", &SGBasis, &SGBasisTotal, SBasis, GBasis, *Parms->SBasisBound, *Parms->GBasisBound, Parms->Embarrassing);
	BuildMatrixBasis("SSBasis", &SSBasis, &SSBasisTotal, SBasis, SBasis, *Parms->SBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	BuildMatrixBasis("DSBasis", &DSBasis, &DSBasisTotal, DBasis, SBasis, *Parms->DBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	
	BuildHieperators(&SGHieperator, &SGAntiHieperator, SGBasis, SGBasisTotal, SBasis, GBasis, Parms);
	BuildHieperators(&SSHieperator, &SSAntiHieperator, SSBasis, SSBasisTotal, SBasis, SBasis, Parms);
	BuildHieperators(&DSHieperator, &DSAntiHieperator, DSBasis, DSBasisTotal, DBasis, SBasis, Parms);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildBasicHamiltonian("DHamiltonian", &BasicDHamiltonian, NULL, XEnergies, VOverlaps, DBasis, Parms->DBasisBound, Parms);
	
	BuildMatrixOperator("SGLiouville", &SGLiouville, &SGLiouvilleMask, BasicSHamiltonian, NULL, SGBasis, SGBasis, *Parms->SBasisBound, *Parms->GBasisBound, *Parms->SBasisBound, *Parms->GBasisBound, Parms->Embarrassing);
	BuildMatrixOperator("SSLiouville", &SSLiouville, &SSLiouvilleMask, BasicSHamiltonian, BasicSHamiltonian, SSBasis, SSBasis, *Parms->SBasisBound, *Parms->SBasisBound, *Parms->SBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	BuildMatrixOperator("DSLiouville", &DSLiouville, &DSLiouvilleMask, BasicDHamiltonian, BasicSHamiltonian, DSBasis, DSBasis, *Parms->DBasisBound, *Parms->SBasisBound, *Parms->DBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	RealTimesRArray(SGLiouville, InvCmToFreq * Parms->TimeStep, IntPow2(SGBasisTotal));
	RealTimesRArray(SSLiouville, InvCmToFreq * Parms->TimeStep, IntPow2(SSBasisTotal));
	RealTimesRArray(DSLiouville, InvCmToFreq * Parms->TimeStep, IntPow2(DSBasisTotal));
	
	BuildDipoles("GSDipoles", &GSDipoles, &GSDipolesMask, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildDipoles("SDDipoles", &SDDipoles, &SDDipolesMask, SBasis, DBasis, *Parms->SBasisBound, *Parms->DBasisBound, VOverlaps, Parms);
	BuildMatrixOperator("GGSGDipoles", &GGSGDipoles, &GGSGDipolesMask, GSDipoles, NULL, GGBasis, SGBasis, *Parms->GBasisBound, *Parms->GBasisBound, *Parms->SBasisBound, *Parms->GBasisBound, Parms->Embarrassing);
	BuildMatrixOperator("SGSSDipoles", &SGSSDipoles, &SGSSDipolesMask, NULL, GSDipoles, SGBasis, SSBasis, *Parms->SBasisBound, *Parms->GBasisBound, *Parms->SBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	BuildMatrixOperator("SSDSDipoles", &SSDSDipoles, &SSDSDipolesMask, SDDipoles, NULL, SSBasis, DSBasis, *Parms->SBasisBound, *Parms->SBasisBound, *Parms->DBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	
	GGCoef = calloc(GGBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType));
	PriSGCoef = calloc(SGBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType));
	SecSGCoef = calloc(SGBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType));
	SSCoef = calloc(SSBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType));
	DSCoef = calloc(DSBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType));
	
	ConjGGCoef = malloc(GGBasisTotal * Parms->HierarchyTotal * sizeof(ComplexType));
	ConjSSCoef = malloc(SSBasisTotal * Parms->HierarchyTotal * sizeof(ComplexType));
	CoefWork = malloc(7 * IntMax(IntMax(GGBasisTotal, SGBasisTotal), IntMax(SSBasisTotal, DSBasisTotal)) * Parms->HierarchyTotal * sizeof(ComplexType));
	
	GGCoef[0] = 1;
	
	Excite(PriSGCoef, GGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, SGBasisTotal, GGBasisTotal, 1);
	
	for (PriTimeCount = 0; PriTimeCount < Parms->ResponseTimeTotal; PriTimeCount++) {
		
		DeExcite(GGCoef, PriSGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, GGBasisTotal, SGBasisTotal, Parms->HierarchyTotal);
		
		for (WaitTimeCount = 0; WaitTimeCount < Parms->WaitTimeTotal; WaitTimeCount++)
			Propagate_HEOM(GGCoef, CoefWork, NULL, NULL, NULL, NULL, GGBasisTotal, Hierarchy, HierarchyCouplings, Parms);
		
		for (HierarchyCount = 0; HierarchyCount < Parms->HierarchyTotal; HierarchyCount++)
			Conjugate(&ConjGGCoef[HierarchyCount * GGBasisTotal], &GGCoef[HierarchyCount * GGBasisTotal], *Parms->GBasisBound);
		
		Excite(SecSGCoef, GGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, SGBasisTotal, GGBasisTotal, Parms->HierarchyTotal);
		
		for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
			DeExcite(GGCoef, SecSGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, GGBasisTotal, SGBasisTotal, 1);
			Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].GBNR = GGCoef[0];
			Propagate_HEOM(SecSGCoef, CoefWork, SGLiouville, SGLiouvilleMask, SGHieperator, SGAntiHieperator, SGBasisTotal, Hierarchy, HierarchyCouplings, Parms);
		}
		
		Excite(SecSGCoef, ConjGGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, SGBasisTotal, GGBasisTotal, Parms->HierarchyTotal);
		for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
			DeExcite(GGCoef, SecSGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, GGBasisTotal, SGBasisTotal, 1);
			Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].GBRR = GGCoef[0];
			Propagate_HEOM(SecSGCoef, CoefWork, SGLiouville, SGLiouvilleMask, SGHieperator, SGAntiHieperator, SGBasisTotal, Hierarchy, HierarchyCouplings, Parms);
		}
		
		Excite(SSCoef, PriSGCoef, PermDipoles, SGSSDipoles, SGSSDipolesMask, SSBasisTotal, SGBasisTotal, Parms->HierarchyTotal);
		for (WaitTimeCount = 0; WaitTimeCount < Parms->WaitTimeTotal; WaitTimeCount++)
			Propagate_HEOM(SSCoef, CoefWork, SSLiouville, SSLiouvilleMask, SSHieperator, SSAntiHieperator, SSBasisTotal, Hierarchy, HierarchyCouplings, Parms);
		
		DeExcite(SecSGCoef, SSCoef, PermDipoles, SGSSDipoles, SGSSDipolesMask, SGBasisTotal, SSBasisTotal, Parms->HierarchyTotal);
		for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
			DeExcite(GGCoef, SecSGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, GGBasisTotal, SGBasisTotal, 1);
			Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].SENR = GGCoef[0];
			Propagate_HEOM(SecSGCoef, CoefWork, SGLiouville, SGLiouvilleMask, SGHieperator, SGAntiHieperator, SGBasisTotal, Hierarchy, HierarchyCouplings, Parms);
		}
		
		for (HierarchyCount = 0; HierarchyCount < Parms->HierarchyTotal; HierarchyCount++)
			Conjugate(&ConjSSCoef[HierarchyCount * SSBasisTotal], &SSCoef[HierarchyCount * SSBasisTotal], *Parms->SBasisBound);
		DeExcite(SecSGCoef, ConjSSCoef, PermDipoles, SGSSDipoles, SGSSDipolesMask, SGBasisTotal, SSBasisTotal, Parms->HierarchyTotal);
		for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
			DeExcite(GGCoef, SecSGCoef, PermDipoles, GGSGDipoles, GGSGDipolesMask, GGBasisTotal, SGBasisTotal, 1);
			Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].SERR = GGCoef[0];
			Propagate_HEOM(SecSGCoef, CoefWork, SGLiouville, SGLiouvilleMask, SGHieperator, SGAntiHieperator, SGBasisTotal, Hierarchy, HierarchyCouplings, Parms);
		}
		
		Excite(DSCoef, SSCoef, PermDipoles, SSDSDipoles, SSDSDipolesMask, DSBasisTotal, SSBasisTotal, Parms->HierarchyTotal);
		for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
			DeExcite(SSCoef, DSCoef, PermDipoles, SSDSDipoles, SSDSDipolesMask, SSBasisTotal, DSBasisTotal, 1);
			Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].EANR = -Trace(SSCoef, *Parms->SBasisBound);
			Propagate_HEOM(DSCoef, CoefWork, DSLiouville, DSLiouvilleMask, DSHieperator, DSAntiHieperator, DSBasisTotal, Hierarchy, HierarchyCouplings, Parms);
		}
		
		Excite(DSCoef, ConjSSCoef, PermDipoles, SSDSDipoles, SSDSDipolesMask, DSBasisTotal, SSBasisTotal, Parms->HierarchyTotal);
		for (SecTimeCount = 0; SecTimeCount < Parms->ResponseTimeTotal; SecTimeCount++) {
			DeExcite(SSCoef, DSCoef, PermDipoles, SSDSDipoles, SSDSDipolesMask, SSBasisTotal, DSBasisTotal, 1);
			Response[PriTimeCount * Parms->ResponseTimeTotal + SecTimeCount].EARR = -Trace(SSCoef, *Parms->SBasisBound);
			Propagate_HEOM(DSCoef, CoefWork, DSLiouville, DSLiouvilleMask, DSHieperator, DSAntiHieperator, DSBasisTotal, Hierarchy, HierarchyCouplings, Parms);
		}
		
		Propagate_HEOM(PriSGCoef, CoefWork, SGLiouville, SGLiouvilleMask, SGHieperator, SGAntiHieperator, SGBasisTotal, Hierarchy, HierarchyCouplings, Parms);
	}
	
	SaveResponse2D("ResponseGBRR", &(*Response).GBRR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseGBNR", &(*Response).GBNR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseSERR", &(*Response).SERR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseSENR", &(*Response).SENR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseEARR", &(*Response).EARR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), true);
	SaveResponse2D("ResponseEANR", &(*Response).EANR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), true);
	
	
	free(Response);
	free(GBasis); free(SBasis); free(DBasis); free(GGBasis); free(SGBasis); free(SSBasis); free(DSBasis);
	free(BasicSHamiltonian); free(BasicDHamiltonian); free(SGLiouville); free(SSLiouville); free(DSLiouville);
	free(GSDipoles); free(SDDipoles); free(GGSGDipoles); free(SGSSDipoles); free(SSDSDipoles);
	free(GGCoef); free(PriSGCoef); free(SecSGCoef); free(SSCoef); free(DSCoef); free(ConjGGCoef); free(ConjSSCoef); free(CoefWork);
	free(Hierarchy); free(HierarchyCouplings); free(SGLiouvilleMask); free(SSLiouvilleMask); free(DSLiouvilleMask);
	free(GSDipolesMask); free(SDDipolesMask); free(GGSGDipolesMask); free(SGSSDipolesMask); free(SSDSDipolesMask);
	free(SGHieperator); free(SGAntiHieperator); free(SSHieperator); free(SSAntiHieperator); free(DSHieperator); free(DSAntiHieperator);
}


static void PopTransfer(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	int TimeCount;
	
	int *VelocityMap;
	
	BasisStruct *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *SEigenValues, *XPopulation, *AllPopulation;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType)), *TrapPopulation = calloc(Parms->WaitTimeTotal, sizeof(RealType)), *Velocity = calloc(Parms->WaitTimeTotal, sizeof(RealType));
	
	ComplexType *SPropagator, *SCoef, *SCoefOld;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(NULL, &SBasis, NULL, NULL, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	VelocityMap = calloc(IntPow2(*Parms->SBasisBound), sizeof(int));
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, VelocityMap, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	SCoef = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	
	XPopulation = calloc(Parms->UnitTotal * Parms->WaitTimeTotal, sizeof(RealType));
	AllPopulation = calloc(IntPow2(Parms->UnitTotal) * Parms->WaitTimeTotal, sizeof(RealType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	if (Parms->Trotter) PrepareTrotter(BasicSHamiltonian, SPropagator, *Parms->SBasisBound, Parms);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(Bath, Parms);
		Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, Parms->Trotter);
        
		FillCArray(SCoef, 0, *Parms->SBasisBound);
        
        int site = Parms->InitialSite;
        printf("Initial excitation on site %d.\n", site);
		SCoef[site * Parms->VTotal] = 1;
		
		for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
			
			CalculateXPopulations(&XPopulation[TimeCount * Parms->UnitTotal], SCoef, SBasis, Parms->SBasisBound);
			CalculateAllPopulations(&AllPopulation[TimeCount * IntPow2(Parms->UnitTotal)], SCoef, SBasis, Parms->SBasisBound, Parms->UnitTotal);
			CalculateVelocity(&Velocity[TimeCount], SCoef, BasicSHamiltonian, VelocityMap, *Parms->SBasisBound);
			if (TimeCount != 0 && Parms->TrapRate) Trap(&TrapPopulation[TimeCount - 1], SCoef, SBasis, Parms->SBasisBound, Parms->TrapRate * Parms->TimeStep, Parms->UnitTotal - 1); 
			UpdateBath(Bath, NULL, NULL, NULL, Parms, false);
			Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, Parms->Trotter);
			Propagate(SCoef, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, Parms->Trotter);
		}
	}	
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesRArray(XPopulation, 1.0 / Parms->SampleTotal, Parms->UnitTotal * Parms->WaitTimeTotal);
	SaveRArray("XPopulations", XPopulation, Parms->WaitTimeTotal, Parms->UnitTotal, TimeGrid, Parms->Embarrassing);
	RealTimesRArray(AllPopulation, 1.0 / Parms->SampleTotal, IntPow2(Parms->UnitTotal) * Parms->WaitTimeTotal);
	SaveRArray("AllPopulations", AllPopulation, Parms->WaitTimeTotal, IntPow2(Parms->UnitTotal), TimeGrid, Parms->Embarrassing);
	RealTimesRArray(Velocity, InvCmToFreq / Parms->SampleTotal, Parms->WaitTimeTotal);
	SaveRArray("XVelocities", Velocity, Parms->WaitTimeTotal, 1, TimeGrid, Parms->Embarrassing);
	
	CumulateRArray(TrapPopulation, Parms->WaitTimeTotal);
	RealTimesRArray(TrapPopulation, 1.0 / Parms->SampleTotal, Parms->WaitTimeTotal);
	SaveRArray("TrapPopulation", TrapPopulation, Parms->WaitTimeTotal, 1, TimeGrid, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(VelocityMap);
	free(Bath);
	free(SBasis);
	free(BasicSHamiltonian); free(SMatrix); free(SEigenValues); free(XPopulation); free(AllPopulation); free(TimeGrid); free(TrapPopulation); free(Velocity);
	free(SPropagator); free(SCoef); free(SCoefOld);
}


static void PopTransfer_SH(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	int TimeCount, Feedback, ExciteState;
	
	BasisStruct *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *OldEigenVectors, *SEigenValues, *Populations;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	ComplexType *SPropagator, *SCoef, *SCoefOld;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(NULL, &SBasis, NULL, NULL, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType)); OldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	SCoef = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	Populations = calloc(Parms->WaitTimeTotal * *Parms->SBasisBound, sizeof(RealType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	
	ExciteState = *Parms->SBasisBound - 1;
	//ExciteState = 0;
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(Bath, Parms);
		Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
		Feedback = ExciteState;
		RealToComplex(SCoef, &SMatrix[ExciteState * *Parms->SBasisBound], *Parms->SBasisBound);
		
		OldEigenVectors[0] = -2;
		
		for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
			
			Populations[TimeCount * *Parms->SBasisBound + Feedback]++;
			UpdateBath(Bath, &SMatrix[Feedback * *Parms->SBasisBound], NULL, SBasis, Parms, false);
			Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
			Propagate(SCoef, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
			UpdateFeedbackState(&Feedback, Bath, OldEigenVectors, SCoef, SMatrix, SEigenValues, SBasis, Parms, FMKetOnly);
		}
	}
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesRArray(Populations, 1.0 / Parms->SampleTotal, *Parms->SBasisBound * Parms->WaitTimeTotal);
	SaveRArray("Populations", Populations, Parms->WaitTimeTotal, *Parms->SBasisBound, TimeGrid, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(Bath);
	free(SBasis);
	free(BasicSHamiltonian); free(SMatrix); free(OldEigenVectors); free(SEigenValues); free(Populations); free(TimeGrid);
	free(SPropagator); free(SCoef); free(SCoefOld);
}


static void PopTransfer_HEOM(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	int TimeCount, BasisCount, HierarchyTotal, SSBasisTotal, XEHBasisBound[5], ExciteUnit = 0;
	
	BasisStruct *SBasis;
	
	MatrixBasisStruct *SSBasis;
	
	RealType *BasicSHamiltonian, *SSLiouville, *BoltzmannFactor, *XPopulation, *AllPopulation;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	ComplexType *SSCoef, *SSCoefWork;
	
	int *Hierarchy, *HierarchyCouplings, *SSLiouvilleMask, *SSHieperator, *SSAntiHieperator;
	
	MonteCarloStruct MC = {0,0,0};
	
		
	BuildHierarchy(&Hierarchy, &HierarchyCouplings, &Parms->HierarchyTotal, Parms->UnitTotal, Parms->HierarchyDepth, Parms->Embarrassing);
		
	BuildBases(NULL, &SBasis, NULL, NULL, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	BuildMatrixBasis("SSBasis", &SSBasis, &SSBasisTotal, SBasis, SBasis, *Parms->SBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	
	BuildHieperators(&SSHieperator, &SSAntiHieperator, SSBasis, SSBasisTotal, SBasis, SBasis, Parms);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildMatrixOperator("SSLiouville", &SSLiouville, &SSLiouvilleMask, BasicSHamiltonian, BasicSHamiltonian, SSBasis, SSBasis, *Parms->SBasisBound, *Parms->SBasisBound, *Parms->SBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	RealTimesRArray(SSLiouville, InvCmToFreq * Parms->TimeStep, IntPow2(SSBasisTotal));
	
	
	SSCoef = calloc(SSBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType)); SSCoefWork = malloc(7 * SSBasisTotal * Parms->HierarchyTotal * sizeof(ComplexType));
	XPopulation = calloc(Parms->UnitTotal * Parms->WaitTimeTotal, sizeof(RealType));
	AllPopulation = calloc(IntPow2(Parms->UnitTotal) * Parms->WaitTimeTotal, sizeof(RealType));
	
	SSCoef[ExciteUnit * (Parms->UnitTotal + 1)] = 1;
	
	CalculateBoltzmannFactors(&BoltzmannFactor, BasicSHamiltonian, Parms->Temperature * KJPMolToInvCm, *Parms->SBasisBound);
	PrintRArray(BoltzmannFactor, 1, *Parms->SBasisBound);
	
	
	for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
		
		if (TimeCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->WaitTimeTotal, Parms->Embarrassing);
		
		//for (BasisCount = 0; BasisCount < *Parms->SBasisBound; BasisCount++)
		//	Population[TimeCount * *Parms->SBasisBound + BasisCount] = creal(SSCoef[BasisCount * (*Parms->SBasisBound + 1)]);
		SSCoefToXPopulations(&XPopulation[TimeCount * Parms->UnitTotal], SSCoef, SBasis, Parms->SBasisBound);
		SSCoefToAllPopulations(&AllPopulation[TimeCount * IntPow2(Parms->UnitTotal)], SSCoef, SBasis, Parms->SBasisBound, Parms->UnitTotal);
		Propagate_HEOM(SSCoef, SSCoefWork, SSLiouville, SSLiouvilleMask, SSHieperator, SSAntiHieperator, SSBasisTotal, Hierarchy, HierarchyCouplings, Parms);
	}
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	SaveRArray("XPopulations", XPopulation, Parms->WaitTimeTotal, Parms->UnitTotal, TimeGrid, Parms->Embarrassing);
	SaveRArray("AllPopulations", AllPopulation, Parms->WaitTimeTotal, IntPow2(Parms->UnitTotal), TimeGrid, Parms->Embarrassing);
	
	
	free(SBasis); free(SSBasis);
	free(BasicSHamiltonian); free(SSLiouville); free(BoltzmannFactor); free(XPopulation); free(AllPopulation); free(TimeGrid);
	free(SSCoef); free(SSCoefWork);
	free(Hierarchy); free(HierarchyCouplings); free(SSLiouvilleMask); free(SSHieperator); free(SSAntiHieperator);
}


static void TempCoh(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	int TimeCount, Sign;
	
	BasisStruct *GBasis, *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *OldEigenVectors, *SEigenValues, *GSDipoles;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	ComplexType *SPropagator, *GCoef, *SCoefKet, *SCoefBra, *SCoefOld;
	
	double complex *Coherence = calloc(Parms->WaitTimeTotal, sizeof(double complex));
	
	int *GSDipolesMask, *SOS;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(&GBasis, &SBasis, NULL, Parms->GBasisBound, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildDipoles("GSDipoles", &GSDipoles, &GSDipolesMask, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildSOS(&SOS, Parms->SOSChar, *Parms->SBasisBound);
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType)); OldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	
	GCoef = malloc(*Parms->GBasisBound * sizeof(ComplexType));
	SCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(Bath, Parms);
		Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
		
		InitGCoef(GCoef, *Parms->GBasisBound);
		ExciteAdiabat(SCoefKet, GCoef, GSDipoles, GSDipolesMask, &SMatrix[SOS[0] * *Parms->SBasisBound], Parms);
		InitGCoef(GCoef, *Parms->GBasisBound);
		ExciteAdiabat(SCoefBra, GCoef, GSDipoles, GSDipolesMask, &SMatrix[SOS[1] * *Parms->SBasisBound], Parms);
		
		Sign = 1; *OldEigenVectors = -2;
		
		for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
			Coherence[TimeCount] += Sign * RCInnerProduct(&SMatrix[SOS[0] * *Parms->SBasisBound], SCoefKet, *Parms->SBasisBound) * conj(RCInnerProduct(&SMatrix[SOS[1] * *Parms->SBasisBound], SCoefBra, *Parms->SBasisBound));
			UpdateBath(Bath, NULL, NULL, NULL, Parms, false);
			Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
			UpdateSign(&Sign, SMatrix, OldEigenVectors, SOS[0], *Parms->SBasisBound);
			UpdateSign(&Sign, SMatrix, OldEigenVectors, SOS[1], *Parms->SBasisBound);
			Propagate(SCoefKet, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, Parms->Trotter);
			Propagate(SCoefBra, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, Parms->Trotter);
			memcpy(OldEigenVectors, SMatrix, IntPow2(*Parms->SBasisBound) * sizeof(RealType));
		}
	}
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesDCArray(Coherence, 1.0 / Parms->SampleTotal, Parms->WaitTimeTotal);
	SaveDCVector("TempCoh", Coherence, TimeGrid, Parms->WaitTimeTotal, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(GBasis); free(SBasis); free(Bath);
	free(BasicSHamiltonian); free(SMatrix); free(SEigenValues); free(GSDipoles); free(TimeGrid);
	free(SPropagator); free(GCoef); free(SCoefKet); free(SCoefBra); free(SCoefOld);
	free(GSDipolesMask); free(SOS); free(Coherence);
}


static void TempCoh_(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const VOverlapStruct VOverlaps)
{
	int TimeCount, Sign;
	
	BasisStruct *GBasis, *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *OldEigenVectors, *SEigenValues, *GSDipoles;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	RealType Factor;
	
	ComplexType *SPropagator, *GCoef, *SCoefKet, *SCoefBra, *SCoefOld;
	
	ComplexType Temp;
	
	double complex *Coherence = calloc(Parms->WaitTimeTotal, sizeof(double complex));
	
	int *GSDipolesMask, *SOS;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(&GBasis, &SBasis, NULL, Parms->GBasisBound, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildDipoles("GSDipoles", &GSDipoles, &GSDipolesMask, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildSOS(&SOS, Parms->SOSChar, *Parms->SBasisBound);
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType)); OldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	
	GCoef = malloc(*Parms->GBasisBound * sizeof(ComplexType));
	SCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(Bath, Parms);
		Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
		
		InitGCoef(GCoef, *Parms->GBasisBound);
		Factor = ExciteAdiabat(SCoefKet, GCoef, GSDipoles, GSDipolesMask, &SMatrix[SOS[0] * *Parms->SBasisBound], Parms);
		InitGCoef(GCoef, *Parms->GBasisBound);
		Factor *= ExciteAdiabat(SCoefBra, GCoef, GSDipoles, GSDipolesMask, &SMatrix[SOS[1] * *Parms->SBasisBound], Parms);
		
		Sign = 1; *OldEigenVectors = -2;
		
		for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
			DeExcite(GCoef, SCoefKet, PermDipoles, GSDipoles, GSDipolesMask, *Parms->GBasisBound, *Parms->SBasisBound, 1);
			Temp = GCoef[0];
			DeExcite(GCoef, SCoefBra, PermDipoles, GSDipoles, GSDipolesMask, *Parms->GBasisBound, *Parms->SBasisBound, 1);
			
			//printf("%f + %fi, %f + %fi\n", creal(GCoef[0]), cimag(GCoef[0]), creal(Temp), cimag(Temp));
			
			Coherence[TimeCount] += Factor * Temp * conj(GCoef[0]);
			
			UpdateBath(Bath, NULL, NULL, NULL, Parms, false);
			Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
			//UpdateSign(&Sign, SMatrix, OldEigenVectors, SOS[0], *Parms->SBasisBound);
			//UpdateSign(&Sign, SMatrix, OldEigenVectors, SOS[1], *Parms->SBasisBound);
			//printf("%d, %d\n", TimeCount, Sign);
			//PrintRArray(SMatrix, *Parms->SBasisBound, *Parms->SBasisBound);
			Propagate(SCoefKet, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, Parms->Trotter);
			Propagate(SCoefBra, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, Parms->Trotter);
			memcpy(OldEigenVectors, SMatrix, IntPow2(*Parms->SBasisBound) * sizeof(RealType));
		}
	}
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesDCArray(Coherence, 1.0 / Parms->SampleTotal, Parms->WaitTimeTotal);
	SaveDCVector("TempCoh", Coherence, TimeGrid, Parms->WaitTimeTotal, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(GBasis); free(SBasis); free(Bath);
	free(BasicSHamiltonian); free(SMatrix); free(SEigenValues); free(GSDipoles); free(TimeGrid);
	free(SPropagator); free(GCoef); free(SCoefKet); free(SCoefBra); free(SCoefOld);
	free(GSDipolesMask); free(SOS); free(Coherence);
}


static void TempCoh_Site(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	int TimeCount;
	
	BasisStruct *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *SEigenValues;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	ComplexType *SPropagator, *SCoefKet, *SCoefBra, *SCoefOld;
	
	double complex *Coherence = calloc(Parms->WaitTimeTotal, sizeof(double complex));
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(NULL, &SBasis, NULL, NULL, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	
	SCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	if (Parms->Trotter) PrepareTrotter(BasicSHamiltonian, SPropagator, *Parms->SBasisBound, Parms);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(Bath, Parms);
		Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, Parms->Trotter);
		
		FillCArray(SCoefKet, 0, *Parms->SBasisBound);
		memcpy(SCoefBra, SCoefKet, *Parms->SBasisBound * sizeof(ComplexType));
		SCoefKet[0] = 1; SCoefBra[1] = 1;
		
		for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
			Coherence[TimeCount] += SCoefKet[0] * conj(SCoefBra[1]);
			UpdateBath(Bath, NULL, NULL, NULL, Parms, false);
			Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, Parms->Trotter);
			Propagate(SCoefKet, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, Parms->Trotter);
			Propagate(SCoefBra, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, Parms->Trotter);
		}
	}
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesDCArray(Coherence, 1.0 / Parms->SampleTotal, Parms->WaitTimeTotal);
	SaveDCVector("TempCoh", Coherence, TimeGrid, Parms->WaitTimeTotal, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(SBasis); free(Bath);
	free(BasicSHamiltonian); free(SMatrix); free(SEigenValues); free(TimeGrid);
	free(SPropagator); free(SCoefKet); free(SCoefBra); free(SCoefOld);
	free(Coherence);
}


static void TempCoh_SH(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	int TimeCount, PriSId = 0, SecSId = 1, FeedbackKet, FeedbackBra;
	
	BasisStruct *SBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *OldEigenVectors, *SEigenValues;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	RealType PreCalc = InvCmToFreq * Parms->TimeStep;
	
	ComplexType *SPropagator, *SCoefKet, *SCoefBra, *SCoefOld;
	
	ComplexType Factor;
	
	double complex *Coherence = calloc(Parms->WaitTimeTotal, sizeof(double complex));
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(NULL, &SBasis, NULL, NULL, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType)); OldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	
	SCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(Bath, Parms);
		Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
		RealToComplex(SCoefKet, &SMatrix[PriSId * *Parms->SBasisBound], *Parms->SBasisBound);
		RealToComplex(SCoefBra, &SMatrix[SecSId * *Parms->SBasisBound], *Parms->SBasisBound);
		
		FeedbackKet = PriSId; FeedbackBra = SecSId; *OldEigenVectors = -2, Factor = 1;
		
		for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
			Coherence[TimeCount] += (FeedbackKet == PriSId) * (FeedbackBra == SecSId) * Factor;
			
			UpdateBath(Bath, &SMatrix[FeedbackKet * *Parms->SBasisBound], &SMatrix[FeedbackBra * *Parms->SBasisBound], SBasis, Parms, true);
			Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
			Propagate(SCoefKet, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
			Propagate(SCoefBra, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
			Factor *= cexp(-I * (SEigenValues[FeedbackKet] - SEigenValues[FeedbackBra]) * PreCalc);
			UpdateFeedbackState(&FeedbackKet, Bath, OldEigenVectors, SCoefKet, SMatrix, SEigenValues, SBasis, Parms, FMKet);
			UpdateFeedbackState(&FeedbackBra, Bath, OldEigenVectors, SCoefBra, SMatrix, SEigenValues, SBasis, Parms, FMBra);
		}
	}
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesDCArray(Coherence, 1.0 / Parms->SampleTotal, Parms->WaitTimeTotal);
	SaveDCVector("TempCoh", Coherence, TimeGrid, Parms->WaitTimeTotal, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(SBasis); free(Bath);
	free(BasicSHamiltonian); free(SMatrix); free(OldEigenVectors); free(SEigenValues); free(TimeGrid);
	free(Coherence);
	free(SPropagator); free(SCoefKet); free(SCoefBra); free(SCoefOld);
}


static void TempCoh_SH_Site(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	int TimeCount, SiteKet = 0, SiteBra = 1, SCountKet, SCountBra, FeedbackKet, FeedbackBra;
	
	BasisStruct *SBasis;
	
	BathStruct *Bath = malloc(Parms->UnitTotal * sizeof(BathStruct)), *StartBath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicSHamiltonian, *SMatrix, *OldEigenVectors, *SEigenValues;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	RealType PreCalc = InvCmToFreq * Parms->TimeStep;
	
	ComplexType *SPropagator, *SCoefKet, *SCoefBra, *SCoefOld;
	
	ComplexType Factor;
	
	double complex *Coherence = calloc(Parms->WaitTimeTotal, sizeof(double complex));
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(NULL, &SBasis, NULL, NULL, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	
	SMatrix = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType)); OldEigenVectors = malloc(IntPow2(*Parms->SBasisBound) * sizeof(RealType));
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	SPropagator = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	
	SCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoefOld = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(StartBath, Parms);
		
		for (SCountKet = 0; SCountKet < *Parms->SBasisBound; SCountKet++)
			for (SCountBra = 0; SCountBra < *Parms->SBasisBound; SCountBra++) {
				
				//if (SCountKet == SCountBra) {
				//	for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
				//		Coherence[TimeCount] +=
				//	}
				//}
				
				memcpy(Bath, StartBath, Parms->UnitTotal * sizeof(BathStruct));
				Initialize(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
				RealToComplex(SCoefKet, &SMatrix[SCountKet * *Parms->SBasisBound], *Parms->SBasisBound);
				RealToComplex(SCoefBra, &SMatrix[SCountBra * *Parms->SBasisBound], *Parms->SBasisBound);
				
				FeedbackKet = SCountKet; FeedbackBra = SCountBra; *OldEigenVectors = -2; 
				Factor = SMatrix[SCountKet * *Parms->SBasisBound + SiteKet] * SMatrix[SCountBra * *Parms->SBasisBound + SiteBra];
				
				for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
					Coherence[TimeCount] += Factor * SMatrix[FeedbackKet * *Parms->SBasisBound + SiteKet] * SMatrix[FeedbackBra * *Parms->SBasisBound + SiteBra];
					
					UpdateBath(Bath, &SMatrix[FeedbackKet * *Parms->SBasisBound], &SMatrix[FeedbackBra * *Parms->SBasisBound], SBasis, Parms, true);
					Update(SMatrix, SEigenValues, SPropagator, Bath, BasicSHamiltonian, SBasis, Parms, *Parms->SBasisBound, false);
					Propagate(SCoefKet, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
					Propagate(SCoefBra, SCoefOld, SMatrix, SPropagator, *Parms->SBasisBound, Parms, false);
					Factor *= cexp(-I * (SEigenValues[FeedbackKet] - SEigenValues[FeedbackBra]) * PreCalc);
					UpdateFeedbackState(&FeedbackKet, Bath, OldEigenVectors, SCoefKet, SMatrix, SEigenValues, SBasis, Parms, FMKet);
					UpdateFeedbackState(&FeedbackBra, Bath, OldEigenVectors, SCoefBra, SMatrix, SEigenValues, SBasis, Parms, FMBra);
				}
			}
	}
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesDCArray(Coherence, 1.0 / Parms->SampleTotal, Parms->WaitTimeTotal);
	SaveDCVector("TempCoh", Coherence, TimeGrid, Parms->WaitTimeTotal, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(SBasis); free(Bath); free(StartBath);
	free(BasicSHamiltonian); free(SMatrix); free(OldEigenVectors); free(SEigenValues); free(TimeGrid);
	free(Coherence);
	free(SPropagator); free(SCoefKet); free(SCoefBra); free(SCoefOld);
}


static void TempCoh_HEOM(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)

// Change BasicSHamiltonian to SMatrix!

{
	int TimeCount, SSBasisTotal;
	
	BasisStruct *SBasis;
	
	MatrixBasisStruct *SSBasis;
	
	RealType *BasicSHamiltonian, *SEigenValues, *SSLiouville;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	ComplexType *SSCoef, *SSCoefWork;
	
	ComplexType *Coherence = malloc(Parms->WaitTimeTotal * sizeof(ComplexType));
	
	int *Hierarchy, *HierarchyCouplings, *SSLiouvilleMask, *SSHieperator, *SSAntiHieperator;
	
	
	MonteCarloStruct MC = {0,0,0};
	
	BuildHierarchy(&Hierarchy, &HierarchyCouplings, &Parms->HierarchyTotal, Parms->UnitTotal, Parms->HierarchyDepth, Parms->Embarrassing);
	
	BuildBases(NULL, &SBasis, NULL, NULL, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildMatrixBasis("SSBasis", &SSBasis, &SSBasisTotal, SBasis, SBasis, *Parms->SBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	
	BuildHieperators(&SSHieperator, &SSAntiHieperator, SSBasis, SSBasisTotal, SBasis, SBasis, Parms);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildMatrixOperator("SSLiouville", &SSLiouville, &SSLiouvilleMask, BasicSHamiltonian, BasicSHamiltonian, SSBasis, SSBasis, *Parms->SBasisBound, *Parms->SBasisBound, *Parms->SBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	RealTimesRArray(SSLiouville, InvCmToFreq * Parms->TimeStep, IntPow2(SSBasisTotal));
	
	SEigenValues = malloc(*Parms->SBasisBound * sizeof(RealType));
	
	SSCoef = calloc(SSBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType)); SSCoefWork = malloc(7 * SSBasisTotal * Parms->HierarchyTotal * sizeof(ComplexType));
	
	
	InitializeDiag(BasicSHamiltonian, NULL, *Parms->SBasisBound, Parms->DAndC);
	DiagonalizeMatrix(BasicSHamiltonian, SEigenValues, *Parms->SBasisBound, Parms->DAndC);
	
	//SSCoef[1] = 1;
	ExciteMatrixAdiabat(SSCoef, BasicSHamiltonian, 0, 1, *Parms->SBasisBound);
	
	for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
		if (TimeCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->WaitTimeTotal, Parms->Embarrassing);
		//Coherence[TimeCount] = SSCoef[1];
		Coherence[TimeCount] = DeExciteMatrixAdiabat(SSCoef, BasicSHamiltonian, 0, 1, *Parms->SBasisBound);
		Propagate_HEOM(SSCoef, SSCoefWork, SSLiouville, SSLiouvilleMask, SSHieperator, SSAntiHieperator, SSBasisTotal, Hierarchy, HierarchyCouplings, Parms);
	}
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesCArray(Coherence, 1.0 / Parms->SampleTotal, Parms->WaitTimeTotal);
	SaveCVector("TempCoh", Coherence, TimeGrid, Parms->WaitTimeTotal, Parms->Embarrassing);
	
	CloseDiag(Parms->DAndC);
	
	
	free(SBasis); free(SSBasis);
	free(BasicSHamiltonian); free(SEigenValues); free(SSLiouville); free(TimeGrid);
	free(SSCoef); free(SSCoefWork); free(Coherence);
	free(Hierarchy); free(HierarchyCouplings); free(SSLiouvilleMask); free(SSHieperator); free(SSAntiHieperator);
}


static void TempCoh_HEOM_Site(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	int TimeCount, SSBasisTotal;
	
	BasisStruct *SBasis;
	
	MatrixBasisStruct *SSBasis;
	
	RealType *BasicSHamiltonian, *SSLiouville;
	
	RealType *TimeGrid = malloc(Parms->WaitTimeTotal * sizeof(RealType));
	
	ComplexType *SSCoef, *SSCoefWork;
	
	ComplexType *Coherence = malloc(Parms->WaitTimeTotal * sizeof(ComplexType));
	
	int *Hierarchy, *HierarchyCouplings, *SSLiouvilleMask, *SSHieperator, *SSAntiHieperator;
	
	
	MonteCarloStruct MC = {0,0,0};
	
	BuildHierarchy(&Hierarchy, &HierarchyCouplings, &Parms->HierarchyTotal, Parms->UnitTotal, Parms->HierarchyDepth, Parms->Embarrassing);
	
	BuildBases(NULL, &SBasis, NULL, NULL, Parms->SBasisBound, NULL, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildMatrixBasis("SSBasis", &SSBasis, &SSBasisTotal, SBasis, SBasis, *Parms->SBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	
	BuildHieperators(&SSHieperator, &SSAntiHieperator, SSBasis, SSBasisTotal, SBasis, SBasis, Parms);
	
	BuildBasicHamiltonian("SHamiltonian", &BasicSHamiltonian, NULL, XEnergies, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildMatrixOperator("SSLiouville", &SSLiouville, &SSLiouvilleMask, BasicSHamiltonian, BasicSHamiltonian, SSBasis, SSBasis, *Parms->SBasisBound, *Parms->SBasisBound, *Parms->SBasisBound, *Parms->SBasisBound, Parms->Embarrassing);
	RealTimesRArray(SSLiouville, InvCmToFreq * Parms->TimeStep, IntPow2(SSBasisTotal));
	
	SSCoef = calloc(SSBasisTotal * Parms->HierarchyTotal, sizeof(ComplexType)); SSCoefWork = malloc(7 * SSBasisTotal * Parms->HierarchyTotal * sizeof(ComplexType));
	
	
	SSCoef[1] = 1;
	
	for (TimeCount = 0; TimeCount < Parms->WaitTimeTotal; TimeCount++) {
		if (TimeCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->WaitTimeTotal, Parms->Embarrassing);
		Coherence[TimeCount] = SSCoef[1];
		Propagate_HEOM(SSCoef, SSCoefWork, SSLiouville, SSLiouvilleMask, SSHieperator, SSAntiHieperator, SSBasisTotal, Hierarchy, HierarchyCouplings, Parms);
	}
	
	CreateGrid(TimeGrid, 0, (Parms->WaitTimeTotal - 1) * Parms->TimeStep, Parms->WaitTimeTotal);
	RealTimesCArray(Coherence, 1.0 / Parms->SampleTotal, Parms->WaitTimeTotal);
	SaveCVector("TempCoh", Coherence, TimeGrid, Parms->WaitTimeTotal, Parms->Embarrassing);
	
	
	free(SBasis); free(SSBasis);
	free(BasicSHamiltonian); free(SSLiouville); free(TimeGrid);
	free(SSCoef); free(SSCoefWork); free(Coherence);
	free(Hierarchy); free(HierarchyCouplings); free(SSLiouvilleMask); free(SSHieperator); free(SSAntiHieperator);
}


static void Annihilation(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps)
{
	int StateCount, BasisCount, UnitCount;
	
	BasisStruct *DBasis;
	
	BathStruct *Bath = calloc(Parms->UnitTotal, sizeof(BathStruct));
	
	RealType *BasicDHamiltonian, *DMatrix, *DEigenValues, *Coefficients;
	
	RealType AnnihilationProbability = 0, AnnProbContribution, BoltzmannFactor, FirstBoltzmannFactor, PartitionFunction;
	
	MonteCarloStruct MC = {0,0,0};
	
	
	BuildBases(NULL, NULL, &DBasis, NULL, NULL, Parms->DBasisBound, Parms->UnitTotal, Parms->EHRadius, Parms->HDomain, Parms->VTotal, Parms->VRadius, Parms->VSlope, Parms->VExtra, Parms->Periodic, Parms->Embarrassing);
	
	BuildBasicHamiltonian("DHamiltonian", &BasicDHamiltonian, NULL, XEnergies, VOverlaps, DBasis, Parms->DBasisBound, Parms);
	
	DMatrix = malloc(IntPow2(*Parms->DBasisBound) * sizeof(RealType));
	DEigenValues = malloc(*Parms->DBasisBound * sizeof(RealType));
	Coefficients = calloc(IntPow2(*Parms->DBasisBound), sizeof(RealType));
	
	
	InitializeDiag(BasicDHamiltonian, NULL, *Parms->DBasisBound, Parms->DAndC);
	
	for (MC.SampleCount = 0; MC.SampleCount < Parms->SampleTotal; MC.SampleCount++) {
		
		if (MC.SampleCount == MC.ProgressSample) PrintProgress(&MC.ProgressPercent, &MC.ProgressSample, Parms->SampleTotal, Parms->Embarrassing);
		
		InitBath(Bath, Parms);
		BuildHamiltonian(DMatrix, BasicDHamiltonian, Bath, DBasis, *Parms->DBasisBound);
		DiagonalizeMatrix(DMatrix, DEigenValues, *Parms->DBasisBound, Parms->DAndC);
		
		for (BasisCount = 0; BasisCount < *Parms->DBasisBound; BasisCount++)
			Coefficients[DBasis[BasisCount].X[0].U * Parms->UnitTotal + DBasis[BasisCount].X[1].U] += DMatrix[BasisCount];
		
		PartitionFunction = 0;
		for (StateCount = 0; StateCount < *Parms->DBasisBound; StateCount++)
			PartitionFunction += (Parms->Temperature == 0) ? (StateCount == 0) : RealExp(-(DEigenValues[StateCount] - DEigenValues[0]) / (Parms->Temperature * KJPMolToInvCm));
		
		for (StateCount = 0; StateCount < *Parms->DBasisBound; StateCount++) {
			BoltzmannFactor = (Parms->Temperature == 0) ? (StateCount == 0) : RealExp(-(DEigenValues[StateCount] - DEigenValues[0]) / (Parms->Temperature * KJPMolToInvCm)) / PartitionFunction;
			for (UnitCount = 0; UnitCount < Parms->UnitTotal; UnitCount++) {
				AnnProbContribution = 0;
				for (BasisCount = 0; BasisCount < *Parms->DBasisBound; BasisCount++)
					if (DBasis[BasisCount].X[0].U == UnitCount || DBasis[BasisCount].X[1].U == UnitCount)
						AnnProbContribution += DMatrix[StateCount * *Parms->DBasisBound + BasisCount] / ((RealType) IntPowN(Separation(DBasis[BasisCount].X[0].U, DBasis[BasisCount].X[1].U, Parms->UnitTotal, Parms->Periodic), 3));
				AnnihilationProbability += BoltzmannFactor * RealPow2(AnnProbContribution * Parms->XCouplingNN);
			}
		}
	}
	
	RealTimesRArray(Coefficients, 1.0 / Parms->SampleTotal, IntPow2(Parms->UnitTotal));
	SaveRArray("Coefficients", Coefficients, Parms->UnitTotal, Parms->UnitTotal, NULL, Parms->Embarrassing);
	
	AnnihilationProbability /= Parms->SampleTotal;
	SaveRArray("AProbability", &AnnihilationProbability, 1, 1, NULL, Parms->Embarrassing);
	
	
	CloseDiag(Parms->DAndC);
	
	free(DBasis); free(Bath);
	free(BasicDHamiltonian); free(DMatrix); free(DEigenValues); free(Coefficients);
}
