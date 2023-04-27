#ifndef ParmsMod
#define ParmsMod

#include "GlobalsMod.h"

static const RealType Pi = 3.14159265;
static const RealType SqrtOfHalf = 0.707106781;
static const RealType LightSpeed = 2.99792458e-2;							// [cm/ps]
static const RealType InvCmToFreq = 2 * 3.14159265 * 2.99792458e-2;			// 2Pi c [cm/ps]
static const RealType FreqToInvCm = 1 / (2 * 3.14159265 * 2.99792458e-2);	// 1/(2Pi c) [ps/cm]
static const RealType KJPMolToInvCm = 1 / (0.39903132 * 2.99792458e-2);		// 1/hc [/cm mol/kJ]
static const RealType InvCmToKJPMol = 0.39903132 * 2.99792458e-2;			// hc [kJ/mol cm]
static const RealType BoltzmannCnst = 8.314510e-3;							// [kJ/mol/K]

enum ModeList { ModeNothing, ModeDebug, ModeDirect1D, ModeResponse1D, ModeResponse2D, ModePopTransfer, ModeTempCoh, ModeAnnihilation };
enum BroadList { Gauss, Lorentz };

typedef struct {
	
	enum ModeList Mode;
    enum BroadList BroadFunc;
	
	int Embarrassing;
	int DAndC;
	int Trotter;
	int SiteRep;
	int PriRandomSeed, SecRandomSeed;
	int SampleTotal;
	int ResponseTimeTotal, WaitTimeTotal;
	int GBasisBound[4], SBasisBound[5], DBasisBound[4]; // XEHBasisBound[5]
	int UnitTotal;
    int InitialSite;
    int MolPerCell;
	int XPointDipoleRadius, Periodic;
	int EHRadius, HDomain;
	int VTotal, VRadius, VExtra;
	int BathDynamics;
	int QFeedback;
	int HierarchyDepth, HierarchyTotal;
	int SpectrumTotal;
	int PulseTimeTotal;
	
	RealType TimeStep;
	RealType TwistAngle;
	RealType XEnergyMean, XEnergyOffset,  XCouplingNN, XLifeTime;
	RealType XHuangRhys, EHuangRhys, HHuangRhys;
	RealType VSlope;
	RealType VEnergy;
	RealType EHEnergyNN, EHEnergyInf, ECouplingNN, HCouplingNN, ECouplingNNRight, ECouplingNNLeft, HCouplingNNRight, HCouplingNNLeft;
	RealType BathReorganization, BathStaticDev, BathDynamicDev, BathTempDev, BathTime;
	RealType BathCoupling, BathMass, BathSpring, BathDamping;
	RealType BathPreCalc[2];
	RealType HierarchyDamping;
	RealType HierarchyPreCalc[3];
	RealType TrapRate;
	RealType Temperature;
	RealType HBroad;
	RealType EnergyMin, EnergyMax;
    RealType CrystalShift, MonomerShift;
	
	char SOSChar[1000];
	char XEnergiesFile[100], XDiagsFile[100], XCouplingsFile[100], PermDipolesFile[100], PulseFile[100];
	
} ParmStruct;

#endif
