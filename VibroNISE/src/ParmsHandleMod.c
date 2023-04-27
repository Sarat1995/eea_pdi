#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "ToolsMod.h"
#include "ParmsHandleMod.h"

#include <string.h>
#include <math.h>


void InitializeParms(ParmStruct *Parms)
{
	PrintMessage("Applying default settings.", "", "", true, Parms->Embarrassing);
	
	Parms->DAndC = true;
	Parms->PriRandomSeed = 206;
	Parms->SecRandomSeed = 1985;
	Parms->SampleTotal = 1;
	Parms->UnitTotal = 1;
	Parms->InitialSite = 0;
	Parms->XPointDipoleRadius = 1;
	Parms->VTotal = 1;
	Parms->VRadius = -1;
	Parms->HBroad = 10; // [/cm]
	Parms->SpectrumTotal = 1;
	Parms->QFeedback = false; // Necessary?
    Parms->MolPerCell = 1;
	
	sprintf(Parms->XEnergiesFile, "|");
	sprintf(Parms->XDiagsFile, "|");
	sprintf(Parms->XCouplingsFile, "|");
	sprintf(Parms->PermDipolesFile, "|");
	sprintf(Parms->PulseFile, "|");
}


void ReadInputFile(char *Args[], ParmStruct *Parms, int *StartPoint)
{
	double ReadDouble;
	int ReadInt;
	char InputFile[100], Tag[100], DisplayNumber[100];
	FILE *FileID;
	
	sprintf(InputFile, "%s", Args[1]);
	PrintMessage("Reading input file '", InputFile, "'.", true, Parms->Embarrassing);
	
	FileID = fopen(InputFile, "r");
	
	if (StartPoint != NULL)
	fseek(FileID, *StartPoint, 0), *StartPoint = -1;
	
	while (!feof(FileID)) {
		fscanf(FileID, "%s\n", Tag);
		
		if (!strcmp(Tag, "Run") && StartPoint != NULL) {
			*StartPoint = ftell(FileID);
			break;
		}
		else if (!strcmp(Tag, "Mode")) {
			fscanf(FileID, "%s\n", Tag);
			if (!strcmp(Tag, "Debug"))
				Parms->Mode = ModeDebug,		PrintMessage(" - Running in debugging mode.", "", "", false, Parms->Embarrassing);
			else if (!strcmp(Tag, "Direct1D"))
				Parms->Mode = ModeDirect1D,		PrintMessage(" - Linear spectra are calculated.", "", "", false, Parms->Embarrassing);
			else if (!strcmp(Tag, "Response1D"))
				Parms->Mode = ModeResponse1D,	PrintMessage(" - Linear response is calculated.", "", "", false, Parms->Embarrassing);
			else if (!strcmp(Tag, "Response2D"))
				Parms->Mode = ModeResponse2D,	PrintMessage(" - 2D response is calculated.", "", "", false, Parms->Embarrassing);
			else if (!strcmp(Tag, "PopTransfer"))
				Parms->Mode = ModePopTransfer,	PrintMessage(" - Population transfer is calculated.", "", "", false, Parms->Embarrassing);
			else if (!strcmp(Tag, "TempCoh"))
				Parms->Mode = ModeTempCoh,	PrintMessage(" - Temporal coherence is calculated.", "", "", false, Parms->Embarrassing);
			else if (!strcmp(Tag, "Annihilation"))
				Parms->Mode = ModeAnnihilation,	PrintMessage(" - Singlet-singlet annihilation is calculated.", "", "", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "BroadFunc")) {
			fscanf(FileID, "%s\n", Tag);
			if (!strcmp(Tag, "Gauss"))
				Parms->BroadFunc = Gauss,		PrintMessage(" - Broadening linear spectra with Gaussians.", "", "", false, Parms->Embarrassing);
			else if (!strcmp(Tag, "Lorentz"))
				Parms->BroadFunc = Lorentz,		PrintMessage(" - Broadening linear spectra with Lorentzians.", "", "", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "NoDAndC")) {
			Parms->DAndC = false;	PrintMessage(" - Divide&Conquer diagonalization algorythm is disabled.", "", "", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "Trotter")) {
			Parms->Trotter = true;	PrintMessage(" - Trotter formula is applied.", "", "", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "UberTrotter")) {
			Parms->Trotter = 2;		PrintMessage(" - UberTrotter formula is applied.", "", "", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "SiteRep")) {
			Parms->SiteRep = true;	PrintMessage(" - Site represenation applied.", "", "", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "RandomSeed")) {
			fscanf(FileID, "%d ", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting primary random number seed to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->PriRandomSeed = ReadInt;
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting secondary random number seed to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->SecRandomSeed = ReadInt;
		}
		else if (!strcmp(Tag, "SOS")) {
			fscanf(FileID, "%s\n", Parms->SOSChar);
			PrintMessage(" - Limiting the sum-over-states to ", Parms->SOSChar, ".", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "SampleTotal")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the number of samples to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->SampleTotal = ReadInt;
		}
		else if (!strcmp(Tag, "TimeStep")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the time step size to ", DisplayNumber, " ps.", false, Parms->Embarrassing);
			Parms->TimeStep = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "ResponseTimeTotal")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting response time limit to ", DisplayNumber, " time steps.", false, Parms->Embarrassing);
			Parms->ResponseTimeTotal = ReadInt;
		}
		else if (!strcmp(Tag, "WaitTimeTotal")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting waiting time to ", DisplayNumber, " time steps.", false, Parms->Embarrassing);
			Parms->WaitTimeTotal = ReadInt;
		}
		else if (!strcmp(Tag, "UnitTotal")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the number of molecular units to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->UnitTotal = ReadInt;
		}
		else if (!strcmp(Tag, "InitialSite")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the initial site to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->InitialSite = ReadInt;
		}
		else if (!strcmp(Tag, "MolPerCell")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the number of molecules per unit cell to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->MolPerCell = ReadInt;
		}
		else if (!strcmp(Tag, "Periodic")) {
			Parms->Periodic = true; PrintMessage(" - Periodic boundaries enabled.", "", "", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "Dipoles")) {
			fscanf(FileID, "%s\n", Parms->PermDipolesFile);
			PrintMessage(" - Associating dipoles with file '", Parms->PermDipolesFile, "'.", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "XEnergies")) {
			fscanf(FileID, "%s\n", Parms->XEnergiesFile);
			PrintMessage(" - Associating quantum energies with file '", Parms->XEnergiesFile, "'.", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "XDiags")) {
			fscanf(FileID, "%s\n", Parms->XDiagsFile);
			PrintMessage(" - Associating quantum transition energies with file '", Parms->XDiagsFile, "'.", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "XCouplings")) {
			fscanf(FileID, "%s\n", Parms->XCouplingsFile);
			PrintMessage(" - Associating quantum couplings with file '", Parms->XCouplingsFile, "'.", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "XPointDipoleRadius")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the radius of point-dipole couplings to ", DisplayNumber, " units.", false, Parms->Embarrassing);
			Parms->XPointDipoleRadius = ReadInt;
		}
		else if (!strcmp(Tag, "XEnergyMean")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the average transition energy to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->XEnergyMean = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "XEnergyOffset")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the transition energy offset to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->XEnergyOffset = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "XCouplingNN")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the nearest-neighbor exciton coupling to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->XCouplingNN = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "XLifeTime")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the exciton life time to ", DisplayNumber, " ps.", false, Parms->Embarrassing);
			Parms->XLifeTime = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "HBroad")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting homogeneous broadening to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->HBroad = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "CT")) {
			PrintMessage(" - Charge transfer is enabled.", "", "", false, Parms->Embarrassing);
			Parms->EHRadius = -1;
		}
		else if (!strcmp(Tag, "EHRadius")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the maximum electron-hole separation to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->EHRadius = ReadInt;
		}
		else if (!strcmp(Tag, "HDomain")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Limit the hole domain to ", DisplayNumber, " sites.", false, Parms->Embarrassing);
			Parms->HDomain = ReadInt;
		}
		else if (!strcmp(Tag, "EHEnergyNN")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the nearest-neighbor electron-hole energy to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->EHEnergyNN = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "EHEnergyInf")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the fully-separated electron-hole energy to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->EHEnergyInf = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "ECouplingNN")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the nearest-neighbor electron coupling to ", DisplayNumber, " wvib.", false, Parms->Embarrassing);
			Parms->ECouplingNN = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "ECouplingNNRight")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the nearest-neighbor electron coupling on right to ", DisplayNumber, " wvib.", false, Parms->Embarrassing);
			Parms->ECouplingNNRight = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "ECouplingNNLeft")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the nearest-neighbor electron coupling on left to ", DisplayNumber, " wvib.", false, Parms->Embarrassing);
			Parms->ECouplingNNLeft = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "HCouplingNN")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the nearest-neighbor hole coupling to ", DisplayNumber, " wvib.", false, Parms->Embarrassing);
			Parms->HCouplingNN = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "HCouplingNNRight")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the nearest-neighbor hole coupling on right to ", DisplayNumber, " wvib.", false, Parms->Embarrassing);
			Parms->HCouplingNNRight = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "HCouplingNNLeft")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the nearest-neighbor hole coupling on left to ", DisplayNumber, " wvib.", false, Parms->Embarrassing);
			Parms->HCouplingNNLeft = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "VMax")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the maximum number of vibrations to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->VTotal = ReadInt + 1;
		}
		else if (!strcmp(Tag, "VRadius")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the vibrational radius to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->VRadius = ReadInt;
		}
		else if (!strcmp(Tag, "VSlope")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the vibrational slope to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->VSlope = ReadDouble;
		}
		else if (!strcmp(Tag, "VExtra")) {
			Parms->VExtra = true; PrintMessage(" - Extra vibrational state is included.", "", "", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "XHuangRhys")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the exciton Huang-Rhys factor to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->XHuangRhys = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "EHuangRhys")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the electron Huang-Rhys factor to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->EHuangRhys = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "HHuangRhys")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the hole Huang-Rhys factor to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->HHuangRhys = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "VEnergy")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the vibrational energy to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->VEnergy = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "HierarchyDepth")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the hierarchy depth to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->HierarchyDepth = ReadInt;
		}
		else if (!strcmp(Tag, "HierarchyDamping")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the hierarchy damping parameter to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->HierarchyDamping = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "QFeedback")) {
			Parms->QFeedback = true; PrintMessage(" - Quantum feedback is enabled.", "", "", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "BathStaticDev")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the static bath energy deviation to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->BathStaticDev = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "BathDynamicDev")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the dynamic bath energy deviation to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->BathDynamicDev = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "BathTime")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the bath correlation time to ", DisplayNumber, " ps.", false, Parms->Embarrassing);
			Parms->BathTime = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "BathReorganization")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the bath reorganization energy to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->BathReorganization = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "BathCoupling")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the bath coupling constant to ", DisplayNumber, " /cm /nm.", false, Parms->Embarrassing);
			Parms->BathCoupling = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "BathMass")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the bath mass to ", DisplayNumber, " u.", false, Parms->Embarrassing);
			Parms->BathMass = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "BathSpring")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the bath spring constant to ", DisplayNumber, " u/ps^2.", false, Parms->Embarrassing);
			Parms->BathSpring = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "BathDamping")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the bath damping to ", DisplayNumber, " /ps.", false, Parms->Embarrassing);
			Parms->BathDamping = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "TrapRate")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the trap rate to ", DisplayNumber, " /ps.", false, Parms->Embarrassing);
			Parms->TrapRate = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "Temperature")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the temperature to ", DisplayNumber, " K.", false, Parms->Embarrassing);
			Parms->Temperature = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "EnergyMin")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the spectral lower bound to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->EnergyMin = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "EnergyMax")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the spectral upper bound to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->EnergyMax = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "SpectrumTotal")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the spectral resolution to ", DisplayNumber, ".", false, Parms->Embarrassing);
			Parms->SpectrumTotal = ReadInt;
		}
		else if (!strcmp(Tag, "Pulse")) {
			fscanf(FileID, "%s\n", Parms->PulseFile);
			PrintMessage(" - Associating laser pulse with file '", Parms->PulseFile, "'.", false, Parms->Embarrassing);
		}
		else if (!strcmp(Tag, "PulseTimeTotal")) {
			fscanf(FileID, "%d\n", &ReadInt);
			IntToString(DisplayNumber, ReadInt);
			PrintMessage(" - Setting the maximum pulse length to ", DisplayNumber, " time steps.", false, Parms->Embarrassing);
			Parms->PulseTimeTotal = ReadInt;
		}
		else if (!strcmp(Tag, "CrystalShift")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the gas-to-crystal shift to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->CrystalShift = (RealType) ReadDouble;
		}
		else if (!strcmp(Tag, "MonomerShift")) {
			fscanf(FileID, "%lf\n", &ReadDouble);
			RealToString(DisplayNumber, ReadDouble);
			PrintMessage(" - Setting the monomeric energy shift to ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
			Parms->MonomerShift = (RealType) ReadDouble;
		}
	}
	fclose(FileID);
}


void DerivedParms(ParmStruct *Parms)
{
	char DisplayNumber[100];
	
	PrintMessage("Derived parameters are calculated.", "", "", true, Parms->Embarrassing);
	
	if (Parms->Embarrassing) {
		Parms->SecRandomSeed += Parms->Embarrassing;
		IntToString(DisplayNumber, Parms->SecRandomSeed);
		PrintMessage("Initializing random number generator with secondary seed ", DisplayNumber, ".", true, Parms->Embarrassing);
	}
	
	Parms->TwistAngle *= Pi / 180;
	
	if (Parms->XLifeTime) {
		Parms->HBroad = 1 / (4 * Pi * LightSpeed * Parms->XLifeTime);
		RealToString(DisplayNumber, Parms->HBroad);
		PrintMessage(" - Homogeneous broadening is ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
	}
	
	if (Parms->EHRadius == -1)	Parms->EHRadius = Parms->UnitTotal;
	
	if (Parms->VRadius == -1) Parms->VRadius = Parms->UnitTotal;
	if (Parms->EHuangRhys && !Parms->HHuangRhys) Parms->HHuangRhys = Parms->EHuangRhys;
	
	if (Parms->Temperature != 0) {
		Parms->Temperature *= BoltzmannCnst;
		RealToString(DisplayNumber, Parms->Temperature * KJPMolToInvCm);
		PrintMessage(" - Thermal energy equals ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
	}
	
	if (Parms->BathMass) {
		
		Parms->BathTempDev = sqrt(2 * Parms->Temperature * Parms->BathDamping * Parms->BathMass / Parms->TimeStep);	// [u nm/ps^2]
		
		Parms->BathDynamicDev = Parms->BathCoupling * sqrt(Parms->Temperature / Parms->BathSpring);
		RealToString(DisplayNumber, Parms->BathDynamicDev);
		PrintMessage(" - Dynamic bath energy deviation equals ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
		
		Parms->BathTime = Parms->BathMass * Parms->BathDamping / Parms->BathSpring;
		RealToString(DisplayNumber, Parms->BathTime);
		PrintMessage(" - Bath correlation time (over-damped limit) equals ", DisplayNumber, " ps.", false, Parms->Embarrassing);
		
		IntToString(DisplayNumber, 10 * Parms->BathDamping * sqrt(Parms->BathMass / Parms->BathSpring));
		PrintMessage(" - Over-damped limit reached for ", DisplayNumber, "%.", false, Parms->Embarrassing);
	}
	
	if (Parms->Temperature != 0 && Parms->BathReorganization == 0) {
		Parms->BathReorganization = RealPow2(Parms->BathDynamicDev) / (2 * Parms->Temperature * KJPMolToInvCm);
		RealToString(DisplayNumber, Parms->BathReorganization);
		PrintMessage(" - Bath reorganization energy equals ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
	}
	
	if (Parms->BathTime) {
		Parms->HierarchyDamping = FreqToInvCm / Parms->BathTime;
		RealToString(DisplayNumber, Parms->HierarchyDamping);
		PrintMessage(" - Hierarchy damping parameter equals ", DisplayNumber, " /cm.", false, Parms->Embarrassing);
		
		Parms->BathPreCalc[0] = exp(-Parms->TimeStep / Parms->BathTime), Parms->BathPreCalc[1] = sqrt(1 - RealPow2(Parms->BathPreCalc[0]));
		Parms->BathDynamics = true;
	}
	
	Parms->HierarchyPreCalc[0] = Parms->HierarchyDamping * InvCmToFreq * Parms->TimeStep;
		
	// Old
	//Parms->HierarchyPreCalc[1] = 2 * Parms->BathReorganization * Parms->Temperature * KJPMolToInvCm * InvCmToFreq * Parms->TimeStep;
	//Parms->HierarchyPreCalc[2] = Parms->BathReorganization * Parms->BathDamping * InvCmToFreq * Parms->TimeStep;
	// New
	//Parms->HierarchyPreCalc[1] = 2 * Parms->BathReorganization * Parms->Temperature * KJPMolToInvCm * RealPow2(InvCmToFreq) * RealPow2(Parms->TimeStep);
	//Parms->HierarchyPreCalc[2] = Parms->BathReorganization * Parms->HierarchyDamping * RealPow2(InvCmToFreq) * RealPow2(Parms->TimeStep);
	// Trial
	Parms->HierarchyPreCalc[1] = Parms->BathReorganization * Parms->Temperature * KJPMolToInvCm * InvCmToFreq * RealPow2(Parms->TimeStep);
	Parms->HierarchyPreCalc[2] = .5 * Parms->BathReorganization * Parms->HierarchyDamping * InvCmToFreq * RealPow2(Parms->TimeStep);

    Parms->ECouplingNN *= Parms->VEnergy;
    Parms->ECouplingNNRight *= Parms->VEnergy;
    Parms->ECouplingNNLeft *= Parms->VEnergy;
    Parms->HCouplingNN *= Parms->VEnergy;
    Parms->HCouplingNNRight *= Parms->VEnergy;
    Parms->HCouplingNNLeft *= Parms->VEnergy;
}
