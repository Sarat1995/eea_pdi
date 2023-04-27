#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "RandomMod.h"
#include "ToolsMod.h"
#include "BasisMod.h"
#include "OverlapsMod.h"
#include "ParmsHandleMod.h"
#include "SubsMod.h"

#include <string.h>
#include <math.h>


typedef struct {
	
	int SampleCount, ProgressPercent, ProgressSample;
	
} MonteCarloStruct;


static void Debug(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const PolStruct *Pols, const VOverlapStruct VOverlaps);
static void Direct1D(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const PolStruct *Pols, const VOverlapStruct VOverlaps);
static void Response1D(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, PolStruct *Pols, const VOverlapStruct VOverlaps);
static void Response1D_SH(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const VOverlapStruct VOverlaps);
static void Response1D_HEOM(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const VOverlapStruct VOverlaps);
static void Response2D_SH(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void Response2D_HEOM(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const VOverlapStruct VOverlaps);
static void PopTransfer(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void PopTransfer_SH(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void PopTransfer_HEOM(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void TempCoh(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void TempCoh_(ParmStruct *Parms, const RealType *XEnergies, const RealType *PermDipoles, const VOverlapStruct VOverlaps);
static void TempCoh_Site(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void TempCoh_SH(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void TempCoh_SH_Site(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void TempCoh_HEOM(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void TempCoh_HEOM_Site(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);
static void Annihilation(ParmStruct *Parms, const RealType *XEnergies, const VOverlapStruct VOverlaps);