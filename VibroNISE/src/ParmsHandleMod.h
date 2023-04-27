#ifndef ParmsHandleMod
#define ParmsHandleMod

void InitializeParms(ParmStruct *Parms);
void ReadInputFile(char *Args[], ParmStruct *Parms, int *StartPoint);
void DerivedParms(ParmStruct *Parms);

#endif