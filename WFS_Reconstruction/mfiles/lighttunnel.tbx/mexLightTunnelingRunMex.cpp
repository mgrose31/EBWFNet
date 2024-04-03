// tempus 2007A

//
// tempus, version 2007A, (c) Copyright 1995-2007, MZA Associates Corporation
// Generated on: Tue May 01 10:25:12 MDT 2007
//

#include "tempus/tempus.h"
#include "tnl/TempusInitializers.h"
#include "wtlib/LaserGridInitializers.h"
#include "recorder/FileSys.h"
#include "mli/mliIO.h"
#ifndef NO_TEMPUS_SMF_MONITOR
#include "tempus/TempusStatusSMF.h"
#endif
#include "mli/mliConv.h"
#include "mli/mliSystem.h"


#include "LightTunneling.h"

const int numParams   = 3;
const int numRecords = 2;
const char* sysRecordItems[numRecords] = {
   "lighttunnel.img  Grid<float>  true 0 0",
   "lighttunnel.fpaImage  Grid<float>  true 0 0"
};

const int numDefRecords = 0;
const char** sysDefRecordItems = NULL;

mxArray* mliLoadRecords()
{
   mxArray* a = mxCreateCellMatrix(1, numRecords);
   for(int ic=0; ic < numRecords; ic++) {
      mxSetCell(a, ic, mxCreateString(sysRecordItems[ic]));
   }
   return a;
}

mxArray* mliLoadDefRecords()
{
   mxArray* a = mxCreateCellMatrix(1, numDefRecords);
   for(int ic=0; ic < numDefRecords; ic++) {
      mxSetCell(a, ic, mxCreateString(sysDefRecordItems[ic]));
   }
   return a;
}

const char* sysParams[numParams] = {
   "FL",
   "z",
   "magnification"
};

const char* sysParamsVals[numParams] = {
   "float  FL  1.0f",
   "float  z  1.0e4f",
   "float  magnification  FL/z"
};

/* MEX Functionhandlers */

mxArray* mliLoadParams()
{
   mxArray* a = mxCreateStructMatrix(1, 1, numParams, sysParams);

   mxSetField(a, 0, sysParams[0], mxCreateString(sysParamsVals[0]));
   mxSetField(a, 0, sysParams[1], mxCreateString(sysParamsVals[1]));
   mxSetField(a, 0, sysParams[2], mxCreateString(sysParamsVals[2]));

   return a;
}

class MLILightTunnelingRunMex
{
public:

   // Run variables
   float  FL;
   float  z;

   // Parameters
   float  magnification;

   // External Inputs
   Object< Grid<float> > ___mliei_reflectedIntensity;
   Object< Grid<float> > ___mliei_relativeIntensity;
   Object< Grid<float> > ___mliei_xDisplacement;
   Object< Grid<float> > ___mliei_yDisplacement;
   Object< Grid<float> > ___mliei_xSpread;
   Object< Grid<float> > ___mliei_ySpread;
   Object< Grid<float> > ___mliei_orientation;

   // Member functions
   MLILightTunnelingRunMex(mxArray* mliparams,
                            float  _FL,
                            float  _z,
                            float  _magnification) :
                            FL(_FL),
                            z(_z),
                            magnification(_magnification),
                            ___mliei_reflectedIntensity(uniformGrid(1.0)),
                            ___mliei_relativeIntensity(uniformGrid(1.0)),
                            ___mliei_xDisplacement(uniformGrid(0.0)),
                            ___mliei_yDisplacement(uniformGrid(0.0)),
                            ___mliei_xSpread(uniformGrid(1.0)),
                            ___mliei_ySpread(uniformGrid(1.0)),
                            ___mliei_orientation(uniformGrid(0.0))   {
   }
   ~MLILightTunnelingRunMex()
   {
   }
   MLISystemHandle* createSystem(mxArray* filename, mxArray* verbose)
   {
      mxArray* mlia = NULL;

      double stopTime = 1.0e3f;

      bool verb = _mliGetVerbose(verbose);
      char *trffile;
#ifndef NO_TEMPUS_SMF_MONITOR
      char *smfile;
      TempusStatusSMFWriter* smfWriter = mliCreateSMF(filename, &trffile, &smfile);
#else
      mliCreateTRF(filename, &trffile)
#endif
      Universe* mliuniverse = new Universe("LightTunneling");
      LightTunneling* mlisystem = new LightTunneling(NULL, "LightTunneling",
                                                     magnification);

      // Connect external inputs
      mlisystem->reflectedIntensity <<= ___mliei_reflectedIntensity;
      mlisystem->relativeIntensity <<= ___mliei_relativeIntensity;
      mlisystem->xDisplacement <<= ___mliei_xDisplacement;
      mlisystem->yDisplacement <<= ___mliei_yDisplacement;
      mlisystem->xSpread <<= ___mliei_xSpread;
      mlisystem->ySpread <<= ___mliei_ySpread;
      mlisystem->orientation <<= ___mliei_orientation;

      MLISystemHandle* handle = new MLISystemHandle(mliuniverse, mlisystem, this);
      handle->SetTrfFileName(trffile);
#ifndef NO_TEMPUS_SMF_MONITOR
      handle->SetSMFFileName(smfile);
      handle->SetSMFWriter(smfWriter);
#endif
      if (handle->recorderFile == NULL) {
         writeParam(handle->pset, "stopTime", DOUBLE, &stopTime);
         writeParam(handle->pset, "FL", FLOAT, &FL);
         writeParam(handle->pset, "z", FLOAT, &z);
         writeParam(handle->pset, "magnification", FLOAT, &magnification);
         handle->OpenRecorderFile();
         if (handle->recorderFile == NULL) {
            if(verb)
               mexPrintf("Error -- Null recorder file may cause problems.\n");
         }
      }
      return(handle);
   }
};

MLILightTunnelingRunMex* mliLightTunnelingRunMexCompParams(mxArray* mliparamsin, mxArray*& mliparamsout, mxArray* mliverbose)
{
   bool verb = _mliGetVerbose(mliverbose);
   mxArray* mlia = NULL;

   double stopTime = 1.0e3f;

   float  FL = 1.0f;
   if(mliparamsin) mlia = mxGetField(mliparamsin, 0, "FL");
   if(!_mliCheckValue(mlia, sysParamsVals[0])) fromMxArray(mlia, FL);

   float  z = 1.0e4f;
   if(mliparamsin) mlia = mxGetField(mliparamsin, 0, "z");
   if(!_mliCheckValue(mlia, sysParamsVals[1])) fromMxArray(mlia, z);

   float  magnification = FL/z;
   if(mliparamsin) mlia = mxGetField(mliparamsin, 0, "magnification");
   if(!_mliCheckValue(mlia, sysParamsVals[2])) fromMxArray(mlia, magnification);

   mliparamsout = mxCreateStructMatrix(1, 1, numParams, sysParams);

   mxSetField(mliparamsout, 0, "FL", toMxArray(FL));
   mxSetField(mliparamsout, 0, "z", toMxArray(z));
   mxSetField(mliparamsout, 0, "magnification", toMxArray(magnification));

   return(new MLILightTunnelingRunMex(mliparamsout,
                                      FL,
                                      z,
                                      magnification));
}
mxArray* mliCompParams(mxArray* paramsin, mxArray* verbose)
{
   mxArray* paramsout;
   MLILightTunnelingRunMex* tmpcls = mliLightTunnelingRunMexCompParams(paramsin, paramsout, verbose);
   delete tmpcls;
   return paramsout;
}

mxArray* mliCreateSystem(mxArray* paramsin, mxArray* filename, mxArray* verbose)
{
   mxArray* paramsout;
   mxArray* runhandle;
   MLILightTunnelingRunMex* tmpcls = mliLightTunnelingRunMexCompParams(paramsin, paramsout, verbose);
   MLISystemHandle* handle = tmpcls->createSystem(filename, verbose);
   runhandle = mliMakePointer(handle);
   return(runhandle);
}

mxArray* mliRecordOutputs(mxArray* mhandle_in,  mxArray* outputs,
                          mxArray* mstatus,  mxArray* matLeast,  mxArray* mexactly,
                          mxArray* verbose)
{
   mxArray *mhandle = mxDuplicateArray(mhandle_in);

   if (!_mliVerifyOutputs(outputs)) {
      mexPrintf("Outputs should be a list of structures.\n");
      return mhandle;
   }
   bool verb = _mliGetVerbose(verbose);
   mxClassID outType = mxGetClassID(outputs);
   int nOutputs = 0;
   char outstr[256] = "";
   mxArray* cname = outputs;
   bool   status      = true;
   float  atLeastRate = 0.0;
   double exactlyRate = 0.0;

   MLISystemHandle* handle = reinterpret_cast<MLISystemHandle*>(mliGetPointer(mhandle));
   LightTunneling* system = reinterpret_cast<LightTunneling*>(handle->system);

   if(outType == mxCELL_CLASS)
      nOutputs = mxGetN(outputs);
   else nOutputs = 1;

   handle->NewRecorders(nOutputs);

   for (int ic = 0; ic < nOutputs; ic++) 
   {
      atLeastRate = mxGetPr(matLeast)[ic];
      exactlyRate = mxGetPr(mexactly)[ic];
      status      = mxGetPr(mstatus) [ic] > 0.0;
      if (outType == mxCELL_CLASS) cname = mxGetCell(outputs, ic); 

      if (mxGetString(cname, outstr, 256) != 0) {
         mexPrintf("Error -- Not a string.\n");
      }
      if(strcmp(outstr, "lighttunnel.img") == 0) {
         if(verb) mexPrintf("Creating recorder %s\n", outstr);

         GridRecorder<float>*  rfMex1 = new GridRecorder<float>(NULL, "rfMex1", "lighttunnel.img", "Grid<float>", "imgp modified with tilt, blur, and intensity", status, atLeastRate, exactlyRate);
         if(handle->recorderFile != NULL)
            rfMex1->dr <<= handle->recorderFile->dr;
         rfMex1->i <<= system->lighttunnel.img;
         handle->AddRecorder(rfMex1);
      }
      if(strcmp(outstr, "lighttunnel.fpaImage") == 0) {
         if(verb) mexPrintf("Creating recorder %s\n", outstr);

         GridRecorder<float>*  rfMex2 = new GridRecorder<float>(NULL, "rfMex2", "lighttunnel.fpaImage", "Grid<float>", "Synthetic camera image formed by flipping & magnifying output to simulate camera geometry", status, atLeastRate, exactlyRate);
         if(handle->recorderFile != NULL)
            rfMex2->dr <<= handle->recorderFile->dr;
         rfMex2->i <<= system->lighttunnel.fpaImage;
         handle->AddRecorder(rfMex2);
      }
      if (verb) mexPrintf("Output %s not found.\n", outstr);
   }

   return mhandle;
}

mxArray* mliStopRecOutputs(mxArray* mhandle_in, mxArray* outputs, mxArray* verbose)
{
   mxArray *mhandle = mxDuplicateArray(mhandle_in);

   if (!_mliVerifyOutputs(outputs)) {
      mexPrintf("Outputs should be a list of structures.\n");
      return mhandle;
   }  
   bool verb = _mliGetVerbose(verbose);
   int       nOutputs = mxGetN(outputs);
   mxClassID outType  = mxGetClassID(outputs);
   char      outstr[256] = "";
   mxArray*  cname = outputs;
   MLISystemHandle* handle = reinterpret_cast<MLISystemHandle*>(mliGetPointer(mhandle));
   for (int ic = 0; ic < nOutputs; ic++) 
   {
      RecordableObject* recobj = NULL;
      if (outType == mxCELL_CLASS) cname = mxGetCell(outputs, ic);
      if (mxGetString(cname, outstr, 256) == 0) 
          recobj = handle->RemoveRecorder(outstr);  
      if (recobj == NULL) {
         if (verb) mexPrintf("Cannot find recorder %s.\n", outstr);
         continue;
      }
      if(strcmp(outstr, "lighttunnel.img") == 0) {
         if(verb) mexPrintf("Detaching recorder %s\n", outstr);

         GridRecorder<float>*  rfMex1= (GridRecorder<float>*)reinterpret_cast< RecordableObject* >(recobj);
           rfMex1->i.detach();
      }
      if(strcmp(outstr, "lighttunnel.fpaImage") == 0) {
         if(verb) mexPrintf("Detaching recorder %s\n", outstr);

         GridRecorder<float>*  rfMex2= (GridRecorder<float>*)reinterpret_cast< RecordableObject* >(recobj);
           rfMex2->i.detach();
      }
   }
   return mhandle;
}

mxArray* mliGetOutput(mxArray* mhandle, mxArray* output, mxArray* verbose)
{
   bool verb = _mliGetVerbose(verbose);
   MLISystemHandle* handle = reinterpret_cast<MLISystemHandle*>(mliGetPointer(mhandle));
   LightTunneling* system = reinterpret_cast<LightTunneling*>(handle->system);
   char outstr[256] = "";
   if (mxGetString(output, outstr, 256) != 0) {
      mexPrintf("Output name should be a string.\n");
   }  
   if(strcmp(outstr, "lighttunnel.img") == 0) {
      return toMxArray(system->lighttunnel.img.value());
   }
   if(strcmp(outstr, "lighttunnel.fpaImage") == 0) {
      return toMxArray(system->lighttunnel.fpaImage.value());
   }
   mexPrintf("Output %s not found.\n", outstr);
   return toMxArray(0);
}
mxArray* mliSetInputs(mxArray* mhandle_in, mxArray* inputs, mxArray* values, mxArray* verbose)
{
   mxArray *mhandle = mxDuplicateArray(mhandle_in);

   bool verb = _mliGetVerbose(verbose);
   mxClassID inType = mxGetClassID(inputs);
   mxClassID valType = mxGetClassID(values);
   int nInputs = 0;
   char instr[256] = "";
   mxArray* cname = inputs;
   mxArray* val = values;
   
   MLISystemHandle* handle = reinterpret_cast<MLISystemHandle*>(mliGetPointer(mhandle));
   MLILightTunnelingRunMex* runset = reinterpret_cast<MLILightTunnelingRunMex*>(handle->parameters);
   
   if(inType == mxCELL_CLASS) nInputs = mxGetN(inputs);
   else nInputs = 1;

   for(int ic = 0; ic < nInputs; ic++)
   {
      if(inType == mxCELL_CLASS) cname = mxGetCell(inputs, ic);
      if(valType == mxCELL_CLASS) val = mxGetCell(values, ic);
      
      if(mxGetString(cname, instr, 256) != 0) mexPrintf("Error -- Not a string.\n");

      if(strcmp(instr, "reflectedIntensity") == 0) {
         if(verb) mexPrintf("Setting input %s\n", instr);
         fromMxArray(val, runset->___mliei_reflectedIntensity);
      }
      else {
         if(verb) mexPrintf("Input %s not found.\n", instr);
      }
   
      if(strcmp(instr, "relativeIntensity") == 0) {
         if(verb) mexPrintf("Setting input %s\n", instr);
         fromMxArray(val, runset->___mliei_relativeIntensity);
      }
      else {
         if(verb) mexPrintf("Input %s not found.\n", instr);
      }
   
      if(strcmp(instr, "xDisplacement") == 0) {
         if(verb) mexPrintf("Setting input %s\n", instr);
         fromMxArray(val, runset->___mliei_xDisplacement);
      }
      else {
         if(verb) mexPrintf("Input %s not found.\n", instr);
      }
   
      if(strcmp(instr, "yDisplacement") == 0) {
         if(verb) mexPrintf("Setting input %s\n", instr);
         fromMxArray(val, runset->___mliei_yDisplacement);
      }
      else {
         if(verb) mexPrintf("Input %s not found.\n", instr);
      }
   
      if(strcmp(instr, "xSpread") == 0) {
         if(verb) mexPrintf("Setting input %s\n", instr);
         fromMxArray(val, runset->___mliei_xSpread);
      }
      else {
         if(verb) mexPrintf("Input %s not found.\n", instr);
      }
   
      if(strcmp(instr, "ySpread") == 0) {
         if(verb) mexPrintf("Setting input %s\n", instr);
         fromMxArray(val, runset->___mliei_ySpread);
      }
      else {
         if(verb) mexPrintf("Input %s not found.\n", instr);
      }
   
      if(strcmp(instr, "orientation") == 0) {
         if(verb) mexPrintf("Setting input %s\n", instr);
         fromMxArray(val, runset->___mliei_orientation);
      }
      else {
         if(verb) mexPrintf("Input %s not found.\n", instr);
      }
   
   }
   return mhandle;
}
