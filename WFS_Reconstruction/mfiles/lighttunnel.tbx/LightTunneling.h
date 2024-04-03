// tempus 2007A
#ifndef LIGHTTUNNELING_SYSTEM_CLASS
#define LIGHTTUNNELING_SYSTEM_CLASS
//
// (c) Copyright 1995-2007, MZA Associates Corporation
//

#include "tempus/tempus.h"


#include "c:\mza\wavetrain\trunk\sensors.lib/LightTunnel.h"

class LightTunneling : public System
{
public:

   // Parameters
   float magnification;

   // Subsystems
   LightTunnel lighttunnel;


   // Inputs
   Input< Grid<float> > reflectedIntensity;
   Input< Grid<float> > relativeIntensity;
   Input< Grid<float> > xDisplacement;
   Input< Grid<float> > yDisplacement;
   Input< Grid<float> > xSpread;
   Input< Grid<float> > ySpread;
   Input< Grid<float> > orientation;

   // Outputs

   LightTunneling(SystemNode *parent, char *name,
                  float _magnification) : 
      System(parent, name),
      magnification(_magnification),
      lighttunnel(this, "lighttunnel",
         1.0f,
         1.0f/magnification),
      reflectedIntensity(this, "reflectedIntensity"),
      relativeIntensity(this, "relativeIntensity"),
      xDisplacement(this, "xDisplacement"),
      yDisplacement(this, "yDisplacement"),
      xSpread(this, "xSpread"),
      ySpread(this, "ySpread"),
      orientation(this, "orientation")
   {
      lighttunnel.imgp <<= reflectedIntensity;
      lighttunnel.p <<= relativeIntensity;
      lighttunnel.rho <<= orientation;
      lighttunnel.by <<= ySpread;
      lighttunnel.bx <<= xSpread;
      lighttunnel.ty <<= yDisplacement;
      lighttunnel.tx <<= xDisplacement;
   }

   // void respondToScheduledEvent (const Event & /* event */);
   // void respondToChangedInputs();
   // void respondToInputWarning (InputBase* /* input */);
   // void respondToOutputRequest (const OutputBase* /* output */);
};

#endif // LIGHTTUNNELING_SYSTEM_CLASS
