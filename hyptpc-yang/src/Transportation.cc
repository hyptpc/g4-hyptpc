/*
  Transportation.cc
*/

#include "Transportation.hh"
#include "G4SystemOfUnits.hh"

const double MinimumStep = 1.0*cm;

G4double Transportation::
AlongStepGetPhysicalInteractionLength( const G4Track &track,
                                       G4double previousStepSize, 
                                       G4double currentMinimumStep,
                                       G4double &currentSafety, 
                                       G4GPILSelection *selection )
{
  
  //if( DoesGlobalFieldExist() &&
  //    currentMinimumStep>MinimumStep ) 
  //  currentMinimumStep=MinimumStep;

  //geant4 10.07-p03
  if( DoesAnyFieldExist() &&
      currentMinimumStep>MinimumStep ) 
    currentMinimumStep=MinimumStep;

  
  return G4Transportation::
    AlongStepGetPhysicalInteractionLength( track, previousStepSize, 
                                           currentMinimumStep, currentSafety, 
                                           selection );
}
