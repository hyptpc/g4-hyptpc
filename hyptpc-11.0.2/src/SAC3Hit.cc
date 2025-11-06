// -*- C++ -*-

#include "SAC3Hit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<SAC3Hit> SAC3HitAllocator;

//_____________________________________________________________________________
SAC3Hit::SAC3Hit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
SAC3Hit::~SAC3Hit()
{
}

//_____________________________________________________________________________
void
SAC3Hit::Draw()
{
}

//_____________________________________________________________________________
void
SAC3Hit::Print()
{
  VHitInfo::Print();
}
