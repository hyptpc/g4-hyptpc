// -*- C++ -*-

#include "CVCHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<CVCHit> CVCHitAllocator;

//_____________________________________________________________________________
CVCHit::CVCHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
CVCHit::~CVCHit()
{
}

//_____________________________________________________________________________
void
CVCHit::Draw()
{
}

//_____________________________________________________________________________
void
CVCHit::Print()
{
  VHitInfo::Print();
}