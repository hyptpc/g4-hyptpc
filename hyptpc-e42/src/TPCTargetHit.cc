// -*- C++ -*-

#include "TPCTargetHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCTargetHit> TPCTargetHitAllocator;

//_____________________________________________________________________________
TPCTargetHit::TPCTargetHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCTargetHit::~TPCTargetHit( void )
{
}

//_____________________________________________________________________________
void
TPCTargetHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCTargetHit::Print( void )
{
  VHitInfo::Print();
}
