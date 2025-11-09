// -*- C++ -*-

#ifndef SFV_HIT_HH
#define SFV_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class SFVHit : public G4VHit, public VHitInfo
{
public:
  SFVHit(const G4String& name, G4Step* step);
  virtual ~SFVHit();

  SFVHit(const SFVHit& right);
  const SFVHit& operator=(const SFVHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
SFVHit::SFVHit(const SFVHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const SFVHit&
SFVHit::operator =(const SFVHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<SFVHit> SFVHitAllocator;

//_____________________________________________________________________________
inline void*
SFVHit::operator new(size_t)
{
  void* aHit = (void*)SFVHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
SFVHit::operator delete(void* aHit)
{
  SFVHitAllocator.FreeSingle((SFVHit*) aHit);
}

#endif