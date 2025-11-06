// -*- C++ -*-

#ifndef SAC3_HIT_HH
#define SAC3_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class SAC3Hit : public G4VHit, public VHitInfo
{
public:
  SAC3Hit(const G4String& name, G4Step* step);
  virtual ~SAC3Hit();

  SAC3Hit(const SAC3Hit& right);
  const SAC3Hit& operator=(const SAC3Hit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
SAC3Hit::SAC3Hit(const SAC3Hit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const SAC3Hit&
SAC3Hit::operator =(const SAC3Hit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<SAC3Hit> SAC3HitAllocator;

//_____________________________________________________________________________
inline void*
SAC3Hit::operator new(size_t)
{
  void* aHit = (void*)SAC3HitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
SAC3Hit::operator delete(void* aHit)
{
  SAC3HitAllocator.FreeSingle((SAC3Hit*) aHit);
}

#endif
