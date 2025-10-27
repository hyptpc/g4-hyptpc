// -*- C++ -*-

#ifndef CVC_HIT_HH
#define CVC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class CVCHit : public G4VHit, public VHitInfo
{
public:
  CVCHit(const G4String& name, G4Step* step);
  virtual ~CVCHit();

  CVCHit(const CVCHit& right);
  const CVCHit& operator=(const CVCHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
CVCHit::CVCHit(const CVCHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const CVCHit&
CVCHit::operator =(const CVCHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<CVCHit> CVCHitAllocator;

//_____________________________________________________________________________
inline void*
CVCHit::operator new(size_t)
{
  void* aHit = (void*)CVCHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
CVCHit::operator delete(void* aHit)
{
  CVCHitAllocator.FreeSingle((CVCHit*) aHit);
}

#endif