// ====================================================================
//   CH1Hit.hh
//
// ====================================================================
#ifndef CH1_HIT_H
#define CH1_HIT_H

#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4SystemOfUnits.hh"

class CH1Hit : public G4VHit {

  private:
    G4ThreeVector xyz;
    G4ThreeVector pxyz;
    G4double tof;
    G4double edep;
    G4int trackID;
    G4int particleID;
    G4int detectorID;
    G4double massSH;
    G4double qqSH;
    G4int parentID;
    G4ThreeVector vtxmome;
    G4ThreeVector vtxposi;
    G4double vtxene;
    G4double flength;
  public:
    CH1Hit();
    //  CH1Hit(G4ThreeVector& axyz, G4double t);
    CH1Hit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t);
    CH1Hit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
        G4int tid, G4int pid, G4int did,G4double mass, G4int qq, G4int parentid,
        G4ThreeVector& vtxmom, G4ThreeVector& vtxpos, G4double vtxene, G4double tlength);

    virtual ~CH1Hit();

    // copy constructor & assignment operator
    CH1Hit(const CH1Hit& right);
    const CH1Hit& operator=(const CH1Hit& right);

    // new/delete operators
    void* operator new(size_t);
    void operator delete(void* aHit);

    const G4ThreeVector& GetVtxPosition() const { return vtxposi; }
    const G4ThreeVector& GetVtxMomentum() const { return vtxmome; }
    G4double GetVtxEnergy() const { return vtxene; }

    // set/get functions
    const G4ThreeVector& GetPosition() const { return xyz; }
    const G4ThreeVector& GetMomentum() const { return pxyz; }
    G4double GetTOF() const { return tof; }

    G4int GetParentID() const { return parentID; }
    G4int GetTrackID() const { return trackID; }
    G4int GetParticleID() const { return particleID; }
    G4int GetDetectorID() const { return detectorID; }
    G4double GetParticleMassID() const { return massSH; }
    G4double GetParticleQqID() const { return qqSH; }
    void SetEdep(G4double aedep) { edep = aedep; }
    G4double GetEdep() const { return edep; }
    G4double GetLength() const { return flength; }

    // methods
    virtual void Draw();
    virtual void Print();
};

// ====================================================================
// inline functions
// ====================================================================
  inline CH1Hit::CH1Hit(const CH1Hit& right)
: G4VHit()
{
  xyz= right.xyz;
  pxyz= right.pxyz;
  tof= right.tof;
  edep = right.edep;
}

  inline const CH1Hit& CH1Hit::operator=
(const CH1Hit& right)
{
  xyz= right.xyz;
  pxyz= right.pxyz;
  tof= right.tof;
  edep = right.edep;
  return *this;
}

// externally instanciated.
typedef G4THitsCollection<CH1Hit> CH1HitsCollection;
extern G4Allocator<CH1Hit> CH1HitAllocator;

inline void* CH1Hit::operator new(size_t)
{
  void* aHit= (void*)CH1HitAllocator.MallocSingle();
  return aHit;
}

inline void CH1Hit::operator delete(void* aHit)
{
  CH1HitAllocator.FreeSingle((CH1Hit*) aHit);
}

#endif
