#ifndef PMT15Hit_hh
#define PMT15Hit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"

class G4AttDef;
class G4AttValue;

class PMT15Hit : public G4VHit
{
public:
  PMT15Hit();
  PMT15Hit(G4ThreeVector& axyz, G4double t);
  PMT15Hit(G4ThreeVector& axyz, G4ThreeVector& wxyz, G4double t, G4int pid, G4double Wavelength, G4int cn);
  virtual ~PMT15Hit();

  //copy constructor & assignment operator
  PMT15Hit(const PMT15Hit& right);
  const PMT15Hit& operator=(const PMT15Hit &right);
  //G4bool operator==(const PMT15Hit &right) const;

  //new/delete operators
  void *operator new(size_t);
  void operator delete(void *aHit);

  /*
  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;
  */

  const G4ThreeVector& GetPosition() const { return xyz; }
  const G4ThreeVector& GetWorldPosition() const { return worldxyz; }
  G4double GetTOF() const { return tof; }
  //void SetEdep(G4double de) {fEdep = de;}
  //void AddEdep(G4double de) {fEdep += de;}
  //G4double GetEdep() const {return fEdep;}
  G4double GetWavelength() const {return wavelengthMP; }

  G4int GetParticleID() const {return particleID;}

  //void SetTime(G4double dt) {fTime = dt;}
  //G4double GetTime() const {return fTime;}

  void IncPhotonCount() {fPhotons++;}
  G4double GetPhotonCount() {return fPhotons;}
  void ClearPhotonCount() {fPhotons=0.0;}

  G4int GetCopyNum() const {return copynum;}




  //void SetPos(G4ThreeVector xyz) {fPos = xyz;}
  //G4ThreeVector GetPos() const {return fPos;}
  
private:
  G4ThreeVector xyz;
  G4ThreeVector worldxyz;
  G4double tof;
  //G4int fId;
  G4int particleID;
  G4double wavelengthMP;
  //G4double fEdep;
  //G4double fTime;
  G4double fPhotons;
  //G4ThreeVector fPos;
  G4double copynum;

};

inline PMT15Hit::PMT15Hit(const PMT15Hit& right)
  : G4VHit()
{
  xyz = right.xyz;
  worldxyz = right.worldxyz;
  tof = right.tof;
  wavelengthMP = right.wavelengthMP;
  copynum = right.copynum;
}

inline const PMT15Hit& PMT15Hit::operator=
(const PMT15Hit& right)
{
  xyz = right.xyz;
  worldxyz = right.worldxyz;
  tof = right.tof;
  wavelengthMP = right.wavelengthMP;
  copynum = right.copynum;
  return *this;
}

typedef G4THitsCollection<PMT15Hit> PMT15HitsCollection;
extern G4Allocator<PMT15Hit> PMT15HitAllocator;


/*
using PMT15HitsCollection = G4THitsCollection<PMT15Hit>;

extern G4ThreadLocal G4Allocator<PMT15Hit>* PMT15HitAllocator;
*/
inline void* PMT15Hit::operator new(size_t)
{
  return static_cast<void*>(PMT15HitAllocator.MallocSingle());
}

inline void PMT15Hit::operator delete(void* aHit)
{
  PMT15HitAllocator.FreeSingle(static_cast<PMT15Hit*> (aHit));
}

#endif
  
