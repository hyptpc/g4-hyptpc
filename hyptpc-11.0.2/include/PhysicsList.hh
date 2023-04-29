// -*- C++ -*-

#ifndef PHYSICS_LIST_H
#define PHYSICS_LIST_H

#include <G4VModularPhysicsList.hh>

//_____________________________________________________________________________
class PhysicsList : public G4VModularPhysicsList
{
public:
  PhysicsList(G4int verbose=0);
  virtual ~PhysicsList();

private:
  G4VPhysicsConstructor*              m_em_physics_list;
  std::vector<G4VPhysicsConstructor*> m_hadron_physics_list;

protected:
  void ConstructBaryons();
  void ConstructBosons();
  void ConstructEM();
  void ConstructGeneral();
  void ConstructHadron();
  void ConstructIons();
  void ConstructLeptons();
  void ConstructMesons();
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  void ConstructShortLived();
  void ConstructStableHyperons();
  virtual void SetCuts();
};

#endif
