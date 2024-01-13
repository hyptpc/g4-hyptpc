// -*- C++ -*-

#include "EventAction.hh"

#include <fstream>

#include <G4RunManager.hh>
#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4VVisManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4VTrajectory.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>

#include "AnaManager.hh"
#include "BACSD.hh"
#include "ConfMan.hh"
#include "DetectorConstruction.hh"
#include "FuncName.hh"
#include "BH2SD.hh"
#include "FTOFSD.hh"
#include "HTOFSD.hh"
#include "BACSD.hh"
#include "TPCSD.hh"
#include "SCHSD.hh"
#include "SDCSD.hh"
#include "TargetSD.hh"
#include "VPSD.hh"
#include "KVCSD.hh"

#define SDCLASS(x) x ## Hit

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
EventAction::EventAction()
  : G4UserEventAction()
{
}

//_____________________________________________________________________________
EventAction::~EventAction()
{
}

//_____________________________________________________________________________
void
EventAction::BeginOfEventAction(const G4Event*)
{
#ifdef DEBUG
  G4cout << FUNC_NAME << G4endl;
#endif
  gAnaMan.BeginOfEventAction();
}

//_____________________________________________________________________________
void
EventAction::EndOfEventAction(const G4Event* anEvent)
{
  G4int eventID = anEvent->GetEventID();
  if (eventID % 100 == 0) {
    G4cout << FUNC_NAME << G4endl
	   << "   Event number = " << eventID << G4endl;
  }
  auto SDManager = G4SDManager::GetSDMpointer();
  auto HCTE = anEvent->GetHCofThisEvent();
  if (!HCTE) return;

  {
    static const auto id = SDManager->GetCollectionID("BH2/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<BH2Hit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("BH2", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("TGT/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<TargetHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("TGT", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("HTOF/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<HTOFHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("HTOF", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("FTOF/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<FTOFHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("FTOF", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("BAC/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<BACHit>*>(HCTE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("BAC", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("KVC/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<KVCHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("KVC", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("VP/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<VPHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("VP", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("TPC/hit");
    if(id > 0){
      auto HC = dynamic_cast<G4THitsCollection<TPCHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("TPC", HC->entries());
    }
  }

#if 0
  auto trajectoryContainer = anEvent->GetTrajectoryContainer();
  if(trajectoryContainer && G4VVisManager::GetConcreteInstance()){
    G4int n_trajectories = trajectoryContainer->entries();
    for(G4int i=0; i<n_trajectories; ++i){
      auto trj = (G4Trajectory*)((*(anEvent->GetTrajectoryContainer()))[i]);
      trj->DrawTrajectory();
    }
  }
#endif

  gAnaMan.EndOfEventAction();
}
