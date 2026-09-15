// -*- C++ -*-

#include "PhysicsList.hh"

#include <iomanip>

#include <CLHEP/Units/SystemOfUnits.h>

#include <globals.hh>
#include <G4ios.hh>

#include <G4DecayPhysics.hh>
#include <G4EmStandardPhysics.hh>
#include <G4EmExtraPhysics.hh>
#include <G4IonPhysics.hh>
#include <G4StoppingPhysics.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4NeutronTrackingCut.hh>
#include <G4HadronPhysicsQGSP_BERT.hh>
#include <G4VPhysicsConstructor.hh>

// IonGas models for TPC P10 region only (standard EM kept elsewhere)
#include <G4BraggIonGasModel.hh>
#include <G4BetheBlochIonGasModel.hh>
#include <G4IonFluctuations.hh>
#include <G4UniversalFluctuation.hh>
#include <G4EmParameters.hh>
#include <G4ProcessManager.hh>
#include <G4RegionStore.hh>
#include <G4VEnergyLossProcess.hh>

#include "ConfMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"

#define G4MT_physicsVector                                                    \
  ((G4VMPLsubInstanceManager.offset[g4vmplInstanceID]).physicsVector)

namespace
{
  const auto& confMan = ConfMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
PhysicsList::PhysicsList(G4int ver)
  : G4VUserPhysicsList()
{
  verboseLevel = ver;
  defaultCutValue = 0.7*CLHEP::mm;
  // confMan.Get<G4double>("DefaultCutValue");

  // EM Physics
  RegisterPhysics(new G4EmStandardPhysics(ver));
  // Synchrotron Radiation & GN Physics
  RegisterPhysics(new G4EmExtraPhysics(ver));

  // Decays
  RegisterPhysics(new G4DecayPhysics(ver));

  // Hadron Elastic scattering
  RegisterPhysics(new G4HadronElasticPhysics(ver));
  // Hadron Physics
  RegisterPhysics(new G4HadronPhysicsQGSP_BERT(ver));

  ///// Others
  // Stopping Physics
  RegisterPhysics(new G4StoppingPhysics(ver));
  // Ion Physics
  RegisterPhysics(new G4IonPhysics(ver));
  // Neutron tracking cut
  RegisterPhysics(new G4NeutronTrackingCut(ver));
}

//_____________________________________________________________________________
void
PhysicsList::ConstructParticle()
{
  for(auto itr = G4MT_physicsVector->cbegin();
      itr != G4MT_physicsVector->cend(); ++itr)
  {
    (*itr)->ConstructParticle();
  }
}

//_____________________________________________________________________________
void
PhysicsList::ConstructProcess()
{

  const G4int pad_configure = gSize.Get("TpcPadConfigure");

  // G4AutoLock l(&constructProcessMutex);
  AddTransportation();

  for(auto itr = G4MT_physicsVector->cbegin();
      itr != G4MT_physicsVector->cend(); ++itr)
  {
    auto name = (*itr)->GetPhysicsName();
    if(name == "G4EmStandard" || name == "G4GammaLeptoNuclearPhys"){
      if(!confMan.Get<G4bool>("EM")) continue;
    }else if(name == "Decay"){
      if(!confMan.Get<G4bool>("DECAY")) continue;
    }else{
      if(!confMan.Get<G4bool>("HADRON")) continue;
    }

    if(verboseLevel > 0)
      G4cout << FUNC_NAME << " Construct " << name << G4endl;

    (*itr)->ConstructProcess();

    if(pad_configure ==4){
      G4EmParameters* emParameters = G4EmParameters::Instance();
      emParameters->SetMinEnergy(10*CLHEP::eV);
      emParameters->SetMaxEnergy(2.*CLHEP::GeV);
      emParameters->SetNumberOfBinsPerDecade(100);
    }
  }

  // After standard EM: add IonGas ionization models only in TpcP10Region.
  if(pad_configure ==4)
    AddIonGasProcess();
}

//_____________________________________________________________________________
void PhysicsList::AddIonGasProcess()
{
  // configure==4: keep standard EM everywhere; overlay IonGas models
  // only inside TpcP10Region. e+/e- stay on standard eIoni.
  G4Region* gasRegion =
    G4RegionStore::GetInstance()->GetRegion("TpcP10Region", false);
  if(!gasRegion){
    G4cerr << FUNC_NAME
           << " TpcP10Region not found; skip IonGas models" << G4endl;
    return;
  }

  auto pIterator = GetParticleIterator();
  pIterator->reset();
  while((*pIterator)())
  {
    G4ParticleDefinition *pDefinition = pIterator->value();
    G4String pName = pDefinition->GetParticleName();
    const G4bool use_iongas =
      pName == "proton" || pName == "anti_proton" ||
      pName == "pi+" || pName == "pi-" ||
      pName == "kaon+" || pName == "kaon-" ||
      pName == "mu+" || pName == "mu-" ||
      pName == "sigma+" || pName == "sigma-" ||
      pName == "anti_sigma+" || pName == "anti_sigma-" ||
      pName == "xi-" || pName == "xi+" ||
      pName == "deuteron" || pName == "triton" ||
      pName == "He3" || pName == "alpha";
    if(!use_iongas)
      continue;

    auto* pmanager = pDefinition->GetProcessManager();
    if(!pmanager) continue;

    G4VProcess* ioni = pmanager->GetProcess("hIoni");
    if(!ioni) ioni = pmanager->GetProcess("muIoni");
    if(!ioni) ioni = pmanager->GetProcess("ionIoni");
    auto* eloss = dynamic_cast<G4VEnergyLossProcess*>(ioni);
    if(!eloss) continue;

    G4BraggIonGasModel *bIgm = new G4BraggIonGasModel();
    G4BetheBlochIonGasModel *bbIgm = new G4BetheBlochIonGasModel();

    bIgm->SetActivationHighEnergyLimit(
      2.*CLHEP::MeV*pDefinition->GetPDGMass()/CLHEP::proton_mass_c2);
    bbIgm->SetActivationLowEnergyLimit(
      2.*CLHEP::MeV*pDefinition->GetPDGMass()/CLHEP::proton_mass_c2);

    eloss->AddEmModel(0, bIgm, new G4IonFluctuations, gasRegion);
    eloss->AddEmModel(0, bbIgm, new G4UniversalFluctuation, gasRegion);

    // suppress delta rays only inside TPC P10 gas
    eloss->ActivateSecondaryBiasing("TpcP10Region", 1e-10, 100*CLHEP::TeV);
  }
}

//_____________________________________________________________________________
void
PhysicsList::SetCuts()
{
  SetCutsWithDefault();
}
