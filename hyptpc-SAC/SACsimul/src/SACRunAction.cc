#include "SACRunAction.hh"
#include "SACAnalysisManager.hh"
#include "SACPrimaryGeneratorAction.hh"
#include "SACDetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <stdlib.h>

SACRunAction::SACRunAction(SACAnalysisManager *analysisManager)
  :G4UserRunAction(), anaMan(analysisManager)
{
  
}
	     

SACRunAction::~SACRunAction()
{
}

void SACRunAction::BeginOfRunAction(const G4Run *aRun)
{
  if(G4VVisManager::GetConcreteInstance()){
    G4UImanager *UI = G4UImanager::GetUIpointer();
  }
  if(anaMan) anaMan->BeginOfRun(aRun);
}


  


void SACRunAction::EndOfRunAction(const G4Run *aRun)
{
  if(G4VVisManager::GetConcreteInstance()){
    G4UImanager *UI = G4UImanager::GetUIpointer();
  }
  if (anaMan) anaMan->EndOfRun(aRun);
}


