#include "PMT15SD.hh"
#include "PMT15Hit.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include "TGraph.h"

PMT15SD::PMT15SD(G4String name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("Pmt15Collection");
  SetQETable();
}

PMT15SD::~PMT15SD()
{}

void PMT15SD::Initialize(G4HCofThisEvent *HCTE)
{
  static int HCID = -1;
  Pmt15Collection = new PMT15HitsCollection(GetName(), collectionName[0]);

  if (HCID<0) HCID = GetCollectionID(0);

  HCTE->AddHitsCollection(HCID, Pmt15Collection);


  

}

G4bool PMT15SD::ProcessHits(G4Step *astep, G4TouchableHistory *ROhist)
{
  
  const G4StepPoint* preStepPoint = astep-> GetPreStepPoint();
  G4Track* atrack = astep->GetTrack();
  G4int pid = astep->GetTrack()->GetDefinition()-> GetPDGEncoding();
  G4ThreeVector worldPos =preStepPoint->GetPosition();
  G4ThreeVector pos = preStepPoint->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(worldPos);

  G4double energy = atrack->GetTotalEnergy();
  G4ThreeVector hitmom = atrack->GetMomentum()/CLHEP::eV;
  G4double E_p = hitmom.mag();

  const G4double h = 6.628e-34;
  const G4double c = 3.0e+8;
  G4double wavelength = ((h*c)/(energy*1.6e-13))*(1e+9); //nm

  G4int copynumber = preStepPoint->GetTouchableHandle()->GetCopyNumber();

  
  
  
  auto hitTime = astep->GetPreStepPoint()->GetGlobalTime();


  atrack->SetTrackStatus(fStopAndKill);

  G4double random = G4UniformRand();

  G4double qe_value = QETable->Eval(E_p);
  if (qe_value < 0 || qe_value > 1) qe_value = 0;

  if(random<qe_value){
    PMT15Hit* ahit = new PMT15Hit(pos,worldPos, hitTime,pid,wavelength, copynumber);
    Pmt15Collection->insert(ahit);
  }
  return true;
      

}


void PMT15SD::EndOfEvent(G4HCofThisEvent *HCTE)
{
  Pmt15Collection->PrintAllHits();

}

void PMT15SD::SetQETable(){

  //Original 400K
  /*
  const int num_qe = 13;
  G4double pmt15_ep1[num_qe] = {1.72,1.74,1.90,1.98,2.17,2.38,2.78,3.29,3.73,4.04,4.44,4.78,5.08};
  G4double pmt15_effi[num_qe] = {0.0,0.0002,0.004,0.015,0.059,0.13,0.22,0.24,0.19,0.12,0.019,0.0047,0.0};
  */

 

  //400S PMT -> UV enhanced
  const int num_qe = 46;
  G4double pmt15_ep1[num_qe] = {1.74,1.77,1.79,1.81,1.82,1.84,1.86,1.88,1.90,1.92,1.94,1.97,1.99,2.02,2.06,2.10,2.16,2.21,2.25,2.28,2.32,2.36,2.41,2.46,2.52,2.60,2.69,2.80,2.90,3.09,3.32,3.64,3.87,3.99,4.15,4.30,4.41,4.56,4.76,4.94,5.13,5.33,5.52,5.76,5.99,6.18};
  G4double pmt15_effi[num_qe] = {0.0002,0.0003,0.0005,0.0007,0.0010,0.0014,0.0020,0.0028,0.0041,0.0057,0.0080,0.0115,0.0154,0.0211,0.0277,0.0373,0.0496,0.0655,0.0779,0.0903,0.1047,0.1193,0.1392,0.1555,0.1760,0.1906,0.2052,0.2169,0.2264,0.2377,0.2405,0.2403,0.2314,0.2243,0.2107,0.1956,0.1838,0.1644,0.1507,0.1356,0.1251,0.1140,0.1078,0.1091,0.1118,0.1167};





  
  

  QETable = new TGraph(num_qe,pmt15_ep1,pmt15_effi);

  
}
  
