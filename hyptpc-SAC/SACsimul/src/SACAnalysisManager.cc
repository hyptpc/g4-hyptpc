#include "SACAnalysisManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "AeroHit.hh"
#include "PMT15Hit.hh"

#include "Randomize.hh"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

#include <string>
#include <sstream>

extern int gCerenkovCounter;
extern double wavecheck[10000];

SACAnalysisManager::SACAnalysisManager(const G4String & histname)
  :outfile(histname), fActive_(true)
{}

SACAnalysisManager::~SACAnalysisManager()
{
  SaveFile();
}

void SACAnalysisManager::SaveFile(void) const
{
  if (fActive_)
    hfile->Write();
}

void SACAnalysisManager::Terminate(void) const
{
  if(fActive_)
    {
      hfile->Write();
      hfile->Close();
    }
}

void SACAnalysisManager::BeginOfRun(const G4Run*)
{

  

  G4SDManager* SDManager = G4SDManager::GetSDMpointer();

  hfile = new TFile(outfile, "RECREATE");
  tree = new TTree("tree","EvtGen tree");

  tree->Branch("event",&event, "event/I");
  tree->Branch("nEvt",&nEvt, "nEvt/I");
  //tree->Branch("evtid",evtid,"evtid[nEvt]/I");
  //tree->Branch("evtpid",&evtpid,"evtpid/I");
  tree->Branch("evtposx",&evtposx,"evtposx/D");
  tree->Branch("evtposy",&evtposy,"evtposy/D");
  //tree->Branch("evtposz",&evtposz,"evtposz/D");
  //tree->Branch("evtpid",evtpid,"evtpid[nEvt]/I");
  tree->Branch("evtnumce",&evtnumce,"evtnumce/I");
  tree->Branch("evtwavelength",evtwavelength,"evtwavelength[evtnumce]/D");


  //Aerogel
  
  tree->Branch("nhAero",&nhAero,"nhAero/I");
  //tree->Branch("aeroangle",aeroangle,"aeroangle[nhAero]/D");
  /*
  tree->Branch("aeropid",aeropid,"aeropid[nhAero]/I");
  tree->Branch("aeroposx",aeroposx,"aeroposx[nhAero]/D");
  tree->Branch("aeroposy",aeroposy,"aeroposy[nhAero]/D");
  tree->Branch("aeroposz",aeroposz,"aeroposz[nhAero]/D");
  tree->Branch("aerotime",aerotime,"aerotime[nhAero]/D");
  */
  

  //PMT15
  tree->Branch("nhPmt15",&nhPmt15,"nhPmt15/I");
  //tree->Branch("mppcmulti"&mppcmulti,"mppcmulti/I");
  //tree->Branch("mppcpid",mppcpid,"mppcpid[nhMppc]/I");
  //tree->Branch("pmt15posx",pmt15posx,"pmt15posx[nhPmt15]/D");
  //tree->Branch("pmt15posy",pmt15posy,"pmt15posy[nhPmt15]/D");
  //tree->Branch("pmt15posz",pmt15posz,"pmt15posz[nhPmt15]/D");
  //tree->Branch("mppctime",mppctime,"mppctime[nhMppc]/D");
  tree->Branch("pmt15wavelength",pmt15wavelength,"pmt15wavelength[nhPmt15]/D");
  tree->Branch("pmt15num",pmt15num,"pmt15num[nhPmt15]/I");


  event = 0;
  nEvt = 0;


}

void SACAnalysisManager::EndOfRun(const G4Run*)
{
  tree->Write();
  hfile->Write();
  hfile->Close();
}

void SACAnalysisManager::BeginOfEvent(const G4Event* anEvent)
{

}

void SACAnalysisManager::EndOfEvent(const G4Event* anEvent)
{

  G4HCofThisEvent* HCTE = anEvent->GetHCofThisEvent();
  if(!HCTE) return;
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  G4int nhaero = 0;
  G4int nhpmt15 = 0;
  G4int multiplicity = 0;

  G4int pdg = anEvent->GetPrimaryVertex(0)->GetPrimary(0)->GetPDGcode();
  evtpid = pdg;
  
  G4ThreeVector pripos = anEvent->GetPrimaryVertex(0)->GetPosition();
  G4double prix = pripos.x();
  G4double priy = pripos.y();
  G4double priz = pripos.z();
  evtposx = prix;
  evtposy = priy;
  evtposz = priz;

  evtnumce = gCerenkovCounter;
  for(int i=0;i<gCerenkovCounter;i++)evtwavelength[i] = wavecheck[i];
  
  

  PMT15HitsCollection *PMT15HC = 0;
  G4int ColIdPMT15 = SDMan->GetCollectionID("Pmt15Collection");
  if(ColIdPMT15>=0)
    {
      PMT15HC=dynamic_cast<PMT15HitsCollection *>(HCTE->GetHC( ColIdPMT15 ));
      if(PMT15HC)
	{
	  nhpmt15 = PMT15HC->entries();
	}
    }

  for(int i=0;i<nhpmt15;i++)
    {
      PMT15Hit* aHit = (*PMT15HC)[i];
      pmt15time[i] = aHit->GetTOF();
      pmt15posx[i] = aHit->GetPosition().x();
      pmt15posy[i] = aHit->GetPosition().y();
      //if(aHit->GetWorldPosition().x()<=0)mppcposz[i]=-1.0;
      //else if(aHit->GetWorldPosition().x()>0)mppcposz[i]=1.0;
      pmt15posz[i] = aHit->GetPosition().z();
      //mppcposz[i] = aHit->GetWorldPosition().y();
      pmt15pid[i] = aHit->GetParticleID();
      pmt15wavelength[i] = aHit->GetWavelength();
      pmt15num[i] = aHit-> GetCopyNum();
      
    }
  nhPmt15 = nhpmt15;


  AeroHitsCollection *AEROHC = 0;
  G4int ColIdAERO = SDMan->GetCollectionID("AeroCollection");
  if(ColIdAERO>=0)
    {
      AEROHC=dynamic_cast<AeroHitsCollection *>(HCTE->GetHC(ColIdAERO));
      if(AEROHC)
	{
	  nhaero = AEROHC->entries();
	}
    }

  for(int i=0;i<nhaero;i++)
    {
      AeroHit* aHit = (*AEROHC)[i];
      aeroangle[i] = aHit->GetAngle();
      aerotime[i] = aHit->GetTOF();
      aeroposx[i] = aHit->GetPosition().x();
      aeroposy[i] = aHit->GetPosition().y();
      aeroposz[i] = aHit->GetPosition().z();
      aeropid[i] = aHit->GetParticleID();

    }
  nhAero = nhaero;

  tree->Fill();
  event++;

  nEvt=0;
  nhaero = 0;
  nhpmt15 = 0;


  
}

void SACAnalysisManager::SetEvtGen(int j,int partnum)
{
  if(nEvt<j) nEvt=j;

}


  
      



