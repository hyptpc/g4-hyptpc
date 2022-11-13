#include "AnalysisManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "TpcHit.hh"
#include "TofHit.hh"
#include "TargetHit.hh"
#include "TargetHit.hh"
//#include "TrackerHit.hh"
#include "Randomize.hh"
#include "FACHit.hh"
#include "BDCHit.hh"
#include "DCHit.hh"
#include "CH1Hit.hh"
#include "FTOFHit.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TString.h"

#include "EvtGenBase/EvtVector4R.hh"

//#include "RootHelper.hh"
#include <string>
#include <sstream>

#include "TF1.h"
#include "TMath.h"

AnalysisManager::AnalysisManager( const G4String & histname )
  :outfile(histname), fActive_(true)
{
}

AnalysisManager::~AnalysisManager()
{
  SaveFile();
}

void AnalysisManager::SaveFile( void ) const
{
  if( fActive_ )
    hfile->Write();
}

void AnalysisManager::Terminate ( void ) const
{
  G4cout << "[Analysis] Terminate()"<< G4endl;
  if ( fActive_ )
    {
      hfile->Write();
      hfile->Close();
    }
}

void AnalysisManager::BeginOfRun(const G4Run*)
{
  G4SDManager* SDManager = G4SDManager::GetSDMpointer();

  hfile = new TFile(outfile, "RECREATE");
  tree = new TTree("tree","EvtGen tree");

  G4cout<<"Output file made: "<< outfile << G4endl;

  tree->Branch("event", &event, "event/I");
  tree->Branch("nEvt", &nEvt, "nEvt/I");
  tree->Branch("evtid", evtid, "evtid[nEvt]/I");
  tree->Branch("evtpid", evtpid, "evtpid[nEvt]/I");
  tree->Branch("evttrid", evttrid, "evttrid[nEvt]/I");
  tree->Branch("evtfm", evtfm, "evtfm[nEvt]/I");
  tree->Branch("evtlm", evtlm, "evtlm[nEvt]/I");
  tree->Branch("evtfd", evtfd, "evtfd[nEvt]/I");
  tree->Branch("evtld", evtld, "evtld[nEvt]/I");
  tree->Branch("evttime", evttime, "evttime[nEvt]/D");
  tree->Branch("evtpx", evtpx, "evtpx[nEvt]/D");
  tree->Branch("evtpy", evtpy, "evtpy[nEvt]/D");
  tree->Branch("evtpz", evtpz, "evtpz[nEvt]/D");
  tree->Branch("evtvx", evtvx, "evtvx[nEvt]/D");
  tree->Branch("evtvy", evtvy, "evtvy[nEvt]/D");
  tree->Branch("evtvz", evtvz, "evtvz[nEvt]/D");
  tree->Branch("evtm", evtm, "evtm[nEvt]/D");

  tree->Branch("nBeam", &nBeam, "nBeam/I");
  tree->Branch("bpx", bpx, "bpx[nBeam]/D");
  tree->Branch("bpy", bpy, "bpy[nBeam]/D");
  tree->Branch("bpz", bpz, "bpz[nBeam]/D");
  tree->Branch("bvx", bvx, "bvx[nBeam]/D");
  tree->Branch("bvy", bvy, "bvy[nBeam]/D");
  tree->Branch("bvz", bvz, "bvz[nBeam]/D");

  tree->Branch("nhTpcPad", &nhTpcPad, "nhTpcPad/I");
  tree->Branch("tpcpad", tpcpad, "tpcpad[nhTpcPad]/I");
  tree->Branch("tpcpadedep", tpcpadedep, "tpcpadedep[nhTpcPad]/D");
  tree->Branch("tpcpadtime", tpcpadtime, "tpcpadtime[nhTpcPad]/D");
  tree->Branch("tpcpadposx", tpcpadposx, "tpcpadposx[nhTpcPad]/D");
  tree->Branch("tpcpadposy", tpcpadposy, "tpcpadposy[nhTpcPad]/D");
  tree->Branch("tpcpadposz", tpcpadposz, "tpcpadposz[nhTpcPad]/D");
  tree->Branch("tpcpadmomx", tpcpadmomx, "tpcpadmomx[nhTpcPad]/D");
  tree->Branch("tpcpadmomy", tpcpadmomy, "tpcpadmomy[nhTpcPad]/D");
  tree->Branch("tpcpadmomz", tpcpadmomz, "tpcpadmomz[nhTpcPad]/D");
  tree->Branch("tpcpadtrid", tpcpadtrid, "tpcpadtrid[nhTpcPad]/I");
  tree->Branch("tpcpadpath", tpcpadpath, "tpcpadpath[nhTpcPad]/D");
  tree->Branch("tpcpadpid", tpcpadpid, "tpcpadpid[nhTpcPad]/I");

  tree->Branch("nhTof", &nhTof, "nhTof/I");
  tree->Branch("tofseg", tofseg, "tofseg[nhTof]/I");
  tree->Branch("tofedep", tofedep, "tofedep[nhTof]/D");
  tree->Branch("toftime", toftime, "toftime[nhTof]/D");
  tree->Branch("tofposx", tofposx, "tofposx[nhTof]/D");
  tree->Branch("tofposy", tofposy, "tofposy[nhTof]/D");
  tree->Branch("tofposz", tofposz, "tofposz[nhTof]/D");
  tree->Branch("tofmomx", tofmomx, "tofmomx[nhTof]/D");
  tree->Branch("tofmomy", tofmomy, "tofmomy[nhTof]/D");
  tree->Branch("tofmomz", tofmomz, "tofmomz[nhTof]/D");
  tree->Branch("toftrid", toftrid, "toftrid[nhTof]/I");
  tree->Branch("tofpath", tofpath, "tofpath[nhTof]/D");
  tree->Branch("tofpid", tofpid, "tofpid[nhTof]/I");
  tree->Branch("tofparentid1", tofparentid1, "tofparentid1[nhTof]/I");
  tree->Branch("tofparentid2", tofparentid2, "tofparentid2[nhTof]/I");
  tree->Branch("tofparentid3", tofparentid3, "tofparentid3[nhTof]/I");
  tree->Branch("tofparentpid1", tofparentpid1, "tofparentpid1[nhTof]/I");
  tree->Branch("tofparentpid2", tofparentpid2, "tofparentpid2[nhTof]/I");
  tree->Branch("tofparentpid3", tofparentpid3, "tofparentpid3[nhTof]/I");

  tree->Branch("nhTarget", &nhTarget, "nhTarget/I");
  tree->Branch("targetseg", targetseg, "targetseg[nhTarget]/I");
  tree->Branch("targetedep", targetedep, "targetedep[nhTarget]/D");
  tree->Branch("targettime", targettime, "targettime[nhTarget]/D");
  tree->Branch("targetposx", targetposx, "targetposx[nhTarget]/D");
  tree->Branch("targetposy", targetposy, "targetposy[nhTarget]/D");
  tree->Branch("targetposz", targetposz, "targetposz[nhTarget]/D");
  tree->Branch("targetmomx", targetmomx, "targetmomx[nhTarget]/D");
  tree->Branch("targetmomy", targetmomy, "targetmomy[nhTarget]/D");
  tree->Branch("targetmomz", targetmomz, "targetmomz[nhTarget]/D");
  tree->Branch("targettrid", targettrid, "targettrid[nhTarget]/I");
  tree->Branch("targetpath", targetpath, "targetpath[nhTarget]/D");
  tree->Branch("targethitdiff", targethitdiff, "targethitdiff[nhTarget]/D");
  tree->Branch("targetpid", targetpid, "targetpid[nhTarget]/I");

  // elastic test //
  tree->Branch("elcosth", &elcosth, "elcosth/D");


  //KURAMA Part//

  // DC //
  tree->Branch("nhDc", &nhDc, "nhDc/I");
  tree->Branch("dcid", dcid, "dcid[nhDc]/I");
  tree->Branch("dclayer", dclayer, "dclayer[nhDc]/I");
  tree->Branch("dcedep", dcedep, "dcedep[nhDc]/D");
  tree->Branch("dctime", dctime, "dctime[nhDc]/D");
  tree->Branch("dcposx", dcposx, "dcposx[nhDc]/D");
  tree->Branch("dcposy", dcposy, "dcposy[nhDc]/D");
  tree->Branch("dcposz", dcposz, "dcposz[nhDc]/D");
  tree->Branch("dcmomx", dcmomx, "dcmomx[nhDc]/D");
  tree->Branch("dcmomy", dcmomy, "dcmomy[nhDc]/D");
  tree->Branch("dcmomz", dcmomz, "dcmomz[nhDc]/D");
  tree->Branch("dctrid", dctrid, "dctrid[nhDc]/I");
  tree->Branch("dcpath", dcpath, "dcpath[nhDc]/D");
  tree->Branch("dcpid", dcpid, "dcpid[nhDc]/I");

  //FAC//                                                                                  
  tree->Branch("nhFac", &nhFac, "nhFac/I");
  tree->Branch("facseg", facseg, "facseg[nhFac]/I");
  tree->Branch("facedep", facedep, "facedep[nhFac]/D");
  tree->Branch("factime", factime, "factime[nhFac]/D");
  tree->Branch("facposx", facposx, "facposx[nhFac]/D");
  tree->Branch("facposy", facposy, "facposy[nhFac]/D");
  tree->Branch("facposz", facposz, "facposz[nhFac]/D");
  tree->Branch("facmomx", facmomx, "facmomx[nhFac]/D");
  tree->Branch("facmomy", facmomy, "facmomy[nhFac]/D");
  tree->Branch("facmomz", facmomz, "facmomz[nhFac]/D");
  tree->Branch("factrid", factrid, "factrid[nhFac]/I");
  tree->Branch("facpath", facpath, "facpath[nhFac]/D");
  tree->Branch("facpid", facpid, "facpid[nhFac]/I");

  //BDC//                                                                      
  tree->Branch("nhBdc", &nhBdc, "nhBdc/I");
  tree->Branch("bdcseg", bdcseg, "bdcseg[nhBdc]/I");
  tree->Branch("bdcedep", bdcedep, "bdcedep[nhBdc]/D");
  tree->Branch("bdctime", bdctime, "bdctime[nhBdc]/D");
  tree->Branch("bdcposx", bdcposx, "bdcposx[nhBdc]/D");
  tree->Branch("bdcposy", bdcposy, "bdcposy[nhBdc]/D");
  tree->Branch("bdcposz", bdcposz, "bdcposz[nhBdc]/D");
  tree->Branch("bdcmomx", bdcmomx, "bdcmomx[nhBdc]/D");
  tree->Branch("bdcmomy", bdcmomy, "bdcmomy[nhBdc]/D");
  tree->Branch("bdcmomz", bdcmomz, "bdcmomz[nhBdc]/D");
  tree->Branch("bdctrid", bdctrid, "bdctrid[nhBdc]/I");
  tree->Branch("bdcpath", bdcpath, "bdcpath[nhBdc]/D");
  tree->Branch("bdcpid", bdcpid, "bdcpid[nhBdc]/I");



  event = 0;
  nEvt = 0;

}


void AnalysisManager::EndOfRun(const G4Run*)
{
  tree->Write();
  hfile->Write();
  hfile->Close();
}

void AnalysisManager::BeginOfEvent(const G4Event*)
{
  //nEvt = 0;
  //G4cout<<"[AnalysisManager] Begin of Event: "<<nEvt<<G4endl;
}

void AnalysisManager::EndOfEvent(const G4Event* anEvent)
{

  G4HCofThisEvent* HCTE = anEvent-> GetHCofThisEvent();
  if(!HCTE) return;
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  G4int nhtpcpad=0, nhtof=0, nhtarget=0, nhfac=0, nhbdc=0, nhdc=0, nhch1 = 0, nhftof=0 ;
  TpcHitsCollection *TpcHC=0;
  TofHitsCollection *TofHC=0;
  TargetHitsCollection *TargetHC=0;
  DCHitsCollection *DCHC = 0;
  FACHitsCollection *FACHC = 0;
  BDCHitsCollection *BDCHC = 0;


  G4int colIdTpc = SDMan->GetCollectionID( "TpcCollection" );
  if(colIdTpc>=0)
    {
      TpcHC=dynamic_cast<TpcHitsCollection *>( HCTE->GetHC( colIdTpc ) );
      if(TpcHC)
	{
	  nhtpcpad=TpcHC->entries();
	}
    }

  G4int colIdTof = SDMan->GetCollectionID( "TofCollection" );
  if(colIdTof>=0)
    {
      TofHC=dynamic_cast<TofHitsCollection *>( HCTE->GetHC( colIdTof ) );
      if(TofHC)
	{
	  nhtof=TofHC->entries();
	}
    }

  G4int colIdTarget = SDMan->GetCollectionID( "TargetCollection" );
  if(colIdTarget>=0)
    {
      TargetHC=dynamic_cast<TargetHitsCollection *>( HCTE->GetHC( colIdTarget ) );
      if(TargetHC)
	{
	  nhtarget=TargetHC->entries();
	}
    }

  G4int colIdDC = SDMan->GetCollectionID( "DcCollection" );
  if(colIdDC>=0)
    {
      DCHC=dynamic_cast<DCHitsCollection *>( HCTE->GetHC( colIdDC ) );
      if(DCHC)
	{
	  nhdc=DCHC->entries();
	}
    }

  G4int colIdFAC = SDMan->GetCollectionID( "FacCollection" );                            
  if(colIdFAC>=0)
    {
      FACHC=dynamic_cast<FACHitsCollection *>( HCTE->GetHC( colIdFAC ) );
      if(FACHC)
        {
          nhfac=FACHC->entries();
        }
    }

  G4int colIdBDC = SDMan->GetCollectionID( "BdcCollection" );                            
  if(colIdBDC>=0)
    {
      BDCHC=dynamic_cast<BDCHitsCollection *>( HCTE->GetHC( colIdBDC ) );
      if(BDCHC)
        {
          nhbdc=BDCHC->entries();
        }
    }



  if(nhtpcpad > 200 -1)
    {
      G4cout<<"[AnalysisManager] Number of TPC Hit > 200"<<G4endl;
      return;
    }
  if(nhtof > 100 -1)
    {
      G4cout<<"[AnalysisManager] Number of Tof Hit > 100"<<G4endl;
      return;
    }
  if(nhtarget > 100 -1)
    {
      G4cout<<"[AnalysisManager] Number of Target Hit > 100"<<G4endl;
      return;
    }

  if(nhdc > 200 -1)
    {
      G4cout<<"[AnalysisManager] Number of DC Hit > 200"<<G4endl;
      return;
    }

  if(nhfac > 100 -1)
    {
      G4cout<<"[AnalysisManager] Number of FAC Hit > 100"<<G4endl;
      return;
    }

  if(nhbdc > 100 -1)
    {
      G4cout<<"[AnalysisManager] Number of BDC Hit > 100"<<G4endl;
      return;
    }


  for( int i=0; i<nhtpcpad; ++i )
    {
      TpcHit *aHit=(*TpcHC)[i];
      tpcpad[i] = aHit->GetPadID();
      tpcpadedep[i] = aHit->GetEdep();
      tpcpadtime[i] = aHit->GetTime();
      tpcpadposx[i] = aHit->GetPos().x();
      tpcpadposy[i] = aHit->GetPos().y();
      tpcpadposz[i] = aHit->GetPos().z();
      tpcpadmomx[i] = aHit->GetMom().x()*0.001;
      tpcpadmomy[i] = aHit->GetMom().y()*0.001;
      tpcpadmomz[i] = aHit->GetMom().z()*0.001;
      tpcpadtrid[i] = aHit->GetTrackNo();
      tpcpadpath[i] = aHit->GetPathLength();
      tpcpadpid[i] = aHit->GetHitParticleID();
    }
  nhTpcPad = nhtpcpad;

  for( int i=0; i<nhtof; ++i )
    {
      TofHit *aHit=(*TofHC)[i];
      tofseg[i] = aHit->GetSegmentID();
      tofedep[i] = aHit->GetEdep();
      toftime[i] = aHit->GetTime();
      tofposx[i] = aHit->GetPos().x();
      tofposy[i] = aHit->GetPos().y();
      tofposz[i] = aHit->GetPos().z();
      tofmomx[i] = aHit->GetMom().x()*0.001;
      tofmomy[i] = aHit->GetMom().y()*0.001;
      tofmomz[i] = aHit->GetMom().z()*0.001;
      toftrid[i] = aHit->GetTrackNo();
      tofpath[i] = aHit->GetPathLength();
      tofpid[i] = aHit->GetHitParticleID();
      tofparentid1[i] = aHit->GetParentID1();
      tofparentid2[i] = aHit->GetParentID2();
      tofparentid3[i] = aHit->GetParentID3();
      tofparentpid1[i] = aHit->GetParentPID1();
      tofparentpid2[i] = aHit->GetParentPID2();
      tofparentpid3[i] = aHit->GetParentPID3();

    }
  nhTof = nhtof;

  for( int i=0; i<nhtarget; ++i )
    {
      TargetHit *aHit=(*TargetHC)[i];
      targetseg[i] = aHit->GetLayerID();
      targetedep[i] = aHit->GetEdep();
      targettime[i] = aHit->GetTime();
      targetposx[i] = aHit->GetPos().x();
      targetposy[i] = aHit->GetPos().y();
      targetposz[i] = aHit->GetPos().z();
      targetmomx[i] = aHit->GetMom().x()*0.001;
      targetmomy[i] = aHit->GetMom().y()*0.001;
      targetmomz[i] = aHit->GetMom().z()*0.001;
      targettrid[i] = aHit->GetTrackNo();
      targetpath[i] = aHit->GetPathLength();
      targethitdiff[i] = aHit->GetHitDiff();
      targetpid[i] = aHit->GetHitParticleID();
    }
  nhTarget = nhtarget;

  for( int i=0; i<nhfac; ++i )
    {
      FACHit *aHit=(*FACHC)[i];
      facseg[i] = aHit->GetDetectorID();
      facedep[i] = aHit->GetEdep();
      factime[i] = aHit->GetTOF();
      facposx[i] = aHit->GetPosition().x();
      facposy[i] = aHit->GetPosition().y();
      facposz[i] = aHit->GetPosition().z();
      facmomx[i] = aHit->GetMomentum().x()*0.001;
      facmomy[i] = aHit->GetMomentum().y()*0.001;
      facmomz[i] = aHit->GetMomentum().z()*0.001;
      factrid[i] = aHit->GetTrackID();
      facpath[i] = aHit->GetLength();
      facpid[i] = aHit->GetParticleID();
    }
  nhFac = nhfac;

  for( int i=0; i<nhbdc; ++i )
    {
      BDCHit *aHit=(*BDCHC)[i];
      bdcseg[i] = aHit->GetDetectorID();
      bdcedep[i] = aHit->GetEdep();
      bdctime[i] = aHit->GetTOF();
      bdcposx[i] = aHit->GetPosition().x();
      bdcposy[i] = aHit->GetPosition().y();
      bdcposz[i] = aHit->GetPosition().z();
      bdcmomx[i] = aHit->GetMomentum().x()*0.001;
      bdcmomy[i] = aHit->GetMomentum().y()*0.001;
      bdcmomz[i] = aHit->GetMomentum().z()*0.001;
      bdctrid[i] = aHit->GetTrackID();
      bdcpath[i] = aHit->GetLength();
      bdcpid[i] = aHit->GetParticleID();
    }
  nhBdc = nhbdc;


  // KURAMA part // 
  for( int i=0; i<nhdc; ++i )
    {
      DCHit *aHit=(*DCHC)[i];

      dclayer[i] = aHit->GetDetectorID();
      if(dclayer[i] < 9)
	{
	  dcid[i] = 0;
	}
      else if(dclayer[i] > 9 && dclayer[i] < 19)
	{
	  dcid[i] = 1;
	}
      else if(dclayer[i] > 19 && dclayer[i] < 29)
	{
	  dcid[i] = 2;
	}

      else if(dclayer[i] == 100)
	{
	  dcid[i] = 100;
	}
      else if(dclayer[i] == 101)
	{
	  dcid[i] = 101;
	}
      else if(dclayer[i] == 102)
	{
	  dcid[i] = 102;
	}
      else if(dclayer[i] == 103)
	{
	  dcid[i] = 103;
	}
      else if(dclayer[i] == 104)
	{
	  dcid[i] = 104;
	}
      else if(dclayer[i] == 105)
	{
	  dcid[i] = 105;
	}
      else 
	{
	  G4cout<< "[AnalysisManager] DC layer is wrong"<<G4endl;
	}
      
      dcedep[i] = aHit->GetEdep();
      dctime[i] = aHit->GetTOF();
      dcposx[i] = aHit->GetPosition().x();
      dcposy[i] = aHit->GetPosition().y();
      dcposz[i] = aHit->GetPosition().z();
      dcmomx[i] = aHit->GetMomentum().x()*0.001;
      dcmomy[i] = aHit->GetMomentum().y()*0.001;
      dcmomz[i] = aHit->GetMomentum().z()*0.001;
      dctrid[i] = aHit->GetTrackID();
      dcpath[i] = aHit->GetLength();
      dcpid[i] = aHit->GetParticleID();
    }
  nhDc = nhdc;

  tree->Fill();
  event++;

  nEvt = 0;
  nhTpcPad = 0;
  nhTof = 0;
  nhTarget = 0;
  nBeam = 0;
  nhFac = 0;
  nhBdc = 0;

  // KURAMA part //
  nhDc = 0;
}

void AnalysisManager::SetEvtGen(int j, int partnum, 
				int jmotherfirst, 
				int jmotherlast, 
				int jdaugfirst, 
				int jdauglast, 
				int tr_id,
				EvtVector4R x,
				EvtVector4R p
				)
{
  //G4cout<<"begin nEvt: "<<nEvt<<" j: "<<j<<G4endl;
  if(nEvt < j) nEvt =j;
  //G4cout<<"nEvt: "<<nEvt<<G4endl;
  evtid[j-1]=j;
  evtpid[j-1] = partnum;
  evttrid[j-1] = tr_id;
  evtfm[j-1] = jmotherfirst;
  evtlm[j-1] = jmotherlast;
  evtfd[j-1] = jdaugfirst;
  evtld[j-1] = jdauglast;
  evttime[j-1] = x.get(0);
  evtpx[j-1] = p.get(1);
  evtpy[j-1] = p.get(2);
  evtpz[j-1] = p.get(3);
  evtvx[j-1] = x.get(1);
  evtvy[j-1] = x.get(2);
  evtvz[j-1] = x.get(3);
  evtm[j-1] = p.mass();
  //G4cout<<"[Analysis Manager] j: "<<j << " ParticleNum: " << partnum << " nEvt: "<<nEvt<<G4endl;
}

void AnalysisManager::SetBeam(int j,
			      G4ThreeVector D,
			      G4ThreeVector P
			      )
{
  nBeam = j;
  bpx[j-1] = P.x();
  bpy[j-1] = P.y();
  bpz[j-1] = P.z();
  bvx[j-1] = D.x();
  bvy[j-1] = D.y();
  bvz[j-1] = D.z();

}

void AnalysisManager::SetCosTheta( double costheta )
{
  elcosth = costheta;
}
