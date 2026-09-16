// -*- C++ -*-

#include "EventAction.hh"

#include <fstream>
#include <map>
#include <vector>

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
#include "TPCPadHelper.hh"
#include "BACSD.hh"
#include "ConfMan.hh"
#include "DetectorConstruction.hh"
#include "FuncName.hh"
#include "BH2SD.hh"
#include "HTOFSD.hh"
#include "BACSD.hh"
#include "TPCSD.hh"
#include "TPCPadSD.hh"
#include "TPCEdepSD.hh"
#include "SCHSD.hh"
#include "SDCSD.hh"
#include "TargetSD.hh"
#include "VPSD.hh"
#include "KVCSD.hh"
#include "CVCSD.hh"
#include "SAC3SD.hh"
#include "SFVSD.hh"
#include "DetSizeMan.hh"

#define SDCLASS(x) x ## Hit

namespace
{
  auto& gAnaMan = AnaManager::GetInstance();
  const auto& gConf = ConfMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
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

  // Legacy KURAMA FTOF (FTOFSD); not E72 IncludeFTOF (CVC/SAC3/SFV).
  // ConstructFTOF() is unused, so this SD is never registered.
#if 0
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
#endif

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

  if (gConf.Get<G4bool>("IncludeFTOF")) {
    {
      static const auto id = SDManager->GetCollectionID("CVC/hit");
      if (id >= 0) {
        auto HC = dynamic_cast<G4THitsCollection<CVCHit>*>(HCTE->GetHC(id));
        for (G4int i=0, n=HC->entries(); i<n; ++i) {
	  gAnaMan.SetHitData((*HC)[i]);
        }
        gAnaMan.SetNhits("CVC", HC->entries());
      }
    }

    {
      static const auto id = SDManager->GetCollectionID("SAC3/hit");
      if (id >= 0) {
        auto HC = dynamic_cast<G4THitsCollection<SAC3Hit>*>(HCTE->GetHC(id));
        for (G4int i=0, n=HC->entries(); i<n; ++i) {
	  gAnaMan.SetHitData((*HC)[i]);
        }
        gAnaMan.SetNhits("SAC3", HC->entries());
      }
    }

    {
      static const auto id = SDManager->GetCollectionID("SFV/hit");
      if (id >= 0) {
        auto HC = dynamic_cast<G4THitsCollection<SFVHit>*>(HCTE->GetHC(id));
        for (G4int i=0, n=HC->entries(); i<n; ++i) {
	  gAnaMan.SetHitData((*HC)[i]);
        }
        gAnaMan.SetNhits("SFV", HC->entries());
      }
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
    if(id >= 0){
      auto HC = dynamic_cast<G4THitsCollection<TPCHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
	gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("TPC", HC->entries());
    }
  }

  {
    static const G4int id = SDManager->GetCollectionID("TPCPad/hit");
    static const G4int id_edep = SDManager->GetCollectionID("TPCEdep/hit");
    if( id >= 0 && id_edep >= 0 ){
      //test
      G4int tid_check;
      auto HC = (G4THitsCollection<TPCPadHit>*)(HCTE->GetHC(id));
      auto HC_edep = (G4THitsCollection<TPCEdepHit>*)(HCTE->GetHC(id_edep));

      G4int nhits= HC -> entries();
      G4int nhits_edep = HC_edep -> entries();

      G4int pidtr[MaxHitsTPC]={0};
      G4int ptidtpc[MaxHitsTPC]={0};
      G4int ptidtpc_pid[MaxHitsTPC]={0};
      G4double pmtpc[MaxHitsTPC]={0};
      G4int qqtpc[MaxHitsTPC]={0};
      G4double pxtpc[MaxHitsTPC]={0};
      G4double pytpc[MaxHitsTPC]={0};
      G4double pztpc[MaxHitsTPC]={0};
      G4double pptpc[MaxHitsTPC]={0};
      G4double vtxxtpc[MaxHitsTPC]={0};
      G4double vtxytpc[MaxHitsTPC]={0};
      G4double vtxztpc[MaxHitsTPC]={0};
      G4double vtxpxtpc[MaxHitsTPC]={0};
      G4double vtxpytpc[MaxHitsTPC]={0};
      G4double vtxpztpc[MaxHitsTPC]={0};
      // G4double vtxpptpc[MaxHitsTPC]={0};
      G4double vtxenetpc[MaxHitsTPC]={0};
      G4double detpc[MaxHitsTPC]={0};
      G4int laytpc[MaxHitsTPC]={0};
      G4double lentpc[MaxHitsTPC]={0};
      G4int nparticle=0;
      std::map<G4int, G4int> tid_to_slot; // tid -> summary slot
      // G4cout << "TPC  " << nhits << G4endl;

      const G4int pad_configure = gSize.Get("TpcPadConfigure");

      std::vector<std::vector<G4ThreeVector>> remain_pos;
      std::vector<std::vector<double>> remain_edep;
      std::vector<std::vector<int>> remain_ilay;
      std::vector<int> remain_pid;
      std::vector<int> remain_tid;
      std::vector<int> remain_parentid;
      std::vector<int> remain_parentpid;

      // Pad lookup built after child→parent merge (pad_configure==4 only).
      std::map<G4int, std::size_t> pad_tid_to_remain;
      std::vector<std::map<G4int, std::vector<std::size_t>>> pad_ilay_to_ks;

      if (pad_configure == 4) {
	// Group Edep hits by tid (non-contiguous tid stays one remain row).
	std::map<G4int, std::size_t> tid_to_remain_idx;
	std::vector<std::vector<double>> open_posx;
	std::vector<std::vector<double>> open_posy;
	std::vector<std::vector<double>> open_posz;
	std::vector<std::vector<double>> open_edep;
	std::vector<G4int> open_ilay;

	auto flush_open_layer = [&](std::size_t r) {
	  if (open_posx[r].empty())
	    return;
	  G4ThreeVector ave_pos(CalculateAverage(open_posx[r]),
				CalculateAverage(open_posy[r]),
				CalculateAverage(open_posz[r]));
	  remain_pos[r].push_back(ave_pos);
	  remain_edep[r].push_back(CalculateSum(open_edep[r]));
	  remain_ilay[r].push_back(open_ilay[r]);
	  open_posx[r].clear();
	  open_posy[r].clear();
	  open_posz[r].clear();
	  open_edep[r].clear();
	};

	for (G4int i = 0; i < nhits_edep; ++i) {
	  G4int pid_edep = (*HC_edep)[i]->GetParticleID();
	  G4int parentid_edep = (*HC_edep)[i]->GetParentID();
	  G4int parentpid_edep = (*HC_edep)[i]->GetParentID_pid();
	  G4double edep_edep = (*HC_edep)[i]->GetEdep();
	  G4int tid_edep = (*HC_edep)[i]->GetTrackID();
	  G4ThreeVector xyz_edep = (*HC_edep)[i]->GetPosition();
	  G4int ilay_edep = (*HC_edep)[i]->GetPadLay();

	  auto it = tid_to_remain_idx.find(tid_edep);
	  if (it == tid_to_remain_idx.end()) {
	    std::size_t r = remain_tid.size();
	    tid_to_remain_idx[tid_edep] = r;
	    remain_pos.emplace_back();
	    remain_edep.emplace_back();
	    remain_ilay.emplace_back();
	    remain_pid.push_back(pid_edep);
	    remain_tid.push_back(tid_edep);
	    remain_parentid.push_back(parentid_edep);
	    remain_parentpid.push_back(parentpid_edep);
	    open_posx.emplace_back();
	    open_posy.emplace_back();
	    open_posz.emplace_back();
	    open_edep.emplace_back();
	    open_ilay.push_back(ilay_edep);
	    open_posx[r].push_back(xyz_edep.x());
	    open_posy[r].push_back(xyz_edep.y());
	    open_posz[r].push_back(xyz_edep.z());
	    open_edep[r].push_back(edep_edep);
	  } else {
	    std::size_t r = it->second;
	    if (!open_posx[r].empty() && open_ilay[r] == ilay_edep) {
	      open_posx[r].push_back(xyz_edep.x());
	      open_posy[r].push_back(xyz_edep.y());
	      open_posz[r].push_back(xyz_edep.z());
	      open_edep[r].push_back(edep_edep);
	    } else {
	      flush_open_layer(r);
	      open_ilay[r] = ilay_edep;
	      open_posx[r].push_back(xyz_edep.x());
	      open_posy[r].push_back(xyz_edep.y());
	      open_posz[r].push_back(xyz_edep.z());
	      open_edep[r].push_back(edep_edep);
	    }
	  }
	}
	for (std::size_t r = 0; r < remain_tid.size(); ++r)
	  flush_open_layer(r);

	// Child → parent edep add (same ilay, dist < 30 mm). Local indexes.
	{
	  std::map<G4int, std::vector<std::size_t>> parent_to_children;
	  for (std::size_t u = 0; u < remain_tid.size(); ++u)
	    parent_to_children[remain_parentid[u]].push_back(u);

	  std::vector<std::map<G4int, std::vector<std::size_t>>> ilay_to_ks(
	    remain_tid.size());
	  for (std::size_t n = 0; n < remain_tid.size(); ++n) {
	    for (std::size_t k = 0; k < remain_ilay[n].size(); ++k)
	      ilay_to_ks[n][remain_ilay[n][k]].push_back(k);
	  }

	  for (std::size_t n = 0; n < remain_tid.size(); ++n) {
	    auto ch_it = parent_to_children.find(remain_tid[n]);
	    if (ch_it == parent_to_children.end())
	      continue;
	    for (std::size_t u : ch_it->second) {
	      for (std::size_t j = 0; j < remain_edep[u].size(); ++j) {
		auto ks_it = ilay_to_ks[n].find(remain_ilay[u][j]);
		if (ks_it == ilay_to_ks[n].end())
		  continue;
		for (std::size_t k : ks_it->second) {
		  if ((remain_pos[n][k] - remain_pos[u][j]).mag() < 30.)
		    remain_edep[n][k] += remain_edep[u][j];
		}
	      }
	    }
	  }
	}

	// Pad → Edep lookup (after merge). Separate local indexes.
	pad_ilay_to_ks.resize(remain_tid.size());
	for (std::size_t j = 0; j < remain_tid.size(); ++j) {
	  pad_tid_to_remain.insert({remain_tid[j], j}); // first j wins
	  for (std::size_t k = 0; k < remain_ilay[j].size(); ++k)
	    pad_ilay_to_ks[j][remain_ilay[j][k]].push_back(k);
	}
      }

      // Primary TrackID → PDG (once per event; parentpid only for primary parents).
      std::map<G4int, G4int> primary_tid_to_pdg;
      for (G4int k = 0; k < anEvent->GetNumberOfPrimaryVertex(); ++k) {
	G4PrimaryVertex* primaryVertex = anEvent->GetPrimaryVertex(k);
	for (G4int j = 0; j < primaryVertex->GetNumberOfParticle(); ++j) {
	  G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary(j);
	  primary_tid_to_pdg[primaryParticle->GetTrackID()] =
	    primaryParticle->GetPDGcode();
	}
      }

      for (G4int i = 0; i < nhits; ++i) {
	G4ThreeVector vtxpos = (*HC)[i]->GetVtxPosition();
	G4ThreeVector vtxmom = (*HC)[i]->GetVtxMomentum();
	G4double vtxene = (*HC)[i]->GetVtxEnergy();
	G4ThreeVector xyz = (*HC)[i]->GetPosition();
	G4ThreeVector mom = (*HC)[i]->GetMomentum();
	G4double tof = (*HC)[i]->GetTOF();
	G4int tid = (*HC)[i]->GetTrackID();
	G4int ptid = (*HC)[i]->GetParentID();
	G4int ptid_pid = (*HC)[i]->GetParentID_pid();
	G4int pid = (*HC)[i]->GetParticleID();
	G4double mass = (*HC)[i]->GetMass();
	G4int charge = (*HC)[i]->GetCharge();
	G4int ilay = (*HC)[i]->GetPadLay();
	G4int parentid = (*HC)[i]->GetParentID();

	G4int parentpid = -9999;
	if (parentid > 0) {
	  auto pp_it = primary_tid_to_pdg.find(parentid);
	  if (pp_it != primary_tid_to_pdg.end())
	    parentpid = pp_it->second;
	}

	G4double tlength = (*HC)[i]->GettLength();
	G4int irow = (*HC)[i]->GetPadRow();
	G4double beta = (*HC)[i]->GetBeta();

	G4double edep;

	if (pad_configure == 3) {
	  edep = (*HC)[i]->GetEdep();
#if 0
	  // Legacy momentum-diff edep (unreachable; OOB on HC[i+1] if re-enabled).
	  // const G4int edep_configure = gSize.Get("TpcEdep");
	  int edep_configure=0;
	  if(edep_configure==0)
	    edep = (*HC)[i]-> GetEdep();
	  else if(edep_configure==1){
	    if(i==0){
	      G4ThreeVector mom_pre = (*HC)[i+1]-> GetMomentum();
	      G4double beta_pre = (*HC)[i+1]-> GetBeta();
	      G4double E_pre = mom_pre.mag() / beta_pre;
	      G4double E_post = mom.mag() / beta;
	      Double_t MeVToeV = TMath::Power(10.,6);
	      edep = (E_post - E_pre) * MeVToeV;
	    }
	    else{
	      G4ThreeVector xyz_pre = (*HC)[i-1]-> GetPosition();
	      if(abs((xyz - xyz_pre).mag())>30){
		G4ThreeVector mom_pre = (*HC)[i+1]-> GetMomentum();
		G4double beta_pre = (*HC)[i+1]-> GetBeta();
		G4double E_pre = mom_pre.mag() / beta_pre;
		G4double E_post = mom.mag() / beta;
		Double_t MeVToeV = TMath::Power(10.,6);
		edep = (E_post - E_pre) * MeVToeV;
	      }
	      else if(tid!=(*HC)[i-1]-> GetTrackID()){
		G4ThreeVector mom_pre = (*HC)[i+1]-> GetMomentum();
		G4double beta_pre = (*HC)[i+1]-> GetBeta();
		G4double E_pre = mom_pre.mag() / beta_pre;
		G4double E_post = mom.mag() / beta;
		Double_t MeVToeV = TMath::Power(10.,6);
		edep = (E_post - E_pre) * MeVToeV;
	      }
	      else{
		G4ThreeVector mom_pre = (*HC)[i-1]-> GetMomentum();
		G4double beta_pre = (*HC)[i-1]-> GetBeta();
		G4double E_pre = mom_pre.mag() / beta_pre;
		G4double E_post = mom.mag() / beta;
		Double_t MeVToeV = TMath::Power(10.,6);
		edep = (E_pre - E_post) * MeVToeV;
		std::cout<<mom_pre.mag()<<" , "<<mom.mag()<<", diff : "<<mom_pre.mag()-mom.mag()<<std::endl;
	      }
	    }
	  }
#endif
	}

	G4bool find_track = false;
	G4bool find_hit = false;
	if (pad_configure == 4) {
	  // Match Pad hit to layer-summed Edep by (tid, ilay) and position.
	  // Pad is a thin-shell crossing; Edep position is the average over the
	  // thick layer — reject if they differ by >= 50 mm (large-angle mismatch).
	  auto tr_it = pad_tid_to_remain.find(tid);
	  if (tr_it != pad_tid_to_remain.end()) {
	    find_track = true;
	    std::size_t j = tr_it->second;
	    auto ks_it = pad_ilay_to_ks[j].find(ilay);
	    if (ks_it != pad_ilay_to_ks[j].end()) {
	      for (std::size_t k : ks_it->second) {
		G4double pos_diff = (remain_pos[j][k] - xyz).mag();
		if (pos_diff < 50.) {
		  edep = remain_edep[j][k];
		  find_hit = true;
		  break;
		}
	      }
	    }
	    if (!find_hit) {
	      G4cout << "[EventAction] no Edep for Pad hit"
		     << " tid=" << tid << " ilay=" << ilay << G4endl;
	    }
	  }
	  if (!find_track) {
	    G4cout << "[EventAction] no Edep track for Pad hit"
		   << " tid=" << tid << G4endl;
	  }

	  if (!find_hit)
	    edep = -9999;
	}






	G4double slength = (*HC)[i]-> GetsLength();
	// Group Pad hits into track summary slots by Geant4 track ID.
	G4int slot = -1;
	auto tid_it = tid_to_slot.find(tid);
	if(tid_it != tid_to_slot.end()){
	  slot = tid_it->second;
	  detpc[slot] = detpc[slot]+edep;
	  if(ilay>-1){
	    laytpc[slot]=laytpc[slot]+1;
	  }
	  lentpc[slot]=tlength;
	}else if(nparticle <= 20){
	  slot = nparticle;
	  qqtpc[slot]=charge;
	  pmtpc[slot]=mass;
	  detpc[slot]=detpc[slot]+edep;
	  if(ilay>-1){
	    laytpc[slot]=laytpc[slot]+1;
	  }
	  pidtr[slot]=pid;
	  pxtpc[slot]=mom[0];
	  pytpc[slot]=mom[1];
	  pztpc[slot]=mom[2];

	  //////////////////////vertex information /////////////////////////
	  vtxpxtpc[slot]=vtxmom[0];
	  vtxpytpc[slot]=vtxmom[1];
	  vtxpztpc[slot]=vtxmom[2];
	  // vtxpptpc[slot]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));

	  vtxxtpc[slot]=vtxpos[0];
	  vtxytpc[slot]=vtxpos[1];
	  vtxztpc[slot]=vtxpos[2];

	  vtxenetpc[slot]=vtxene;

	  pptpc[slot]=mom.mag();
	  lentpc[slot]=tlength;
	  ptidtpc[slot]=ptid;
	  ptidtpc_pid[slot]=ptid_pid;
	  tid_to_slot[tid]=slot;
	  nparticle=nparticle+1;
	}

	if(slot >= 0 && ilay>-1){ //-->  -1 : TPC, layer is from 0 to 38. 2012.10.30
	  // Dead = center-frame pads; Noise = abnormal waveform (E72).
	  // TPCDropNoisePad (default true): GetDeadCon == IsDead || (drop_noise && Noise).
	  static const G4bool drop_noise =
	    gConf.GetOrDefault<G4bool>("TPCDropNoisePad", true);
	  if (TPCPadHelper::GetDeadCon((*HC)[i]->GetPadLay(),
				       (*HC)[i]->GetPadRow(), drop_noise)) {
	    continue;
	  }
	  static const G4int restype = gConf.Get<G4int>("ResType");
	  switch(restype){
	  case 0:
	    gAnaMan.SetCounterDataSimple( slot,tof, xyz, mom, tid, pid, ilay,
					  irow, beta, edep/CLHEP::MeV, parentid, parentpid, charge, tlength,slength );
	    break;
	  case 1:
	    gAnaMan.SetCounterDataExp( slot,tof, xyz, mom, tid, pid, ilay,
				       irow, beta, edep/CLHEP::MeV, parentid, parentpid, charge, tlength,slength );
	    break;
	  default:
	    gAnaMan.SetCounterDataSimple( slot,tof, xyz, mom, tid, pid, ilay,
					  irow, beta, edep/CLHEP::MeV, parentid, parentpid, charge, tlength,slength );
	    G4cout<<"TPC Resolution type is not determined. Now using constant resolution"<<G4endl;
	    break;
	  }

	}
      }

      for(G4int i=0;i<nparticle;i++){
	gAnaMan.SetTPCData(i, pidtr[i], ptidtpc[i], ptidtpc_pid[i], pxtpc[i],pytpc[i],pztpc[i],pptpc[i], qqtpc[i], pmtpc[i], detpc[i], lentpc[i],laytpc[i],
			   vtxpxtpc[i],vtxpytpc[i],vtxpztpc[i],
			   vtxxtpc[i],vtxytpc[i],vtxztpc[i], vtxenetpc[i]);
      }


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

G4double
EventAction::CalculateAverage(const std::vector<double>& pos)
{
  if(pos.size()==0)return -1000;

  G4double pos_sum = 0;
  for(size_t i = 0; i < pos.size(); ++i){
    pos_sum+=pos[i];
  }

  return pos_sum / pos.size();
}

G4double
EventAction::CalculateSum(const std::vector<double>& edep)
{
  G4double edep_sum = 0;

  for(size_t i = 0; i < edep.size(); ++i){
    edep_sum+=edep[i];
  }

  return edep_sum;
}
