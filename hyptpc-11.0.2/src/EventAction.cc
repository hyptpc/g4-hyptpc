// -*- C++ -*-

#include "EventAction.hh"

#include <fstream>
#include <map>

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


      const G4int pad_configure =  gSize.Get("TpcPadConfigure");
      std::vector<G4ThreeVector> temp_pos;
      std::vector<double> temp_edep;
      std::vector<int> temp_ilay;

      G4int temp_pid;
      G4int temp_tid;
      G4int temp_parentid;
      G4int temp_parentpid;

      std::vector<std::vector<G4ThreeVector>> remain_pos;
      std::vector<std::vector<double>> remain_edep;
      std::vector<std::vector<int>> remain_ilay;
      std::vector<int> remain_pid;
      std::vector<int> remain_tid;
      std::vector<int> remain_parentid;
      std::vector<int> remain_parentpid;


      std::vector<double> cal_posx;
      std::vector<double> cal_posy;
      std::vector<double> cal_posz;
      std::vector<double> cal_edep;

      if(pad_configure==4){
	//Edep Hit sort
	for(int i=0; i<nhits_edep;i++){
	  G4int pid_edep = (*HC_edep)[i]-> GetParticleID();
	  G4int parentid_edep = (*HC_edep)[i]-> GetParentID();
	  G4int parentpid_edep = (*HC_edep)[i]-> GetParentID_pid();
	  G4double edep_edep = (*HC_edep)[i]-> GetEdep();
	  G4int tid_edep = (*HC_edep)[i]-> GetTrackID();
	  G4ThreeVector xyz_edep = (*HC_edep)[i]-> GetPosition();
	  G4int ilay_edep = (*HC_edep)[i]-> GetPadLay();

	  if(i==0){
	    cal_posx.push_back(xyz_edep.x());
	    cal_posy.push_back(xyz_edep.y());
	    cal_posz.push_back(xyz_edep.z());
	    cal_edep.push_back(edep_edep);

	    temp_pid = pid_edep;
	    temp_tid = tid_edep;
	    temp_parentid = parentid_edep;
	    temp_parentpid = parentpid_edep;

	    temp_ilay.push_back(ilay_edep);

	    if(nhits_edep==1){
	      G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
	      temp_pos.push_back(ave_pos);
	      temp_edep.push_back(CalculateSum(cal_edep));

	      //put track info into remain
	      remain_pos.push_back(temp_pos);
	      remain_edep.push_back(temp_edep);
	      remain_ilay.push_back(temp_ilay);
	      remain_pid.push_back(temp_pid);
	      remain_tid.push_back(temp_tid);
	      remain_parentid.push_back(temp_parentid);
	      remain_parentpid.push_back(temp_parentpid);
	    }

	  }

	  else if(i>0){
	    if(tid_edep == temp_tid){
	      if(ilay_edep == temp_ilay.back()){
		cal_posx.push_back(xyz_edep.x());
		cal_posy.push_back(xyz_edep.y());
		cal_posz.push_back(xyz_edep.z());
		cal_edep.push_back(edep_edep);
	      }

	      else if(ilay_edep != temp_ilay.back()){
		//put previous layer info into temp
		G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
		temp_pos.push_back(ave_pos);
		temp_edep.push_back(CalculateSum(cal_edep));

		//clear cal
		cal_posx.clear();
		cal_posy.clear();
		cal_posz.clear();
		cal_edep.clear();

		//put new hit info to cal
		cal_posx.push_back(xyz_edep.x());
		cal_posy.push_back(xyz_edep.y());
		cal_posz.push_back(xyz_edep.z());
		cal_edep.push_back(edep_edep);

		temp_ilay.push_back(ilay_edep);
	      }

	      if(i == nhits_edep-1){
		//End of the edep hits
		//Finalize previous track!!!
		//put current layer info into temp (last layer of previous track)
		G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
		temp_pos.push_back(ave_pos);
		temp_edep.push_back(CalculateSum(cal_edep));

		//put track info into remain
		remain_pos.push_back(temp_pos);
		remain_edep.push_back(temp_edep);
		remain_ilay.push_back(temp_ilay);
		remain_pid.push_back(temp_pid);
		remain_tid.push_back(temp_tid);
		remain_parentid.push_back(temp_parentid);
		remain_parentpid.push_back(temp_parentpid);
	      }

	    }
	    else if(tid_edep != temp_tid){
	      //Finalize previous track!!!
	      //put previous layer info into temp (last layer of previous track)
	      G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
	      temp_pos.push_back(ave_pos);
	      temp_edep.push_back(CalculateSum(cal_edep));

	      //put track info into remain
	      remain_pos.push_back(temp_pos);
	      remain_edep.push_back(temp_edep);
	      remain_ilay.push_back(temp_ilay);
	      remain_pid.push_back(temp_pid);
	      remain_tid.push_back(temp_tid);
	      remain_parentid.push_back(temp_parentid);
	      remain_parentpid.push_back(temp_parentpid);

	      //clear cal

	      cal_posx.clear();
	      cal_posy.clear();
	      cal_posz.clear();
	      cal_edep.clear();

	      //clear temp

	      temp_pos.clear();
	      temp_edep.clear();
	      temp_ilay.clear();

	      temp_pid = -9999;
	      temp_tid = -9999;
	      temp_parentid = -9999;
	      temp_parentpid = -9999;

	      //put new hit of new track to temp and cal
	      cal_posx.push_back(xyz_edep.x());
	      cal_posy.push_back(xyz_edep.y());
	      cal_posz.push_back(xyz_edep.z());
	      cal_edep.push_back(edep_edep);

	      temp_pid = pid_edep;
	      temp_tid = tid_edep;
	      temp_parentid = parentid_edep;
	      temp_parentpid = parentpid_edep;

	      temp_ilay.push_back(ilay_edep);



	      if(i == nhits_edep-1){
		//End of the edep hits
		//Finalize previous track!!!
		//put current layer info into temp (last layer of previous track)
		G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
		temp_pos.push_back(ave_pos);
		temp_edep.push_back(CalculateSum(cal_edep));

		//put track info into remain
		remain_pos.push_back(temp_pos);
		remain_edep.push_back(temp_edep);
		remain_ilay.push_back(temp_ilay);
		remain_pid.push_back(temp_pid);
		remain_tid.push_back(temp_tid);
		remain_parentid.push_back(temp_parentid);
		remain_parentpid.push_back(temp_parentpid);
	      }
	    }
	  }
	}


	for(int n=0;n<remain_tid.size();n++){
	  for(int u=0;u<remain_tid.size();u++){
	    if(remain_tid[n]==remain_parentid[u]){
	      for(int j=0;j<remain_edep[u].size();j++){
		for(int k=0;k<remain_edep[n].size();k++){
		  if(remain_ilay[n][k]==remain_ilay[u][j] && (remain_pos[n][k]-remain_pos[u][j]).mag()<30){
		    remain_edep[n][k]+=remain_edep[u][j];
		  }

		}
	      }
	    }
	  }
	}




      }


      for( G4int i=0; i<nhits; ++i ){
	G4ThreeVector vtxpos = (*HC)[i]-> GetVtxPosition();
	G4ThreeVector vtxmom = (*HC)[i]-> GetVtxMomentum();
	G4double vtxene =(*HC)[i]-> GetVtxEnergy();
	G4ThreeVector xyz = (*HC)[i]-> GetPosition();
	G4ThreeVector mom = (*HC)[i]-> GetMomentum();
	G4double tof= (*HC)[i]-> GetTOF();
	G4int tid = (*HC)[i]-> GetTrackID();
	G4int ptid = (*HC)[i]-> GetParentID();
	G4int ptid_pid = (*HC)[i]-> GetParentID_pid();
	G4int pid = (*HC)[i]-> GetParticleID();
	G4double mass = (*HC)[i]-> GetMass();
	G4int charge = (*HC)[i]-> GetCharge();
	// std::cout<<"pid="<<pid<<", mass="<<mass<<", charge="<<charge<<std::endl;
	// getchar();
	G4int ilay = (*HC)[i]-> GetPadLay();
	//      G4double mass = (*HC)[i]-> GetPDGMass(); //mass(GeV)
	G4int parentid = (*HC)[i]-> GetParentID();
	//G4int parentpid = (*HC)[i]-> GetParentID_pid();
	//Get Parent pid
	G4int parentpid = -9999;
	if(parentid>0){
	  const G4Track* parentTrack = nullptr;

	  for (G4int k = 0; k < anEvent->GetNumberOfPrimaryVertex(); k++) {
	    G4PrimaryVertex* primaryVertex = anEvent->GetPrimaryVertex(k);
	    for (G4int j = 0; j < primaryVertex->GetNumberOfParticle(); j++) {
	      G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary(j);
	      if (primaryParticle->GetTrackID() == parentid) {
		parentpid = primaryParticle->GetPDGcode();
		break;
	      }
	    }
	  }
	}

	G4double tlength = (*HC)[i]-> GettLength();
	G4int irow=(*HC)[i]-> GetPadRow();
	G4double beta = (*HC)[i]-> GetBeta();

	//test edep info save
	G4double edep;

	//test end

	if(pad_configure==3){
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
	if(pad_configure == 4){
	  // Match Pad hit to layer-summed Edep by (tid, ilay) and position.
	  // Pad is a thin-shell crossing; Edep position is the average over the
	  // thick layer — reject if they differ by >= 50 mm (large-angle mismatch).
	  for(int j=0;j<remain_tid.size();j++){
	    find_hit = false;
	    if(tid == remain_tid[j]){
	      find_track = true;
	      for(int k=0;k<remain_ilay[j].size();k++){
		if(ilay == remain_ilay[j][k]){
		  G4double pos_diff = (remain_pos[j][k] - xyz).mag();
		  if(pos_diff < 50.){
		    edep = remain_edep[j][k];
		    find_hit = true;
		    break;
		  }
		}
	      }
	      if(!find_hit){
		G4cout << "[EventAction] no Edep for Pad hit"
		       << " tid=" << tid << " ilay=" << ilay << G4endl;
	      }
	      break;
	    }
	  }
	  if(!find_track){
	    G4cout << "[EventAction] no Edep track for Pad hit"
		   << " tid=" << tid << G4endl;
	  }

	  if(!find_hit)edep = -9999;
	  //std::cout<<"Layer : "<<ilay<<", Cal : "<<(*HC)[i]->GetEdep()<<", Exp : "<<edep<<std::endl;


	  int hit_check = 0;
	  if(remain_ilay.size()>0){
	    for(int f=0;f<remain_ilay.size();++f){
	      if(remain_ilay[f].size()>0){
		hit_check+=remain_ilay[f].size();
		for(int k=0;k<remain_ilay[f].size();++k){
		  //std::cout<<"layer : "<<remain_ilay[f][k]<<std::endl;
		}
	      }
	    }
	  }
	  //std::cout<<"nhit Cal : "<<nhits<<"nhit Edep : "<<hit_check<<std::endl;
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
