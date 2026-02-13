// -*- C++ -*-

#include "BeamMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TFile.h>
#include <TTree.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
namespace
{
  int ntBeam, evnum, runnum;
  double xout[5];
  double yout[5];
  double uout[5];
  double vout[5];
  double pBeam[5];
  double qBeam[5];
  double m2Beam[5];
  int trigpat[32];
  int trigflag[32];

  double fermi_mass;
  double fermi_mom_x;
  double fermi_mom_y;
  double fermi_mom_z;
  double fermi_weight;




  int ntK18, ntKurama;
  TTreeReaderValue<int> *ntTPCK18;
  TTreeReaderValue<int> *ntTPCKurama;

  TTreeReaderValue<vector<int> > *isgoodTPCK18 = nullptr;
  TTreeReaderValue<vector<int> > *isgoodTPCKurama = nullptr;

  TTreeReaderValue<vector<int> > *inside = nullptr;
  TTreeReaderValue<vector<int> > *kflagTPCKurama = nullptr;
  TTreeReaderValue<vector<double> > *pHS = nullptr;
  TTreeReaderValue<vector<double> > *xbTPC = nullptr;
  TTreeReaderValue<vector<double> > *ybTPC = nullptr;
  TTreeReaderValue<vector<double> > *ubTPC = nullptr;
  TTreeReaderValue<vector<double> > *vbTPC = nullptr;

  TTreeReaderValue<vector<double> > *xsTPC = nullptr;
  TTreeReaderValue<vector<double> > *ysTPC = nullptr;
  TTreeReaderValue<vector<double> > *usTPC = nullptr;
  TTreeReaderValue<vector<double> > *vsTPC = nullptr;

  TTreeReaderValue<vector<double> > *vtxTPC = nullptr;
  TTreeReaderValue<vector<double> > *vtyTPC = nullptr;
  TTreeReaderValue<vector<double> > *vtzTPC = nullptr;

  TTreeReaderValue<vector<double> > *xtgtHS = nullptr;
  TTreeReaderValue<vector<double> > *ytgtHS = nullptr;
  TTreeReaderValue<vector<double> > *ztgtHS = nullptr;
  TTreeReaderValue<vector<double> > *utgtHS = nullptr;
  TTreeReaderValue<vector<double> > *vtgtHS = nullptr;
  TTreeReaderValue<vector<double> > *pTPCKurama = nullptr;
  TTreeReaderValue<vector<double> > *qTPCKurama = nullptr;
  TTreeReaderValue<vector<double> > *m2TPCKurama = nullptr;
  TTreeReaderValue<vector<double> > *MissMass = nullptr;
  TTreeReaderValue<vector<double> > *utgtTPCKurama = nullptr;
  TTreeReaderValue<vector<double> > *vtgtTPCKurama = nullptr;
  TTreeReaderValue<vector<double> > *xtgtKurama = nullptr;
  TTreeReaderValue<vector<double> > *ytgtKurama = nullptr;

  // For kurama data
  TTreeReaderValue<vector<vector<double> > > *xvpHS = nullptr;
  TTreeReaderValue<vector<vector<double> > > *yvpHS = nullptr;
  TTreeReaderValue<vector<vector<double> > > *zvpHS = nullptr;
  TTreeReaderValue<vector<vector<double> > > *xvpKurama = nullptr;
  TTreeReaderValue<vector<vector<double> > > *yvpKurama = nullptr;
  TTreeReaderValue<vector<vector<double> > > *zvpKurama = nullptr;

  TTreeReaderValue<vector<double> > *xoutK18 = nullptr;
  TTreeReaderValue<vector<double> > *youtK18 = nullptr;
  TTreeReaderValue<vector<double> > *uoutK18 = nullptr;
  TTreeReaderValue<vector<double> > *voutK18 = nullptr;
  TTreeReaderValue<vector<double> > *p_3rd = nullptr;
  TTreeReaderValue<vector<vector<double> > > *layerK18 = nullptr;
  TTreeReaderValue<vector<vector<double> > > *wireK18 = nullptr;
  TTreeReaderValue<vector<vector<double> > > *localhitposK18 = nullptr;

  TTreeReaderValue<vector<double> > *xout_ = nullptr;
  TTreeReaderValue<vector<double> > *yout_ = nullptr;
  TTreeReaderValue<vector<double> > *zout_ = nullptr;
  TTreeReaderValue<vector<double> > *pxout = nullptr;
  TTreeReaderValue<vector<double> > *pyout = nullptr;
  TTreeReaderValue<vector<double> > *pzout = nullptr;
  TTreeReaderValue<vector<vector<double> > > *layer = nullptr;
  TTreeReaderValue<vector<vector<double> > > *wire = nullptr;
  TTreeReaderValue<vector<vector<double> > > *localhitpos = nullptr;
  TTreeReaderValue<vector<double> > *m2Kurama = nullptr;
  TTreeReaderValue<vector<int> > *Kflag = nullptr;


  TTreeReaderValue<vector<double> > *KmMom_x = nullptr;
  TTreeReaderValue<vector<double> > *KmMom_y = nullptr;
  TTreeReaderValue<vector<double> > *KmMom_z = nullptr;
  TTreeReaderValue<vector<double> > *KpMom_x = nullptr;
  TTreeReaderValue<vector<double> > *KpMom_y = nullptr;
  TTreeReaderValue<vector<double> > *KpMom_z = nullptr;
  

  TTreeReaderValue<int> *XiResidualsMultiplicity = nullptr;
  TTreeReaderValue<double> *KFXiPval = nullptr;
  TTreeReaderValue<double> *KFXiPx = nullptr;
  TTreeReaderValue<double> *KFXiPy = nullptr;
  TTreeReaderValue<double> *KFXiPz = nullptr;

  TTreeReaderValue<bool> *XiFlight = nullptr;
  TTreeReaderValue<bool> *Xiflag = nullptr;
  TTreeReaderValue<bool> *Lflag = nullptr;
  TTreeReaderValue<bool> *LLflag = nullptr;

  TTreeReaderValue<double> *xiprodvtx_x = nullptr;
  TTreeReaderValue<double> *xiprodvtx_y = nullptr;
  TTreeReaderValue<double> *xiprodvtx_z = nullptr;
  TTreeReaderValue<double> *xiprodmom_x = nullptr;
  TTreeReaderValue<double> *xiprodmom_y = nullptr;
  TTreeReaderValue<double> *xiprodmom_z = nullptr;

  TTreeReaderValue<double> *KFXiMom_x = nullptr;
  TTreeReaderValue<double> *KFXiMom_y = nullptr;
  TTreeReaderValue<double> *KFXiMom_z = nullptr;
  TTreeReaderValue<double> *GFXiDecayVtx_x = nullptr;
  TTreeReaderValue<double> *GFXiDecayVtx_y = nullptr;
  TTreeReaderValue<double> *GFXiDecayVtx_z = nullptr;

  TTreeReaderValue<double> *GFLambdaDecayVtx_x = nullptr;
  TTreeReaderValue<double> *GFLambdaDecayVtx_y = nullptr;
  TTreeReaderValue<double> *GFLambdaDecayVtx_z = nullptr;
  

  TTreeReaderValue<vector<double> > *KFDecaysMom_x = nullptr;
  TTreeReaderValue<vector<double> > *KFDecaysMom_y = nullptr;
  TTreeReaderValue<vector<double> > *KFDecaysMom_z = nullptr;
  
  TTreeReaderValue<vector<double> > *KFXiDecaysMom_x = nullptr;
  TTreeReaderValue<vector<double> > *KFXiDecaysMom_y = nullptr;
  TTreeReaderValue<vector<double> > *KFXiDecaysMom_z = nullptr;
  TTreeReaderValue<double> *KFXiProductionVtx_x = nullptr;
  TTreeReaderValue<double> *KFXiProductionVtx_y = nullptr;
  TTreeReaderValue<double> *KFXiProductionVtx_z = nullptr;
  TTreeReaderValue<double> *KFXiProductionVtxMom_x = nullptr;
  TTreeReaderValue<double> *KFXiProductionVtxMom_y = nullptr;
  TTreeReaderValue<double> *KFXiProductionVtxMom_z = nullptr;
}
G4double
BeamInfo::GetX(G4double offset) const
{
  return x + std::tan(-1. * u * CLHEP::mrad) * offset;
}

//_____________________________________________________________________________
G4double
BeamInfo::GetY(G4double offset) const
{
  return y + std::tan(-1. * v * CLHEP::mrad) * offset;
}

//_____________________________________________________________________________
G4int BeamInfo::GetTrigPat(G4int flag) const
{
  return trigpat[flag];
}

//_____________________________________________________________________________
void BeamInfo::Print(void) const
{
  PrintHelper helper(4, std::ios::fixed, G4cout);
  const G4int w = 8;
  G4cout << "   "
	 << "x=" << std::setw(w) << x << " "
	 << "y=" << std::setw(w) << y << " "
	 << "u=" << std::setw(w) << u << " "
	 << "v=" << std::setw(w) << v << " "
	 << "p=(" << std::setw(w) << p.x() << ", "
	 << std::setw(w) << p.y() << ", "
	 << std::setw(w) << p.z() << ")" << G4endl;
}

//_____________________________________________________________________________
BeamMan::BeamMan(void)
  : m_is_ready(false),
    m_file_name(),
    m_file(),
    m_param_array(),
    m_n_param(),
    m_is_vi(false),
    m_primary_z(0.)
{
}

//_____________________________________________________________________________
BeamMan::~BeamMan(void)
{
}

//_____________________________________________________________________________
G4bool
BeamMan::Initialize(void)
{
  const auto &gConf = ConfMan::GetInstance();
  const auto &gGeom = DCGeomMan::GetInstance();
  const G4double p0 = gConf.Get<G4double>("BeamMom");
	
  if (m_file_name.empty())
    return true;

	m_param_array.clear();
	m_mm_array.clear();
	m_kmkpl_array.clear();
	G4int generator = gConf.Get<G4int>("Generator");
	if(generator == 25) return true;
	m_is_vi = (gConf.Get<G4int>("Generator") == 10);
	if (abs(generator) == 135 or abs(generator) == 493 or abs(generator) == 938){
//		m_is_k18 = 1;
		G4cout << "Generating K18 Beam" << G4endl;
	}
	if (abs(generator) == 100){
		m_is_kurama = 1;
	}
	if (generator == 181321 or generator == 1321938 or generator == 1321139 or generator == -1321139){
		m_is_missmassXi = 1;
	}
	if (generator == 181530){
		m_is_missmassXi1530 = 1;
	}
	if (generator == 1001321 or generator == 1001321938){
		m_is_TPCXi = 1;
	}
	if (generator == 1101321){
		m_is_HMCXi = 1;
	}
	if (generator == 25){
		m_is_KpUniform = 1;
	}
	if (generator == 26){
		m_is_CosKp = 1;
	}
	if (generator == 1811161116){
		m_is_LL_BE = 1;
		G4cout<<"LL"<<G4endl;
	}
	m_primary_z = gGeom.GetLocalZ("Vertex");
	//  m_target_z = -143;
	// G4cout<<"BeamMan"<<G4endl;
	m_target_z = gGeom.GetLocalZ("SHSTarget");
	if (!m_is_vi)
		m_primary_z -= 1318.9 * CLHEP::mm; // from VO
	TTree *tree = nullptr;
	TTreeReader *reader = nullptr;
	m_file = new TFile(m_file_name);
	if (m_is_k18)
	{
		tree = (TTree *)m_file->Get("k18track");
	}
	else if (m_is_kurama)
	{
		tree = (TTree *)m_file->Get("kurama");
	}
	else if (m_is_missmassXi or m_is_missmassXi1530 or m_is_TPCXi or m_is_LL_BE or m_is_HMCXi)
	{
		tree = (TTree *)m_file->Get("tpc");
	}
	else if (m_is_KpUniform)
	{
		auto h = (TH2D *)m_file->Get("K18HSTgtProfile");
		HitProfile = (TH2D *)h->Clone("K18HSTgtProfile");
		delete h;
		G4cout << "HitProfile is set:: " << HitProfile->GetEntries() << G4endl;
		m_is_ready = true;
	}
	else if(m_is_CosKp){
		auto h = (TH1D *)m_file->Get("CarbonCosKp");
		m_hCosKp = (TH1D *)h->Clone("CarbonCosKp");
		delete h;
		G4cout << "CosKp is set:: " << m_hCosKp->GetEntries() << G4endl;
		m_is_ready = true;
	}
	else
	{
		tree = (TTree *)m_file->Get("tree");
	}

  if (!m_file->IsOpen() || !tree)
    return m_is_ready;
  G4cout << tree->GetEntries() << G4endl;
  BeamInfo beam;
  if (m_is_k18)
    {
      tree->SetBranchAddress("ntK18", &ntBeam);
      tree->SetBranchAddress("evnum", &evnum);
      tree->SetBranchAddress("runnum", &runnum);
      tree->SetBranchAddress("xout", xout);
      tree->SetBranchAddress("yout", yout);
      tree->SetBranchAddress("uout", uout);
      tree->SetBranchAddress("vout", vout);
      tree->SetBranchAddress("pHS", pBeam);
      tree->SetBranchAddress("trigpat", trigpat);
      tree->SetBranchAddress("trigflag", trigflag);
    }
  else if (m_is_kurama)
    {
      tree->SetBranchAddress("ntKurama", &ntBeam);
      tree->SetBranchAddress("evnum", &evnum);
      tree->SetBranchAddress("runnum", &runnum);
      tree->SetBranchAddress("xtgtKurama", xout);
      tree->SetBranchAddress("ytgtKurama", yout);
      tree->SetBranchAddress("utgtKurama", uout);
      tree->SetBranchAddress("vtgtKurama", vout);
      tree->SetBranchAddress("pKurama", pBeam);
      tree->SetBranchAddress("qKurama", qBeam);
      tree->SetBranchAddress("m2", m2Beam);
      tree->SetBranchAddress("trigpat", trigpat);
    }
#if 0
	else if (m_is_missmassXi or m_is_missmassXi1530){
		reader = new TTreeReader("tpc",m_file);
//		tree->SetBranchAddress( "evnum",&evnum);
//		tree->SetBranchAddress( "runnum",&runnum);
		ntTPCK18 = new TTreeReaderValue<int>(*reader,"ntK18");
		isgoodTPCK18 = new TTreeReaderValue<vector<int>>(*reader,"isgoodTPCK18");
		inside = new TTreeReaderValue<vector<int>>(*reader,"inside");
		MissMass = new TTreeReaderValue<vector<double>>(*reader,"MissMassCorrDETPC");
		pHS = new TTreeReaderValue<vector<double>>(*reader,"pK18");
		xbTPC = new TTreeReaderValue<vector<double>>(*reader,"xbTPC");
		ybTPC = new TTreeReaderValue<vector<double>>(*reader,"ybTPC");
		ubTPC = new TTreeReaderValue<vector<double>>(*reader,"ubTPC");
		vbTPC = new TTreeReaderValue<vector<double>>(*reader,"vbTPC");
		vtxTPC = new TTreeReaderValue<vector<double>>(*reader,"vtxTPC");
		vtyTPC = new TTreeReaderValue<vector<double>>(*reader,"vtyTPC");
		vtzTPC = new TTreeReaderValue<vector<double>>(*reader,"vtzTPC");



    //		tree->SetBranchAddress("ntKurama",&ntKurama);
    ntTPCKurama = new TTreeReaderValue<int>(*reader,"ntKurama");
    isgoodTPCKurama = new TTreeReaderValue<vector<int>>(*reader,"isgoodTPCKurama");
    pTPCKurama = new TTreeReaderValue<vector<double>>(*reader,"pCorrDETPC");
    qTPCKurama = new TTreeReaderValue<vector<double>>(*reader,"qTPCKurama");
    xsTPC = new TTreeReaderValue<vector<double>>(*reader,"xsTPC");
    ysTPC = new TTreeReaderValue<vector<double>>(*reader,"ysTPC");
    usTPC = new TTreeReaderValue<vector<double>>(*reader,"usTPC");
    vsTPC = new TTreeReaderValue<vector<double>>(*reader,"vsTPC");
    m2TPCKurama = new TTreeReaderValue<vector<double>>(*reader,"m2TPCKurama");

    // For kurama data
    p_3rd = new TTreeReaderValue<vector<double>>(*reader,"p_3rd");
    xoutK18 = new TTreeReaderValue<vector<double>>(*reader,"xoutK18");
    youtK18 = new TTreeReaderValue<vector<double>>(*reader,"youtK18");
    uoutK18 = new TTreeReaderValue<vector<double>>(*reader,"uoutK18");
    voutK18 = new TTreeReaderValue<vector<double>>(*reader,"voutK18");
    xvpHS = new TTreeReaderValue<vector<vector<double>>>(*reader,"xvpHS");
    yvpHS = new TTreeReaderValue<vector<vector<double>>>(*reader,"yvpHS");
    zvpHS = new TTreeReaderValue<vector<vector<double>>>(*reader,"zvpHS");
    xtgtHS = new TTreeReaderValue<vector<double>>(*reader,"xtgtHS");
    ytgtHS = new TTreeReaderValue<vector<double>>(*reader,"ytgtHS");
    ztgtHS = new TTreeReaderValue<vector<double>>(*reader,"ztgtHS");
    layerK18 = new TTreeReaderValue<vector<vector<double>>>(*reader,"layerK18");
    wireK18 = new TTreeReaderValue<vector<vector<double>>>(*reader,"wireK18");
    localhitposK18 = new TTreeReaderValue<vector<vector<double>>>(*reader,"localhitposK18");

    xvpKurama = new TTreeReaderValue<vector<vector<double>>>(*reader,"xvpKurama");
    yvpKurama = new TTreeReaderValue<vector<vector<double>>>(*reader,"yvpKurama");
    zvpKurama = new TTreeReaderValue<vector<vector<double>>>(*reader,"zvpKurama");
    xtgtKurama = new TTreeReaderValue<vector<double>>(*reader,"xtgtKurama");
    ytgtKurama = new TTreeReaderValue<vector<double>>(*reader,"ytgtKurama");
    xout_ = new TTreeReaderValue<vector<double>>(*reader,"xout");
    yout_ = new TTreeReaderValue<vector<double>>(*reader,"yout");
    zout_ = new TTreeReaderValue<vector<double>>(*reader,"zout");
    pxout = new TTreeReaderValue<vector<double>>(*reader,"pxout");
    pyout = new TTreeReaderValue<vector<double>>(*reader,"pyout");
    pzout = new TTreeReaderValue<vector<double>>(*reader,"pzout");
    layer = new TTreeReaderValue<vector<vector<double>>>(*reader,"layer");
    wire = new TTreeReaderValue<vector<vector<double>>>(*reader,"wire");
    localhitpos = new TTreeReaderValue<vector<vector<double>>>(*reader,"localhitpos");

  }
  else if(m_is_TPCXi){
#else
	else if (m_is_TPCXi or m_is_missmassXi or m_is_missmassXi1530 or m_is_HMCXi){
#endif
		tree->SetBranchAddress("evnum", &evnum);
		tree->SetBranchAddress("runnum", &runnum);
	    reader = new TTreeReader("tpc", m_file);

	    ntTPCK18 = new TTreeReaderValue<int>(*reader, "ntK18");
	    pHS = new TTreeReaderValue<vector<double> >(*reader, "pK18");
	    xbTPC = new TTreeReaderValue<vector<double> >(*reader, "xtgtK18");
	    ybTPC = new TTreeReaderValue<vector<double> >(*reader, "ytgtK18");
	    ubTPC = new TTreeReaderValue<vector<double> >(*reader, "utgtK18");
	    vbTPC = new TTreeReaderValue<vector<double> >(*reader, "vtgtK18");

		ntTPCKurama = new TTreeReaderValue<int>(*reader, "ntKurama");
		isgoodTPCKurama = new TTreeReaderValue<vector<int> >(*reader, "isgoodTPCKurama");
		pTPCKurama = new TTreeReaderValue<vector<double> >(*reader, "pCorrDETPC");
		qTPCKurama = new TTreeReaderValue<vector<double> >(*reader, "qTPCKurama");
		xsTPC = new TTreeReaderValue<vector<double> >(*reader, "xtgtTPCKurama");
		ysTPC = new TTreeReaderValue<vector<double> >(*reader, "ytgtTPCKurama");
		usTPC = new TTreeReaderValue<vector<double> >(*reader, "utgtTPCKurama");
		vsTPC = new TTreeReaderValue<vector<double> >(*reader, "vtgtTPCKurama");
		m2TPCKurama = new TTreeReaderValue<vector<double> >(*reader, "m2TPCKurama");
		inside = new TTreeReaderValue<vector<int> >(*reader, "insideTPC");
		kflagTPCKurama = new TTreeReaderValue<vector<int> >(*reader, "kflagTPCKurama");
		MissMass = new TTreeReaderValue<vector<double> >(*reader, "MissMassCorrDETPC");
		vtxTPC = new TTreeReaderValue<vector<double> >(*reader, "vtxTPC");
		vtyTPC = new TTreeReaderValue<vector<double> >(*reader, "vtyTPC");
		vtzTPC = new TTreeReaderValue<vector<double> >(*reader, "vtzTPC");
		
		KmMom_x = new TTreeReaderValue<vector<double> >(*reader, "KmMom_x");
		KmMom_y = new TTreeReaderValue<vector<double> >(*reader, "KmMom_y");
		KmMom_z = new TTreeReaderValue<vector<double> >(*reader, "KmMom_z");
		KpMom_x = new TTreeReaderValue<vector<double> >(*reader, "KpMom_x");
		KpMom_y = new TTreeReaderValue<vector<double> >(*reader, "KpMom_y");
		KpMom_z = new TTreeReaderValue<vector<double> >(*reader, "KpMom_z");


	    Xiflag = new TTreeReaderValue<bool>(*reader, "Xiflag");

	    KFXiProductionVtx_x = new TTreeReaderValue<double>(*reader, "KFXiProductionVtx_x");
	    KFXiProductionVtx_y = new TTreeReaderValue<double>(*reader, "KFXiProductionVtx_y");
	    KFXiProductionVtx_z = new TTreeReaderValue<double>(*reader, "KFXiProductionVtx_z");
	    KFXiProductionVtxMom_x = new TTreeReaderValue<double>(*reader, "KFXiProductionVtxMom_x");
	    KFXiProductionVtxMom_y = new TTreeReaderValue<double>(*reader, "KFXiProductionVtxMom_y");
	    KFXiProductionVtxMom_z = new TTreeReaderValue<double>(*reader, "KFXiProductionVtxMom_z");

		KFXiPval = new TTreeReaderValue<double>(*reader, "KFXiPval");
		XiResidualsMultiplicity = new TTreeReaderValue<int>(*reader, "XiResidualsMultiplicity");
		GFXiDecayVtx_x = new TTreeReaderValue<double >(*reader, "GFXiDecayVtx_x");
		GFXiDecayVtx_y = new TTreeReaderValue<double >(*reader, "GFXiDecayVtx_y");
		GFXiDecayVtx_z = new TTreeReaderValue<double >(*reader, "GFXiDecayVtx_z");

		GFLambdaDecayVtx_x = new TTreeReaderValue<double >(*reader, "GFLambdaDecayVtx_x");
		GFLambdaDecayVtx_y = new TTreeReaderValue<double >(*reader, "GFLambdaDecayVtx_y");
		GFLambdaDecayVtx_z = new TTreeReaderValue<double >(*reader, "GFLambdaDecayVtx_z");

	  	KFXiDecaysMom_x = new TTreeReaderValue<vector<double> >(*reader, "KFXiDecaysMom_x");
	  	KFXiDecaysMom_y = new TTreeReaderValue<vector<double> >(*reader, "KFXiDecaysMom_y");
	  	KFXiDecaysMom_z = new TTreeReaderValue<vector<double> >(*reader, "KFXiDecaysMom_z");
	  }
	else if (m_is_LL_BE)
	  {
	    reader = new TTreeReader("tpc", m_file);

		ntTPCK18 = new TTreeReaderValue<int>(*reader, "ntK18");
		pHS = new TTreeReaderValue<vector<double> >(*reader, "pK18");
		xbTPC = new TTreeReaderValue<vector<double> >(*reader, "xtgtK18");
		ybTPC = new TTreeReaderValue<vector<double> >(*reader, "ytgtK18");
		ubTPC = new TTreeReaderValue<vector<double> >(*reader, "utgtK18");
		vbTPC = new TTreeReaderValue<vector<double> >(*reader, "vtgtK18");
		
		ntTPCKurama = new TTreeReaderValue<int>(*reader, "ntKurama");
		isgoodTPCKurama = new TTreeReaderValue<vector<int> >(*reader, "isgoodTPCKurama");
		pTPCKurama = new TTreeReaderValue<vector<double> >(*reader, "pCorrDETPC");
		qTPCKurama = new TTreeReaderValue<vector<double> >(*reader, "qTPCKurama");
		xsTPC = new TTreeReaderValue<vector<double> >(*reader, "xtgtTPCKurama");
		ysTPC = new TTreeReaderValue<vector<double> >(*reader, "ytgtTPCKurama");
		usTPC = new TTreeReaderValue<vector<double> >(*reader, "utgtTPCKurama");
		vsTPC = new TTreeReaderValue<vector<double> >(*reader, "vtgtTPCKurama");
		Kflag = new TTreeReaderValue<vector<int> >(*reader, "Kflag");

		vtxTPC = new TTreeReaderValue<vector<double> >(*reader, "vtxTPC");
		vtyTPC = new TTreeReaderValue<vector<double> >(*reader, "vtyTPC");
	
		Lflag = new TTreeReaderValue<bool>(*reader, "Lflag");
		LLflag = new TTreeReaderValue<bool>(*reader, "LLflag");
		Xiflag = new TTreeReaderValue<bool>(*reader, "Xiflag");
		KFDecaysMom_x = new TTreeReaderValue<vector<double> >(*reader, "KFDecaysMom_x");	
		KFDecaysMom_y = new TTreeReaderValue<vector<double> >(*reader, "KFDecaysMom_y");	
		KFDecaysMom_z = new TTreeReaderValue<vector<double> >(*reader, "KFDecaysMom_z");	
	/*	
		xvpHS = new TTreeReaderValue<vector<vector<double> > >(*reader, "xvpHS");
		yvpHS = new TTreeReaderValue<vector<vector<double> > >(*reader, "yvpHS");
		zvpHS = new TTreeReaderValue<vector<vector<double> > >(*reader, "zvpHS");
		xtgtHS = new TTreeReaderValue<vector<double> >(*reader, "xtgtHS");
		ytgtHS = new TTreeReaderValue<vector<double> >(*reader, "ytgtHS");
		ztgtHS = new TTreeReaderValue<vector<double> >(*reader, "ztgtHS");

		xvpKurama = new TTreeReaderValue<vector<vector<double> > >(*reader, "xvpKurama");
		yvpKurama = new TTreeReaderValue<vector<vector<double> > >(*reader, "yvpKurama");
		zvpKurama = new TTreeReaderValue<vector<vector<double> > >(*reader, "zvpKurama");
		xtgtKurama = new TTreeReaderValue<vector<double> >(*reader, "xtgtKurama");
		ytgtKurama = new TTreeReaderValue<vector<double> >(*reader, "ytgtKurama");	
	*/
	}
	else
	  {
	    tree->SetBranchAddress("x", &beam.x);
	    tree->SetBranchAddress("y", &beam.y);
	    tree->SetBranchAddress("u", &beam.u);
	    tree->SetBranchAddress("v", &beam.v);
	    tree->SetBranchAddress("p", &beam.dp);
	  }
	//for (Long64_t i = 0, n = tree->GetEntries(); i < n; ++i)
	int scale_down = 1;//For test purpose fast loading
	std::cout << "BeamEvents = " << tree->GetEntries() << std::endl;
	for (Long64_t i = 0, n = tree->GetEntries()/scale_down; i < n; ++i)
	{
		tree->GetEntry(i);
		if (i % 10000 == 0)
			G4cout << Form("Event %lld/%lld", i, tree->GetEntries()) << G4endl;
		if (reader)
			reader->Next();
		if (m_is_k18 or m_is_kurama)
		{
			beam.x = 0;
			beam.y = 0;
			beam.z = 0;
			beam.u = 0;
			beam.v = 0;
			beam.p.set(0, 0, 0);
			beam.evnum = evnum;
			beam.runnum = runnum;
			beam.ntBeam = ntBeam;
			if (m_is_k18 and (trigflag[14] < 0 or trigflag[23] < 0))
				continue;
			for (int it = 0; it < ntBeam; ++it)
			{
				beam.x = xout[0];
				beam.y = yout[0];
				beam.u = uout[0]; // u,v definition = dxdz,dydz, not mrad.
				beam.v = yout[0];
				double pz = pBeam[0] / sqrt(uout[0] * uout[0] + vout[0] * vout[0] + 1);
				beam.p.set(pz * uout[0], pz * vout[0], pz);
				if (m_is_kurama)
				{
					beam.z = m_target_z;
					beam.m2 = m2Beam[0];
					beam.q = qBeam[0];
				}
				else
				{
					beam.z = m_primary_z; // VO
				}
			}
			for (int itrg = 0; itrg < 32; ++itrg)
			{
				beam.trigpat[itrg] = trigpat[itrg];
			}
			m_nBeam++;
			m_param_array.push_back(beam);
		}
		else if (m_is_missmassXi)
		{
			int in = (*kflagTPCKurama)->at(0);
			double mm = (*MissMass)->at(0);
			int ntKurama = **ntTPCKurama;
			int ntK18 = **ntTPCK18;
			if(! in) continue;
			if(!**Xiflag) continue;

			for (int itk18 = 0; itk18 < ntK18; ++itk18)
			{
				if (ntK18 != 1 or ntKurama != 1)
					continue;
				if ((*ubTPC)->size() == 0)
					continue;
				double ub = (*ubTPC)->at(itk18);
				double vb = (*vbTPC)->at(itk18);
				double nb = hypot(hypot(1, ub), vb);
				double pb = (*pHS)->at(itk18);
				double pzb = pb / nb;
				G4ThreeVector TVKm(pzb * ub, pzb * vb, pzb);
				for (int itkurama = 0; itkurama < ntKurama; ++itkurama)
				{
					if ((*pTPCKurama)->size() == 0)
						continue;
					if ((*pTPCKurama)->at(itkurama) == 0)
						continue;
					double qKp = (*qTPCKurama)->at(itkurama);
					double m2Kp = (*m2TPCKurama)->at(itkurama);
					double pKp = (*pTPCKurama)->at(itkurama);

			double us = (*usTPC)->at(itkurama);
			double vs = (*vsTPC)->at(itkurama);
			double ns = hypot(hypot(1, us), vs);
			double pzs = pKp / ns;
			G4ThreeVector TVKp(pzs * us, pzs * vs, pzs);
			if ((*vtzTPC)->size() == 0)
			  continue;
			if ((*vtzTPC)->at(0) == 0)
			  continue;
			if (pKp < 1.4 and qKp > 0 and m2Kp > 0.14 and m2Kp < 0.34 and in and abs(mm - 1.321) < 0.13)
			  {
			    MMVertex MMVert;
			    MMVert.x = ((*vtxTPC)->at(0));
			    MMVert.y = ((*vtyTPC)->at(0));
			    if (abs(MMVert.x) > 15 or abs(MMVert.y) > 10)
			      continue;
			    MMVert.Moms.push_back(TVKm);
			    MMVert.Moms.push_back(TVKp);
			    MMVert.Moms.push_back(TVKm - TVKp);

						//					vector<double>& xoutk18 = *(xoutK18->Get());
#if 0	
						MMVert.p_3rd = *(p_3rd->Get());
						MMVert.xoutK18 = *(xoutK18->Get());
						MMVert.youtK18 = *(youtK18->Get());
						MMVert.uoutK18 = *(uoutK18->Get());
						MMVert.voutK18 = *(voutK18->Get());
						MMVert.xvpHS = *(xvpHS->Get());
						MMVert.yvpHS = *(yvpHS->Get());
						MMVert.zvpHS = *(zvpHS->Get());
						MMVert.xtgtHS = *(xtgtHS->Get());
						MMVert.ytgtHS = *(ytgtHS->Get());
						MMVert.ztgtHS = *(ztgtHS->Get());
						MMVert.layerK18 = *(layerK18->Get());
						MMVert.wireK18 = *(wireK18->Get());
						MMVert.localhitposK18 = *(localhitposK18->Get());

						MMVert.xvpKurama = *(xvpKurama->Get());
						MMVert.yvpKurama = *(yvpKurama->Get());
						MMVert.zvpKurama = *(zvpKurama->Get());
						MMVert.xtgtKurama = *(xtgtKurama->Get());
						MMVert.ytgtKurama = *(ytgtKurama->Get());
						MMVert.xoutK18 = *(xoutK18->Get());
						MMVert.youtK18 = *(youtK18->Get());
						MMVert.uoutK18 = *(uoutK18->Get());
						MMVert.voutK18 = *(voutK18->Get());
						MMVert.xout = *(xout_->Get());
						MMVert.yout = *(yout_->Get());
						MMVert.zout = *(zout_->Get());
						MMVert.pxout = *(pxout->Get());
						MMVert.pyout = *(pyout->Get());
						MMVert.pzout = *(pzout->Get());
						MMVert.layer = *(layer->Get());
						MMVert.wire = *(wire->Get());
						MMVert.localhitpos = *(localhitpos->Get());
#endif
						m_mm_array.push_back(MMVert);
					}
				}
			}
		}
		else if (m_is_missmassXi1530){
			int in = (*kflagTPCKurama)->at(0);
			double mm = (*MissMass)->at(0);
			int ntKurama = **ntTPCKurama;
			int ntK18 = **ntTPCK18;
			for (int itk18 = 0; itk18 < ntK18; ++itk18)
			{
				if (ntK18 != 1 or ntKurama != 1)
					continue;
				if ((*ubTPC)->size() == 0)
					continue;
				double ub = (*ubTPC)->at(itk18);
				double vb = (*vbTPC)->at(itk18);
				double nb = hypot(hypot(1, ub), vb);
				double pb = (*pHS)->at(itk18);
				double pzb = pb / nb;
				G4ThreeVector TVKm(pzb * ub, pzb * vb, pzb);
				for (int itkurama = 0; itkurama < ntKurama; ++itkurama)
				{
					if ((*pTPCKurama)->size() == 0)
						continue;
					if ((*pTPCKurama)->at(itkurama) == 0)
						continue;
					double qKp = (*qTPCKurama)->at(itkurama);
					double m2Kp = (*m2TPCKurama)->at(itkurama);
					double pKp = (*pTPCKurama)->at(itkurama);

					double us = (*usTPC)->at(itkurama);
					double vs = (*vsTPC)->at(itkurama);
					double ns = hypot(hypot(1, us), vs);
					double pzs = pKp / ns;
					G4ThreeVector TVKp(pzs * us, pzs * vs, pzs);
					if ((*vtzTPC)->size() == 0)
						continue;
					if ((*vtzTPC)->at(0) == 0)
						continue;
					if (in and abs(mm - 1.535) < 0.05)
					{
						MMVertex MMVert;
						MMVert.x = ((*vtxTPC)->at(0));
						MMVert.y = ((*vtyTPC)->at(0));
						if (abs(MMVert.x) > 15 or abs(MMVert.y) > 10)
							continue;
						MMVert.Moms.push_back(TVKm);
						MMVert.Moms.push_back(TVKp);
						MMVert.Moms.push_back(TVKm - TVKp);
/*
						//					vector<double>& xoutk18 = *(xoutK18->Get());
						MMVert.p_3rd = *(p_3rd->Get());
						MMVert.xoutK18 = *(xoutK18->Get());
						MMVert.youtK18 = *(youtK18->Get());
						MMVert.uoutK18 = *(uoutK18->Get());
						MMVert.voutK18 = *(voutK18->Get());
						MMVert.xvpHS = *(xvpHS->Get());
						MMVert.yvpHS = *(yvpHS->Get());
						MMVert.zvpHS = *(zvpHS->Get());
						MMVert.xtgtHS = *(xtgtHS->Get());
						MMVert.ytgtHS = *(ytgtHS->Get());
						MMVert.ztgtHS = *(ztgtHS->Get());
						MMVert.layerK18 = *(layerK18->Get());
						MMVert.wireK18 = *(wireK18->Get());
						MMVert.localhitposK18 = *(localhitposK18->Get());

			    MMVert.xvpKurama = *(xvpKurama->Get());
			    MMVert.yvpKurama = *(yvpKurama->Get());
			    MMVert.zvpKurama = *(zvpKurama->Get());
			    MMVert.xtgtKurama = *(xtgtKurama->Get());
			    MMVert.ytgtKurama = *(ytgtKurama->Get());

						MMVert.xoutK18 = *(xoutK18->Get());
						MMVert.youtK18 = *(youtK18->Get());
						MMVert.uoutK18 = *(uoutK18->Get());
						MMVert.voutK18 = *(voutK18->Get());
						MMVert.xout = *(xout_->Get());
						MMVert.yout = *(yout_->Get());
						MMVert.zout = *(zout_->Get());
						MMVert.pxout = *(pxout->Get());
						MMVert.pyout = *(pyout->Get());
						MMVert.pzout = *(pzout->Get());
						MMVert.layer = *(layer->Get());
						MMVert.wire = *(wire->Get());
						MMVert.localhitpos = *(localhitpos->Get());
*/
						m_mm_array.push_back(MMVert);
					}
				}
			}
		}
		else if (m_is_TPCXi){
			if ((*MissMass)->size() == 0) continue;
			double mm = (*MissMass)->at(0);
			if (abs(mm - 1.321) > 0.13) continue;
			if((*kflagTPCKurama)->at(0) == 0) continue;
			if(!**Xiflag) continue;
			if(isnan(**KFXiProductionVtx_z) or **KFXiProductionVtx_z==0) continue;
			if(isnan(**KFXiProductionVtxMom_z) or **KFXiProductionVtxMom_z==0) continue;
			if(**KFXiPval <0.01) continue;
			if(**XiResidualsMultiplicity > 0) continue;
			if(abs(**KFXiProductionVtx_z+143)>10 or abs(**KFXiProductionVtx_x)> 15 or abs(**KFXiProductionVtx_y)> 10) continue;
			// For Accurate Xi Production momentum
			G4ThreeVector TVKm((*KmMom_x)->at(0), (*KmMom_y)->at(0), (*KmMom_z)->at(0));
			G4ThreeVector TVKp((*KpMom_x)->at(0), (*KpMom_y)->at(0), (*KpMom_z)->at(0));	
			MMVertex MMVert;
			double PxXi = **KFXiProductionVtxMom_x;
			double PyXi = **KFXiProductionVtxMom_y;
			double PzXi = **KFXiProductionVtxMom_z;
			MMVert.x = **KFXiProductionVtx_x;
			MMVert.y = **KFXiProductionVtx_y;
			MMVert.z = **KFXiProductionVtx_z;
			MMVert.Moms.push_back(TVKm);
			MMVert.Moms.push_back(TVKp);
			G4ThreeVector TVXi(PxXi, PyXi, PzXi);
			MMVert.Moms.push_back(TVXi);
			m_mm_array.push_back(MMVert);
	      }
	    else if (m_is_LL_BE){
	      KmKpL KmKpLVert;
	      if (**ntTPCK18 != 1 or **ntTPCKurama != 1) continue;
	      if ((*xbTPC)->size() == 0) continue;
	      if ((*pTPCKurama)->at(0) == 0) continue;
	      if ((*qTPCKurama)->at(0) <0) continue;
	      if((*Kflag)->at(0) == 0) continue;
	      if(!**Lflag or **LLflag or **Xiflag) continue;
	      double vx = (*xbTPC)->at(0);
	      double vy = (*ybTPC)->at(0);
	      if(abs(vx)>15) continue;
	      if(abs(vy)>10) continue;
	      double ub = (*ubTPC)->at(0);
	      double vb = (*vbTPC)->at(0);
	      double nb = hypot(hypot(1, ub), vb);
	      double pb = (*pHS)->at(0);
	      double pzb = pb / nb;
	      G4ThreeVector TVKm(pzb * ub, pzb * vb, pzb);
	      double us = (*usTPC)->at(0);
	      double vs = (*vsTPC)->at(0);
	      double ns = hypot(hypot(1, us), vs);
	      double ps = (*pTPCKurama)->at(0);
	      double pzs = ps / ns;
	      G4ThreeVector TVKp(pzs * us, pzs * vs, pzs);
	      double KFpx = (*KFDecaysMom_x)->at(0);
	      double KFpy = (*KFDecaysMom_y)->at(0);
	      double KFpz = (*KFDecaysMom_z)->at(0);
	      double KFpix = (*KFDecaysMom_x)->at(1);
	      double KFpiy = (*KFDecaysMom_y)->at(1);
	      double KFpiz = (*KFDecaysMom_z)->at(1);
	      G4ThreeVector TVL(KFpx+KFpix, KFpy+KFpiy, KFpz+KFpiz);

	      KmKpLVert.x = vx;
	      KmKpLVert.y = vy;
	      KmKpLVert.z = -143;
	      KmKpLVert.Moms.push_back(TVKm);
	      KmKpLVert.Moms.push_back(TVKp);
	      KmKpLVert.Moms.push_back(TVL);
	      m_kmkpl_array.push_back(KmKpLVert);
	    }
	 	else if (m_is_HMCXi){
			if ((*MissMass)->size() == 0)
				continue;
			double mm = (*MissMass)->at(0);
			if (abs(mm - 1.321) > 0.1)//loose Xi selection, 
				continue;
			if (!(**Xiflag))
				continue;
			if(isnan(**KFXiProductionVtx_z) or **KFXiProductionVtx_z==0)
				continue;
			if(isnan(**KFXiProductionVtxMom_z) or **KFXiProductionVtxMom_z==0)
				continue;
			HMCXi HMCXiVert;
			HMCXiVert.runnum = runnum;
			HMCXiVert.evnum = evnum;
			HMCXiVert.x = **KFXiProductionVtx_x;
			HMCXiVert.y = **KFXiProductionVtx_y;
			HMCXiVert.z = **KFXiProductionVtx_z;
			G4ThreeVector TVKm((*KmMom_x)->at(0), (*KmMom_y)->at(0), (*KmMom_z)->at(0));
			G4ThreeVector TVKp((*KpMom_x)->at(0), (*KpMom_y)->at(0), (*KpMom_z)->at(0));
			HMCXiVert.KMoms.push_back(TVKm);
			HMCXiVert.KMoms.push_back(TVKp);
			
			G4ThreeVector TVP((*KFXiDecaysMom_x)->at(0), (*KFXiDecaysMom_y)->at(0), (*KFXiDecaysMom_z)->at(0));
			G4ThreeVector TVPi1((*KFXiDecaysMom_x)->at(1), (*KFXiDecaysMom_y)->at(1), (*KFXiDecaysMom_z)->at(1));
			G4ThreeVector TVPi2((*KFXiDecaysMom_x)->at(2), (*KFXiDecaysMom_y)->at(2), (*KFXiDecaysMom_z)->at(2));
			G4ThreeVector TVL = TVP + TVPi1;
			G4ThreeVector TVXi = TVL + TVPi2;
			HMCXiVert.BaryonMoms.push_back(TVXi);
			HMCXiVert.BaryonMoms.push_back(TVL);
			HMCXiVert.BaryonMoms.push_back(TVP);
			G4ThreeVector VtxL((**GFLambdaDecayVtx_x),(**GFLambdaDecayVtx_y),(**GFLambdaDecayVtx_z));
			G4ThreeVector VtxXi((**GFXiDecayVtx_x),(**GFXiDecayVtx_y),(**GFXiDecayVtx_z));
			HMCXiVert.Vertices.push_back(VtxXi);
			HMCXiVert.Vertices.push_back(VtxL);
			m_hmcxi_array.push_back(HMCXiVert);
		}   
		else
	      {
		beam.x *= -1. * CLHEP::cm;							  // -cm -> mm
		beam.y *= -1. * CLHEP::cm;							  // -cm -> mm
		G4double dxdz = std::tan(-1. * beam.u * CLHEP::mrad); // -mrad -> tan
		G4double dydz = std::tan(-1. * beam.v * CLHEP::mrad); // -mrad -> tan
		G4double pp = p0 * (1. + beam.dp * CLHEP::perCent);	  // dp/p[%] -> GeV/c
		G4double pz = pp / std::sqrt(dxdz * dxdz + dydz * dydz + 1.);
		beam.x += dxdz * m_primary_z;
		beam.y += dydz * m_primary_z;
		beam.z = m_primary_z;
		beam.p.set(pz * dxdz, pz * dydz, pz);
		m_param_array.push_back(beam);
	      }
	  }
	G4cout << "Closing files" << G4endl;
	m_file->Close();
	G4cout << "File Closed" << G4endl;
	
	
	
	if(m_is_TPCXi or m_is_missmassXi or m_is_missmassXi1530){
		G4cout<<"Number of seeds: "<<m_mm_array.size()<<G4endl;
	}
	if(m_is_HMCXi){
		G4cout<<"Number of seeds: "<<m_hmcxi_array.size()<<G4endl;
		G4cout<<"Filling dummy HMCXi Vertices"<<G4endl;
		for(int i=0;i<1e5;++i){
			HMCXi dummy;
			dummy.runnum = 0;
			dummy.evnum = 0;
			dummy.x = 0;
			dummy.y = 0;
			dummy.z = 0;
			G4ThreeVector TV(0,0,0);
			dummy.KMoms.push_back(TV);
			dummy.KMoms.push_back(TV);
			dummy.BaryonMoms.push_back(TV);
			dummy.BaryonMoms.push_back(TV);
			dummy.BaryonMoms.push_back(TV);
			dummy.Vertices.push_back(TV);
			dummy.Vertices.push_back(TV);
			m_hmcxi_array.push_back(dummy);
		}
	}
	if(m_is_LL_BE){
	  G4cout<<"Number of seeds: "<<m_kmkpl_array.size()<<G4endl;
	}
	m_n_param = m_param_array.size();
	m_is_ready = true;
	

  if (m_accidental_name.empty()){
  }
  else{
	m_accidental_file = new TFile(m_accidental_name);
	TTree* tree_acc = (TTree*)m_accidental_file->Get("k18track");
	tree_acc->SetBranchAddress("ntK18", &ntBeam);
	tree_acc->SetBranchAddress("xtgtHS", xout);
	tree_acc->SetBranchAddress("ytgtHS", yout);
	tree_acc->SetBranchAddress("utgtHS", uout);
	tree_acc->SetBranchAddress("vtgtHS", vout);
	tree_acc->SetBranchAddress("pHS", pBeam);
	int entries_acc = tree_acc->GetEntries();
	G4cout << "Accidental events = " << entries_acc << G4endl;
	BeamInfo accidental;
	double prop = abs(m_target_z - (-250))+ 20;
	for( int iev=0;iev < entries_acc/scale_down;++iev){
		tree_acc->GetEntry(iev);
		if(iev % 100000 == 0) G4cout<<"Reading Accidental Events "<<iev<<" / "<<entries_acc<<G4endl;
		if(ntBeam != 1) continue;
		double mom = pBeam[0];
		double radi = mom * 1000 / 0.3; // mm
		double dth = prop / radi;	
		double u = uout[0] + dth,v = vout[0];
		accidental.u = u;// bending due to B field
		accidental.x = xout[0] - (u + uout[0])/2 *prop;// at the enterance of the TPC	
		double y_rand = G4RandFlat::shoot(-300,300);
		accidental.y = yout[0] + y_rand;
		accidental.z = m_target_z - prop; // just outside the TPC. accidental should not hit the HS Magnet.
		accidental.v = vout[0];
		G4ThreeVector p(pBeam[0] * accidental.u, pBeam[0] * vout[0], pBeam[0]/sqrt(u*u + v*v + 1));
		accidental.p = p;
		m_accidental_array.push_back(accidental);
	}
	m_n_accidental = m_accidental_array.size();	
  }
  if (m_fermi_name.empty()){
	}
	else{
		m_fermi_file = new TFile(m_fermi_name);
		TTree* tree_fermi = (TTree*)m_fermi_file->Get("tree");
		tree_fermi -> SetBranchAddress("mm_mass_fermi",&fermi_mass);
		tree_fermi -> SetBranchAddress("mm_mom_beam_x",&fermi_mom_x);
		tree_fermi -> SetBranchAddress("mm_mom_beam_y",&fermi_mom_y);
		tree_fermi -> SetBranchAddress("mm_mom_beam_z",&fermi_mom_z);
		tree_fermi -> SetBranchAddress("weight",&fermi_weight);
		int entries_fermi = tree_fermi->GetEntries();
		BeamInfo fermimom;
		G4cout<< "Fermi Momentum Seeds = "<< entries_fermi <<G4endl;
		double weight_sum = 0;
		for( int iev=0;iev<entries_fermi;++iev){
			tree_fermi->GetEntry(iev);
			if(iev % 100000 == 0) G4cout<<"Reading Fermi Seeds "<<iev<<" / "<<entries_fermi<<G4endl;
			G4ThreeVector fermi_mom(fermi_mom_x,fermi_mom_y,fermi_mom_z);
			if(fermi_mom.mag()>0.5) continue;
			weight_sum += fermi_weight;
			fermimom.u = fermi_mass;	
			fermimom.x = fermi_mom_x;	
			fermimom.y = fermi_mom_y;	
			fermimom.z = fermi_mom_z;	
			fermimom.weight = fermi_weight;	
			fermimom.weight_cum = weight_sum;	
			m_fermi_array.push_back(fermimom);
		}
		m_n_fermi = m_fermi_array.size();
		m_weight_fermi = weight_sum;
		for(auto f : m_fermi_array){
			f.weight_tot = weight_sum;
		}
	}
  
	const int PolarizedDecay = gConf.Get<G4int>("PolarizedDecay");
	const TString PolarizationFile = gConf.Get<G4String>("PolarizationFile");
	G4cout<<"PolarizedDecay: "<<PolarizedDecay<<G4endl;
	G4cout<<"PolarizationFile: "<<PolarizationFile<<G4endl;
	if(PolarizationFile == ""){
		const double pol = gConf.Get<G4double>("XiPolarization");
		cth_0s.push_back(1.0);
		cth_1s.push_back(-1.0);
		sqrts_0s.push_back(0);
		sqrts_1s.push_back(999);
		pols.push_back(0.0);
	}
	else{
		TFile* pol_file = new TFile(PolarizationFile);
		int conf,AngBin;
		double cth_0,cth_1,x_0,x_1,PXi;
		TTree* pol_tree = (TTree*)pol_file->Get("tree");
		pol_tree->SetBranchAddress("conf",&conf);
		pol_tree->SetBranchAddress("AngBin",&AngBin);
		pol_tree->SetBranchAddress("cth_0",&cth_0);
		pol_tree->SetBranchAddress("cth_1",&cth_1);
		pol_tree->SetBranchAddress("x_0",&x_0);
		pol_tree->SetBranchAddress("x_1",&x_1);
		pol_tree->SetBranchAddress("PXi",&PXi);
		int nentries = pol_tree->GetEntries();
		int nc = 0;
		for(int i=0;i<nentries;++i){
			pol_tree->GetEntry(i);
			if(conf != 0) continue;
			if(AngBin == -1) continue;
			cth_0s.push_back(cth_0);
			cth_1s.push_back(cth_1);
			sqrts_0s.push_back(x_0);
			sqrts_1s.push_back(x_1);
			pols.push_back(PXi);
			nc++;
		}
		G4cout<<"Loaded "<<nc<<" polarization bins"<<G4endl;
	}
	G4cout<<"Loading Polarization Tables"<<G4endl;
	LoadPolarizationTables();
//	PrintPolarizationTables();
	return true;
	
  }

  //_____________________________________________________________________________
  G4bool
    BeamMan::Initialize(const G4String &filename, const G4String &accidental_name)
  {
    m_file_name = filename;
    m_accidental_name = accidental_name;
    return Initialize();
  }
  
	G4bool
    BeamMan::Initialize(const G4String &filename, const G4String &accidental_name, const G4String &fermi_name)
  {
    m_file_name = filename;
    m_accidental_name = accidental_name;
    m_fermi_name = fermi_name;
    return Initialize();
  }

  //_____________________________________________________________________________
  const BeamInfo &
    BeamMan::Get(void) const
  {
    return m_param_array.at(G4RandFlat::shootInt(m_n_param));
  }

  const BeamInfo &
    BeamMan::Get(G4int iev) const
  {
    auto b = m_param_array.at(iev);
    return m_param_array.at(iev % m_param_array.size());
  }

  const MMVertex &
    BeamMan::GetVertex(void) const
  {
    int nev = m_mm_array.size();
    return m_mm_array.at(G4RandFlat::shootInt(nev));
  }

  const MMVertex &
    BeamMan::GetVertex(G4int iev) const
  {
    int nev = m_mm_array.size();
    return m_mm_array.at(iev % nev);
  }

  const KmKpL
    BeamMan::GetKmKpL(void) const
  {
    int nev = m_kmkpl_array.size();
    int ev = G4RandFlat::shootInt(nev);
    auto l = m_kmkpl_array.at(ev);
    return l;
  }

  //_____________________________________________________________________________
  const HMCXi&
    BeamMan::GetHMCXi(G4int iev,G4int nmulti) const
  {
	G4int ev = iev / nmulti;
	int nev = m_hmcxi_array.size();
	while(ev >= nev){
		ev -= nev;
	}
	return m_hmcxi_array.at(ev);
  }
  void BeamMan::Print(void) const
  {
    PrintHelper helper(4, std::ios::fixed, G4cout);
    const G4int w = 8;

    G4cout << FUNC_NAME << G4endl;
    for (const auto &b : m_param_array)
      {
	G4cout << "   "
	       << "x=" << std::setw(w) << b.x << " "
	       << "y=" << std::setw(w) << b.y << " "
	       << "u=" << std::setw(w) << b.u << " "
	       << "v=" << std::setw(w) << b.v << " "
	       << "p=(" << std::setw(w) << b.p.x() << ", "
	       << std::setw(w) << b.p.y() << ", "
	       << std::setw(w) << b.p.z() << ")" << G4endl;
      }
    G4cout << "   nparam = " << m_param_array.size() << G4endl;
  }
  void BeamMan::GetHitProfile(double &x, double &y) const
  {
    HitProfile->GetRandom2(x, y);
  }

  const BeamInfo &
    BeamMan::GetAccidental(void) const
  {
    return m_accidental_array.at(G4RandFlat::shootInt(m_n_accidental));
  }
  
	const BeamInfo &
    BeamMan::GetFermiMom(void) const
  {
    return m_fermi_array.at(G4RandFlat::shootInt(m_n_fermi));
  }
	const BeamInfo &
    BeamMan::GetWeightedFermiMom(void) const
  {
	G4double r = G4RandFlat::shoot(0., m_weight_fermi);
	int i = 0;  
	for( int iev = 0; iev < m_n_fermi; ++iev){
	  if(r < m_fermi_array[iev].weight_cum){
	    i = iev;
	    break;
	  }
	}
	return m_fermi_array.at(i);
  }
  void BeamMan::LoadPolarizationTables(){
	G4cout<<"sqrts_data size: "<<sqrts_0s.size()<<G4endl;
	G4cout<<"cth_data size: "<<cth_0s.size()<<G4endl;
	for(int ix=0;ix<sqrts_0s.size();ix++){
		double x0,x1;
		x0 = sqrts_0s[ix];
		x1 = sqrts_1s[ix];
		double x_tmp = x0;
		if(x0 > x1){
			x0 = x1;
			x1 = x_tmp;
		}
		bool dupl = 0;
		for(int jx=0;jx<sorted_sqrts_0.size();jx++){
			if(abs(x0 - sorted_sqrts_0[jx])<1e-6 and abs(x1 - sorted_sqrts_1[jx])<1e-6){
				dupl = 1;
				break;
			}
		}
		if(!dupl){
			sorted_sqrts_0.push_back(x0);
			sorted_sqrts_1.push_back(x1);
		}
	}
	G4cout<<"Number of sqrts bins: "<<sorted_sqrts_0.size()<<G4endl;
	for(int is=0;is<sorted_sqrts_0.size();is++){
		G4cout<<std::setw(8)<<sorted_sqrts_0[is]<<"-"<<std::setw(8)<<sorted_sqrts_1[is]<<G4endl;
	}
	nbin_sqrts = sorted_sqrts_0.size() + 2;// for underflow and overflow
	for(int ic=0;ic<cth_0s.size();ic++){
		double cth0,cth1;
		cth0 = cth_0s[ic];
		cth1 = cth_1s[ic];
		double cth_tmp = cth0;
		if(cth0 > cth1){
			cth0 = cth1;
			cth1 = cth_tmp;
		}
		bool dupl = 0;
		for(int jc=0;jc<sorted_cth_0.size();jc++){
			if(abs(cth0 - sorted_cth_0[jc])<1e-6 and abs(cth1 - sorted_cth_1[jc])<1e-6){
				dupl = 1;
				break;
			}
		}
		if(!dupl){
			sorted_cth_0.push_back(cth0);
			sorted_cth_1.push_back(cth1);
		}
	}
	nbin_cth = sorted_cth_0.size() + 2;// for underflow and overflow
	G4cout<<"Number of cth bins: "<<sorted_cth_0.size()<<G4endl;
	for(int ic=0;ic<sorted_cth_0.size();ic++){
		G4cout<<std::setw(8)<<sorted_cth_0[ic]<<"-"<<std::setw(8)<<sorted_cth_1[ic]<<G4endl;
	}
	for(int is=0;is<nbin_sqrts;++is){
		vector<double> pol_cth;
		for(int ic=0;ic<nbin_cth;++ic){
			pol_cth.push_back(0.0);
		}
		pol_table.push_back(pol_cth);
	}
	for(int iev = 0; iev < pols.size(); ++iev){
		int cth_bin = GetPolaCThBin((cth_0s[iev] + cth_1s[iev])/2);
		int sqrt_bin = GetPolaSqrtsBin((sqrts_0s[iev] + sqrts_1s[iev])/2);
		pol_table[cth_bin][sqrt_bin] = pols[iev];
	}
	//Fill underflow and overflow bins
	for(int is=1;is<nbin_sqrts-1;++is){
		pol_table[0][is] = pol_table[1][is];
		pol_table[nbin_cth-1][is] = pol_table[nbin_cth-2][is];
	}
	for(int ic=0;ic<nbin_cth;++ic){
		pol_table[ic][0] = pol_table[ic][1];
		pol_table[ic][nbin_sqrts-1] = pol_table[ic][nbin_sqrts-2];
	}
  }

	int BeamMan::GetPolaCThBin(double cth){
		int ibin = -1;
		for(int i=0;i<nbin_cth-2;++i){
			if(cth >= sorted_cth_0[i] and cth < sorted_cth_1[i]){
				ibin = i+1;
				break;
			}
		}
		if(ibin == -1){
			if(cth < sorted_cth_0[0]) ibin = 0;
			else if(cth >= sorted_cth_1[nbin_cth-3]) ibin = nbin_cth -1;
		}
		return ibin;
	}
	int BeamMan::GetPolaSqrtsBin(double sqrt){
		int ibin = -1;
		for(int i=0;i<nbin_sqrts-2;++i){
			if(sqrt >= sorted_sqrts_0[i] and sqrt < sorted_sqrts_1[i]){
				ibin = i+1;
				break;
			}
		}
		if(ibin == -1){
			if(sqrt < sorted_sqrts_0[0]) ibin = 0;
			else if(sqrt >= sorted_sqrts_1[nbin_sqrts-3]) ibin = nbin_sqrts -1;
		}
		return ibin;
	}
	double BeamMan::GetPolaCth(int bin){
		if(bin == 0) return sorted_cth_0[0];
		if(bin == nbin_cth-1) return sorted_cth_1[nbin_cth-3];
		return (sorted_cth_0[bin-1] + sorted_cth_1[bin-1])/2;
	}
	double BeamMan::GetPolaSqrts(int bin){
		if(bin == 0) return sorted_sqrts_0[0];
		if(bin == nbin_sqrts-1) return sorted_sqrts_1[nbin_sqrts-3];
		return (sorted_sqrts_0[bin-1] + sorted_sqrts_1[bin-1])/2;
	}
	void BeamMan::PrintPolarizationTables(){
		G4cout<<"Xi Polarization Table:"<<G4endl;
		for(int is=0;is<nbin_sqrts-2;++is){
			G4cout<<sorted_sqrts_0[is]<<"-"<<sorted_sqrts_1[is]<<"         ";
		}
		G4cout<<G4endl;
		for(int ic=0;ic<nbin_cth-2;++ic){
			G4cout<<sorted_cth_0[ic]<<"-"<<sorted_cth_1[ic]<<" ";
			for(int is=0;is<nbin_sqrts-2;++is){
				G4cout<<std::setw(8)<<pol_table[ic+1][is+1]<<" ";
			}
			G4cout<<G4endl;
		}
	}
	G4double BeamMan::GetXiPolarization(double cth,double sqrt){
		int cth_bin = GetPolaCThBin(cth);
		int sqrt_bin = GetPolaSqrtsBin(sqrt);
		double cth0 = GetPolaCth(cth_bin);
		double sqrt0 = GetPolaSqrts(sqrt_bin);
		double p_cth_left,p_cth_right;
		double p_sqrt_left,p_sqrt_right;
		double cth0_left,cth0_right;
		double sqrt0_left,sqrt0_right;

		return pol_table[cth_bin][sqrt_bin];
	}
