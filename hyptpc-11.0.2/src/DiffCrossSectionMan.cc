// -*- C++ -*-

#include "DiffCrossSectionMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TSpline.h>

//_____________________________________________________________________________
DiffCrossSectionMan::DiffCrossSectionMan()
  : m_is_ready(false),
    m_file_name(),
    m_file(),
    m_n_spline(),
    m_spline_container(),
    m_range_min(),
    m_range_max()
{
}

//_____________________________________________________________________________
DiffCrossSectionMan::~DiffCrossSectionMan()
{
}

//_____________________________________________________________________________
G4bool
DiffCrossSectionMan::Initialize()
{
  if (m_file_name.empty()) return true;

  TFile* file = TFile::Open(m_file_name.c_str(), "READ");
  if (!file || file->IsZombie() || !file->IsOpen()) {
    G4cerr << "[DCS] Failed to open file: " << m_file_name << G4endl;
    return false;
  }

  for (auto* p : m_spline_container) delete p;
  m_spline_container.clear();

  for (int i = 0; ; ++i) {
    G4String name = "A" + G4String(std::to_string(i)) + "_spline";

    TSpline3* sp = nullptr;
    file->GetObject(name.c_str(), sp);
    if (!sp) {
      if (i == 0) {
        G4cerr << "[DCS] No TSpline3 found (expected " << name << ") in "
               << m_file_name << G4endl;
        file->Close();
        return false;
      }
      break;
    }

    TSpline3* sp_clone = static_cast<TSpline3*>(sp->Clone());
    if (!sp_clone) {
      G4cerr << "[DCS] Clone failed for: " << name << G4endl;
      file->Close();
      return false;
    }
    m_spline_container.push_back(sp_clone);
  }

  if (m_spline_container.empty()) {
    G4cerr << "[DCS] No TSpline3 cloned from: " << m_file_name << G4endl;
    file->Close();
    return false;
  }

  m_range_min = m_spline_container.front()->GetXmin();
  m_range_max = m_spline_container.front()->GetXmax();

  file->Close();

  m_n_spline = static_cast<G4int>(m_spline_container.size());
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
DiffCrossSectionMan::Initialize(const G4String& filename)
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
G4double
DiffCrossSectionMan::GetCoeff(G4int order, G4double mom_kaon) const
{
  if (order < 0 || order >= m_n_spline || !m_spline_container[order]) {
    G4cerr << "Error: Invalid spline index: " << order << G4endl;
    return 0.0;
  }
  
  if (mom_kaon < m_range_min) {
    return m_spline_container[order]->Eval(m_range_min);
  } else if (mom_kaon > m_range_max) {
    return m_spline_container[order]->Eval(m_range_max);
  } else {
    return m_spline_container[order]->Eval(mom_kaon);
  }
}
