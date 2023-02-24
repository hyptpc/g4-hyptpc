// -*- C++ -*-

#ifndef BEAM_MAN_HH
#define BEAM_MAN_HH

#include <string>
#include <map>
#include <vector>

#include <globals.hh>
#include <G4ThreeVector.hh>

#include <Rtypes.h>

class TFile;

//_____________________________________________________________________________
struct BeamInfo
{
  G4double      x; // [mm]
  G4double      y; // [mm]
  G4double      u; // [mrad]
  G4double      v; // [mrad]
  G4double      dp; // [%]
  G4ThreeVector p; // [GeV/c]
  G4double      z; // [mm]
  G4double GetX(G4double offset=0.) const;
  G4double GetY(G4double offset=0.) const;
  void     Print() const;
};

//_____________________________________________________________________________
class BeamMan
{
public:
  static G4String ClassName();
  static BeamMan& GetInstance();
  ~BeamMan();

private:
  BeamMan();
  BeamMan(const BeamMan& );
  BeamMan& operator =(const BeamMan&);

private:
  typedef std::vector<BeamInfo> ParamArray;
  G4bool        m_is_ready;
  G4String      m_file_name;
  TFile*        m_file;
  ParamArray    m_param_array;
  G4int         m_n_param;
  G4bool        m_is_vi; // true:VI or false:VO
  G4double      m_primary_z; // from VI or VO
  G4ThreeVector m_vi_pos;

public:
  const BeamInfo&      Get() const;
  G4double             GetPrimaryZ() const { return m_primary_z; }
  const G4ThreeVector& GetVIPosition() const { return m_vi_pos; }
  G4bool               Initialize();
  G4bool               Initialize(const G4String& filename);
  G4bool               IsReady() const { return m_is_ready; }
  void                 Print() const;
  void                 SetPrimaryZ(G4double z){ m_primary_z = z; }
  void                 SetVIPosition(G4ThreeVector pos){ m_vi_pos = pos; }
};

//_____________________________________________________________________________
inline G4String
BeamMan::ClassName()
{
  static G4String s_name("BeamMan");
  return s_name;
}

//_____________________________________________________________________________
inline BeamMan&
BeamMan::GetInstance()
{
  static BeamMan s_instance;
  return s_instance;
}

#endif
