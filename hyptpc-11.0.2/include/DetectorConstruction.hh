// -*- C++ -*-

#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <map>

#include <G4VUserDetectorConstruction.hh>
#include <G4RotationMatrix.hh>
#include <G4String.hh>

class G4Element;
class G4Material;
class G4LogicalVolume;
class G4PVPlacement;
class SDCSD;
class MagneticField;

//_____________________________________________________________________________
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  static G4String ClassName();
  DetectorConstruction();
  ~DetectorConstruction();

private:
  G4int                           m_experiment;
  std::map<G4String, G4Element*>  m_element_map;
  std::map<G4String, G4Material*> m_material_map;
  G4LogicalVolume*                m_world_lv;
  G4LogicalVolume*                m_tpc_lv;
  G4double                        m_rotation_angle;
  G4RotationMatrix*               m_rotation_matrix;
  SDCSD*                          m_sdc_sd;
  MagneticField*                  m_field;

private:
  // Materials
  void ConstructElements();
  void ConstructMaterials();
  // Detectors
  void ConstructBAC();
  void ConstructBH2();
  void ConstructFTOF();
  void ConstructHTOF();
  void ConstructHypTPC();
  void ConstructKVC();
  void ConstructShsMagnet();
  void ConstructTarget();
  // For K1.8
  void ConstructAreaTent();
  void ConstructBC3();
  void ConstructBC4();
  void ConstructK18BeamlineSpectrometer();
  void ConstructKuramaMagnet();
  void ConstructLAC();
  void ConstructSCH();
  void ConstructSDC1();
  void ConstructSDC2();
  void ConstructSDC3();
  void ConstructSDC4();
  void ConstructWC();

private:
  virtual G4VPhysicalVolume* Construct();
};

//_____________________________________________________________________________
inline G4String
DetectorConstruction::ClassName()
{
  static G4String s_name("DetectorConstruction");
  return s_name;
}

#endif
