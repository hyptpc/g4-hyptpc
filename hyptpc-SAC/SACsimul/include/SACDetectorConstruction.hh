#ifndef SACDetectorConstruction_hh
#define SACDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4PVPlacement;
class G4VisAttributes;


class SACDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  SACDetectorConstruction();

  virtual ~SACDetectorConstruction();

  virtual void ConstructSDandField();

  virtual G4VPhysicalVolume* Construct();

private:

  G4VPhysicalVolume* physWorld;
  G4VPhysicalVolume* physDetect;
  G4LogicalVolume* AeroLW;
  G4LogicalVolume* FrameLW;
  G4LogicalVolume* PMT15LW;
  G4LogicalVolume* ReflectLW;
  G4LogicalVolume* Reflect1LW;


  std::vector<G4VisAttributes*> fVisAttributes;
  
};

#endif
