#ifndef BACDetectorConstruction_110_3_hh
#define BACDetectorConstruction_110_3_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4PVPlacement;
class G4VisAttributes;


class BACDetectorConstruction_110_3 : public G4VUserDetectorConstruction
{
public:
  BACDetectorConstruction_110_3();
  virtual ~BACDetectorConstruction_110_3();

  virtual void ConstructSDandField();

  virtual G4VPhysicalVolume* Construct();

private:


  
  G4VPhysicalVolume* physWorld;
  G4VPhysicalVolume* physDetect;
  G4LogicalVolume* Aero1LW;
  G4LogicalVolume* Aero2LW;
  G4LogicalVolume* Aero3LW;
  G4LogicalVolume* MPPCLW;
  G4LogicalVolume* FrameLW;
  


  G4LogicalVolume* ReflectLW;
  G4LogicalVolume* Reflect1LW;


  std::vector<G4VisAttributes*> fVisAttributes;


  
};

#endif
