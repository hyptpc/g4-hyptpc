#ifndef BACDetectorConstruction_167_3_hh
#define BACDetectorConstruction_167_3_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4PVPlacement;
class G4VisAttributes;


class BACDetectorConstruction_167_3 : public G4VUserDetectorConstruction
{
public:
  BACDetectorConstruction_167_3();
  virtual ~BACDetectorConstruction_167_3();

  virtual void ConstructSDandField();

  virtual G4VPhysicalVolume* Construct();

private:


  
  G4VPhysicalVolume* physWorld;
  G4VPhysicalVolume* physDetect;
  G4LogicalVolume* AeroLW;
  G4LogicalVolume* Aero_airLW;

  G4LogicalVolume* MPPCLW;
  G4LogicalVolume* ReflectLW;



  std::vector<G4VisAttributes*> fVisAttributes;


  
};

#endif
