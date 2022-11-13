/*
  UNIField.hh
*/

#ifndef UNIField_h 
#define UNIField_h 1

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "SimpleFieldElement.hh"

#include <vector>


class UNIField : public G4MagneticField
{
public:
  UNIField();
  explicit UNIField( const std::string &FieldMapName, double scaleFactor1=1.0, double scaleFactor2=1.0);
  ~UNIField();

private:
  UNIField( const UNIField & );
  UNIField & operator = ( const UNIField & );

  double scaleFactor1_;
  double scaleFactor2_;
  G4ThreeVector pos1_;
  G4ThreeVector pos2_;
  G4RotationMatrix rot1_;
  G4RotationMatrix rot2_;

public:
  void GetFieldValue( const double Point[4], double *Bfield ) const;
  void cleanupSimpleElementList( void );
  void AddSimpleElement( SimpleFieldElement *elem );
  void AddPosRot(int name, G4ThreeVector pos, G4RotationMatrix rot);

private:
  
};

#endif

