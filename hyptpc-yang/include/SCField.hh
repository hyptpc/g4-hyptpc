/*
  SCField.hh
*/

#ifndef SCField_h 
#define SCField_h 1

#include "SCFieldMap.hh"
#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include <vector>

class SimpleFieldElement;
//class G4ThreeVector;
//class G4RotationMatrix;

class SCField : public G4MagneticField
{
public:
  SCField();
  explicit SCField( const std::string &FieldMapName, double scaleFactor=1.0 );
  ~SCField();

private:
  SCField( const SCField & );
  SCField & operator = ( const SCField & );

public:
  void GetFieldValue( const double Point[4], double *Bfield ) const;
  void cleanupSimpleElementList( void );
  void AddSimpleElement( SimpleFieldElement *elem );
  void AddFieldMap( const std::string &FieldMapName, 
		    double scaleFactor=1.0 );
  void AddPosRot(G4ThreeVector pos, G4RotationMatrix rot);

private:
  //SCFieldMap fMap;

  typedef std::vector <SimpleFieldElement *> FMContainer;
  typedef std::vector <SimpleFieldElement *>
  ::const_iterator FMIterator;
  
  FMContainer elemList_;

  typedef std::vector <SCFieldMap *> fMapContainer;
  fMapContainer fMapCont;
  typedef std::vector <G4ThreeVector> TVVec; 
  typedef std::vector <G4RotationMatrix> RMVec;
  TVVec TVCont;
  RMVec RMCont;
  
};

#endif

