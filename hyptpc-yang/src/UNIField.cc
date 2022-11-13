/*
  UNIField.cc
*/

#include "UNIField.hh"
#include "SimpleFieldElement.hh"
#include "G4SystemOfUnits.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

UNIField::UNIField( const std::string &FieldMapName, double scaleFactor1, double scaleFactor2 )
//  : fMap( FieldMapName.c_str(),scaleFactor )
{
  scaleFactor1_ = scaleFactor1;
  scaleFactor2_ = scaleFactor2;
}

UNIField::~UNIField()
{
}

void UNIField::AddPosRot(int name, G4ThreeVector pos, G4RotationMatrix rot)
{
  if( name == 1)
    {
      pos1_ = pos;
      rot1_ = rot;
    }
  else if( name == 2)
    {
      pos2_ = pos;
      rot2_ = rot;
    }

}

void UNIField::GetFieldValue( const double Point[4], 
			      double *Bfield ) const
{
  double X[3];
  X[0]=Point[0]/mm; X[1]=Point[1]/mm; X[2]=Point[2]/mm;
  G4ThreeVector pos(X[0], X[1], X[2]);
  G4RotationMatrix inverse_rotation_matrix = rot1_.inverse(); // rot1_ = rot2_ 
  G4ThreeVector inverse_rotation_pos = inverse_rotation_matrix*pos; //inverse rotation

  //KURAMA
  double kurama_p[3];
  //kurama_p[0] = 150.;
  //kurama_p[1] = 0.;
  //kurama_p[2] = 1620.;

  kurama_p[0] = pos2_.x(); 
  kurama_p[1] = pos2_.y(); 
  kurama_p[2] = pos2_.z(); 

  //SC
  double sc_p[3];

  sc_p[0] = pos1_.x(); 
  sc_p[1] = pos1_.y(); 
  sc_p[2] = pos1_.z(); 


  if( fabs(inverse_rotation_pos.x() - sc_p[0]) < 310 && fabs(inverse_rotation_pos.y() - sc_p[1]) < 300 && fabs(inverse_rotation_pos.z() - sc_p[2]) <310 )
    {
      Bfield[0] = 0.;
      Bfield[1] = -1. * scaleFactor1_;
      Bfield[2] = 0.;
    }

  else if( fabs(inverse_rotation_pos.x() - kurama_p[0])<500. && fabs(inverse_rotation_pos.y()-kurama_p[1])< 400. && fabs(inverse_rotation_pos.z()-kurama_p[2]) < 400.)
    {
      Bfield[0] = 0.;
      Bfield[1] = -0.7 * scaleFactor2_;
      Bfield[2] = 0.;
    }
  else
    {
      Bfield[0] = 0.;
      Bfield[1] = 0.;
      Bfield[2] = 0.;
    }

  Bfield[0] *= tesla;
  Bfield[1] *= tesla;
  Bfield[2] *= tesla;

  //G4cout<<"[UNIField]BField: "<<Bfield[0]<<" "<<Bfield[1]<<" "<<Bfield[2]<<G4endl;                                                                       
  if( X[2] > -100000. && X[2] < 100000.0 )
    {
      //G4cout<<"[UNIField]test"<<G4endl;
      //G4cout<<"[UNIField]position: "<<X[0]<<" "<<X[1]<<" "<<X[2]<<G4endl;
      //G4cout<<"[UNIField]position1: "<<pointx1.x()/mm <<" "<< pointx1.y()/mm<<" "<<pointx1.z()/mm<<G4endl;
      //G4cout<<"[UNIField]position2: "<<pointx2.x()/mm <<" "<< pointx2.y()/mm<<" "<<pointx2.z()/mm<<G4endl;
      //G4cout<<"[UNIField]bfield: "<<Bfield[0]<<" "<<Bfield[1]<<" "<<Bfield[2]<<G4endl;
    }                                                                                                                                                     

#if 0
  if(Bfield[1]/tesla > -0.5 && Bfield[1]/tesla < -0.4)
    {
      G4cout << "X=(" << X[0] << "," << X[1] << "," << X[2] << ") "
	     << "B=(" << Bfield[0]/tesla << "," << Bfield[1]/tesla
	     << "," << Bfield[2]/tesla << ")" << G4endl;
    }
#endif
}


void UNIField::cleanupSimpleElementList( void )
{
}

void UNIField::AddSimpleElement( SimpleFieldElement *elem )
{
}
