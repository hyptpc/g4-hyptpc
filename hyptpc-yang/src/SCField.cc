/*
  SCField.cc
*/

#include "SCField.hh"
#include "SimpleFieldElement.hh"
#include "G4SystemOfUnits.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

SCField::SCField( const std::string &FieldMapName, double scaleFactor )
//  : fMap( FieldMapName.c_str(),scaleFactor )
{
  fMapCont.clear();
  TVCont.clear();
  RMCont.clear();

  SCFieldMap *fMap = new SCFieldMap( FieldMapName.c_str(),scaleFactor );
  fMap->Initialize();
  fMapCont.push_back(fMap);
  //TVCont.push_back(pos);
  //RMCont.push_back(rot);
}

void SCField::AddFieldMap(const std::string &FieldMapName, double scaleFactor)
{
  G4cout<<"Adding one more magnetic field map <KURAMA MAP>"<<G4endl;
  SCFieldMap *fMap = new SCFieldMap( FieldMapName.c_str(),scaleFactor );
  fMap->Initialize();
  fMapCont.push_back(fMap);
  //TVCont.push_back(pos);
  //RMCont.push_back(rot);
}

void SCField::AddPosRot(G4ThreeVector pos, G4RotationMatrix rot)
{
  TVCont.push_back(pos);
  RMCont.push_back(rot);
}

SCField::~SCField()
{
}

void SCField::GetFieldValue( const double Point[4], 
			      double *Bfield ) const
{
  double X[3];
  X[0]=Point[0]/mm; X[1]=Point[1]/mm; X[2]=Point[2]/mm;
  G4ThreeVector pos(X[0], X[1], X[2]);

  int fieldmap_size = fMapCont.size();
  int pos_size = TVCont.size();
  int rotation_matrix_size = RMCont.size();
  if( fieldmap_size != pos_size || pos_size != rotation_matrix_size || fieldmap_size != rotation_matrix_size )
    {
      G4cout<<"[SCField] Field Maps not matched"<<G4endl;
      //break;
    }

  //G4cout<<"[SCField]test: "<< elemList_.size()<<G4endl;

  double bfield[3]={0.};
  //double bfield1[3]={0.};
  //double bfield2[3]={0.};
  //G4ThreeVector pointx1, pointx2;

  for( int i=0 ; i < fieldmap_size; i++)
    {
      SCFieldMap *fMap = fMapCont[i];
      G4ThreeVector pos_shift = TVCont[i];
      G4RotationMatrix rotation_matrix = RMCont[i];
      G4RotationMatrix inverse_rotation_matrix = rotation_matrix.inverse(); //inverse rotation matrix
      G4ThreeVector pos_unrotate = inverse_rotation_matrix*pos; //inverse rotation
      G4ThreeVector point_x(pos_unrotate.x() - pos_shift.x()/mm, pos_unrotate.y() - pos_shift.y()/mm, pos_unrotate.z()-pos_shift.z()/mm); //go to center position
      double XX[3];
      XX[0]=point_x.x();
      XX[1]=point_x.y();
      XX[2]=point_x.z();

      /*
      G4cout<<"[SCField]x: "<<X[0]<<" y: "<<X[1]<<" z: "<<X[2] << G4endl;
      G4cout<<"[SCField]pos: "<< pos << G4endl;
      G4cout<<"[SCField]after unrotation: "<< pos_unrotate << G4endl;
      G4cout<<"[SCField]after shift: "<< point_x << G4endl;
      */

      if( fMap->GetFieldValue( XX, Bfield ) )
	{
	  bfield[0] += Bfield[0];
	  bfield[1] += Bfield[1];
	  bfield[2] += Bfield[2];
	  /*
	  if(i==0)
	    {
	      bfield1[0] = Bfield[0]*tesla;
	      bfield1[1] = Bfield[1]*tesla;
	      bfield1[2] = Bfield[2]*tesla;
	      pointx1 = point_x;
	    }
	  if(i==1)
	    {
	      bfield2[0] = Bfield[0]*tesla;
	      bfield2[1] = Bfield[1]*tesla;
	      bfield2[2] = Bfield[2]*tesla;
	      pointx2 = point_x;
	    }
	  */
	}
      else
	{
	  G4cout<<"[SCField]I cannot find the field map"<<G4endl;
	  //bfield[0] = bfield[1] = bfield[2]  = 0.0;
	  continue;
	}
    }

  Bfield[0] = bfield[0];
  Bfield[1] = bfield[1];
  Bfield[2] = bfield[2];

  Bfield[0] *= tesla;
  Bfield[1] *= tesla;
  Bfield[2] *= tesla;

  
  /*
  if( fMapCont[0]->GetFieldValue( X, Bfield ) )
    {
      Bfield[0] *= tesla;
      Bfield[1] *= tesla;
      Bfield[2] *= tesla;
    }
  else
    {
      Bfield[0]=Bfield[1]=Bfield[2]=0.0;
    }
  */

  //G4cout<<"[SCField]BField: "<<Bfield[0]<<" "<<Bfield[1]<<" "<<Bfield[2]<<G4endl;
  /*
  if( fabs(X[0]-1000.0) < 100.1 && fabs(X[1]-0.0) < 20.0 )
    {
      G4cout<<"[SCField]test"<<G4endl;
      G4cout<<"[SCField]position: "<<X[0]<<" "<<X[1]<<" "<<X[2]<<G4endl;
      //G4cout<<"[SCField]position1: "<<pointx1.x()/mm <<" "<< pointx1.y()/mm<<" "<<pointx1.z()/mm<<G4endl;
      //G4cout<<"[SCField]position2: "<<pointx2.x()/mm <<" "<< pointx2.y()/mm<<" "<<pointx2.z()/mm<<G4endl;
      G4cout<<"[SCField]bfield: "<<bfield[0]<<" "<<bfield[1]<<" "<<bfield[2]<<G4endl;
      //G4cout<<"[SCField]bfield1: "<<bfield1[0]*1000.<<" "<<bfield1[1]*1000.<<" "<<bfield1[2]*1000.<<G4endl;
      //G4cout<<"[SCField]bfield2: "<<bfield2[0]*1000.<<" "<<bfield2[1]*1000.<<" "<<bfield2[2]*1000.<<G4endl;
    }
  */

  G4ThreeVector gPos( Point[0], Point[1], Point[2] );
  G4ThreeVector B( 0., 0., 0. );
  FMIterator end=elemList_.end();
  for( FMIterator itr=elemList_.begin(); itr!=end; ++itr )
    {
      if( (*itr)->ExistMagneticField() )
	B += (*itr)->GetMagneticField( gPos );
      std::cout<<"FMIterator!!!"<<std::endl;
    }
  
  Bfield[0] += B.x(); Bfield[1] += B.y(); Bfield[2] += B.z();

#if 0
  if(Bfield[1]/tesla > -0.5 && Bfield[1]/tesla < -0.4)
    {
      G4cout << "X=(" << X[0] << "," << X[1] << "," << X[2] << ") "
	     << "B=(" << Bfield[0]/tesla << "," << Bfield[1]/tesla
	     << "," << Bfield[2]/tesla << ")" << G4endl;
    }
#endif
}


void SCField::cleanupSimpleElementList( void )
{
  elemList_.clear();
}

void SCField::AddSimpleElement( SimpleFieldElement *elem )
{
  elemList_.push_back( elem );
}
