// -*- C++ -*-

#ifndef VIS_MANAGER_HH
#define VIS_MANAGER_HH

#include <G4VisManager.hh>

//_____________________________________________________________________________
class VisManager : public G4VisManager
{
public:
  VisManager();
  virtual ~VisManager();

private:
  virtual void RegisterGraphicsSystems();
};

#endif
