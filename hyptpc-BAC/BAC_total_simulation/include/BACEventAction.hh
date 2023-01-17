#ifndef BACEVENTACTION_HH
#define BACEventAction_hh 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "TTree.h"
#include "TFile.h"

#include <vector>
#include <array>


class BACAnalysisManager;


class BACEventAction : public G4UserEventAction
{
public:
  BACEventAction(BACAnalysisManager *analysisManager=0);
  virtual ~BACEventAction();
private:
  BACEventAction(const BACEventAction &);
  BACEventAction & operator = (const BACEventAction &);

public:
  
  virtual void BeginOfEventAction(const G4Event *anEvent);
  virtual void EndOfEventAction(const G4Event *anEvent);

protected:
  BACAnalysisManager *anaMan;



  
};

#endif

 
  
