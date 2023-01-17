#ifndef PMT15SD_h
#define PMT15SD_h 1

#include "G4VSensitiveDetector.hh"
#include "PMT15Hit.hh"
#include "TGraph.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class PMT15SD : public G4VSensitiveDetector
{
public:
  PMT15SD(G4String name);
  //PMT15SD(G4String name, ParamMan*);
  virtual ~PMT15SD();

  virtual void Initialize( G4HCofThisEvent *HCTE );
  virtual G4bool ProcessHits(G4Step* astep, G4TouchableHistory *ROhist);
  virtual void EndOfEvent(G4HCofThisEvent *HCTE);

  void SetQETable();
  TGraph* QETable;
  

private:
  PMT15HitsCollection *Pmt15Collection;

};

#endif
