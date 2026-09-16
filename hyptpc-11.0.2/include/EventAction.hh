// -*- C++ -*-

#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH

#include <vector>

#include <G4String.hh>
#include <G4Types.hh>
#include <G4UserEventAction.hh>

class G4Event;
class G4RunAction;

//_____________________________________________________________________________
class EventAction : public G4UserEventAction
{
public:
  static G4String ClassName();
  EventAction();
  virtual ~EventAction();
  virtual void BeginOfEventAction(const G4Event* anEvent);
  virtual void EndOfEventAction(const G4Event* anEvent);
private:
  G4double CalculateAverage(const std::vector<G4double>& pos);
  G4double CalculateSum(const std::vector<G4double>& edep);
};

//_____________________________________________________________________________
inline G4String
EventAction::ClassName()
{
  static G4String s_name("EventAction");
  return s_name;
}

#endif
