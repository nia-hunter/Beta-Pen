#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include <vector>
#include "g4root.hh"
#include <fstream>
//can change to g4csv or g4xml

//====================================
class HistoManager
{
public:HistoManager();
  virtual ~HistoManager();

  void Book();
  void Save();

  void FillHisto(G4int id, G4double e, G4double weight = 1.0);
  void FillNtuple(G4double Edep,G4int Scint, G4int SiPM, G4double tot);
  void FillTimeAndLoc(std::vector<G4double> x, std::vector<G4double> y, std::vector<G4double>t1 ,std::vector<G4double> t2);
  
private:
  G4bool fFactoryOn;
    std::ofstream fout;
};
#endif
