#ifndef G4LUTTREE_H
#define G4LUTTREE_H

#include "G4EventTree.h"

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <TVector3.h>
#include <TClonesArray.h>

// Forward declarations
class PHCompositeNode;
class PHG4HitContainer;
class TFile;
class TTree;
class PHParameters;

class G4LUTTree : public SubsysReco
{
 public:
  //! constructor
  G4LUTTree(const std::string &name = "G4LUTTree", const std::string &filename = "G4LUTTree.root");

  //! destructor
  ~G4LUTTree() override {}

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! hit processing method
  int process_hit(PHG4HitContainer *hits, const std::string &dName, int detid, int &nhits, TVector3 dir_vec);

  //! end of run method
  int End(PHCompositeNode *) override;

  void AddNode(const std::string &name, const int detid = 0);
  int evt_num = 0;

 protected:
  std::string _filename;
  std::set<std::string> _node_postfix;
  std::map<std::string, int> _detid;

  TTree *fLutTree;
  TFile *outfile;
  TClonesArray *fLut[10];

};

#endif
