#include "G4LUTTree.h"

#include "PrtHit.h"
#include "PrtLutNode.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>  // for PHG4HitContainer, PHG4Hit...
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phparameter/PHParameters.h>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <cmath>     // for atan2, sqrt
#include <cstring>   // for strcmp
#include <iostream>  // for ostringstream, operator<<
#include <sstream>
#include <utility>  // for pair

G4LUTTree::G4LUTTree(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , _filename(filename)
  , outfile(nullptr)
{
}

int G4LUTTree::Init(PHCompositeNode *)
{
  outfile = new TFile(_filename.c_str(), "RECREATE");
  
  //-------- LUT --------
  
  fLut = new TClonesArray("PrtLutNode");
  fLutTree = new TTree("prtlut","Look-up table for the geometrical reconstruction.");
  fLutTree->Branch("LUT",&fLut,256000,2); 
  Int_t Nnodes = 100000;
    
  TClonesArray &fLuta = *fLut; 
  for (Long64_t n=0; n<Nnodes; n++) 
    {
      new((fLuta)[n]) PrtLutNode(n);
    }    

  return 0;
}

int G4LUTTree::process_event(PHCompositeNode *topNode)
{
  // get the primary particle which did this to us....
  PHG4TruthInfoContainer *truthInfoList = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  const PHG4TruthInfoContainer::Range primRange = truthInfoList->GetPrimaryParticleRange();
  double px = primRange.first->second->get_px();
  double py = primRange.first->second->get_py();
  double pz = primRange.first->second->get_pz();
  double p = sqrt(px * px + py * py + pz * pz);
  double pt = sqrt(px * px + py * py);
  double phi = atan2(py, px);
  double theta = atan2(pt, pz);
  double pid = primRange.first->second->get_pid();
  TVector3 dir_vec(px, py, pz);
  
  int nhits = 0;

  std::ostringstream nodename;
  std::set<std::string>::const_iterator iter;

  for (iter = _node_postfix.begin(); iter != _node_postfix.end(); ++iter)
  {
    int detid = (_detid.find(*iter))->second;
    nodename.str("");
    nodename << "G4HIT_" << *iter;
    PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());

    if (!strcmp("G4HIT_DIRC", nodename.str().c_str()))  // DIRC
    {
      process_hit(hits, "G4HIT_DIRC", detid, nhits, dir_vec);
    }
  }

  evt_num++;

  return 0;
}

int G4LUTTree::End(PHCompositeNode *topNode)
{
  outfile->cd();
  fLutTree->Fill();
  fLutTree->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  return 0;
}

void G4LUTTree::AddNode(const std::string &name, const int detid)
{
  _node_postfix.insert(name);
  _detid[name] = detid;
  return;
}

int G4LUTTree::process_hit(PHG4HitContainer *hits, const std::string &dName, int detid, int &nhits, TVector3 dir_vec)
{
  if (hits)
  {
    PHG4HitContainer::ConstRange hit_range_0 = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter_0 = hit_range_0.first; hit_iter_0 != hit_range_0.second; hit_iter_0++)
    {
      PrtHit *dirc_hit = dynamic_cast<PrtHit *>(hit_iter_0->second);
      
      int id = 256*dirc_hit->GetMcpId() + dirc_hit->GetPixelId();
      ((PrtLutNode*)(fLut->At(id)))->
	AddEntry(id, dir_vec, dirc_hit->GetPathInPrizm(),
		 dirc_hit->GetNreflectionsInPrizm(),
		 dirc_hit->GetLeadTime(), dirc_hit->GetGlobalPos(), dirc_hit->GetDigiPos());	

      nhits++;
    }

  }

  return 0;
}
