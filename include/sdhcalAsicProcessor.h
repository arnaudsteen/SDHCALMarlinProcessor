#ifndef sdhcalAsicProcessor_h
#define sdhcalAsicProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <EVENT/CalorimeterHit.h>
#include <vector>
#include <map>
#include <limits>

#include "CaloObject/CaloGeom.h"
#include "CaloObject/CaloHit.h"
#include "CaloObject/Asic.h"
#include "Algorithm/Cluster.h"
#include "Algorithm/Tracking.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/InteractionFinder.h"
#include "Algorithm/Efficiency.h"
#include "Algorithm/AsicKeyFinder.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>

using namespace lcio ;
using namespace marlin ;

const int sdhcal_asic_table[48]={
  4,3,2,1,
  5,6,7,8,
  12,11,10,9,
  13,14,15,16,
  20,19,18,17,
  21,22,23,24,
  28,27,26,25,
  29,30,31,32,
  36,35,34,33,
  37,38,39,40,
  44,43,42,41,
  45,46,47,48
};

class sdhcalAsicProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new sdhcalAsicProcessor ; }
  
  
  sdhcalAsicProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

  void AlgorithmRegistrationParameters(); 
  void LayerProperties(std::vector<caloobject::CaloCluster2D*> &clusters);
  void clearVec();
  void DoTracking();
  inline int findDifID(int key){return _difList.at(key/1000*3+2-key%1000%12/4);}
  inline int findAsicID(int key){return sdhcal_asic_table[4*(key%1000/12) + key%1000%12%4];}

 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

 private:
  std::map<int,std::vector<caloobject::CaloHit*> > hitMap;
  
  /*--------------------Global parameters--------------------*/
  int _nActiveLayers;
  int numElements;
  LCCollection * col;
  int _nAsicX;
  int _nAsicY;
  std::vector<int> _difList;
  std::vector<float> edges; //vector to recover geometry parameters
  CLHEP::Hep3Vector posShift;
  /*------------------------------------------------------------------------------*/

  /*--------------------Algorithms list to initialise--------------------*/
  algorithm::Cluster *algo_Cluster;
  algorithm::ClusteringHelper *algo_ClusteringHelper;
  algorithm::Tracking *algo_Tracking;
  algorithm::InteractionFinder *algo_InteractionFinder;
  algorithm::Efficiency *algo_Efficiency;
  algorithm::AsicKeyFinder *algo_AsicKeyFinder;
  /*------------------------------------------------------------------------------*/
  
  /*--------------------Algorithms setting parameter structure--------------------*/
  algorithm::clusterParameterSetting m_ClusterParameterSetting; 
  algorithm::ClusteringHelperParameterSetting m_ClusteringHelperParameterSetting; 
  algorithm::TrackingParameterSetting m_TrackingParameterSetting; 
  algorithm::InteractionFinderParameterSetting m_InteractionFinderParameterSetting; 
  algorithm::EfficiencyParameterSetting m_EfficiencyParameterSetting; 
  algorithm::AsicKeyFinderParameterSetting m_AsicKeyFinderParameterSetting; 
  /*------------------------------------------------------------------------------*/
  
  /*--------------------CaloObject setting parameter structure--------------------*/
  caloobject::GeomParameterSetting m_CaloGeomSetting;
  /*------------------------------------------------------------------------------*/

  /*--------------------CaloObject list to initialise--------------------*/
  std::vector<caloobject::CaloLayer*> layers;
  std::map<int,caloobject::SDHCAL_Asic*> asics;
  /*---------------------------------------------------------------------*/
  
  /*--------------------Root output object--------------------*/
  std::string outputRootName;
  TFile *file; 
  TTree* tree; 

  TH1D *ntrack;
  TH1D *effGlobal1;
  TH1D *effGlobal2;
  TH1D *effGlobal3;
  TH1D *mulGlobal;
  TH1D *keyList;
  
  TH2D *mul2D;
  TH2D *eff2D_thr1;
  TH2D *eff2D_thr2;
  TH2D *eff2D_thr3;
  TH2D *trackPosition;

  int _layerID;
  int _difID;
  int _asicID;
  float _efficiency1;
  float _efficiency2;
  float _efficiency3;
  float _efficiency1_error;
  float _efficiency2_error;
  float _efficiency3_error;
  float _multiplicity;
  float _multiplicity_error;
  int _ntrack;


} ;

#endif
