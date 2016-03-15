#ifndef sdhcalEfficiencyProcessor_h
#define sdhcalEfficiencyProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <EVENT/CalorimeterHit.h>
#include <vector>
#include <map>

#include "CaloObject/CaloHit.h"
#include "Algorithm/Cluster.h"
#include "Algorithm/Tracking.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/InteractionFinder.h"
#include "Algorithm/Efficiency.h"

#include <TFile.h>
#include <TTree.h>

using namespace lcio ;
using namespace marlin ;

class sdhcalEfficiencyProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new sdhcalEfficiencyProcessor ; }
  
  
  sdhcalEfficiencyProcessor() ;
  
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
  void LayerProperties(std::vector<caloobject::CaloCluster*> &clusters);
  void clearVec();
  void DoTracking();
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
  /*------------------------------------------------------------------------------*/

  /*--------------------Algorithms list to initialise--------------------*/
  algorithm::Cluster *algo_Cluster;
  algorithm::ClusteringHelper *algo_ClusteringHelper;
  algorithm::Tracking *algo_Tracking;
  algorithm::InteractionFinder *algo_InteractionFinder;
  algorithm::Efficiency *algo_Efficiency;
  /*------------------------------------------------------------------------------*/
  
  /*--------------------Algorithms setting parameter structure--------------------*/
   algorithm::clusterParameterSetting m_ClusterParameterSetting; 
   algorithm::ClusteringHelperParameterSetting m_ClusteringHelperParameterSetting; 
   algorithm::TrackingParameterSetting m_TrackingParameterSetting; 
   algorithm::InteractionFinderParameterSetting m_InteractionFinderParameterSetting; 
   algorithm::EfficiencyParameterSetting m_EfficiencyParameterSetting; 
  /*------------------------------------------------------------------------------*/
  
  /*--------------------CaloObject setting parameter structure--------------------*/
   caloobject::LayerParameterSetting m_LayerParameterSetting;
  /*------------------------------------------------------------------------------*/

  /*--------------------CaloObject list to initialise--------------------*/
   std::vector<caloobject::CaloLayer*> layers;
  /*---------------------------------------------------------------------*/
  
  /*--------------------Root output object--------------------*/
  std::string outputRootName;
   TFile *file; 
   TTree* tree; 

   float _eventChi2; 
   float _transverseRatio; 
   float _trackCosTheta; 
   float _trackX0; 
   float _trackY0; 
   float _trackZ0; 

   std::vector<double> _efficiency; 
   std::vector<double> _effEnergy; 
   std::vector<double> _multiplicity; 
   std::vector<double> _chi2; 

} ;

#endif
