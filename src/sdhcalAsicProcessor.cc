#include "sdhcalAsicProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
using namespace lcio ;
using namespace marlin ;

sdhcalAsicProcessor asdhcalAsicProcessor ;

sdhcalAsicProcessor::sdhcalAsicProcessor() : Processor("sdhcalAsicProcessor") {

  // modify processor description
  _description = "sdhcalAsicProcessor calculates a SDHCAL efficiency and pad multiplcity for each SHDCAL layer" ;
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("HCALBarrel"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "CollectionName" , 
			    "HCAL Collection Names"  ,
			    _hcalCollections  ,
			    hcalCollections);

  registerProcessorParameter( "RootFileName" ,
			      "File name for the root output",
			      outputRootName,
			      std::string("toto.root") ); 

  registerProcessorParameter( "NActiveLayers" ,
			      "Number of active layers",
			      _nActiveLayers,
			      int(48) ); 

  registerProcessorParameter( "N_ASIC" ,
			      "Number of ASIC per layer in x direction",
			      _nAsicX,
			      int(12) ); 

  registerProcessorParameter( "N_ASIC" ,
			      "Number of ASIC per layer in y direction",
			      _nAsicY,
			      int(12) ); 

  int difTab[]={ 181,94,30, 174,175,176, 158,142,141, 129,118,119, 164,152,151,  74,61,75,
		 156,111,110, 102,177,103,  133,136,134,  128,120,121,  65,64,58,  148,72,73,
		 78,79,60,  44,43,113,  243,242,241,   186,127,154,  147,70,71,   47,139,140,
		 143,77,76,   159,91,36,   179,178,183,  41,42,67,  137,46,138,  131,173,144,
		 189,184,160,  172,167,171,  146,135,145,  185,170,180,  187,188,190,  169,165,166,
		 155,57,50,  153,108,25,   51,56,109,   107,150,116,  126,124,49,  117,149,115,
		 48,45,114,   98,93,40,   92,97,100,  62,106,132,  101,35,99,  122,123,130,
		 163,161,162,  104,29,112,  59,53,54,  96,90,27,  95,8,5,  63,87,18 };
  std::vector<int> difVec(difTab, difTab + sizeof(difTab) / sizeof(int) );

  registerProcessorParameter( "DifList" ,
    			      "Vector of dif number",
    			      _difList,
    			      difVec ); 


  std::vector<float> vec;
  std::vector<float> _posShift;
  vec.push_back(499.584);
  vec.push_back(499.584);
  vec.push_back(0);
  // registerProcessorParameter( "PositionShift" ,
  // 			      "3 Vector to shift to have the right (0,0,0) position",
  // 			      _posShift,
  // 			      vec );
  posShift=CLHEP::Hep3Vector( 0.0,
			      0.0,
			      0.0 );

  AlgorithmRegistrationParameters();
}


void sdhcalAsicProcessor::AlgorithmRegistrationParameters()
{
  /*------------algorithm::Cluster------------*/
  registerProcessorParameter( "Cluster::MaxTransversalCellID" ,
 			      "Maximum difference between two hits cellID (0 and 1) to build a cluster",
 			      m_ClusterParameterSetting.maxTransversal,
 			      (int) 1 ); 

  registerProcessorParameter( "Cluster::MaxLongitudinalCellID" ,
 			      "Maximum difference between two hits cellID (2) to build a cluster",
 			      m_ClusterParameterSetting.maxLongitudinal,
 			      (int) 0 ); 

  registerProcessorParameter( "Cluster::UseDistanceInsteadCellID" ,
 			      "Boolean to know if clustering algorithms uses distance instead of cellID to cluster hits together",
 			      m_ClusterParameterSetting.useDistanceInsteadCellID,
 			      (bool) false ); 

  registerProcessorParameter( "Cluster::MaxTransversalDistance" ,
 			      "Maximum transversal distance (in mm) between two hits to gathered them in one cluster",
 			      m_ClusterParameterSetting.maxTransversalDistance,
 			      (float) 11.0 ); 

  registerProcessorParameter( "Cluster::MaxLongitudinalDistance" ,
 			      "Maximum longitudinal distance (in mm) between two hits to gathered them in one cluster",
 			      m_ClusterParameterSetting.maxLongitudinalDistance,
 			      (float) 27.0 ); 

  /*------------algorithm::ClusteringHelper------------*/
  registerProcessorParameter( "ClusteringHelper::LongitudinalDistanceForIsolation" ,
    			      "Minimum longitudinal distance (in mm) between one hits and its neighbours to decide if it is isolated",
    			      m_ClusteringHelperParameterSetting.longitudinalDistance,
    			      (float) 100.0 ); 

  registerProcessorParameter( "ClusteringHelper::TransversalDistanceDistanceForIsolation" ,
    			      "Minimum transversal distance (in mm) between one hits and its neighbours to decide if it is isolated",
    			      m_ClusteringHelperParameterSetting.transversalDistance,
    			      (float) 200.0 ); 

  /*------------algorithm::Tracking-----------*/
  registerProcessorParameter( "Tracking::ChiSquareLimit" ,
    			      "Maximum value of tracking fit chi2 to construct a track",
    			      m_TrackingParameterSetting.chiSquareLimit,
    			      (float) 100.0 ); 

  registerProcessorParameter( "Tracking::MaxTransverseRatio" ,
    			      "Maximum value of transverse ratio to construct a track",
    			      m_TrackingParameterSetting.maxTransverseRatio,
    			      (float) 0.05 ); 

  registerProcessorParameter( "Tracking::CosThetaLimit" ,
    			      "Minimum value of cos(Theta) to accept the track",
    			      m_TrackingParameterSetting.cosThetaLimit,
    			      (float) 0.0 ); 

  registerProcessorParameter( "Tracking::PrintDebug" ,
    			      "Boolean to know if debug if printed",
    			      m_TrackingParameterSetting.printDebug,
    			      (bool) false ); 

  /*------------algorithm::Efficiency-----------*/
  registerProcessorParameter( "Efficiency::MaxRadius" ,
    			      "Maximum distance parameter to find a hit to consider the layer as efficient",
    			      m_EfficiencyParameterSetting.maxRadius,
    			      (float) 25.0 ); 

  registerProcessorParameter( "Efficiency::SDHCALReadout" ,
    			      "Boolean to know if the detector used the semi digital readout",
    			      m_EfficiencyParameterSetting.semiDigitalReadout,
    			      (bool) true ); 

  registerProcessorParameter( "Efficiency::PrintDebug" ,
    			      "If true, Efficiency algorithm will print some debug information",
    			      m_EfficiencyParameterSetting.printDebug,
    			      (bool) false ); 

  m_EfficiencyParameterSetting.trackingParams=m_TrackingParameterSetting;

  /*------------algorithm::InteractionFinder-----------*/
  registerProcessorParameter( "InteractionFinder::MinSize" ,
    			      "Minimum cluster size for to define an interaction point",
    			      m_InteractionFinderParameterSetting.minSize,
    			      (int) 4 ); 

  registerProcessorParameter( "InteractionFinder::MaxRadius" ,
    			      "Maximum transversal distance to look for clusters",
    			      m_InteractionFinderParameterSetting.maxRadius,
    			      (float) 50.0 ); 

  registerProcessorParameter( "InteractionFinder::MaxDepth" ,
    			      "Maximum depth (number of layers) to look for clusters",
    			      m_InteractionFinderParameterSetting.maxDepth,
    			      (int) 4 ); 

  registerProcessorParameter( "InteractionFinder::MinNumberOfCluster" ,
    			      "Minimum number of found clusters (big enough) after the interaction point",
    			      m_InteractionFinderParameterSetting.minNumberOfCluster,
    			      (int) 3 ); 

  registerProcessorParameter( "InteractionFinder::UseAnalogEnergy" ,
    			      "Boolean to know if interaction finder algo should use cluster energy of cluster number of hits",
    			      m_InteractionFinderParameterSetting.useAnalogEnergy,
    			      (bool) false ); 

  registerProcessorParameter( "InteractionFinder::PrintDebug" ,
    			      "Boolean to know if debug if printed",
    			      m_InteractionFinderParameterSetting.printDebug,
    			      (bool) false ); 
  /*------------caloobject::CaloGeom------------*/
  registerProcessorParameter( "Geometry::NLayers" ,
 			      "Number of layers",
 			      m_CaloGeomSetting.nLayers,
 			      (int) 48 ); 
  registerProcessorParameter( "Geometry::NPixelsPerLayer" ,
 			      "Number of pixels per layer (assume square geometry)",
 			      m_CaloGeomSetting.nPixelsPerLayer,
 			      (int) 96 ); 
  registerProcessorParameter( "Geometry::PixelSize" ,
 			      "Pixel size (assume square pixels)",
 			      m_CaloGeomSetting.pixelSize,
 			      (float) 10.408 ); 

  std::vector<float> vec;
  vec.push_back(0.0);
  vec.push_back(1000.0);
  registerProcessorParameter( "Geometry::DetectorTransverseSize" ,
     			      "Define the detector transverse size used by efficiency algorithm (vector size must be 2 or 4; if 2 -> first value is min, second value is max; if 4 -> two first values define x edges , two last values define y edges) ",
     			      edges,
     			      vec ); 
  if( edges.size()==2 ){
     m_CaloGeomSetting.xmin=edges[0];
     m_CaloGeomSetting.ymin=edges[0];
     m_CaloGeomSetting.xmax=edges[1];
     m_CaloGeomSetting.ymax=edges[1];
   } 
   else if( edges.size()==4 ){
     m_CaloGeomSetting.xmin=edges[0];
     m_CaloGeomSetting.xmax=edges[1];
     m_CaloGeomSetting.ymin=edges[2];
     m_CaloGeomSetting.ymax=edges[3];
   }
  else{
    std::cout << "WARING : Wrong number of values in paramater Geometry::DetectorTransverseSize => will use default value -500.0, +500.0" << std::endl;
  }
  /*--------------------------------------------*/

  /*------------algorithm::AsicKeyFinder-----------*/
  std::vector<int> asicKeyFactor;
  asicKeyFactor.push_back(1);
  asicKeyFactor.push_back(12);
  asicKeyFactor.push_back(1000);
  registerProcessorParameter( "AsicKeyFinder::KeyFactor" ,
			      "Define the factor to build asic keys; default for sdhcal : 1000*K + 12*J + I",
			      m_AsicKeyFinderParameterSetting.keyFactors,
			      asicKeyFactor ); 
  
  registerProcessorParameter( "AsicKeyFinder::NPadX" ,
      			      "Number of pads in x direction per layer",
      			      m_AsicKeyFinderParameterSetting.nPadX,
      			      (int) 96 ); 

  registerProcessorParameter( "AsicKeyFinder::NPadY" ,
      			      "Number of pads in x direction per layer",
      			      m_AsicKeyFinderParameterSetting.nPadY,
      			      (int) 96 ); 

  registerProcessorParameter( "AsicKeyFinder::AsicNPad" ,
      			      "number of pads in x or y direction per asic (assuming a square)",
      			      m_AsicKeyFinderParameterSetting.asicNPad,
      			      (int) 8 ); 

  registerProcessorParameter( "AsicKeyFinder::LayerGap" ,
      			      "Gap size (in mm) between 2 layers",
      			      m_AsicKeyFinderParameterSetting.layerGap,
      			      (float) 26.131 ); 

  registerProcessorParameter( "AsicKeyFinder::PadSize" ,
      			      "Size of one pad in mm",
      			      m_AsicKeyFinderParameterSetting.padSize,
      			      (float) 10.408 ); 

  registerProcessorParameter( "AsicKeyFinder::PrintDebug" ,
      			      "Boolean to know if debug if printed",
      			      m_AsicKeyFinderParameterSetting.printDebug,
      			      (bool) false ); 
}

void sdhcalAsicProcessor::init()
{ 
  printParameters() ;

  file = new TFile(outputRootName.c_str(),"RECREATE");
  
  tree = (TTree*)file->Get("tree");
  if(!tree){
    std::cout << "tree creation" << std::endl; 
    tree = new TTree("tree","Shower variables");
  }
  tree->Branch("LayerID",&_layerID);
  tree->Branch("DifID",&_difID);
  tree->Branch("AsicID",&_asicID);
  tree->Branch("Efficiency1",&_efficiency1);
  tree->Branch("Efficiency2",&_efficiency2);
  tree->Branch("Efficiency3",&_efficiency3);
  tree->Branch("Efficiency1_Error",&_efficiency1_error);
  tree->Branch("Efficiency2_Error",&_efficiency2_error);
  tree->Branch("Efficiency3_Error",&_efficiency3_error);
  tree->Branch("Multiplicity",&_multiplicity);
  tree->Branch("Multiplicity_Error",&_multiplicity_error);
  tree->Branch("Ntrack",&_ntrack);
  tree->Branch("AsicPosition",&_asicPosition,"AsicPosition[3]/F");

  ntrack=new TH1D("ntrack","ntrack",10000,0,10000);
  effGlobal1=new TH1D("effGlobal1","effGlobal1",100,0,1);
  effGlobal2=new TH1D("effGlobal2","effGlobal2",100,0,1);
  effGlobal3=new TH1D("effGlobal3","effGlobal3",100,0,1);
  mulGlobal=new TH1D("mulGlobal","mulGlobal",100,0,4);
  keyList=new TH1D("keyList","keyList",50000,0,50000);
  
  mul2D=new TH2D("mul2D","mul2D",12,0,12,12,0,12);
  eff2D_thr1=new TH2D("eff2D_thr1","eff2D_thr1",12,0,12,12,0,12);
  eff2D_thr2=new TH2D("eff2D_thr2","eff2D_thr2",12,0,12,12,0,12);
  eff2D_thr3=new TH2D("eff2D_thr3","eff2D_thr3",12,0,12,12,0,12);
  trackPosition=new TH2D("trackPosition","trackPosition",1000,0,1000,1000,0,1000);

  _nRun = 0 ;
  _nEvt = 0 ;
  _goodTrackCounter = 0 ;

  /*--------------------Algorithms initialisation--------------------*/
  algo_Cluster=new algorithm::Cluster();
  algo_Cluster->SetClusterParameterSetting(m_ClusterParameterSetting);

  algo_ClusteringHelper=new algorithm::ClusteringHelper();
  algo_ClusteringHelper->SetClusteringHelperParameterSetting(m_ClusteringHelperParameterSetting);
  
  algo_Tracking=new algorithm::Tracking();
  algo_Tracking->SetTrackingParameterSetting(m_TrackingParameterSetting);

  algo_InteractionFinder=new algorithm::InteractionFinder();
  algo_InteractionFinder->SetInteractionFinderParameterSetting(m_InteractionFinderParameterSetting);

  m_EfficiencyParameterSetting.geometry=m_CaloGeomSetting;
  algo_Efficiency=new algorithm::Efficiency();
  algo_Efficiency->SetEfficiencyParameterSetting(m_EfficiencyParameterSetting);

  algo_AsicKeyFinder=new algorithm::AsicKeyFinder();
  algo_AsicKeyFinder->SetAsicKeyFinderParameterSetting(m_AsicKeyFinderParameterSetting);

  for(int k=0; k<_nActiveLayers; k++){
    caloobject::CaloLayer* aLayer=new caloobject::CaloLayer(k);
    //aLayer->setLayerParameterSetting(m_LayerParameterSetting);
    layers.push_back(aLayer);
    for( int i=0; i<_nAsicX; i++){
      for( int j=0; j<_nAsicY; j++){
	int key=algo_AsicKeyFinder->BuildAsicKey(i,j,k);
	caloobject::SDHCAL_Asic* asic=new caloobject::SDHCAL_Asic(key);
	int difNum = findDifID(key);
	int asicNum = findAsicID(key);
	asic->setASIC_ID(asicNum,difNum);
	std::pair<int,caloobject::SDHCAL_Asic*> myPair(key,asic);
	asics.insert(myPair);
      }
    }
  }
  std::cout << "asics.size() = " << asics.size() << std::endl;
}

void sdhcalAsicProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void sdhcalAsicProcessor::DoTracking()
{
  std::vector<caloobject::CaloCluster2D*> clusters;
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it){
    algo_Cluster->Run(it->second,clusters);
  }
  std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer);
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if(algo_ClusteringHelper->IsIsolatedCluster(*it,clusters)){
      delete *it; 
      clusters.erase(it); 
      it--;
    }
  }
  caloobject::CaloTrack* track=NULL;
  algo_Tracking->Run(clusters,track);
  if( NULL != track ){
    _goodTrackCounter++;
    algo_InteractionFinder->Run(clusters,track->getTrackParameters());
    if( algo_InteractionFinder->FindInteraction()==false )
      LayerProperties(clusters);
  }
  file->cd();
  
  delete track;
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    delete (*it);
}

void sdhcalAsicProcessor::LayerProperties(std::vector<caloobject::CaloCluster2D*> &clusters)
{
  int trackBegin= (*clusters.begin())->getLayerID();
  int trackEnd=(*(clusters.end()-1))->getLayerID();
  if(trackBegin==1) trackBegin=0;
  if(trackEnd==_nActiveLayers-2) trackEnd=_nActiveLayers-1;
    
  for(int K=trackBegin; K<=trackEnd; K++){
    layers.at(K)->Reset();
    algo_Efficiency->Run(layers.at(K),clusters);

    if( layers.at(K)->getNTracks()!=0 ){
      int key = algo_AsicKeyFinder->FindAsicKey(algo_Efficiency->getExpectedPosition());
      if( asics.find(key)==asics.end() ){
	std::cout << "Problem :: wrong key : " << key << "\t at vec = " << algo_Efficiency->getExpectedPosition() << std::endl;
	continue;
      }
      asics[ key ]->Update( layers.at(K) );
      trackPosition->Fill(algo_Efficiency->getExpectedPosition().x(),algo_Efficiency->getExpectedPosition().y());
    }
  }
}

void sdhcalAsicProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  //std::string initString;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");

  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      //initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      numElements = col->getNumberOfElements();
      //      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	CLHEP::Hep3Vector vec(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
	int cellID[]={IDdecoder(hit)["I"],IDdecoder(hit)["J"],IDdecoder(hit)["K-1"]};
	caloobject::CaloHit *aHit=new caloobject::CaloHit(cellID,vec,hit->getEnergy(),hit->getTime(),posShift);
	hitMap[cellID[2]].push_back(aHit);
	if( hit->getPosition()[0]<0 ||
	    hit->getPosition()[1]<0 ||
	    hit->getPosition()[2]<0 ){
	  std::cout << "WARNING : hit at " 
		    << hit->getPosition()[0] << ",\t" 
		    << hit->getPosition()[1] << ",\t" 
		    << hit->getPosition()[2] << std::endl;
	  getchar();
	}
      }
      DoTracking();
      clearVec();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}

void sdhcalAsicProcessor::clearVec()
{
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    for( std::vector<caloobject::CaloHit*>::iterator jt=(it->second).begin(); jt!=(it->second).end(); ++jt)
      delete *(jt);

  hitMap.clear();
}


void sdhcalAsicProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void sdhcalAsicProcessor::end(){ 

  file->cd();
  for(std::map<int,caloobject::SDHCAL_Asic*>::iterator it=asics.begin(); it!=asics.end(); it++){
    if(it->second->getAsicNtrack()>0){
      _layerID=it->second->getAsicKey()/1000;
      _difID=it->second->getAsicDifNum();
      _asicID=it->second->getAsicNum();
      _ntrack=it->second->getAsicNtrack();
      _efficiency1=it->second->getAsicEfficiency();
      _efficiency2=it->second->getAsicEfficiency2();
      _efficiency3=it->second->getAsicEfficiency3();
      _efficiency1_error=sqrt(_efficiency1*(1-_efficiency1)/_ntrack);
      _efficiency2_error=sqrt(_efficiency2*(1-_efficiency2)/_ntrack);
      _efficiency3_error=sqrt(_efficiency3*(1-_efficiency3)/_ntrack);

      ntrack->Fill(it->second->getAsicNtrack());
      effGlobal1->Fill(it->second->getAsicEfficiency());
      effGlobal2->Fill(it->second->getAsicEfficiency2());
      effGlobal3->Fill(it->second->getAsicEfficiency3());
      eff2D_thr1->Fill(it->second->getPosition().x(),it->second->getPosition().y(),it->second->getAsicEfficiency()/(float)_nActiveLayers);
      eff2D_thr2->Fill(it->second->getPosition().x(),it->second->getPosition().y(),it->second->getAsicEfficiency2()/(float)_nActiveLayers);
      eff2D_thr3->Fill(it->second->getPosition().x(),it->second->getPosition().y(),it->second->getAsicEfficiency3()/(float)_nActiveLayers);
      
      _asicPosition[0]=it->second->getPosition().x();
      _asicPosition[1]=it->second->getPosition().y();
      _asicPosition[2]=it->second->getPosition().z();
      
      if(it->second->getAsicEfficiency()>0.0){
	//std::cout << it->second->getAsicMultiplicity() << std::endl;
   	_multiplicity=it->second->getAsicMultiplicity();
   	_multiplicity_error=it->second->getAsicRMSMultiplicity();
	mulGlobal->Fill(it->second->getAsicMultiplicity());
	mul2D->Fill(it->second->getPosition().x(),it->second->getPosition().y(),it->second->getAsicMultiplicity()/(float)_nActiveLayers);
      }
      else{
   	_multiplicity=0;
   	_multiplicity_error=0;
      }
      tree->Fill();
    }
  }

  delete algo_Cluster;
  delete algo_ClusteringHelper;
  delete algo_Tracking;
  delete algo_InteractionFinder;
  delete algo_Efficiency;
  delete algo_AsicKeyFinder;

  for(std::vector<caloobject::CaloLayer*>::iterator it=layers.begin(); it!=layers.end(); ++it)
    delete (*it);
  layers.clear();
  
  for(std::map<int,caloobject::SDHCAL_Asic*>::iterator it=asics.begin(); it!=asics.end(); ++it)
    delete it->second;
  asics.clear();

  file->Write();
  file->Close();

  std::cout << "_goodTrackCounter " << _goodTrackCounter << std::endl;
}
