/*
 *
* This is a part of CTPPS offline software.
* Author:
*   Fabrizio Ferro (ferro@ge.infn.it)
*   Enrico Robutti (robutti@ge.infn.it)
*   Fabio Ravera   (fabio.ravera@cern.ch)
*
*/
#ifndef RecoCTPPS_PixelLocal_RPixRoadFinder_H
#define RecoCTPPS_PixelLocal_RPixRoadFinder_H

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
//#include "DataFormats/Common/interface/DetSetVectorNew.h"
//#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "DataFormats/CTPPSReco/interface/CTPPSPixelCluster.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"
#include "RecoCTPPS/PixelLocal/interface/RPixClusterToHit.h" 
#include "RecoCTPPS/PixelLocal/interface/RPixDetPatternFinder.h"

#include "FWCore/Framework/interface/ESWatcher.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/VeryForwardRPTopology/interface/RPTopology.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/Records/interface/VeryForwardMisalignedGeometryRecord.h"


#include <vector>
#include <set>



class RPixRoadFinder : public RPixDetPatternFinder{

  public:
    explicit RPixRoadFinder(const edm::ParameterSet& param);
    ~RPixRoadFinder() override;
    void findPattern() override;

  private:
    edm::ParameterSet param_;
    int verbosity_;
    double roadRadius_;
    unsigned int minRoadSize_;
    unsigned int maxRoadSize_;
    void run(const edm::DetSetVector<CTPPSPixelRecHit> &input, const CTPPSGeometry & geometry, std::vector<Road> &roads);

};


#endif
