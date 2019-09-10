// -*- C++ -*-
//
// Package:    ​RecoBTag/​SecondaryVertex
// Class:      DeepNNTagInfoProducer
//
/**\class DeepNNTagInfoProducer DeepNNTagInfoProducer.cc ​RecoBTag/DeepFlavour/plugins/DeepNNTagInfoProducer.cc
 *
 * Description: EDProducer that produces collection of DeepNNTagInfos
 *
 * Implementation:
 *    A collection of CandIPTagInfo and CandSecondaryVertexTagInfo and a CombinedSVComputer ESHandle is taken as input and a collection of DeepNNTagInfos
 *    is produced as output.
 */
//
// Original Author:  Mauro Verzetti (U. Rochester)
// Added templated functionality : Praveen C Tiwari & Jyothsna Komaragiri (IISc)
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include <map>

////////////////////////////////////////////////////////////
//
// TemplatedDeepNNTagInfoProducer
//
// class declaration
//
template<typename IPTag, typename SVTag>
class TemplatedDeepNNTagInfoProducer : public edm::stream::EDProducer<> {
public:
        explicit TemplatedDeepNNTagInfoProducer(const edm::ParameterSet&);
        ~TemplatedDeepNNTagInfoProducer() override;

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
        void beginStream(edm::StreamID) override {}
        void produce(edm::Event&, const edm::EventSetup&) override;
        void endStream() override {}

   // ----------member data ---------------------------
        const edm::EDGetTokenT<std::vector<SVTag> > svSrc_;
        CombinedSVComputer computer_;

        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
template<typename IPTag, typename SVTag>
TemplatedDeepNNTagInfoProducer<IPTag,SVTag>::TemplatedDeepNNTagInfoProducer(const edm::ParameterSet& iConfig) :
  svSrc_( consumes<std::vector<SVTag> >(iConfig.getParameter<edm::InputTag>("svTagInfos")) ),
  computer_(iConfig.getParameter<edm::ParameterSet>("computer"))
{
  produces<std::vector<reco::ShallowTagInfo> >();
  vtxToken_ = consumes<reco::VertexCollection>(
                    edm::InputTag("offlineSlimmedPrimaryVertices4D"));
  puToken_ = consumes<std::vector<PileupSummaryInfo >>(edm::InputTag("slimmedAddPileupInfo"));

}

template<typename IPTag, typename SVTag>
TemplatedDeepNNTagInfoProducer<IPTag,SVTag>::~TemplatedDeepNNTagInfoProducer()
{

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
template<typename IPTag, typename SVTag>
void TemplatedDeepNNTagInfoProducer<IPTag,SVTag>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get input TagInfos
  edm::Handle<std::vector<SVTag> > svTagInfos;
  iEvent.getByToken(svSrc_, svTagInfos);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  edm::Handle<std::vector <PileupSummaryInfo> > PUInfo;
  iEvent.getByToken(puToken_, PUInfo);

  // create the output collection
  auto tagInfos = std::make_unique<std::vector<reco::ShallowTagInfo> >();

  auto npv_0_z = vertices->at(0).z();

  float PUrho = 0;
  std::vector<PileupSummaryInfo>::const_iterator ipu;
  for (ipu = PUInfo->begin(); ipu != PUInfo->end(); ++ipu) {
      if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0

      for (unsigned int i=0; i<ipu->getPU_zpositions().size(); ++i) {
          auto PU_z = (ipu->getPU_zpositions())[i];
          if ( std::abs(PU_z - npv_0_z) < 1) PUrho++;
      }
  }	
  PUrho /= 20.;

  // loop over TagInfos
  for(auto iterTI = svTagInfos->begin(); iterTI != svTagInfos->end(); ++iterTI) {
    // get TagInfos
    const SVTag & svTagInfo = *(iterTI);
    const IPTag & ipInfo = *(iterTI->trackIPTagInfoRef().get());
    
    reco::TaggingVariableList vars = computer_(ipInfo, svTagInfo);
    std::vector<float> tagValList = vars.getList(reco::btau::trackEtaRel,false);
    vars.insert(reco::btau::jetNTracksEtaRel, tagValList.size());
    tagValList = vars.getList(reco::btau::trackSip2dSig,false);
    vars.insert(reco::btau::jetNSelectedTracks, tagValList.size());
    vars.finalize(); //fix the TaggingVariableList, nothing should be added/removed
    
    //If not SV found set it to 0, not to non-existent
    if(!vars.checkTag(reco::btau::jetNSecondaryVertices))
      vars.insert(reco::btau::jetNSecondaryVertices, 0);
    if(!vars.checkTag(reco::btau::vertexNTracks))
      vars.insert(reco::btau::vertexNTracks, 0);
    
    vars.insert(reco::btau::PUrho, PUrho);

    tagInfos->emplace_back(
			   vars,
			   svTagInfo.jet()
			   );    
  }

  // put the output in the event
  iEvent.put( std::move(tagInfos) );
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<typename IPTag, typename SVTag>
void TemplatedDeepNNTagInfoProducer<IPTag,SVTag>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
//For PFJets
typedef TemplatedDeepNNTagInfoProducer<reco::CandIPTagInfo, reco::CandSecondaryVertexTagInfo> DeepNNTagInfoProducer;
//For CaloJets
typedef TemplatedDeepNNTagInfoProducer<reco::TrackIPTagInfo, reco::SecondaryVertexTagInfo> TrackDeepNNTagInfoProducer;
DEFINE_FWK_MODULE(DeepNNTagInfoProducer);
DEFINE_FWK_MODULE(TrackDeepNNTagInfoProducer);
