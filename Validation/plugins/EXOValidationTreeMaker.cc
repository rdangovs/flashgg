// -----------------------
// By Y.Haddad & L.Corpe  12/2014
// Modified by E.Scott    04/2016
//
// Jet validation analyzer: tree maker for the jet studies in flashgg
// -----------------------

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"

#include "flashgg/DataFormats/interface/VertexCandidateMap.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"


#include "TTree.h"
#include "TMatrix.h"
#include "TVector.h"
#include "TLorentzVector.h"

// **********************************************************************
#ifndef FLASHgg_VertexCandidateMap_h
#define FLASHgg_VertexCandidateMap_h
namespace flashgg {
    // typedef std::map<edm::Ptr<reco::Vertex>,edm::PtrVector<pat::PackedCandidate> > VertexCandidateMap;
    typedef std::pair<edm::Ptr<reco::Vertex>, edm::Ptr<pat::PackedCandidate> > VertexCandidatePair;
    typedef std::vector<VertexCandidatePair> VertexCandidateMap;

    struct compare_by_vtx {
        bool operator()( const VertexCandidatePair &left, const VertexCandidatePair &right )
        {
            return ( left.first < right.first );
        }
    };

    struct compare_with_vtx {
        bool operator()( const VertexCandidatePair &left, const edm::Ptr<reco::Vertex> &right )
        {
            return ( left.first < right );
        }
        bool operator()( const edm::Ptr<reco::Vertex> &left, const VertexCandidatePair &right )
        {
            return( left < right.first );
        }
    };

    struct compare_by_cand {
        bool operator()( const VertexCandidatePair &left, const VertexCandidatePair &right )
        {
            return ( left.second < right.second );
        }
    };

    struct compare_with_cand {
        bool operator()( const VertexCandidatePair &left, const edm::Ptr<pat::PackedCandidate> &right )
        {
            return ( left.second < right );
        }
        bool operator()( const edm::Ptr<pat::PackedCandidate> &left, const VertexCandidatePair &right )
        {
            return( left < right.second );
        }
    };
}
#endif



struct diphotonInfo {
    int   eventID;
    float mgg;
    float leadPt           ;
    float subLeadPt         ;
    float leadR9           ;
    float subLeadR9         ;
    float leadEtaSC        ;
    float subLeadEtaSC     ;
    float leadPhiSC        ;
    float subLeadPhiSC     ;

    void init()
    {
        mgg          = -999;
        leadPt       = -999;
        subLeadPt     = -999;
        leadR9       = -999;
        subLeadR9     = -999;
        leadEtaSC    = -999;
        subLeadEtaSC = -999;
        leadPhiSC    = -999;
        subLeadPhiSC = -999;
    }
};


// **********************************************************************

using namespace std;
using namespace edm;
using namespace flashgg;

// **********************************************************************

class EXOValidationTreeMaker : public edm::EDAnalyzer
{
public:
    explicit EXOValidationTreeMaker( const edm::ParameterSet & );
    ~EXOValidationTreeMaker();

    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


private:

    edm::Service<TFileService> fs_;

    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    // Additional methods
    void initEventStructure();


    //bool GenRecoMatching(reco::GenJet genjet, const PtrVector<flashgg::Jet>& RecoJets){}
    //EDGetTokenT< VertexCandidateMap > vertexCandidateMapTokenDz_;
    //EDGetTokenT< VertexCandidateMap > vertexCandidateMapTokenAOD_;

    //EDGetTokenT< edm::View<reco::GenParticle> >          genPartToken_;
    //EDGetTokenT< edm::View<reco::GenJet> >               genJetToken_;
    //EDGetTokenT< edm::View<flashgg::Jet> >               jetDzToken_;
    //std::vector<edm::EDGetTokenT<View<flashgg::Jet> > >  tokenJets_;
    //std::vector<edm::InputTag>                           inputTagJets_;
    EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
    //EDGetTokenT< View<reco::Vertex> >                    vertexToken_;
    //EDGetTokenT< VertexCandidateMap > vertexCandidateMapToken_;

    //edm::InputTag 					qgVariablesInputTag;
    //edm::EDGetTokenT<edm::ValueMap<float>> 		qgToken;

    //typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;


    TTree     *diphotonTree;

    diphotonInfo   eInfo;
    Int_t       event_number;
    //std::string jetCollectionName;

    //GenPartInfo jGenPhotonInfo;

    //bool        usePUJetID;
    //bool        photonJetVeto;
    //bool        homeGenJetMatching_;
    //bool        ZeroVertexOnly_;
    bool        debug_;
    //bool        useVBFTagPhotonMatching_;

};

EXOValidationTreeMaker::EXOValidationTreeMaker( const edm::ParameterSet &iConfig ):
    //genPartToken_( consumes<View<reco::GenParticle> >( iConfig.getUntrackedParameter<InputTag> ( "GenParticleTag", InputTag( "prunedGenParticles" ) ) ) ),
    //genJetToken_( consumes<View<reco::GenJet> >( iConfig.getUntrackedParameter<InputTag> ( "GenJetTag", InputTag( "slimmedGenJets" ) ) ) ),
    //jetDzToken_   ( consumes<View<flashgg::Jet> >( iConfig.getParameter<InputTag>( "JetTagDz" ) ) ),
    //inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
    diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    //vertexToken_( consumes<View<reco::Vertex> >( iConfig.getUntrackedParameter<InputTag> ( "VertexTag", InputTag( "offlineSlimmedPrimaryVertices" ) ) ) ),
    //vertexCandidateMapToken_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTag" ) ) ),

    //qgVariablesInputTag( iConfig.getParameter<edm::InputTag>( "qgVariablesInputTag" ) ),

    //usePUJetID( iConfig.getUntrackedParameter<bool>( "UsePUJetID"   , false ) ),
    //photonJetVeto( iConfig.getUntrackedParameter<bool>( "PhotonJetVeto", true ) ),
    //homeGenJetMatching_( iConfig.getUntrackedParameter<bool>( "homeGenJetMatching", false ) ),
    //ZeroVertexOnly_( iConfig.getUntrackedParameter<bool>( "ZeroVertexOnly", false ) ),
    debug_( iConfig.getUntrackedParameter<bool>( "debug", false ) )
    //useVBFTagPhotonMatching_( iConfig.getUntrackedParameter<bool>( "useVBFTagPhotonMatching", false ) )

{
    //for( uint i = 0; i < inputTagJets_.size(); i++ ) {
    //    auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
    //    tokenJets_.push_back(token); 
    //}

    event_number = 0;
    //jetCollectionName = iConfig.getParameter<string>( "StringTag" );
    //qgToken	= consumes<edm::ValueMap<float>>( edm::InputTag( qgVariablesInputTag.label(), "qgLikelihood" ) );

}

EXOValidationTreeMaker::~EXOValidationTreeMaker()
{
    event_number = 0;
}


void
EXOValidationTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    if( debug_ ) {
        std::cout << "\e[0;31m";
        std::cout << setw( 6 )  << "========================= "            << std::endl;
        std::cout << setw( 12 ) << "Event" << setw( 12 ) << event_number   << std::endl;
        std::cout << setw( 6 )  << "------------------------- "            << std::endl;
        std::cout << "\e[0m" << std::endl;
    }


    //Handle<View<reco::Vertex> > vtxs;
    //iEvent.getByToken( vertexToken_, vtxs );
    //const PtrVector<reco::Vertex>& vtxs = primaryVertices->ptrVector();

    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    iEvent.getByToken( diPhotonToken_, diPhotons );
    initEventStructure();
    size_t diPhotonsSize = diPhotons->size();
    for( unsigned int diphoIndex = 0; diphoIndex < diPhotonsSize; diphoIndex++ ) {
    eInfo.eventID    = event_number;
    eInfo.mgg    = diPhotons->ptrAt(diphoIndex)->mass();
    eInfo.leadPt            = diPhotons->ptrAt(diphoIndex)->leadingPhoton()->pt();
    eInfo.subLeadPt          = diPhotons->ptrAt(diphoIndex)->subLeadingPhoton()->pt();
    eInfo.leadR9            = diPhotons->ptrAt(diphoIndex)->leadingPhoton()->r9();
    eInfo.subLeadR9          = diPhotons->ptrAt(diphoIndex)->subLeadingPhoton()->r9();
    eInfo.leadEtaSC         = diPhotons->ptrAt(diphoIndex)->leadingPhoton()->superCluster()->eta();
    eInfo.subLeadEtaSC      = diPhotons->ptrAt(diphoIndex)->subLeadingPhoton()->superCluster()->eta();
    eInfo.leadPhiSC         = diPhotons->ptrAt(diphoIndex)->leadingPhoton()->superCluster()->phi();
    eInfo.subLeadPhiSC      = diPhotons->ptrAt(diphoIndex)->subLeadingPhoton()->superCluster()->phi();
    std::cout << "DEBUG event : "<< eInfo.eventID << " diphoton " << diphoIndex << " mgg " << eInfo.mgg << std::endl;
    
    diphotonTree->Fill();
    }

    event_number++;
}

void
EXOValidationTreeMaker::beginJob()
{
    // +++ trees
    std::string type( "diphotonTree_" );

    diphotonTree = fs_->make<TTree>( type.c_str(), "tree");
    diphotonTree->Branch( "eventID"    , &eInfo.eventID       , Form("%s/I","eventID"));
    diphotonTree->Branch( "mgg          " , &eInfo.mgg           , Form("%s/F","mgg"));
    diphotonTree->Branch( "leadPt       " , &eInfo.leadPt        , Form("%s/F","leadPt"));
    diphotonTree->Branch( "subLeadPt    " , &eInfo.subLeadPt      , Form("%s/F","subLeadPt"));
    diphotonTree->Branch( "leadR9       " , &eInfo.leadR9        , Form("%s/F","leadR9"));
    diphotonTree->Branch( "subLeadR9    " , &eInfo.subLeadR9      , Form("%s/F","subLeadR9"));
    diphotonTree->Branch( "leadEtaSC    " , &eInfo.leadEtaSC     , Form("%s/F","leadEtaSC"));
    diphotonTree->Branch( "subLeadEtaSC " , &eInfo.subLeadEtaSC  , Form("%s/F","subLeadEtaSC"));
    diphotonTree->Branch( "leadPhiSC    " , &eInfo.leadPhiSC     , Form("%s/F","leadPhiSC"));
    diphotonTree->Branch( "subLeadPhiSC " , &eInfo.subLeadPhiSC  , Form("%s/F","subLeadPhiSC"));
}

void EXOValidationTreeMaker::endJob()
{

}

void EXOValidationTreeMaker::initEventStructure()
{

}


void EXOValidationTreeMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    // The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}

typedef EXOValidationTreeMaker FlashggEXOValidationTreeMaker;
DEFINE_FWK_MODULE( FlashggEXOValidationTreeMaker );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
