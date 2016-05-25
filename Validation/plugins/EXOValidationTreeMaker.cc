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
    int   category;
    float leadChargedHadronIso;
    float leadPfPhoIso03;
    float leadFull5x5_sigmaIetaIeta;
    float leadFull5x5_r9;
    float leadHadronicOverEm;
    float subLeadChargedHadronIso;
    float subLeadPfPhoIso03;
    float subLeadFull5x5_sigmaIetaIeta;
    float subLeadFull5x5_r9;
    float subLeadHadronicOverEm;
    int   leadIsSaturated;
    int   subLeadIsSaturated;
    int   leadPassElectronVeto;
    int   subLeadPassElectronVeto;

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
        category     = -999;
        leadChargedHadronIso = -999;
        leadPfPhoIso03 = -999;
        leadFull5x5_sigmaIetaIeta = -999;
        leadFull5x5_r9 = -999;
        leadHadronicOverEm = -999;
        subLeadChargedHadronIso = -999;
        subLeadPfPhoIso03 = -999;
        subLeadFull5x5_sigmaIetaIeta = -999;
        subLeadFull5x5_r9 = -999;
        subLeadHadronicOverEm = -999;
        leadIsSaturated = -999;
        subLeadIsSaturated = -999;
        leadPassElectronVeto = -999;
        subLeadPassElectronVeto = -999;
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
    bool passPhotonIDCuts(const flashgg::Photon* pho, const double rho); 
    float  correctIsoGam(const flashgg::Photon* pho, const double rho);


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
    EDGetTokenT<double> rhoToken_;
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
    rhoToken_( consumes<double>( iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
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

    /*
    if( debug_ ) {
        std::cout << "\e[0;31m";
        std::cout << setw( 6 )  << "========================= "            << std::endl;
        std::cout << setw( 12 ) << "Event" << setw( 12 ) << event_number   << std::endl;
        std::cout << setw( 6 )  << "------------------------- "            << std::endl;
        std::cout << "\e[0m" << std::endl;
    }
    */


    //Handle<View<reco::Vertex> > vtxs;
    //iEvent.getByToken( vertexToken_, vtxs );
    //const PtrVector<reco::Vertex>& vtxs = primaryVertices->ptrVector();
    
    //Set cat boundaries
    float boundaryEB(1.4442);
    float boundaryEELo(1.566), boundaryEEHi(2.5);
    int leadCat(-1), subLeadCat(-1);
    

    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    iEvent.getByToken( diPhotonToken_, diPhotons );
    initEventStructure();
    size_t diPhotonsSize = diPhotons->size();

    Handle<double> rhoHandle; 
    iEvent.getByToken( rhoToken_, rhoHandle );
    const double rhoFixedGrd = *( rhoHandle.product() );

        
        //Pick highest pt diphoton
        float maxDiphoPt(0.0);
        int  maxDiphoIndex(-1);
        for( unsigned int diphoIndex = 0; diphoIndex < diPhotonsSize; diphoIndex++ ) {
            
            if (!passPhotonIDCuts(diPhotons->ptrAt(diphoIndex)->leadingPhoton(),rhoFixedGrd)){ 
                 //std::cout << "DEBUG LC diphoton" << diphoIndex << " lead photon failed photonID cuts" << std::endl;
                continue;}
            else {
                 //std::cout << "DEBUG LC diphoton" << diphoIndex << " lead photon PASSED  photonID cuts" << std::endl;
            }
            if (!passPhotonIDCuts(diPhotons->ptrAt(diphoIndex)->subLeadingPhoton(),rhoFixedGrd)){ 
                
                 //std::cout << "DEBUG LC diphoton" << diphoIndex << " sublead photon failed photonID cuts" << std::endl;
                continue;
                }
            else {
                 //std::cout << "DEBUG LC diphoton" << diphoIndex << " sublead photon PASSED  photonID cuts" << std::endl;
            }


            float tempPt = diPhotons->ptrAt(diphoIndex)->leadingPhoton()->pt()+diPhotons->ptrAt(diphoIndex)->subLeadingPhoton()->pt();
            if (tempPt > maxDiphoPt){
                maxDiphoPt = tempPt;
                maxDiphoIndex = diphoIndex;
            }
        }
        
        if (diPhotonsSize > 0 &&  maxDiphoIndex>-1){
        //Fill struct
        eInfo.eventID    = event_number;
        eInfo.mgg    = diPhotons->ptrAt(maxDiphoIndex)->mass();
        eInfo.leadPt            = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->pt();
        eInfo.subLeadPt          = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->pt();
        eInfo.leadR9            = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->r9();
        eInfo.subLeadR9          = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->r9();
        eInfo.leadEtaSC         = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->superCluster()->eta();
        eInfo.subLeadEtaSC      = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->superCluster()->eta();
        eInfo.leadPhiSC         = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->superCluster()->phi();
        eInfo.subLeadPhiSC      = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->superCluster()->phi();
        eInfo.leadChargedHadronIso = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->egChargedHadronIso();
        //eInfo.leadChargedHadronIso = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->chargedHadronIso();
        eInfo.leadPfPhoIso03 = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->egPhotonIso();
        //eInfo.leadPfPhoIso03 = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->pfPhoIso03();
        eInfo.leadFull5x5_sigmaIetaIeta = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->full5x5_sigmaIetaIeta();
        eInfo.leadFull5x5_r9 = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->full5x5_r9();
        eInfo.leadHadronicOverEm = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->hadTowOverEm();
        //eInfo.leadHadronicOverEm = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->hadronicOverEm();
        //eInfo.subLeadChargedHadronIso = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->chargedHadronIso();
        eInfo.subLeadChargedHadronIso = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->egChargedHadronIso();
        eInfo.subLeadPfPhoIso03 = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->egPhotonIso();
        //eInfo.subLeadPfPhoIso03 = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->pfPhoIso03();
        eInfo.subLeadFull5x5_sigmaIetaIeta = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->full5x5_sigmaIetaIeta();
        eInfo.subLeadFull5x5_r9 = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->full5x5_r9();
        eInfo.subLeadHadronicOverEm = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->hadTowOverEm();
        //eInfo.subLeadHadronicOverEm = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->hadronicOverEm();
        eInfo.leadIsSaturated = (int)(diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated) 
                                  && !diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
        eInfo.subLeadIsSaturated = (int)(diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated) 
                                     && !diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
        eInfo.leadPassElectronVeto = (int)(diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->passElectronVeto());
        eInfo.subLeadPassElectronVeto = (int)(diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->passElectronVeto());

        if (fabs(diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->eta()) < boundaryEB){
            leadCat = 0;
        }else if (fabs(diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->eta()) > boundaryEELo
                    && fabs(diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->eta()) < boundaryEEHi){
            leadCat = 1;
        }
        if (fabs(diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->eta()) < boundaryEB){
            subLeadCat = 0;
        }else if (fabs(diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->eta()) > boundaryEELo
                    && fabs(diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->eta()) < boundaryEEHi){
            subLeadCat = 1;
        }
        if (leadCat == 0 && subLeadCat == 0){
            eInfo.category = 0;
        }else if ((leadCat == 0 && subLeadCat == 1) || (leadCat == 1 && subLeadCat == 0)){
            eInfo.category = 1;
        }

        //std::cout << "DEBUG event : " << eInfo.eventID << " diphoton " << maxDiphoIndex << " mgg " << eInfo.mgg << " cat " << eInfo.category;

        //Preselection
        float ptCut(75),mggCutEBEB(230),mggCutEBEE(320);
        bool passesPtCut = eInfo.leadPt > ptCut && eInfo.subLeadPt > ptCut;
        bool passesLeadEtaSCCut = fabs(eInfo.leadEtaSC) < boundaryEEHi && !(fabs(eInfo.leadEtaSC) > boundaryEB && fabs(eInfo.leadEtaSC) < boundaryEELo);
        bool passesSubLeadEtaSCCut = fabs(eInfo.subLeadEtaSC) < boundaryEEHi && !(fabs(eInfo.subLeadEtaSC) > boundaryEB && fabs(eInfo.subLeadEtaSC) < boundaryEELo);
        bool passesOneBarrelEtaSCCut = fabs(eInfo.leadEtaSC) < boundaryEB || fabs(eInfo.subLeadEtaSC) < boundaryEB;
        bool passesMassCut(false);
        if (eInfo.category == 0){
            passesMassCut = eInfo.mgg > mggCutEBEB;
        }else if (eInfo.category == 1){
            passesMassCut = eInfo.mgg > mggCutEBEE;
        }

        if (passesPtCut && passesLeadEtaSCCut && passesSubLeadEtaSCCut && passesOneBarrelEtaSCCut && passesMassCut){
            //std::cout << " PASS" << endl;
            diphotonTree->Fill();
            std::cout << "CANDIDATE PHOTON HAS INDEX " << maxDiphoIndex << std::endl;
            std::cout << " mgg is iPhotons->ptrAt(maxDiphoIndex)->mass() " << diPhotons->ptrAt(maxDiphoIndex)->mass()  <<std::endl;
            //std::cout << " CORRECTED mgg is iPhotons->ptrAt(maxDiphoIndex)->mass() " << diPhotonsSystematics->ptrAt(maxDiphoIndex)->mass()  <<std::endl;
        }/*else{
            std::cout << " FAIL" << endl;
        }*/
    }/*else{
        std::cout << "Event has no diphotons" << std::endl;
    }*/
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
    diphotonTree->Branch( "category " , &eInfo.category  , Form("%s/I","category"));
    diphotonTree->Branch( "leadChargedHadronIso " , &eInfo.leadChargedHadronIso  , Form("%s/F","leadChargedHadronIso"));
    diphotonTree->Branch( "leadPfPhoIso03 " , &eInfo.leadPfPhoIso03  , Form("%s/F","leadPfPhoIso03"));
    diphotonTree->Branch( "leadFull5x5_sigmaIetaIeta " , &eInfo.leadFull5x5_sigmaIetaIeta  , Form("%s/F","leadFull5x5_sigmaIetaIeta"));
    diphotonTree->Branch( "leadFull5x5_r9 " , &eInfo.leadFull5x5_r9  , Form("%s/F","leadFull5x5_r9"));
    diphotonTree->Branch( "leadHadronicOverEm " , &eInfo.leadHadronicOverEm  , Form("%s/F","leadHadronicOverEm"));
    diphotonTree->Branch( "subLeadChargedHadronIso " , &eInfo.subLeadChargedHadronIso  , Form("%s/F","subLeadChargedHadronIso"));
    diphotonTree->Branch( "subLeadPfPhoIso03 " , &eInfo.subLeadPfPhoIso03  , Form("%s/F","subLeadPfPhoIso03"));
    diphotonTree->Branch( "subLeadFull5x5_sigmaIetaIeta " , &eInfo.subLeadFull5x5_sigmaIetaIeta  , Form("%s/F","subLeadFull5x5_sigmaIetaIeta"));
    diphotonTree->Branch( "subLeadFull5x5_r9 " , &eInfo.subLeadFull5x5_r9  , Form("%s/F","subLeadFull5x5_r9"));
    diphotonTree->Branch( "subLeadHadronicOverEm " , &eInfo.subLeadHadronicOverEm  , Form("%s/F","subLeadHadronicOverEm"));
    diphotonTree->Branch( "leadIsSaturated " , &eInfo.leadIsSaturated  , Form("%s/I","leadIsSaturated"));
    diphotonTree->Branch( "subLeadIsSaturated " , &eInfo.subLeadIsSaturated  , Form("%s/I","subLeadIsSaturated"));
    diphotonTree->Branch( "leadPassElectronVeto " , &eInfo.leadPassElectronVeto  , Form("%s/I","leadPassElectronVeto"));
    diphotonTree->Branch( "subLeadPassElectronVeto " , &eInfo.subLeadPassElectronVeto  , Form("%s/I","subLeadPassElectronVeto"));

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

float EXOValidationTreeMaker::correctIsoGam(const flashgg::Photon* pho, const double rho){
    float eta = pho->superCluster()->eta();
    float isoGam =  pho->egPhotonIso();
    //float isoGam =  pho->pfPhoIso03();
    float pt = pho->pt();
    float alpha =2.5;
    float A =-1.;
    float kappa =-1.;
        if (fabs(eta) < 0.9 ){ A=0.17 ; kappa =4.5e-3 ;}
        if (fabs(eta) >= 0.9 && fabs(eta)<1.4442 ){ A=0.14 ; kappa =4.5e-3;}
        if (fabs(eta) >= 1.566 && fabs(eta)<2.0 ){ A=0.11 ; kappa =3e-3;}
        if (fabs(eta) >= 2.0 && fabs(eta)<2.2 ){ A=0.14 ; kappa =3e-3;}
        if (fabs(eta) >= 2.2 && fabs(eta)<2.5 ){ A=0.22 ; kappa =3e-3;}
    
    float corrIsoGam = alpha + isoGam - rho*A  - kappa *pt;
    //std::cout << "DEBUG cpriginal isoGam " << isoGam << " corrected isoGHam " << corrIsoGam << " : alpha " << alpha << " A  " << A << " kappa " << kappa <<  " rho " << rho << std::endl;
    return corrIsoGam;


}


bool EXOValidationTreeMaker::passPhotonIDCuts(const flashgg::Photon* pho, const double rho){
    float eta = pho->superCluster()->eta();
    int saturated = int(pho->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated));
    int weird = int(pho->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
    //float isoCh = pho->chargedHadronIso();
    float isoCh = pho->egChargedHadronIso();
    //float isoGam =  pho->pfPhoIso03();
    //float hoe = pho->hadronicOverEm() ;
    float hoe = pho->hadTowOverEm() ;
    float sieie = pho->full5x5_sigmaIetaIeta(); 
    int eleVeto = pho->passElectronVeto();
    bool pass=0;

    float correctedIsoGam = correctIsoGam(pho, rho);
    if (eleVeto){
        if (fabs(eta) < 1.4442 && !saturated){
             if (isoCh <5 && correctedIsoGam < 2.75 && hoe <0.05 && sieie<0.0105){
                pass=1;
             }
        }
        if (fabs(eta) < 1.4442 && saturated && !weird){
             if (isoCh <5 && correctedIsoGam < 2.75 && hoe <0.05 && sieie<0.0112){
                pass=1;
             }
        }
        if (fabs(eta) > 1.566 && !saturated){
             if (isoCh <5 && correctedIsoGam < 2.0 && hoe <0.05 && sieie<0.028){
                pass=1;
             }
        }
        if (fabs(eta) < 1.566 && saturated && !weird){
             if (isoCh <5 && correctedIsoGam < 2.0 && hoe <0.05 && sieie<0.03){
                pass=1;
             }
        }

    }
    return pass;
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
