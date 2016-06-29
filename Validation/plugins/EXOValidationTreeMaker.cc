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
#include "flashgg/DataFormats/interface/Electron.h"
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
    float weight;
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
    int   refinedCategory; //0 for central EBEB, 1 for the remaining EBEB, 2 for the EBEE
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
    float nConv;
    float pullConv;
   
    //cosThetaStar
    float cosThetaStar; 
   
    //Extra objects
    int   jetMultiplicity;
    int   jetMultiplicity_EGT20;
    int   jetMultiplicity_EGT30;
    int   jetMultiplicity_EGT40;

    //Dijet
    float dijetLeadPt;
    float dijetSubleadPt;
    float dijetLeadEta;
    float dijetSubleadEta;
    float dijetMass;
    float dijetDeltaEta;
    float dijetZeppenfeld;
    float dijetDeltaPhi_jj;
    float dijetDeltaPhi_ggjj;
    
    //Electrons
    float electronMultiplicity_EGT35;
    float electronMultiplicity_EGT75;

    //Dielectrons
    float dielecLeadPt;
    float dielecSubleadPt;
    float dielecLeadEta;
    float dielecSubleadEta;
    float dielecMass;
    float dielecDeltaEta;
    float dielecZeppenfeld;
    float dielecDeltaPhi_ee;
    float dielecDeltaPhi_ggee;

    void init()
    {
        mgg          = -999;
        weight       = -999;
        leadPt       = -999;
        subLeadPt     = -999;
        leadR9       = -999;
        subLeadR9     = -999;
        leadEtaSC    = -999;
        subLeadEtaSC = -999;
        leadPhiSC    = -999;
        subLeadPhiSC = -999;
        category     = -999;
        refinedCategory = -999;
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
        nConv = -999;
        pullConv = -999;
        
        cosThetaStar = -999; 
        
        jetMultiplicity = 0;
        jetMultiplicity_EGT20 = 0;
        jetMultiplicity_EGT30 = 0;
        jetMultiplicity_EGT40 = 0;

        dijetLeadPt = -999;
        dijetSubleadPt = -999;
        dijetLeadEta = -999;
        dijetSubleadEta = -999;
        dijetMass = -999;
        dijetDeltaEta = -999;
        dijetZeppenfeld = -999;
        dijetDeltaPhi_jj = -999;
        dijetDeltaPhi_ggjj = -999;
        
        electronMultiplicity_EGT35 = 0;;
        electronMultiplicity_EGT75 = 0;;

        dielecLeadPt = -999;
        dielecSubleadPt = -999;
        dielecLeadEta = -999;
        dielecSubleadEta = -999;
        dielecMass = -999;
        dielecDeltaEta = -999;
        dielecZeppenfeld = -999;
        dielecDeltaPhi_ee = -999;
        dielecDeltaPhi_ggee = -999;

    }

    void printExtraObjects(){

        cout << "Jet Info" << endl;
        cout << "Multiplciity:" << endl;
        cout << setw(12) << "All" << setw(12) << ">20GeV" << setw(12) << ">30GeV" << setw(12) << ">40GeV" << endl;
        cout << setw(12) << jetMultiplicity << setw(12) << jetMultiplicity_EGT20 << setw(12) << jetMultiplicity_EGT30 << setw(12) << jetMultiplicity_EGT40 << endl;
        cout << "Dijet Kinematics" << endl;
        cout << setw(12) << "Lead Pt" << setw(12) << "Sublead Pt" << setw(12) << "Lead Eta" << setw(12) << "Sublead Eta" << endl;
        cout << setw(12) << dijetLeadPt << setw(12) << dijetSubleadPt << setw(12) << dijetLeadEta << setw(12) << dijetSubleadEta << endl;
        cout << setw(12) << "Mass" << setw(12) << "dEta" << setw(12) << "Zepp" << setw(12) << "dPhi_jj" << setw(12) << "dPhi_ggjj" << endl;
        cout << setw(12) << dijetMass << setw(12) << dijetDeltaEta << setw(12) << dijetZeppenfeld << setw(12) << dijetDeltaPhi_jj << setw(12) << dijetDeltaPhi_ggjj << endl;
 
        cout << "Electron Info" << endl;
        cout << "Multiplciity:" << endl;
        cout << setw(12) << ">35GeV" << setw(12) << ">75GeV" << endl;
        cout << setw(12) << electronMultiplicity_EGT35 << setw(12) << electronMultiplicity_EGT75 << endl;
        cout << "Electron Kinematics" << endl;
        cout << setw(12) << "Lead Pt" << setw(12) << "Sublead Pt" << setw(12) << "Lead Eta" << setw(12) << "Sublead Eta" << endl;
        cout << setw(12) << dielecLeadPt << setw(12) << dielecSubleadPt << setw(12) << dielecLeadEta << setw(12) << dielecSubleadEta << endl;
        cout << setw(12) << "Mass" << setw(12) << "dEta" << setw(12) << "Zepp" << setw(12) << "dPhi_ee" << setw(12) << "dPhi_ggee" << endl;
        cout << setw(12) << dielecMass << setw(12) << dielecDeltaEta << setw(12) << dielecZeppenfeld << setw(12) << dielecDeltaPhi_ee << setw(12) << dielecDeltaPhi_ggee << endl;

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

    std::vector<edm::EDGetTokenT<View<flashgg::Jet> > >  tokenJets_;

    std::vector<edm::InputTag>                           inputTagJets_;
    EDGetTokenT< edm::View<flashgg::Electron> > electronToken_;
    EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
    EDGetTokenT<double> rhoToken_;

    typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

    TTree     *diphotonTree;

    diphotonInfo   eInfo;
    Int_t       event_number;
    bool        debug_;

    void FillDielectronInfo(Handle<View<flashgg::Electron>> &electrons,Ptr<flashgg::DiPhotonCandidate> diphoton);
    void FillElectronInfo(Handle<View<flashgg::Electron>> &electrons,Ptr<flashgg::DiPhotonCandidate> diphoton);
    void FillJetInfo(edm::Handle<edm::View<flashgg::Jet> > &jets,Ptr<flashgg::DiPhotonCandidate> diphoton);
    void FillDijetInfo(edm::Handle<edm::View<flashgg::Jet> > &jets,Ptr<flashgg::DiPhotonCandidate> diphoton);
    void FillDiphotonInfo(Ptr<flashgg::DiPhotonCandidate> diphoton);

};

EXOValidationTreeMaker::EXOValidationTreeMaker( const edm::ParameterSet &iConfig ):

    inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
    electronToken_( consumes<View<flashgg::Electron>>(iConfig.getParameter<edm::InputTag>("ElectronTag"))),
    diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    rhoToken_( consumes<double>( iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
    debug_( iConfig.getUntrackedParameter<bool>( "debug", false ) )
{
    for( uint i = 0; i < inputTagJets_.size(); i++ ) {
        auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
        tokenJets_.push_back(token); 
    }

    event_number = 0;
}

EXOValidationTreeMaker::~EXOValidationTreeMaker()
{
    event_number = 0;
}


void
EXOValidationTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    eInfo.init();
    

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
        
 //       cout << "DEBUG: removed the ID cuts to create a reasonable tree for practice." << endl; 
        /*       
        if (!passPhotonIDCuts(diPhotons->ptrAt(diphoIndex)->leadingPhoton(),rhoFixedGrd)){ 
            //std::cout << "DEBUG LC diphoton" << diphoIndex << " lead photon failed photonID cuts" << std::endl;
            continue;
        }
        if (!passPhotonIDCuts(diPhotons->ptrAt(diphoIndex)->subLeadingPhoton(),rhoFixedGrd)){ 
            //std::cout << "DEBUG LC diphoton" << diphoIndex << " sublead photon failed photonID cuts" << std::endl;
            continue;
        }*/

        float tempPt = diPhotons->ptrAt(diphoIndex)->leadingPhoton()->pt()+diPhotons->ptrAt(diphoIndex)->subLeadingPhoton()->pt();
        if (tempPt > maxDiphoPt){
            maxDiphoPt = tempPt;
            maxDiphoIndex = diphoIndex;
        }

    }
    

    if (diPhotonsSize > 0 &&  maxDiphoIndex>-1){
        
        //Category selection
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

        eInfo.refinedCategory = eInfo.category + 1; 
        if (fabs(diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->eta()) < 1 && fabs(diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->eta())  < 1) 
            eInfo.refinedCategory = 0; 
        
        cout << eInfo.refinedCategory << endl; 

        //Preselection
        float ptCut(75),mggCutEBEB(230),mggCutEBEE(320);

        float diphotonLeadPt            = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->pt();
        float diphotonSubLeadPt         = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->pt();
        float diphotonLeadEtaSC         = diPhotons->ptrAt(maxDiphoIndex)->leadingPhoton()->superCluster()->eta();
        float diphotonSubLeadEtaSC      = diPhotons->ptrAt(maxDiphoIndex)->subLeadingPhoton()->superCluster()->eta();

        bool passesPtCut = diphotonLeadPt > ptCut && diphotonSubLeadPt > ptCut;
        bool passesLeadEtaSCCut = fabs(diphotonLeadEtaSC) < boundaryEEHi && !(fabs(diphotonLeadEtaSC) > boundaryEB && fabs(diphotonLeadEtaSC) < boundaryEELo);
        bool passesSubLeadEtaSCCut = fabs(diphotonSubLeadEtaSC) < boundaryEEHi && !(fabs(diphotonSubLeadEtaSC) > boundaryEB && fabs(diphotonSubLeadEtaSC) < boundaryEELo);
        bool passesOneBarrelEtaSCCut = fabs(diphotonLeadEtaSC) < boundaryEB || fabs(diphotonSubLeadEtaSC) < boundaryEB;
        bool passesMassCut(false);
        if (eInfo.category == 0){
            passesMassCut = diPhotons->ptrAt(maxDiphoIndex)->mass() > mggCutEBEB;
        }else if (eInfo.category == 1){
            passesMassCut = diPhotons->ptrAt(maxDiphoIndex)->mass() > mggCutEBEE;
        }
        
        //debug 
 //       std::cout << passesPtCut << " " << passesLeadEtaSCCut << " " << passesSubLeadEtaSCCut << " " << passesOneBarrelEtaSCCut << " " << passesMassCut << std::endl;   

        if (passesPtCut && passesLeadEtaSCCut && passesSubLeadEtaSCCut && passesOneBarrelEtaSCCut && passesMassCut){

            Ptr<flashgg::DiPhotonCandidate> diphoton = diPhotons->ptrAt(maxDiphoIndex);
            FillDiphotonInfo(diphoton);
            
            //std::cout << "CANDIDATE PHOTON HAS INDEX " << maxDiphoIndex << std::endl;
            //std::cout << " mgg is iPhotons->ptrAt(maxDiphoIndex)->mass() " << diPhotons->ptrAt(maxDiphoIndex)->mass()  <<std::endl;

            //EXTRA OBJECTS
            
            //JETS
            unsigned jetCollectionIndex = diPhotons->ptrAt(maxDiphoIndex)->jetCollectionIndex();
            JetCollectionVector Jets(inputTagJets_.size());
            for (size_t i=0;i<inputTagJets_.size();i++){
                iEvent.getByToken(tokenJets_[i],Jets[i]);
            }
            edm::Handle<edm::View<flashgg::Jet>> jets = Jets[jetCollectionIndex];
            //Filling jet info
            FillJetInfo(jets,diphoton);
            //Dijet info
            FillDijetInfo(jets,diphoton);

            //ELECTRONS
            Handle<View<flashgg::Electron>> electrons;
            iEvent.getByToken(electronToken_,electrons);
            //electrons
            FillElectronInfo(electrons,diphoton);
            //dielectrons
            FillDielectronInfo(electrons,diphoton);

            diphotonTree->Fill();

            cout << "Diphoton passes selection" << endl;
            eInfo.printExtraObjects();

        }else{
            cout << "Diphoton fails selection" << endl;
        }


    }
    event_number++;
}

void EXOValidationTreeMaker::FillDielectronInfo(Handle<View<flashgg::Electron>> &electrons,Ptr<flashgg::DiPhotonCandidate> diphoton){

    float leadPt(0), subleadPt(0);
    int leadIndex(-1);
    int subleadIndex(-1);

     for (unsigned i=0;i<electrons->size();i++){
        Ptr<flashgg::Electron> electron = electrons->ptrAt(i);

        float dR_leadDP = deltaR(electron->eta(),electron->phi(),diphoton->leadingPhoton()->eta(),diphoton->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(electron->eta(),electron->phi(),diphoton->subLeadingPhoton()->eta(),diphoton->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;

        if (electron->pt() > leadPt){
            subleadPt = leadPt; 
            leadPt = electron->pt();
            subleadIndex = leadIndex;
            leadIndex = i;
        }else if (electron->pt() > subleadPt){
            subleadPt = electron->pt();
            subleadIndex = i;
        }
    }

    if (leadIndex != -1 && subleadIndex != -1){
        
        Ptr<flashgg::Electron> leadElectron = electrons->ptrAt(leadIndex);;
        Ptr<flashgg::Electron> subleadElectron = electrons->ptrAt(subleadIndex);;

        eInfo.dielecLeadPt = leadElectron->pt();
        eInfo.dielecSubleadPt = subleadElectron->pt();
        eInfo.dielecLeadEta = leadElectron->eta();
        eInfo.dielecSubleadEta = subleadElectron->eta();
        eInfo.dielecMass = (leadElectron->p4() + subleadElectron->p4()).M();
        eInfo.dielecDeltaEta = fabs(leadElectron->eta() - subleadElectron->eta());
        eInfo.dielecZeppenfeld = fabs(diphoton->eta() - 0.5*(leadElectron->eta()+subleadElectron->eta()));
        eInfo.dielecDeltaPhi_ee = fabs(deltaPhi(leadElectron->phi() , subleadElectron->phi()));
        eInfo.dielecDeltaPhi_ggee = fabs(deltaPhi(diphoton->phi() , (leadElectron->p4()+subleadElectron->p4()).Phi()));

    }

}

void EXOValidationTreeMaker::FillElectronInfo(Handle<View<flashgg::Electron>> &electrons,Ptr<flashgg::DiPhotonCandidate> diphoton){
    unsigned EGT_35_count = 0;
    unsigned EGT_75_count = 0;

    for (unsigned i=0;i<electrons->size();i++){
        Ptr<flashgg::Electron> electron = electrons->ptrAt(i);

        float dR_leadDP = deltaR(electron->eta(),electron->phi(),diphoton->leadingPhoton()->eta(),diphoton->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(electron->eta(),electron->phi(),diphoton->subLeadingPhoton()->eta(),diphoton->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;

        if (electron->pt() > 35) EGT_35_count++;
        if (electron->pt() > 75) EGT_75_count++;
    }

    eInfo.electronMultiplicity_EGT35 = EGT_35_count; 
    eInfo.electronMultiplicity_EGT75 = EGT_75_count; 

}

void EXOValidationTreeMaker::FillJetInfo(edm::Handle<edm::View<flashgg::Jet> > &jets,Ptr<flashgg::DiPhotonCandidate> diphoton){

    unsigned jet_count = 0;
    unsigned EGT_20_count = 0;
    unsigned EGT_30_count = 0;
    unsigned EGT_40_count = 0;

    for (size_t i=0;i<jets->size();i++){
        Ptr<flashgg::Jet> jet = jets->ptrAt(i);

        float dR_leadDP = deltaR(jet->eta(),jet->phi(),diphoton->leadingPhoton()->eta(),diphoton->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(jet->eta(),jet->phi(),diphoton->subLeadingPhoton()->eta(),diphoton->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) {
            continue;
        }

        jet_count++;
        if (jet->pt() > 20) EGT_20_count++;
        if (jet->pt() > 30) EGT_30_count++;
        if (jet->pt() > 40) EGT_40_count++;
    }

    eInfo.jetMultiplicity = jet_count;
    eInfo.jetMultiplicity_EGT20 = EGT_20_count;
    eInfo.jetMultiplicity_EGT30 = EGT_30_count;
    eInfo.jetMultiplicity_EGT40 = EGT_40_count;

}

void EXOValidationTreeMaker::FillDijetInfo(edm::Handle<edm::View<flashgg::Jet> > &jets,Ptr<flashgg::DiPhotonCandidate> diphoton){

    //Get lead and sublead jets
    float leadPt(0), subleadPt(0);
    int leadIndex(-1);
    int subleadIndex(-1);
    for (unsigned i=0;i<jets->size();i++){
        Ptr<flashgg::Jet> jet = jets->ptrAt(i);

        float dR_leadDP = deltaR(jet->eta(),jet->phi(),diphoton->leadingPhoton()->eta(),diphoton->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(jet->eta(),jet->phi(),diphoton->subLeadingPhoton()->eta(),diphoton->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;

        if (jet->pt() > leadPt){
            subleadPt = leadPt; 
            leadPt = jet->pt();
            subleadIndex = leadIndex;
            leadIndex = i;
        }else if (jet->pt() > subleadPt){
            subleadPt = jet->pt();
            subleadIndex = i;
        }
    }

    if (leadIndex != -1 && subleadIndex != -1){

        Ptr<flashgg::Jet> leadJet = jets->ptrAt(leadIndex);
        Ptr<flashgg::Jet> subleadJet = jets->ptrAt(subleadIndex);

        eInfo.dijetLeadPt = leadJet->pt();
        eInfo.dijetSubleadPt = subleadJet->pt();
        eInfo.dijetLeadEta = leadJet->eta();
        eInfo.dijetSubleadEta = subleadJet->eta();
        eInfo.dijetMass = (leadJet->p4() + subleadJet->p4()).M();
        eInfo.dijetDeltaEta = fabs(leadJet->eta() - subleadJet->eta());
        eInfo.dijetZeppenfeld = fabs(diphoton->eta() - 0.5*(leadJet->eta()+subleadJet->eta()));
        eInfo.dijetDeltaPhi_jj = fabs(deltaPhi(leadJet->phi() , subleadJet->phi()));
        eInfo.dijetDeltaPhi_ggjj = fabs(deltaPhi(diphoton->phi() , (leadJet->p4()+subleadJet->p4()).Phi()));

        std::cout << "Lead Jet Pt = " << eInfo.dijetLeadPt << std::endl;
        std::cout << "Sublead Jet Pt = " << eInfo.dijetSubleadPt << std::endl;
        std::cout << "Lead Jet Eta = " << eInfo.dijetLeadEta << std::endl;
        std::cout << "Sublead Jet Eta = " << eInfo.dijetSubleadEta << std::endl;
        std::cout << "Dijet Mass = " << eInfo.dijetMass << std::endl;
        std::cout << "Dijet Delta Eta = " << eInfo.dijetDeltaEta << std::endl;
        std::cout << "Dijet Zeppenfeld = " << eInfo.dijetZeppenfeld << std::endl;
        std::cout << "Dijet Delta Phi jj = " << eInfo.dijetDeltaPhi_jj << std::endl;
        std::cout << "Dijet Delta Phit ggjj = " << eInfo.dijetDeltaPhi_ggjj << std::endl;
    }

}


void EXOValidationTreeMaker::FillDiphotonInfo(Ptr<flashgg::DiPhotonCandidate> diphoton){

    eInfo.eventID = event_number;
    eInfo.mgg = diphoton->mass();
    eInfo.leadPt            = diphoton->leadingPhoton()->pt();
    eInfo.subLeadPt          = diphoton->subLeadingPhoton()->pt();
    eInfo.leadR9            = diphoton->leadingPhoton()->r9();
    eInfo.subLeadR9          = diphoton->subLeadingPhoton()->r9();
    eInfo.leadEtaSC         = diphoton->leadingPhoton()->superCluster()->eta();
    eInfo.subLeadEtaSC      = diphoton->subLeadingPhoton()->superCluster()->eta();
    eInfo.leadPhiSC         = diphoton->leadingPhoton()->superCluster()->phi();
    eInfo.subLeadPhiSC      = diphoton->subLeadingPhoton()->superCluster()->phi();
    eInfo.leadChargedHadronIso = diphoton->leadingPhoton()->egChargedHadronIso();
    eInfo.leadPfPhoIso03 = diphoton->leadingPhoton()->egPhotonIso();
    eInfo.leadFull5x5_sigmaIetaIeta = diphoton->leadingPhoton()->full5x5_sigmaIetaIeta();
    eInfo.leadFull5x5_r9 = diphoton->leadingPhoton()->full5x5_r9();
    eInfo.leadHadronicOverEm = diphoton->leadingPhoton()->hadTowOverEm();
    eInfo.subLeadChargedHadronIso = diphoton->subLeadingPhoton()->egChargedHadronIso();
    eInfo.subLeadPfPhoIso03 = diphoton->subLeadingPhoton()->egPhotonIso();
    eInfo.subLeadFull5x5_sigmaIetaIeta = diphoton->subLeadingPhoton()->full5x5_sigmaIetaIeta();
    eInfo.subLeadFull5x5_r9 = diphoton->subLeadingPhoton()->full5x5_r9();
    eInfo.subLeadHadronicOverEm = diphoton->subLeadingPhoton()->hadTowOverEm();
    eInfo.leadIsSaturated = (int)(diphoton->leadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated) 
                              && !diphoton->leadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
    eInfo.subLeadIsSaturated = (int)(diphoton->subLeadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated) 
                                 && !diphoton->subLeadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
    eInfo.leadPassElectronVeto = (int)(diphoton->leadingPhoton()->passElectronVeto());
    eInfo.subLeadPassElectronVeto = (int)(diphoton->subLeadingPhoton()->passElectronVeto());
    eInfo.nConv = diphoton->nConv();
    eInfo.pullConv = diphoton->pullConv();
    
    eInfo.cosThetaStar = (float) (2 * (diphoton->subLeadingPhoton()->energy() * diphoton->leadingPhoton()->pz() - diphoton->leadingPhoton()->energy() * diphoton->subLeadingPhoton()->pz()) / (diphoton->mass() * TMath::Sqrt (diphoton->mass() * diphoton->mass() + diphoton->leadingPhoton()->pt() * diphoton->leadingPhoton()->pt()))); 

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
    diphotonTree->Branch( "refinedCategory ", &eInfo.refinedCategory, Form("%s/I", "refinedCategory")); //addition of refined category
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
    diphotonTree->Branch( "nConv " , &eInfo.nConv  , Form("%s/F","nConv"));
    diphotonTree->Branch( "pullConv " , &eInfo.pullConv  , Form("%s/F","pullConv"));
    
    diphotonTree->Branch( "cosThetaStar" , &eInfo.cosThetaStar          , Form("%s/F","cosThetaStar"));

    diphotonTree->Branch( "jetMultiplicity " , &eInfo.jetMultiplicity  , Form("%s/I","jetMultiplicity"));
    diphotonTree->Branch( "jetMultiplicity_EGT20 " , &eInfo.jetMultiplicity_EGT20  , Form("%s/I","jetMultiplicity_EGT20"));
    diphotonTree->Branch( "jetMultiplicity_EGT30 " , &eInfo.jetMultiplicity_EGT30  , Form("%s/I","jetMultiplicity_EGT30"));
    diphotonTree->Branch( "jetMultiplicity_EGT40 " , &eInfo.jetMultiplicity_EGT40  , Form("%s/I","jetMultiplicity_EGT40"));
    diphotonTree->Branch( "dijetLeadPt          " , &eInfo.dijetLeadPt           , Form("%s/F","dijetLeadPt"));
    diphotonTree->Branch( "dijetSubleadPt          " , &eInfo.dijetSubleadPt           , Form("%s/F","dijetSubleadPt"));
    diphotonTree->Branch( "dijetLeadEta          " , &eInfo.dijetLeadEta           , Form("%s/F","dijetLeadEta"));
    diphotonTree->Branch( "dijetSubleadEta          " , &eInfo.dijetSubleadEta           , Form("%s/F","dijetSubleadEta"));
    diphotonTree->Branch( "dijetMass          " , &eInfo.dijetMass           , Form("%s/F","dijetMass"));
    diphotonTree->Branch( "dijetDeltaEta          " , &eInfo.dijetDeltaEta           , Form("%s/F","dijetDeltaEta"));
    diphotonTree->Branch( "dijetZeppenfeld          " , &eInfo.dijetZeppenfeld           , Form("%s/F","dijetZeppenfeld"));
    diphotonTree->Branch( "dijetDeltaPhi_jj          " , &eInfo.dijetDeltaPhi_jj           , Form("%s/F","dijetDeltaPhi_jj"));
    diphotonTree->Branch( "dijetDeltaPhi_ggjj          " , &eInfo.dijetDeltaPhi_ggjj           , Form("%s/F","dijetDeltaPhi_ggjj"));
    diphotonTree->Branch( "electronMultiplicity_EGT35 " , &eInfo.electronMultiplicity_EGT35  , Form("%s/I","electronMultiplicity_EGT35"));
    diphotonTree->Branch( "electronMultiplicity_EGT75 " , &eInfo.electronMultiplicity_EGT75  , Form("%s/I","electronMultiplicity_EGT75"));
    diphotonTree->Branch( "dielecLeadPt          " , &eInfo.dielecLeadPt           , Form("%s/F","dielecLeadPt"));
    diphotonTree->Branch( "dielecSubleadPt          " , &eInfo.dielecSubleadPt           , Form("%s/F","dielecSubleadPt"));
    diphotonTree->Branch( "dielecLeadEta          " , &eInfo.dielecLeadEta           , Form("%s/F","dielecLeadEta"));
    diphotonTree->Branch( "dielecSubleadEta          " , &eInfo.dielecSubleadEta           , Form("%s/F","dielecSubleadEta"));
    diphotonTree->Branch( "dielecMass          " , &eInfo.dielecMass           , Form("%s/F","dielecMass"));
    diphotonTree->Branch( "dielecDeltaEta          " , &eInfo.dielecDeltaEta           , Form("%s/F","dielecDeltaEta"));
    diphotonTree->Branch( "dielecZeppenfeld          " , &eInfo.dielecZeppenfeld           , Form("%s/F","dielecZeppenfeld"));
    diphotonTree->Branch( "dielecDeltaPhi_ee          " , &eInfo.dielecDeltaPhi_ee           , Form("%s/F","dielecDeltaPhi_ee"));
    diphotonTree->Branch( "dielecDeltaPhi_ggee          " , &eInfo.dielecDeltaPhi_ggee           , Form("%s/F","dielecDeltaPhi_ggee"));
    
    diphotonTree->Branch( "justTest", &eInfo.dielecDeltaPhi_ggee, Form ("%s/F", "dielecDeltaPhi_ggee")); 
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
