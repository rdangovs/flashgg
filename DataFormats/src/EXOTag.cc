#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/EXOTag.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

using namespace flashgg;

EXOTag::EXOTag() {}

EXOTag::~EXOTag() {}

EXOTag::EXOTag( edm::Ptr<DiPhotonCandidate> &diphoton, edm::Handle<edm::View<flashgg::Jet>> &jets, 
                edm::Handle<edm::View<flashgg::Electron>> &electrons, unsigned eventNumber){

    eventNumber_=eventNumber;

    diphoton_=diphoton;
    hasDiphoton_ = !diphoton.isNull();

    jets_=jets;
    setHasJets();
    setDijet();

    electrons_=electrons;
    setHasElectrons();
    setDielectron();

    if (!hasDiphoton_){
        hasDijet_ = false;
        hasElectrons_ = false;
    }

}

void EXOTag::setDijet(){

    float leadPt(0), subleadPt(0);
    int leadIndex(-1), subleadIndex(-1);

    for (unsigned i=0;i<jets_->size();i++){
        edm::Ptr<flashgg::Jet> jet = jets_->ptrAt(i);

        float dR_leadDP = deltaR(jet->eta(),jet->phi(),diphoton_->leadingPhoton()->eta(),diphoton_->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(jet->eta(),jet->phi(),diphoton_->subLeadingPhoton()->eta(),diphoton_->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;

        if (jet->pt() > leadPt){
            subleadPt = leadPt;
            subleadIndex = leadIndex;
            leadPt = jet->pt();
            leadIndex = i;
        }else if (jet->pt() > subleadPt){
            subleadPt = jet->pt();
            subleadIndex = i;
        }
    }

    if (leadIndex == -1 || subleadIndex == -1){
        hasDijet_=false;
    }else{
        hasDijet_=true;
        dijet_.first = jets_->ptrAt(leadIndex);
        dijet_.second = jets_->ptrAt(subleadIndex);
    }

}

void EXOTag::setDielectron(){

    float leadPt(0), subleadPt(0);
    int leadIndex(-1), subleadIndex(-1);

    for (unsigned i=0;i<electrons_->size();i++){
        edm::Ptr<flashgg::Electron> electron = electrons_->ptrAt(i);

        float dR_leadDP = deltaR(electron->eta(),electron->phi(),diphoton_->leadingPhoton()->eta(),diphoton_->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(electron->eta(),electron->phi(),diphoton_->subLeadingPhoton()->eta(),diphoton_->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;

        if (electron->pt() > leadPt){
            subleadPt = leadPt;
            subleadIndex = leadIndex;
            leadPt = electron->pt();
            leadIndex = i;
        }else if (electron->pt() > subleadPt){
            subleadPt = electron->pt();
            subleadIndex = i;
        }
    }

    if (leadIndex == -1 || subleadIndex == -1){
        hasDielectron_=false;
    }else{
        hasDielectron_=true;
        dielectron_.first = electrons_->ptrAt(leadIndex);
        dielectron_.second = electrons_->ptrAt(subleadIndex);
    }

}

void EXOTag::setHasJets(){
    if (EXOTag::countJetsOverPT(0.0) == 0){
        hasJets_ = false;
    }else{
        hasJets_ = true;
    }
}

void EXOTag::setHasElectrons(){
    if (EXOTag::countElectronsOverPT(0.0) == 0){
        hasElectrons_ = false;
    }else{
        hasElectrons_ = true;
    }
}

const unsigned EXOTag::getEventNumber(){ return eventNumber_; }

const float EXOTag::getDiphotonLeadPt(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->pt() : -999.; }

const float EXOTag::getDiphotonSubleadPt(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->pt() : -999.; }

const float EXOTag::getDiphotonLeadEta(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->eta() : -999.; }

const float EXOTag::getDiphotonSubleadEta(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->eta() : -999.; }

const float EXOTag::getDiphotonMass(){ return hasDiphoton_ ? diphoton_->mass() : -999.; }

const float EXOTag::getDiphotonLeadR9(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->r9() : -999.; }

const float EXOTag::getDiphotonSubleadR9(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->r9() : -999.; }

const float EXOTag::getDiphotonLeadEtaSC(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->superCluster()->eta() : -999.; }

const float EXOTag::getDiphotonSubleadEtaSC(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->superCluster()->eta() : -999.; }

const float EXOTag::getDiphotonLeadPhiSC(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->superCluster()->phi() : -999.; }

const float EXOTag::getDiphotonSubleadPhiSC(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->superCluster()->phi() : -999.; }

const int EXOTag::getDiphotonCategory(){

    if (!hasDiphoton_){
        return -999;
    }else{
        
        float boundaryEB(1.4442);
        float boundaryEELo(1.566), boundaryEEHi(2.5);
        int leadCat(-1), subLeadCat(-1);

        if (fabs(diphoton_->leadingPhoton()->eta()) < boundaryEB){
            leadCat = 0;
        }else if (fabs(diphoton_->leadingPhoton()->eta()) > boundaryEELo
                    && fabs(diphoton_->subLeadingPhoton()->eta()) < boundaryEEHi){
            leadCat = 1;
        }
        if (fabs(diphoton_->subLeadingPhoton()->eta()) < boundaryEB){
            subLeadCat = 0;
        }else if (fabs(diphoton_->subLeadingPhoton()->eta()) > boundaryEELo
                    && fabs(diphoton_->subLeadingPhoton()->eta()) < boundaryEEHi){
            subLeadCat = 1;
        }
        if (leadCat == 0 && subLeadCat == 0){
            return 0;
        }else if ((leadCat == 0 && subLeadCat == 1) || (leadCat == 1 && subLeadCat == 0)){
            return 1;
        }else{
            return -999;
        }

    }
}

const float EXOTag::getDiphotonLeadCHI(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->egChargedHadronIso() : -999.; }

const float EXOTag::getDiphotonSubleadCHI(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->egChargedHadronIso() : -999.; }

const float EXOTag::getDiphotonLeadEGPhoIso(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->egPhotonIso() : -999.; }

const float EXOTag::getDiphotonSubleadEGPhoIso(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->egPhotonIso() : -999.; }

const float EXOTag::getDiphotonLeadF5x5SigmaIetaIeta(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->full5x5_sigmaIetaIeta() : -999.; }

const float EXOTag::getDiphotonSubleadF5x5SigmaIetaIeta(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->full5x5_sigmaIetaIeta() : -999.; }

const float EXOTag::getDiphotonLeadFull5x5R9(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->full5x5_r9() : -999.; }

const float EXOTag::getDiphotonSubleadFull5x5R9(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->full5x5_r9() : -999.; }

const float EXOTag::getDiphotonLeadHadronicOverEM(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->hadTowOverEm() : -999.; }

const float EXOTag::getDiphotonSubleadHadronicOverEM(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->hadTowOverEm() : -999.; }

const int EXOTag::getDiphotonLeadIsSaturated(){
    if (!hasDiphoton_){
        return -999.;
    }else{
        return (int)(diphoton_->leadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated) 
                                              && !diphoton_->leadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
    }
}

const int EXOTag::getDiphotonSubleadIsSaturated(){
    if (!hasDiphoton_){
        return -999.;
    }else{
        return (int)(diphoton_->subLeadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated) 
                                              && !diphoton_->subLeadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
    }
}

const int EXOTag::getDiphotonLeadPassElectronVeto(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->passElectronVeto() : -999.; }

const int EXOTag::getDiphotonSubleadPassElectronVeto(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->passElectronVeto() : -999.; }

const int EXOTag::getDiphotonNConv(){ return hasDiphoton_ ? diphoton_->nConv() : -999; }

const float EXOTag::getDiphotonPullConv(){ return hasDiphoton_ ? diphoton_->pullConv() : -999; }

const int EXOTag::countJetsOverPT(float ptCut){
    
    unsigned count(0);
    for (unsigned i=0;i<jets_->size();i++){
        edm::Ptr<flashgg::Jet> jet = jets_->ptrAt(i);

        float dR_leadDP = deltaR(jet->eta(),jet->phi(),diphoton_->leadingPhoton()->eta(),diphoton_->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(jet->eta(),jet->phi(),diphoton_->subLeadingPhoton()->eta(),diphoton_->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;
        if (jet->pt() > ptCut) count++;
    }
    return count;
}

const int EXOTag::countElectronsOverPT(float ptCut){

    unsigned count(0);
    for (unsigned i=0;i<electrons_->size();i++){
        edm::Ptr<flashgg::Electron> electron = electrons_->ptrAt(i);

        float dR_leadDP = deltaR(electron->eta(),electron->phi(),diphoton_->leadingPhoton()->eta(),diphoton_->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(electron->eta(),electron->phi(),diphoton_->subLeadingPhoton()->eta(),diphoton_->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;
        if (electron->pt()) count++;
    }
    return count;
}

       
        
const int EXOTag::getJetMultiplicities_All(){return countJetsOverPT(0.0);}
const int EXOTag::getJetMultiplicities_EGT20(){return countJetsOverPT(20);}
const int EXOTag::getJetMultiplicities_EGT30(){return countJetsOverPT(30);}
const int EXOTag::getJetMultiplicities_EGT40(){return countJetsOverPT(40);}

const float EXOTag::getDijetLeadPt(){ return hasDijet_ ? dijet_.first->pt() : -999.; }
const float EXOTag::getDijetSubleadPt(){ return hasDijet_ ? dijet_.second->pt() : -999.; }
const float EXOTag::getDijetLeadEta(){ return hasDijet_ ? dijet_.first->eta() : -999.; }
const float EXOTag::getDijetSubleadEta(){ return hasDijet_ ? dijet_.second->eta() : -999.; }
const float EXOTag::getDijetMass(){ return hasDijet_ ? (dijet_.first->p4()+dijet_.second->p4()).M() : -999.; }
const float EXOTag::getDijetDeltaEta(){ return hasDijet_ ? fabs(dijet_.first->eta()-dijet_.second->eta()) : -999.; }
const float EXOTag::getDijetZeppenfeld(){ return hasDijet_ ? fabs(diphoton_->eta() - 0.5*(dijet_.first->eta()-dijet_.second->eta())) : -999.; }
const float EXOTag::getDijetDeltaPhi_jj(){ return hasDijet_ ? fabs(deltaPhi(dijet_.first->phi() , dijet_.second->phi())) : -999.; } 
const float EXOTag::getDijetDeltaPhi_ggjj(){ return hasDijet_ ? fabs(deltaPhi(diphoton_->phi() , (dijet_.first->p4()+dijet_.second->p4()).Phi())) : -999; }


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

