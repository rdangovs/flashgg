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
    hasDiphoton_ = !diphoton->isNull();

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

    unsigned count(0);
    for (unsigned i=0;i<jets_->size();i++){
        edm::Ptr<flashgg::Jet> jet = jets_->ptrAt(i);

        float dR_leadDP = deltaR(jet->eta(),jet->phi(),diphoton_->leadingPhoton()->eta(),diphoton_->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(jet->eta(),jet->phi(),diphoton_->subLeadingPhoton()->eta(),diphoton_->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;
        count++;
    }
    hasJets_ = !(count==0);
}

void EXOTag::setHasElectrons(){

    unsigned count(0);
    for (unsigned i=0;i<electrons_->size();i++){
        edm::Ptr<flashgg::Electron> electron = electrons_->ptrAt(i);

        float dR_leadDP = deltaR(electron->eta(),electron->phi(),diphoton_->leadingPhoton()->eta(),diphoton_->leadingPhoton()->phi());
        float dR_subLeadDP = deltaR(electron->eta(),electron->phi(),diphoton_->subLeadingPhoton()->eta(),diphoton_->subLeadingPhoton()->phi());
        if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;
        count++;
    }
    hasElectrons_ = !(count==0);
}

const unsigned EXOTag::getEventNumber(){ return eventNumber_; }

const float EXOTag::getDiphotonLeadPt(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->pt() : -999.; }

const float EXOTag::getDiphotonSubleadPt(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->pt() : -999.; }

const float EXOTag::getDiphotonLeadEta(){ return hasDiphoton_ ? diphoton_->leadingPhoton()->eta() : -999.; }

const float EXOTag::getDiphotonSubleadEta(){ return hasDiphoton_ ? diphoton_->subLeadingPhoton()->eta() : -999.; }

const float EXOTag::getDiphotonMass(){ return hasDiphoton_ ? diphoton_->mass() : -999.; }

const float EXOTag::getDiphotonLeadR9(){ return hasDiphoton_ ? diphoton->leadingPhoton()->r9() : -999.; }

const float EXOTag::getDiphotonSubleadR9(){ return hasDiphoton_ ? diphoton->subleadingPhoton()->r9() : -999.; }

const float EXOTag::getDiphotonLeadEtaSC(){ return hasDiphoton_ ? diphoton->leadingPhoton()->superCluster()->eta() : -999.; }

const float EXOTag::getDiphotonSubleadEtaSC(){ return hasDiphoton_ ? diphoton->subleadingPhoton()->superCluster()->eta() : -999.; }

const float EXOTag::getDiphotonLeadPhiSC(){ return hasDiphoton_ ? diphoton->leadingPhoton()->superCluster()->phi() : -999.; }

const float EXOTag::getDiphotonSubleadPhiSC(){ return hasDiphoton_ ? diphoton->subleadingPhoton()->superCluster()->phi() : -999.; }

const unsigned EXOTag::getDiphotonCategory(){

    if (!hasDiphoton_){
        return -999.;
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
        }

    }
}












// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

