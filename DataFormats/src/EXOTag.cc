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

    hasDiphoton_ = true;   //Placeholder

    jets_=jets;
    electrons_=electrons;

    setHasJets();
    setHasElectrons();

    setDijet();
    setDielectron();

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
            leadPt = jet->pt();
            subleadIndex = leadIndex;
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
            leadPt = electron->pt();
            subleadIndex = leadIndex;
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

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

