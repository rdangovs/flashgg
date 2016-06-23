#ifndef FLASHgg_EXOTag_h
#define FLASHgg_EXOTag_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "DataFormats/Math/interface/deltaPhi.h"

namespace flashgg {

    class EXOTag: public DiPhotonTagBase
    {
    public:
        EXOTag();
        ~EXOTag();

        EXOTag *clone() const { return ( new EXOTag( *this ) ); }

        EXOTag( edm::Ptr<DiPhotonCandidate> &diphoton, edm::Handle<edm::View<flashgg::Jet>> &jets, 
                edm::Handle<edm::View<flashgg::Electron>> &electrons, unsigned eventNumber);

        const unsigned getEventNumber();
        const float getDiphotonLeadPt(); 
        const float getDiphotonSubleadPt(); 
        const float getDiphotonLeadEta(); 
        const float getDiphotonSubleadEta(); 
        const float getDiphotonMass();

    private:
        unsigned eventNumber_;

        edm::Ptr<DiPhotonCandidate> diphoton_;
        edm::Handle<edm::View<flashgg::Jet>> jets_;
        edm::Handle<edm::View<flashgg::Electron>> electrons_;

        typedef std::pair<edm::Ptr<flashgg::Jet>,edm::Ptr<flashgg::Jet>> Dijet;
        typedef std::pair<edm::Ptr<flashgg::Electron>,edm::Ptr<flashgg::Electron>> Dielectron;

        bool hasDiphoton_;
        bool hasJets_;
        bool hasElectrons_;

        Dijet dijet_;
        Dielectron dielectron_;
        bool hasDijet_;
        bool hasDielectron_;

        void setDijet();
        void setDielectron();
        void setHasJets();
        void setHasElectrons();

    };

}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

