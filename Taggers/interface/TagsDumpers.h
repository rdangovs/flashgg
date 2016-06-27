#ifndef flashgg_TagsDumpers_h
#define flashgg_TagsDumpers_h

#include "flashgg/DataFormats/interface/UntaggedTag.h"
#include "flashgg/DataFormats/interface/VBFTag.h"
#include "flashgg/DataFormats/interface/TTHLeptonicTag.h"
#include "flashgg/DataFormats/interface/TTHHadronicTag.h"
#include "flashgg/DataFormats/interface/VHLooseTag.h"
#include "flashgg/DataFormats/interface/VHTightTag.h"
#include "flashgg/DataFormats/interface/VHHadronicTag.h"
#include "flashgg/DataFormats/interface/ZPlusJetTag.h"
#include "flashgg/DataFormats/interface/EXOTag.h"

#include "flashgg/Taggers/interface/CollectionDumper.h"

namespace flashgg {
    typedef CollectionDumper<std::vector<UntaggedTag>,
            UntaggedTag,
            CutBasedClassifier<UntaggedTag> > CutBasedUntaggedTagDumper;
    typedef CollectionDumper<std::vector<VBFTag>,
            VBFTag,
            CutBasedClassifier<VBFTag> > CutBasedVBFTagDumper;
    typedef CollectionDumper<std::vector<TTHLeptonicTag>,
            TTHLeptonicTag,
            CutBasedClassifier<TTHLeptonicTag> > CutBasedTTHLeptonicTagDumper;
    typedef CollectionDumper<std::vector<TTHHadronicTag>,
            TTHHadronicTag,
            CutBasedClassifier<TTHHadronicTag> > CutBasedTTHHadronicTagDumper;
    typedef CollectionDumper<std::vector<VHLooseTag>,
            VHLooseTag,
            CutBasedClassifier<VHLooseTag> > CutBasedVHLooseTagDumper;
    typedef CollectionDumper<std::vector<VHTightTag>,
            VHTightTag,
            CutBasedClassifier<VHTightTag> > CutBasedVHTightTagDumper;
    typedef CollectionDumper<std::vector<VHHadronicTag>,
            VHHadronicTag,
            CutBasedClassifier<VHHadronicTag> > CutBasedVHHadronicTagDumper;
    typedef CollectionDumper<std::vector<ZPlusJetTag>,
            ZPlusJetTag,
            CutBasedClassifier<ZPlusJetTag> > CutBasedZPlusJetTagDumper;
    typedef CollectionDumper<std::vector<EXOTag>,
            EXOTag,
            CutBasedClassifier<EXOTag> > CutBasedEXOTagDumper;
}

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
