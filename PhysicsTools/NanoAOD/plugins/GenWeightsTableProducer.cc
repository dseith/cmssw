#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/NanoAOD/interface/MergeableCounterTable.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/transform.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>

namespace {
    ///  ---- Cache object for running sums of weights ----
    struct Counter {
        Counter() : 
            num(0), sumw(0), sumw2(0), sumPDF(), sumScale(), sumRwgt(), sumNamed(), sumPS() {}

        // the counters
        long long num;
        long double sumw;
        long double sumw2;
        std::vector<long double> sumPDF, sumScale, sumRwgt, sumNamed, sumPS;

        void clear() { 
            num = 0; sumw = 0; sumw2 = 0;
            sumPDF.clear(); sumScale.clear(); sumRwgt.clear(); sumNamed.clear(), sumPS.clear();
        }

        // inc the counters
        void incGenOnly(double w) { 
            num++; sumw += w; sumw2 += (w*w); 
        }

        void incPSOnly(double w0, const std::vector<double> & wPS) {
            if (!wPS.empty()) {
                if (sumPS.empty()) sumPS.resize(wPS.size(), 0);
                for (unsigned int i = 0, n = wPS.size(); i < n; ++i) sumPS[i] += (w0 * wPS[i]);
            }
        }

        void incLHE(double w0, const std::vector<double> & wScale, const std::vector<double> & wPDF, const std::vector<double> & wRwgt, const std::vector<double> & wNamed, const std::vector<double> & wPS) {
            // add up weights
            incGenOnly(w0);
            // then add up variations
            if (!wScale.empty()) {
                if (sumScale.empty()) sumScale.resize(wScale.size(), 0);
                for (unsigned int i = 0, n = wScale.size(); i < n; ++i) sumScale[i] += (w0 * wScale[i]);
            }
            if (!wPDF.empty()) {
                if (sumPDF.empty()) sumPDF.resize(wPDF.size(), 0);
                for (unsigned int i = 0, n = wPDF.size(); i < n; ++i) sumPDF[i] += (w0 * wPDF[i]);
            }
            if (!wRwgt.empty()) {
                if (sumRwgt.empty()) sumRwgt.resize(wRwgt.size(), 0);
                for (unsigned int i = 0, n = wRwgt.size(); i < n; ++i) sumRwgt[i] += (w0 * wRwgt[i]);
            }
            if (!wNamed.empty()) {
                if (sumNamed.empty()) sumNamed.resize(wNamed.size(), 0);
                for (unsigned int i = 0, n = wNamed.size(); i < n; ++i) sumNamed[i] += (w0 * wNamed[i]);
            }
            incPSOnly(w0, wPS);
        }

        void merge(const Counter & other) { 
            num += other.num; sumw += other.sumw; sumw2 += other.sumw2; 
            if (sumScale.empty() && !other.sumScale.empty()) sumScale.resize(other.sumScale.size(),0);
            if (sumPDF.empty() && !other.sumPDF.empty()) sumPDF.resize(other.sumPDF.size(),0);
            if (sumRwgt.empty() && !other.sumRwgt.empty()) sumRwgt.resize(other.sumRwgt.size(),0);
            if (sumNamed.empty() && !other.sumNamed.empty()) sumNamed.resize(other.sumNamed.size(),0);
            if (sumPS.empty() && !other.sumPS.empty()) sumPS.resize(other.sumPS.size(),0);
            if (!other.sumScale.empty()) for (unsigned int i = 0, n = sumScale.size(); i < n; ++i) sumScale[i] += other.sumScale[i];
            if (!other.sumPDF.empty()) for (unsigned int i = 0, n = sumPDF.size(); i < n; ++i) sumPDF[i] += other.sumPDF[i];
            if (!other.sumRwgt.empty()) for (unsigned int i = 0, n = sumRwgt.size(); i < n; ++i) sumRwgt[i] += other.sumRwgt[i];
            if (!other.sumNamed.empty()) for (unsigned int i = 0, n = sumNamed.size(); i < n; ++i) sumNamed[i] += other.sumNamed[i];
            if (!other.sumPS.empty()) for (unsigned int i = 0, n = sumPS.size(); i < n; ++i) sumPS[i] += other.sumPS[i];
        }
    };

    ///  ---- RunCache object for dynamic choice of LHE IDs ----
    struct DynamicWeightChoice {
        // choice of LHE weights
        // ---- scale ----
        std::vector<std::string> scaleWeightIDs; 
        std::string scaleWeightsDoc;
        // ---- pdf ----
        std::vector<std::string> pdfWeightIDs; 
        std::string pdfWeightsDoc;
        // ---- rwgt ----
        std::vector<std::string> rwgtIDs;
        std::string rwgtWeightDoc;
    };

    float stof_fortrancomp(const std::string &str) {
        std::string::size_type match = str.find("d");
        if (match != std::string::npos) {
            std::string pre  = str.substr(0,match);
            std::string post = str.substr(match+1);
            return std::stof(pre) * std::pow(10.0f, std::stof(post));
        } else {
            std::cout << str << std::endl;
            return std::stof(str);
        }
    }
    ///  -------------- temporary objects --------------
    struct ScaleVarWeight {
        std::string wid, label;
        std::pair<float,float> scales;
        ScaleVarWeight(const std::string & id, const std::string & text, const std::string & muR, const std::string & muF) :
            wid(id), label(text), scales(stof_fortrancomp(muR), stof_fortrancomp(muF)) {}
        bool operator<(const ScaleVarWeight & other) { return (scales == other.scales ? wid < other.wid : scales < other.scales); }
    };
    struct PDFSetWeights {
        std::vector<std::string> wids;
        std::pair<unsigned int,unsigned int> lhaIDs;
        PDFSetWeights(const std::string & wid, unsigned int lhaID) : wids(1,wid), lhaIDs(lhaID,lhaID) {}
        bool operator<(const PDFSetWeights & other) const { return lhaIDs < other.lhaIDs; }
        void add(const std::string & wid, unsigned int lhaID) {
            wids.push_back(wid);
            lhaIDs.second = lhaID;
        }
        bool maybe_add(const std::string & wid, unsigned int lhaID) {
            if (lhaID == lhaIDs.second+1) {
                lhaIDs.second++;
                wids.push_back(wid);
                return true;
            } else {
                return false;
            }
        }
    };
}

class GenWeightsTableProducer : public edm::global::EDProducer<edm::StreamCache<Counter>, edm::RunCache<DynamicWeightChoice>, edm::RunSummaryCache<Counter>, edm::EndRunProducer> {
    public:
        GenWeightsTableProducer( edm::ParameterSet const & params ) :
            genTag_(consumes<GenEventInfoProduct>(params.getParameter<edm::InputTag>("genEvent"))),
            lheLabel_(params.getParameter<std::vector<edm::InputTag>>("lheInfo")),
            lheTag_(edm::vector_transform(lheLabel_, [this](const edm::InputTag & tag) { return mayConsume<LHEEventProduct>(tag); })),
            lheRunTag_(edm::vector_transform(lheLabel_, [this](const edm::InputTag & tag) { return mayConsume<LHERunInfoProduct, edm::InRun>(tag); })),
            namedWeightIDs_(params.getParameter<std::vector<std::string>>("namedWeightIDs")),
            namedWeightLabels_(params.getParameter<std::vector<std::string>>("namedWeightLabels")),
            lheWeightPrecision_(params.getParameter<int32_t>("lheWeightPrecision")),
            maxPdfWeights_(params.getParameter<uint32_t>("maxPdfWeights")),
            debug_(params.getUntrackedParameter<bool>("debug",false)), debugRun_(debug_.load()),
            hasIssuedWarning_(false)
        {
            produces<nanoaod::FlatTable>();
            produces<nanoaod::FlatTable>("LHEScale");
            produces<nanoaod::FlatTable>("LHEPdf");
            produces<nanoaod::FlatTable>("LHEReweighting");
            produces<nanoaod::FlatTable>("LHENamed");
            produces<nanoaod::FlatTable>("PS");
            produces<nanoaod::FlatTable>("EFTWeights");
            produces<nanoaod::MergeableCounterTable,edm::Transition::EndRun>();
            if (namedWeightIDs_.size() != namedWeightLabels_.size()) {
                throw cms::Exception("Configuration", "Size mismatch between namedWeightIDs & namedWeightLabels");
            }
            for (const edm::ParameterSet & pdfps : params.getParameter<std::vector<edm::ParameterSet>>("preferredPDFs")) {
                const std::string & name = pdfps.getParameter<std::string>("name");
                uint32_t lhaid = pdfps.getParameter<uint32_t>("lhaid");
                preferredPDFLHAIDs_.push_back(lhaid);
                lhaNameToID_[name] = lhaid;
                lhaNameToID_[name+".LHgrid"] = lhaid;
            }
        }

        ~GenWeightsTableProducer() override {}

        void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override {
            // get my counter for weights
            Counter * counter = streamCache(id);

            // generator information (always available)
            edm::Handle<GenEventInfoProduct> genInfo;
            iEvent.getByToken(genTag_, genInfo);
            double weight = genInfo->weight();

            // table for gen info, always available
            auto out = std::make_unique<nanoaod::FlatTable>(1, "genWeight", true);
            out->setDoc("generator weight");
            out->addColumnValue<float>("", weight, "generator weight", nanoaod::FlatTable::FloatColumn);
            iEvent.put(std::move(out));

            // tables for LHE weights, may not be filled
            std::unique_ptr<nanoaod::FlatTable> lheScaleTab, lhePdfTab, lheRwgtTab, lheNamedTab, eftWeightTab;
            std::unique_ptr<nanoaod::FlatTable> genPSTab;

            edm::Handle<LHEEventProduct> lheInfo;
            for (const auto & lheTag: lheTag_) {
                iEvent.getByToken(lheTag, lheInfo);
                if (lheInfo.isValid()) {
                    break;
                }
            }
            if (lheInfo.isValid()) {
                // get the dynamic choice of weights
                const DynamicWeightChoice * weightChoice = runCache(iEvent.getRun().index());
                // go fill tables
                fillLHEWeightTables(counter, weightChoice, weight, *lheInfo, *genInfo, lheScaleTab, lhePdfTab, lheRwgtTab, lheNamedTab, genPSTab, eftWeightTab);
            } else {
                // Still try to add the PS weights
                fillOnlyPSWeightTable(counter, weight, *genInfo, genPSTab);
                // make dummy values 
                lheScaleTab.reset(new nanoaod::FlatTable(1, "LHEScaleWeights", true));
                lhePdfTab.reset(new nanoaod::FlatTable(1, "LHEPdfWeights", true));
                lheRwgtTab.reset(new nanoaod::FlatTable(1, "LHEReweightingWeights", true));
                lheNamedTab.reset(new nanoaod::FlatTable(1, "LHENamedWeights", true));
                eftWeightTab.reset(new nanoaod::FlatTable(1, "EFTWeights", true));
                if (!hasIssuedWarning_.exchange(true)) {
                    edm::LogWarning("LHETablesProducer") << "No LHEEventProduct, so there will be no LHE Tables\n";
                }
            }

            iEvent.put(std::move(lheScaleTab), "LHEScale");
            iEvent.put(std::move(lhePdfTab), "LHEPdf");
            iEvent.put(std::move(lheRwgtTab), "LHEReweighting");
            iEvent.put(std::move(lheNamedTab), "LHENamed");
            iEvent.put(std::move(genPSTab), "PS");
            iEvent.put(std::move(eftWeightTab), "EFTWeights");
            
        }

        void fillLHEWeightTables(
                Counter * counter,
                const DynamicWeightChoice * weightChoice,
                double genWeight,
                const LHEEventProduct & lheProd, 
                const GenEventInfoProduct & genProd,
                std::unique_ptr<nanoaod::FlatTable> & outScale, 
                std::unique_ptr<nanoaod::FlatTable> & outPdf,
                std::unique_ptr<nanoaod::FlatTable> & outRwgt,
                std::unique_ptr<nanoaod::FlatTable> & outNamed,
                std::unique_ptr<nanoaod::FlatTable> & outPS,
                std::unique_ptr<nanoaod::FlatTable> & outEFT
                ) const
        {
            bool lheDebug = debug_.exchange(false); // make sure only the first thread dumps out this (even if may still be mixed up with other output, but nevermind)

            const std::vector<std::string> & scaleWeightIDs = weightChoice->scaleWeightIDs;
            const std::vector<std::string> & pdfWeightIDs   = weightChoice->pdfWeightIDs;
            const std::vector<std::string> & rwgtWeightIDs  = weightChoice->rwgtIDs;

            // std::cout << "====================" << std::endl;
            // for(auto const& w : scaleWeightIDs)
            //     std::cout << w << " "; 
            // std::cout << std::endl;
            // for(auto const& w : pdfWeightIDs)
            //     std::cout << w << " "; 
            // std::cout << std::endl;
            // for(auto const& w : rwgtWeightIDs)
            //     std::cout << w << " "; 
            // std::cout << std::endl;

            // std::cout << "====================" << std::endl;

            double w0 = lheProd.originalXWGTUP();

            std::vector<double> wScale(scaleWeightIDs.size(), 1), wPDF(pdfWeightIDs.size(), 1), wRwgt(rwgtWeightIDs.size(), 1), wNamed(namedWeightIDs_.size(), 1);
            for (auto & weight : lheProd.weights()) {
                if (lheDebug) printf("Weight  %+9.5f   rel %+9.5f   for id %s\n", weight.wgt, weight.wgt/w0,  weight.id.c_str());
                // now we do it slowly, can be optimized
                auto mScale = std::find(scaleWeightIDs.begin(), scaleWeightIDs.end(), weight.id);
                if (mScale != scaleWeightIDs.end()) wScale[mScale-scaleWeightIDs.begin()] = weight.wgt/w0;

                auto mPDF = std::find(pdfWeightIDs.begin(), pdfWeightIDs.end(), weight.id);
                if (mPDF != pdfWeightIDs.end()) wPDF[mPDF-pdfWeightIDs.begin()] = weight.wgt/w0;

                auto mRwgt = std::find(rwgtWeightIDs.begin(), rwgtWeightIDs.end(), weight.id);
                if (mRwgt != rwgtWeightIDs.end()) wRwgt[mRwgt-rwgtWeightIDs.begin()] = weight.wgt/w0;

                auto mNamed = std::find(namedWeightIDs_.begin(), namedWeightIDs_.end(), weight.id);
                if (mNamed != namedWeightIDs_.end()) wNamed[mNamed-namedWeightIDs_.begin()] = weight.wgt/w0;
            } 

            int vectorSize = (genProd.weights().size() == 14 || genProd.weights().size() == 46) ? 4 : 1;
            std::vector<double> wPS(vectorSize, 1);
            if (vectorSize > 1 ) {
                for (unsigned int i=6; i<10; i++){
                    wPS[i-6] = (genProd.weights()[i])/w0;
                }
            }
            outPS.reset(new nanoaod::FlatTable(wPS.size(), "PSWeight", false));
            outPS->addColumn<float>("", wPS, vectorSize > 1 ? "PS weights (w_var / w_nominal); [0] is ISR=0.5 FSR=1; [1] is ISR=1 FSR=0.5; [2] is ISR=2 FSR=1; [3] is ISR=1 FSR=2 " : "dummy PS weight (1.0) ", nanoaod::FlatTable::FloatColumn, lheWeightPrecision_);

            outScale.reset(new nanoaod::FlatTable(wScale.size(), "LHEScaleWeight", false));
            outScale->addColumn<float>("", wScale, weightChoice->scaleWeightsDoc, nanoaod::FlatTable::FloatColumn, lheWeightPrecision_); 

            outPdf.reset(new nanoaod::FlatTable(wPDF.size(), "LHEPdfWeight", false));
            outPdf->addColumn<float>("", wPDF, weightChoice->pdfWeightsDoc, nanoaod::FlatTable::FloatColumn, lheWeightPrecision_);

            outRwgt.reset(new nanoaod::FlatTable(wRwgt.size(), "LHEReweightingWeight", false));
            outRwgt->addColumn<float>("", wRwgt, weightChoice->rwgtWeightDoc, nanoaod::FlatTable::FloatColumn, lheWeightPrecision_);

            outNamed.reset(new nanoaod::FlatTable(1, "LHEWeight", true));
            outNamed->addColumnValue<float>("originalXWGTUP", lheProd.originalXWGTUP(), "Nominal event weight in the LHE file", nanoaod::FlatTable::FloatColumn);
            for (unsigned int i = 0, n = wNamed.size(); i < n; ++i) {
                outNamed->addColumnValue<float>(namedWeightLabels_[i], wNamed[i], "LHE weight for id "+namedWeightIDs_[i]+", relative to nominal", nanoaod::FlatTable::FloatColumn, lheWeightPrecision_);
            }

            std::vector<std::string> eft_weight_names = {
                    "SM",
                    "ctW_1p87_ctWI_-0p24_cptb_-19p7_cptbI_16p44",
                    "ctW_1p76_ctWI_0p33_cptb_6p86_cptbI_-16p64",
                    "ctW_1p07_ctWI_-1p05_cptb_-18p77_cptbI_11p55",
                    "ctW_-0p62_ctWI_0p49_cptb_4p63_cptbI_-14p06",
                    "ctW_-1p27_ctWI_-1p54_cptb_-19p42_cptbI_-0p53",
                    "ctW_1p86_ctWI_-1p74_cptb_1p64_cptbI_-1p36",
                    "ctW_0p41_ctWI_-1p64_cptb_3p16_cptbI_-9p22",
                    "ctW_0p23_ctWI_0p58_cptb_-0p76_cptbI_-5p79",
                    "ctW_-1p0_ctWI_1p73_cptb_-1p86_cptbI_1p21",
                    "ctW_-1p92_ctWI_0p03_cptb_-19p77_cptbI_-14p25",
                    "ctW_-0p11_ctWI_-0p49_cptb_-17p83_cptbI_3p5",
                    "ctW_-1p34_ctWI_0p23_cptb_-14p23_cptbI_17p49",
                    "ctW_1p08_ctWI_1p83_cptb_-14p35_cptbI_-7p78",
                    "ctW_-1p84_ctWI_-0p89_cptb_12p26_cptbI_-12p91",
                    "ctW_-1p38_ctWI_1p82_cptb_-13p82_cptbI_13p36",
                    "ctW_-1p84_ctWI_-0p46_cptb_-6p02_cptbI_-6p33",
                    "ctW_1p27_ctWI_-0p1_cptb_11p32_cptbI_-1p17",
                    "ctW_1p27_ctWI_1p53_cptb_-2p42_cptbI_11p24",
                    "ctW_1p26_ctWI_-0p82_cptb_-15p04_cptbI_-12p58",
                    "ctW_-0p26_ctWI_-1p52_cptb_1p19_cptbI_13p18",
                    "ctW_-0p06_ctWI_1p27_cptb_6p26_cptbI_5p64",
                    "ctW_-0p62_ctWI_0p81_cptb_12p4_cptbI_-13p71",
                    "ctW_1p63_ctWI_-0p92_cptb_-13p81_cptbI_13p62",
                    "ctW_0p88_ctWI_1p17_cptb_-2p13_cptbI_-17p17",
                    "ctW_-0p42_ctWI_-1p81_cptb_-8p55_cptbI_-18p48",
                    "ctW_0p03_ctWI_-1p65_cptb_17p31_cptbI_7p98",
                    "ctW_-0p73_ctWI_1p78_cptb_-17p35_cptbI_-9p7",
                    "ctW_-1p71_ctWI_-0p29_cptb_-11p92_cptbI_-4p13",
                    "ctW_0p82_ctWI_1p55_cptb_0p02_cptbI_12p74",
                    "ctW_-0p55_ctWI_1p43_cptb_0p6_cptbI_8p1",
                    "ctW_-1p3_ctWI_0p34_cptb_-7p93_cptbI_12p51",
                    "ctW_0p14_ctWI_-0p0_cptb_10p93_cptbI_1p96",
                    "ctW_-0p66_ctWI_-1p48_cptb_4p98_cptbI_17p01",
                    "ctW_1p37_ctWI_-1p72_cptb_-6p99_cptbI_-19p93",
                    "ctW_0p7_ctWI_0p55_cptb_10p32_cptbI_-14p07",
                    "ctW_-1p13_ctWI_-0p27_cptb_9p47_cptbI_-11p72",
                    "ctW_1p29_ctWI_-0p48_cptb_14p97_cptbI_18p4",
                    "ctW_0p15_ctWI_1p68_cptb_-3p12_cptbI_7p41",
                    "ctW_1p39_ctWI_1p36_cptb_-16p29_cptbI_-9p59",
                    "ctW_-0p36_ctWI_1p43_cptb_-8p93_cptbI_-15p45",
                    "ctW_-0p5_ctWI_-1p14_cptb_8p65_cptbI_3p5",
                    "ctW_-1p79_ctWI_1p51_cptb_-2p38_cptbI_11p08",
                    "ctW_-1p34_ctWI_-0p77_cptb_-17p07_cptbI_-4p2",
                    "ctW_-0p28_ctWI_0p73_cptb_-0p82_cptbI_-3p98",
                    "ctW_-0p09_ctWI_-0p87_cptb_-5p22_cptbI_4p02",
                    "ctW_-1p18_ctWI_0p08_cptb_17p19_cptbI_-1p7",
                    "ctW_0p67_ctWI_0p64_cptb_12p05_cptbI_9p02",
                    "ctW_0p64_ctWI_-1p59_cptb_-5p33_cptbI_11p76",
                    "ctW_-1p74_ctWI_-1p44_cptb_16p21_cptbI_3p52",
                    "ctW_1p92_ctWI_1p01_cptb_19p83_cptbI_-8p06",
                    "ctW_0p74_ctWI_-0p68_cptb_17p84_cptbI_-8p25",
                    "ctW_-1p22_ctWI_1p26_cptb_-8p84_cptbI_1p36",
                    "ctW_0p41_ctWI_-0p64_cptb_-6p63_cptbI_-7p16",
                    "ctW_-0p93_ctWI_1p04_cptb_3p57_cptbI_1p93",
                    "ctW_1p3_ctWI_-0p97_cptb_-14p52_cptbI_-18p86",
                    "ctW_-1p85_ctWI_1p96_cptb_18p95_cptbI_6p66",
                    "ctW_1p96_ctWI_0p06_cptb_-2p62_cptbI_0p61",
                    "ctW_-1p72_ctWI_1p04_cptb_-7p79_cptbI_7p92",
                    "ctW_1p64_ctWI_-1p12_cptb_0p18_cptbI_6p13",
                    "ctW_1p24_ctWI_0p73_cptb_-6p76_cptbI_2p57",
                    "ctW_0p125_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p25_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p375_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p5_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p625_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p75_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p875_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_1p0_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_1p125_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_1p25_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_1p375_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_1p5_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_1p625_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_1p75_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_1p875_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_2p0_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_2p125_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_2p25_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_2p375_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_2p5_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-0p125_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-0p25_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-0p375_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-0p5_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-0p625_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-0p75_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-0p875_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-1p0_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-1p125_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-1p25_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-1p375_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-1p5_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-1p625_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-1p75_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-1p875_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-2p0_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-2p125_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-2p25_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-2p375_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_-2p5_ctWI_0p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p125_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p25_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p375_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p5_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p625_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p75_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p875_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_1p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_1p125_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_1p25_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_1p375_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_1p5_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_1p625_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_1p75_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_1p875_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_2p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_2p125_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_2p25_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_2p375_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_2p5_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-0p125_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-0p25_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-0p375_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-0p5_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-0p625_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-0p75_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-0p875_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-1p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-1p125_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-1p25_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-1p375_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-1p5_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-1p625_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-1p75_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-1p875_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-2p0_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-2p125_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-2p25_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-2p375_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_-2p5_cptb_0p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_1p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_2p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_3p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_5p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_6p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_7p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_8p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_10p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_11p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_12p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_13p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_15p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_16p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_17p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_18p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_20p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_21p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_22p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_23p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_25p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-1p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-2p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-3p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-5p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-6p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-7p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-8p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-10p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-11p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-12p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-13p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-15p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-16p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-17p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-18p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-20p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-21p25_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-22p5_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-23p75_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_-25p0_cptbI_0p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_1p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_2p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_3p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_5p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_6p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_7p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_8p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_10p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_11p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_12p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_13p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_15p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_16p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_17p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_18p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_20p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_21p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_22p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_23p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_25p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-1p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-2p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-3p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-5p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-6p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-7p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-8p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-10p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-11p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-12p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-13p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-15p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-16p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-17p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-18p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-20p0",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-21p25",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-22p5",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-23p75",
                    "ctW_0p0_ctWI_0p0_cptb_0p0_cptbI_-25p0",
                   };


            std::vector<double> wEFT(eft_weight_names.size(), 1);
            std::string description = "(";

            auto lheweights = lheProd.weights();
            for(unsigned int index=0; index<eft_weight_names.size(); ++index){
                auto name = eft_weight_names[index];
                auto pos = std::find_if(lheweights.begin(), lheweights.end(), [&name](const gen::WeightsInfo x){return x.id == name;});
                auto weight = *pos;
                wEFT[index] = weight.wgt / w0;

                description += name + " [" + std::to_string(index) + "]; " ;
            }

            description += ")";

            
            outEFT.reset(new nanoaod::FlatTable(wEFT.size(), "EFTWeights", false));
            outEFT->addColumn<float>("", wEFT, "EFT Weights " + description, nanoaod::FlatTable::FloatColumn, lheWeightPrecision_);
            
            counter->incLHE(genWeight, wScale, wPDF, wRwgt, wNamed, wPS);
        }

        void fillOnlyPSWeightTable(
                Counter * counter,
                double genWeight,
                const GenEventInfoProduct & genProd,
                std::unique_ptr<nanoaod::FlatTable> & outPS ) const
        {
            int vectorSize = (genProd.weights().size() == 14 || genProd.weights().size() == 46) ? 4 : 1;


            std::vector<double> wPS(vectorSize, 1);
            if (vectorSize > 1 ){
                for (unsigned int i=6; i<10; i++){
                    wPS[i-6] = (genProd.weights()[i])/genWeight;
                }
            }

            outPS.reset(new nanoaod::FlatTable(wPS.size(), "PSWeight", false));
            outPS->addColumn<float>("", wPS, vectorSize > 1 ? "PS weights (w_var / w_nominal); [0] is ISR=0.5 FSR=1; [1] is ISR=1 FSR=0.5; [2] is ISR=2 FSR=1; [3] is ISR=1 FSR=2 " : "dummy PS weight (1.0) " , nanoaod::FlatTable::FloatColumn, lheWeightPrecision_);

            counter->incGenOnly(genWeight);
            counter->incPSOnly(genWeight,wPS);
        }



        // create an empty counter
        std::shared_ptr<DynamicWeightChoice> globalBeginRun(edm::Run const& iRun, edm::EventSetup const&) const override {
            edm::Handle<LHERunInfoProduct> lheInfo;

            bool lheDebug = debugRun_.exchange(false); // make sure only the first thread dumps out this (even if may still be mixed up with other output, but nevermind)
            // lheDebug=true;
            auto weightChoice = std::make_shared<DynamicWeightChoice>();

            // getByToken throws since we're not in the endRun (see https://github.com/cms-sw/cmssw/pull/18499)
            //if (iRun.getByToken(lheRunTag_, lheInfo)) {
            for (const auto & lheLabel: lheLabel_) {
                iRun.getByLabel(lheLabel, lheInfo);
                if (lheInfo.isValid()) {
                    break;
                }
            }

            if (lheInfo.isValid()) {
                std::vector<ScaleVarWeight> scaleVariationIDs;
                std::vector<PDFSetWeights>  pdfSetWeightIDs;
                std::vector<std::string>    lheReweighingIDs;
                
                std::regex weightgroupmg26x("<weightgroup\\s+(?:name|type)=\"(.*)\"\\s+combine=\"(.*)\"\\s*>");
                std::regex weightgroup("<weightgroup\\s+combine=\"(.*)\"\\s+(?:name|type)=\"(.*)\"\\s*>");
                std::regex weightgroupRwgt("<weightgroup\\s+(?:name|type)=\"(.*)\"\\s*>");
                std::regex endweightgroup("</weightgroup>");
                std::regex scalewmg26x("<weight\\s+(?:.*\\s+)?id=\"(\\d+)\"\\s*(?:lhapdf=\\d+|dyn=\\s*-?\\d+)?\\s*((?:[mM][uU][rR]|renscfact)=\"(\\S+)\"\\s+(?:[mM][uU][Ff]|facscfact)=\"(\\S+)\")(\\s+.*)?</weight>");
                // std::regex scalew("<weight\\s+(?:.*\\s+)?id=\"(\\d+)\">\\s*(?:lhapdf=\\d+|dyn=\\s*-?\\d+)?\\s*((?:mu[rR]|renscfact)=(\\S+)\\s+(?:mu[Ff]|facscfact)=(\\S+)(\\s+.*)?)</weight>");
                std::regex scalew("<weight\\s*MUF=\"(\\S+)\"\\s*MUR=\"(\\S+)\"\\s*PDF=\"(\\d+)\"\\s*id=\"(\\d+)\">\\s*MUR=(?:\\S+)\\s*MUF=(?:\\S+)\\s*</weight>");
                std::regex pdfw("<weight\\s+id=\"(\\d+)\">\\s*(?:PDF set|lhapdf|PDF|pdfset)\\s*=\\s*(\\d+)\\s*(?:\\s.*)?</weight>");
                // std::regex pdfwOld("<weight\\s+(?:.*\\s+)?id=\"(\\d+)\">\\s*Member \\s*(\\d+)\\s*(?:.*)</weight>");
                // std::regex pdfwOld("<weight\\s+(?:.*\\s+)?id=\"(\\d+)\">\\s*Member \\s*(\\d+)\\s*(?:.*)</weight>");
                std::regex pdfwOld("<weight\\s+MUF=\"(?:\\S+)\"\\s+MUR=\"(?:\\S+)\"\\s*PDF=\\s*\"(\\d+)\"\\s*id=\"(\\d+)\"\\s*>\\s*(?:PDF=(\\d+)\\s*MemberID=(\\d+))?\\s*(?:\\s.*)?</weight>");

                // std::regex pdfwOld("<weight\\s+id=\"(\\d+)\"\\s*MUF=\"(?:\\S+)\"\\s*MUR=\"(?:\\S+)\"\\s*(?:PDF set|lhapdf|PDF|pdfset)\\s*=\\s*\"(\\d+)\"\\s*>\\s*(?:PDF=(\\d+)\\s*MemberID=(\\d+))?\\s*(?:\\s.*)?</weight>");
                std::regex pdfwmg26x("<weight\\s+id=\"(\\d+)\"\\s*MUR=\"(?:\\S+)\"\\s*MUF=\"(?:\\S+)\"\\s*(?:PDF set|lhapdf|PDF|pdfset)\\s*=\\s*\"(\\d+)\"\\s*>\\s*(?:PDF=(\\d+)\\s*MemberID=(\\d+))?\\s*(?:\\s.*)?</weight>");
                std::regex rwgt("<weight\\s+id=\"(.+)\">(.+)?(</weight>)?");
                std::smatch groups;
                for (auto iter=lheInfo->headers_begin(), end = lheInfo->headers_end(); iter != end; ++iter) {
                    if (iter->tag() != "initrwgt") {
                        if (lheDebug) std::cout << "Skipping LHE header with tag" << iter->tag() << std::endl;
                        continue;
                    }
                    if (lheDebug) std::cout << "Found LHE header with tag" << iter->tag() << std::endl;
                    std::vector<std::string>  lines = iter->lines();
                    bool missed_weightgroup=false; //Needed because in some of the samples ( produced with MG26X ) a small part of the header info is ordered incorrectly
                    bool ismg26x=false;
                    for (unsigned int iLine = 0, nLines = lines.size(); iLine < nLines; ++iLine) { //First start looping through the lines to see which weightgroup pattern is matched
                        boost::replace_all(lines[iLine],"&lt;", "<");
                        boost::replace_all(lines[iLine],"&gt;", ">");
                        if(std::regex_search(lines[iLine],groups,weightgroupmg26x)){
                            ismg26x=true;
                        }
                    }
                    for (unsigned int iLine = 0, nLines = lines.size(); iLine < nLines; ++iLine) {
                        if (lheDebug) std::cout << lines[iLine];
                        if (std::regex_search(lines[iLine], groups, ismg26x ? weightgroupmg26x : weightgroup) ) {
                            std::string groupname = groups.str(2);
                            if (ismg26x) groupname = groups.str(1);
                            if (lheDebug) std::cout << ">>> Looks like the beginning of a weight group for '" << groupname << "'" << std::endl;
                            if (groupname.find("scale_variation") == 0 || groupname == "Central scale variation") {
                                if (lheDebug) std::cout << ">>> Looks like scale variation for theory uncertainties" << std::endl;
                                for ( ++iLine; iLine < nLines; ++iLine) {
                                    if (lheDebug) std::cout << "    " << lines[iLine];
                                    if (std::regex_search(lines[iLine], groups, ismg26x ? scalewmg26x : scalew)) {
                                        if (lheDebug) std::cout << "    >>> Scale weight " << groups[1].str() << " for " << groups[3].str() << " , " << groups[4].str() << " , " << groups[5].str() << std::endl;
                                        scaleVariationIDs.emplace_back(groups.str(1), groups.str(2), groups.str(3), groups.str(4));
                                    } else if (std::regex_search(lines[iLine], endweightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the end of a weight group" << std::endl;
                                        if (!missed_weightgroup){
                                            break;
                                        } else missed_weightgroup=false;
                                    } else if (std::regex_search(lines[iLine], ismg26x ? weightgroupmg26x : weightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the beginning of a new weight group, I will assume I missed the end of the group." << std::endl;
                                        if (ismg26x) missed_weightgroup=true;
                                        --iLine; // rewind by one, and go back to the outer loop
                                        break;
                                    }
                                }
                            } else if (groupname == "PDF_variation" || groupname.find("PDF_variation ") == 0) {
                                if (lheDebug) std::cout << ">>> Looks like a new-style block of PDF weights for one or more pdfs" << std::endl;
                                for ( ++iLine; iLine < nLines; ++iLine) {
                                    if (lheDebug) std::cout << "    " << lines[iLine];
                                    if (std::regex_search(lines[iLine], groups, pdfw)) {
                                        unsigned int lhaID = std::stoi(groups.str(2));
                                        if (lheDebug) std::cout << "    >>> PDF weight " << groups.str(1) << " for " << groups.str(2) << " = " << lhaID << std::endl;
                                        if (pdfSetWeightIDs.empty() || ! pdfSetWeightIDs.back().maybe_add(groups.str(1),lhaID)) {
                                            pdfSetWeightIDs.emplace_back(groups.str(1),lhaID);
                                        }
                                    } else if (std::regex_search(lines[iLine], endweightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the end of a weight group" << std::endl;
                                        if (!missed_weightgroup){ 
                                            break;
                                        } else missed_weightgroup=false;
                                    } else if (std::regex_search(lines[iLine], ismg26x ? weightgroupmg26x : weightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the beginning of a new weight group, I will assume I missed the end of the group." << std::endl;
                                        if (ismg26x) missed_weightgroup=true;
                                        --iLine; // rewind by one, and go back to the outer loop
                                        break;
                                    }
                                }
                            } else if (groupname == "PDF_variation1" || groupname == "PDF_variation2") { 
                                if (lheDebug) std::cout << ">>> Looks like a new-style block of PDF weights for multiple pdfs" << std::endl;
                                unsigned int lastid = 0;
                                for ( ++iLine; iLine < nLines; ++iLine) {
                                    if (lheDebug) std::cout << "    " << lines[iLine];
                                    if (std::regex_search(lines[iLine], groups, pdfw)) {
                                        unsigned int id = std::stoi(groups.str(1));
                                        unsigned int lhaID = std::stoi(groups.str(2));
                                        if (lheDebug) std::cout << "    >>> PDF weight " << groups.str(1) << " for " << groups.str(2) << " = " << lhaID << std::endl;
                                        if (id != (lastid+1) || pdfSetWeightIDs.empty()) {
                                            pdfSetWeightIDs.emplace_back(groups.str(1),lhaID);
                                        } else {
                                            pdfSetWeightIDs.back().add(groups.str(1),lhaID);
                                        }
                                        lastid = id;
                                    } else if (std::regex_search(lines[iLine], endweightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the end of a weight group" << std::endl;
                                        if(!missed_weightgroup) {
                                            break;
                                        } else missed_weightgroup=false;
                                    } else if (std::regex_search(lines[iLine], ismg26x ? weightgroupmg26x : weightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the beginning of a new weight group, I will assume I missed the end of the group." << std::endl;
                                        if (ismg26x) missed_weightgroup=true;
                                        --iLine; // rewind by one, and go back to the outer loop
                                        break;
                                    }
                                }
                            } else if (lhaNameToID_.find(groupname) != lhaNameToID_.end()) {
                                if (lheDebug) std::cout << ">>> Looks like an old-style PDF weight for an individual pdf" << std::endl;
                                unsigned int firstLhaID = lhaNameToID_.find(groupname)->second;
                                bool first = true;
                                for ( ++iLine; iLine < nLines; ++iLine) {
                                    if (lheDebug) std::cout << "    " << lines[iLine];
                                    if (std::regex_search(lines[iLine], groups, ismg26x ? pdfwmg26x : pdfwOld)) {
                                        unsigned int member = 0;
                                        if (ismg26x==0){
                                            member = std::stoi(groups.str(2));
                                        } else {
                                            if (!groups.str(4).empty()){
                                                member = std::stoi(groups.str(4));
                                             }
                                        }
                                        unsigned int lhaID = member+firstLhaID;
                                        if (lheDebug) std::cout << "    >>> PDF weight " << groups.str(1) << " for " << member << " = " << lhaID << std::endl;
                                        //if (member == 0) continue; // let's keep also the central value for now
                                        if (first) {
                                            pdfSetWeightIDs.emplace_back(groups.str(2),lhaID);
                                            first = false;
                                        } else {
                                            pdfSetWeightIDs.back().add(groups.str(2),lhaID);
                                        }
                                    } else if (std::regex_search(lines[iLine], endweightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the end of a weight group" << std::endl;
                                        if (!missed_weightgroup) {
                                            break;
                                        } else missed_weightgroup=false;
                                    } else if (std::regex_search(lines[iLine], ismg26x ? weightgroupmg26x : weightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the beginning of a new weight group, I will assume I missed the end of the group." << std::endl;
                                        if (ismg26x) missed_weightgroup=true;
                                        --iLine; // rewind by one, and go back to the outer loop
                                        break;
                                    }
                                }
                            } else {
                                for ( ++iLine; iLine < nLines; ++iLine) {
                                    if (lheDebug) std::cout << "    " << lines[iLine];
                                    if (std::regex_search(lines[iLine], groups, endweightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the end of a weight group" << std::endl;
                                        if (!missed_weightgroup){
                                            break;
                                        } else missed_weightgroup=false;
                                    } else if (std::regex_search(lines[iLine], ismg26x ? weightgroupmg26x : weightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the beginning of a new weight group, I will assume I missed the end of the group." << std::endl;
                                        if (ismg26x) missed_weightgroup=true;
                                        --iLine; // rewind by one, and go back to the outer loop
                                        break;
                                    }
                                }
                            }
                        } else if(std::regex_search(lines[iLine], groups, weightgroupRwgt) ) {
                            std::string groupname = groups.str(1);
                            if (groupname == "mg_reweighting") {
                                if (lheDebug) std::cout << ">>> Looks like a LHE weights for reweighting" << std::endl;
                                for ( ++iLine; iLine < nLines; ++iLine) {
                                    if (lheDebug) std::cout << "    " << lines[iLine];
                                    if (std::regex_search(lines[iLine], groups, rwgt)) {
                                        std::string rwgtID = groups.str(1);
                                        if (lheDebug) std::cout << "    >>> LHE reweighting weight: " << rwgtID << std::endl;
                                        if (std::find(lheReweighingIDs.begin(), lheReweighingIDs.end(), rwgtID) == lheReweighingIDs.end()) {
                                            // we're only interested in the beggining of the block
                                            lheReweighingIDs.emplace_back(rwgtID);
                                        }
                                    } else if (std::regex_search(lines[iLine], endweightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the end of a weight group" << std::endl;
                                        if (!missed_weightgroup){
                                            break;
                                        } else missed_weightgroup=false;
                                    } else if (std::regex_search(lines[iLine], ismg26x ? weightgroupmg26x : weightgroup)) {
                                        if (lheDebug) std::cout << ">>> Looks like the beginning of a new weight group, I will assume I missed the end of the group." << std::endl;
                                        if (ismg26x) missed_weightgroup=true;
                                        --iLine; // rewind by one, and go back to the outer loop
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    //std::cout << "============= END [ " << iter->tag() << " ] ============ \n\n" << std::endl;

                    // ----- SCALE VARIATIONS -----
                    // std::cout << scaleVariationIDs.size() << std::endl;
                    // if(scaleVariationIDs.size() == 0){
                    //     std::vector<float> ids = {2, 7, 12, 17, 26, 31, 36, 41};
                    //     std::vector<float> muFs = {0.5, 1.0, 2.0, 0.5, 2.0, 0.5, 1.0, 2.0};
                    //     std::vector<float> muRs = {0.5, 0.5, 0.5, 1.0, 1.0, 2.0, 2.0, 2.0};
                    //     for(unsigned int idx=0; idx<ids.size(); idx++){
                    //         auto id = ids[idx];
                    //         auto muF = muFs[idx];
                    //         auto muR = muRs[idx];
                    //         std::cout << id << " " << muF << " " << muR << std::endl;
                    //         scaleVariationIDs.emplace_back(std::to_string(id), "label", std::to_string(muR), std::to_string(muF));
                    //     }
                    // }
                    std::sort(scaleVariationIDs.begin(), scaleVariationIDs.end());
                    if (lheDebug) std::cout << "Found " << scaleVariationIDs.size() << " scale variations: " << std::endl;
                    std::stringstream scaleDoc; scaleDoc << "LHE scale variation weights (w_var / w_nominal); ";
                    for (unsigned int isw = 0, nsw = scaleVariationIDs.size(); isw < nsw; ++isw) {
                        const auto & sw = scaleVariationIDs[isw];
                        if (isw) scaleDoc << "; ";
                        scaleDoc << "[" << isw << "] is " << sw.label;
                        weightChoice->scaleWeightIDs.push_back(sw.wid);
                        if (lheDebug) printf("    id %s: scales ren = % .2f  fact = % .2f  text = %s\n", sw.wid.c_str(), sw.scales.first, sw.scales.second, sw.label.c_str());
                    }
                    if (!scaleVariationIDs.empty()) weightChoice->scaleWeightsDoc = scaleDoc.str();

                    // ------ PDF VARIATIONS (take the preferred one) -----
                    // if(pdfSetWeightIDs.size() == 0){
                    //     bool first = true;
                    //     for(unsigned int i=0; i<103; ++i){

                    //         if(first){
                    //             pdfSetWeightIDs.emplace_back(std::to_string(i+46), 306000+i);
                    //             first=false;
                    //         }
                    //         else{
                    //             pdfSetWeightIDs.back().add(std::to_string(i+46), 306000+i);
                    //         }
                    //     }
                    // }
                    if (lheDebug) {
                        std::cout << "Found " << pdfSetWeightIDs.size() << " PDF set errors: " << std::endl;
                        for (const auto & pw : pdfSetWeightIDs) printf("lhaIDs %6d - %6d (%3lu weights: %s, ... )\n", pw.lhaIDs.first, pw.lhaIDs.second, pw.wids.size(), pw.wids.front().c_str());
                    }

                    // ------ LHE REWEIGHTING -------
                    if (lheDebug) {
                        std::cout << "Found " << lheReweighingIDs.size() << " reweighting weights" << std::endl;
                    }
                    std::copy(lheReweighingIDs.begin(), lheReweighingIDs.end(), std::back_inserter(weightChoice->rwgtIDs));
                    
                    std::stringstream pdfDoc; pdfDoc << "LHE pdf variation weights (w_var / w_nominal) for LHA IDs ";
                    bool found = false;
                    for (uint32_t lhaid : preferredPDFLHAIDs_) {
                        for (const auto & pw : pdfSetWeightIDs) {
                            if (pw.lhaIDs.first != lhaid && pw.lhaIDs.first != (lhaid+1)) continue; // sometimes the first weight is not saved if that PDF is the nominal one for the sample
                            if (pw.wids.size() == 1) continue; // only consider error sets
                            pdfDoc << pw.lhaIDs.first << " - " << pw.lhaIDs.second;
                            weightChoice->pdfWeightIDs = pw.wids;
                            if (maxPdfWeights_ < pw.wids.size()) {
                                weightChoice->pdfWeightIDs.resize(maxPdfWeights_); // drop some replicas
                                pdfDoc << ", truncated to the first " << maxPdfWeights_ << " replicas";
                            } 
                            weightChoice->pdfWeightsDoc = pdfDoc.str(); 
                            found = true; break;
                        }
                        if (found) break;
                    }
                }
            }
            return weightChoice; 
        }


        // create an empty counter
        std::unique_ptr<Counter> beginStream(edm::StreamID) const override { 
            return std::make_unique<Counter>(); 
        }
        // inizialize to zero at begin run
        void streamBeginRun(edm::StreamID id, edm::Run const&, edm::EventSetup const&) const override { 
            streamCache(id)->clear(); 
        }
        // create an empty counter
        std::shared_ptr<Counter> globalBeginRunSummary(edm::Run const&, edm::EventSetup const&) const override { 
            return std::make_shared<Counter>(); 
        }
        // add this stream to the summary
        void streamEndRunSummary(edm::StreamID id, edm::Run const&, edm::EventSetup const&, Counter* runCounter) const override { 
            runCounter->merge(*streamCache(id)); 
        }
        // nothing to do per se
        void globalEndRunSummary(edm::Run const&, edm::EventSetup const&, Counter* runCounter) const override { 
        }
        // write the total to the run 
        void globalEndRunProduce(edm::Run& iRun, edm::EventSetup const&, Counter const* runCounter) const override {
            auto out = std::make_unique<nanoaod::MergeableCounterTable>();
            out->addInt("genEventCount", "event count", runCounter->num);
            out->addFloat("genEventSumw", "sum of gen weights", runCounter->sumw);
            out->addFloat("genEventSumw2", "sum of gen (weight^2)", runCounter->sumw2);

            double norm = runCounter->sumw ? 1.0/runCounter->sumw : 1;
            auto sumScales = runCounter->sumScale; for (auto & val : sumScales) val *= norm;
            out->addVFloat("LHEScaleSumw", "Sum of genEventWeight * LHEScaleWeight[i], divided by genEventSumw", sumScales);
            auto sumPDFs = runCounter->sumPDF; for (auto & val : sumPDFs) val *= norm;
            out->addVFloat("LHEPdfSumw", "Sum of genEventWeight * LHEPdfWeight[i], divided by genEventSumw", sumPDFs);
            if (!runCounter->sumRwgt.empty()) {
                auto sumRwgts = runCounter->sumRwgt; for (auto & val : sumRwgts) val *= norm;
                out->addVFloat("LHEReweightingSumw", "Sum of genEventWeight * LHEReweightingWeight[i], divided by genEventSumw", sumRwgts);
            }
            if (!runCounter->sumNamed.empty()) { // it could be empty if there's no LHE info in the sample
                for (unsigned int i = 0, n = namedWeightLabels_.size(); i < n; ++i) {
                    out->addFloat("LHESumw_"+namedWeightLabels_[i], "Sum of genEventWeight * LHEWeight_"+namedWeightLabels_[i]+", divided by genEventSumw", runCounter->sumNamed[i] * norm);
                }
            }
            iRun.put(std::move(out));
        }
        // nothing to do here
        void globalEndRun(edm::Run const&, edm::EventSetup const&) const override { }

        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
            edm::ParameterSetDescription desc;
            desc.add<edm::InputTag>("genEvent", edm::InputTag("generator"))->setComment("tag for the GenEventInfoProduct, to get the main weight");
            desc.add<std::vector<edm::InputTag>>("lheInfo", std::vector<edm::InputTag>{{"externalLHEProducer"},{"source"}})->setComment("tag(s) for the LHE information (LHEEventProduct and LHERunInfoProduct)");

            edm::ParameterSetDescription prefpdf;
            prefpdf.add<std::string>("name");
            prefpdf.add<uint32_t>("lhaid");
            desc.addVPSet("preferredPDFs", prefpdf, std::vector<edm::ParameterSet>())->setComment("LHA PDF Ids of the preferred PDF sets, in order of preference (the first matching one will be used)");
            desc.add<std::vector<std::string>>("namedWeightIDs")->setComment("set of LHA weight IDs for named LHE weights");
            desc.add<std::vector<std::string>>("namedWeightLabels")->setComment("output names for the namedWeightIDs (in the same order)");
            desc.add<int32_t>("lheWeightPrecision")->setComment("Number of bits in the mantissa for LHE weights");
            desc.add<uint32_t>("maxPdfWeights")->setComment("Maximum number of PDF weights to save (to crop NN replicas)");
            desc.addOptionalUntracked<bool>("debug")->setComment("dump out all LHE information for one event");
            descriptions.add("genWeightsTable", desc);
        }


    protected:
        const edm::EDGetTokenT<GenEventInfoProduct> genTag_;
        const std::vector<edm::InputTag> lheLabel_;
        const std::vector<edm::EDGetTokenT<LHEEventProduct>> lheTag_;
        const std::vector<edm::EDGetTokenT<LHERunInfoProduct>> lheRunTag_;

        std::vector<uint32_t> preferredPDFLHAIDs_;
        std::unordered_map<std::string,uint32_t> lhaNameToID_;
        std::vector<std::string> namedWeightIDs_;
        std::vector<std::string> namedWeightLabels_;
        int lheWeightPrecision_;
        unsigned int maxPdfWeights_;

        mutable std::atomic<bool> debug_, debugRun_, hasIssuedWarning_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenWeightsTableProducer);

