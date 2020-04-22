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
            lheDebug = false;

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

            std::vector<std::string> eft_weight_names_small = {
                        "sm",
                        "ctw_-5p0",
                        "ctw_-4p75",
                        "ctw_-4p5",
                        "ctw_-4p25",
                        "ctw_-4p0",
                        "ctw_-3p75",
                        "ctw_-3p5",
                        "ctw_-3p25",
                        "ctw_-3p0",
                        "ctw_-2p75",
                        "ctw_-2p5",
                        "ctw_-2p25",
                        "ctw_-2p0",
                        "ctw_-1p75",
                        "ctw_-1p5",
                        "ctw_-1p25",
                        "ctw_-1p0",
                        "ctw_-0p75",
                        "ctw_-0p5",
                        "ctw_-0p25",
                        "ctw_0p0",
                        "ctw_0p25",
                        "ctw_0p5",
                        "ctw_0p75",
                        "ctw_1p0",
                        "ctw_1p25",
                        "ctw_1p5",
                        "ctw_1p75",
                        "ctw_2p0",
                        "ctw_2p25",
                        "ctw_2p5",
                        "ctw_2p75",
                        "ctw_3p0",
                        "ctw_3p25",
                        "ctw_3p5",
                        "ctw_3p75",
                        "ctw_4p0",
                        "ctw_4p25",
                        "ctw_4p5",
                        "ctw_4p75",
                        "ctw_5p0",
                        "ctwi_-5p0",
                        "ctwi_-4p75",
                        "ctwi_-4p5",
                        "ctwi_-4p25",
                        "ctwi_-4p0",
                        "ctwi_-3p75",
                        "ctwi_-3p5",
                        "ctwi_-3p25",
                        "ctwi_-3p0",
                        "ctwi_-2p75",
                        "ctwi_-2p5",
                        "ctwi_-2p25",
                        "ctwi_-2p0",
                        "ctwi_-1p75",
                        "ctwi_-1p5",
                        "ctwi_-1p25",
                        "ctwi_-1p0",
                        "ctwi_-0p75",
                        "ctwi_-0p5",
                        "ctwi_-0p25",
                        "ctwi_0p0",
                        "ctwi_0p25",
                        "ctwi_0p5",
                        "ctwi_0p75",
                        "ctwi_1p0",
                        "ctwi_1p25",
                        "ctwi_1p5",
                        "ctwi_1p75",
                        "ctwi_2p0",
                        "ctwi_2p25",
                        "ctwi_2p5",
                        "ctwi_2p75",
                        "ctwi_3p0",
                        "ctwi_3p25",
                        "ctwi_3p5",
                        "ctwi_3p75",
                        "ctwi_4p0",
                        "ctwi_4p25",
                        "ctwi_4p5",
                        "ctwi_4p75",
                        "ctwi_5p0",
                        "cpq3_-5p0",
                        "cpq3_-4p75",
                        "cpq3_-4p5",
                        "cpq3_-4p25",
                        "cpq3_-4p0",
                        "cpq3_-3p75",
                        "cpq3_-3p5",
                        "cpq3_-3p25",
                        "cpq3_-3p0",
                        "cpq3_-2p75",
                        "cpq3_-2p5",
                        "cpq3_-2p25",
                        "cpq3_-2p0",
                        "cpq3_-1p75",
                        "cpq3_-1p5",
                        "cpq3_-1p25",
                        "cpq3_-1p0",
                        "cpq3_-0p75",
                        "cpq3_-0p5",
                        "cpq3_-0p25",
                        "cpq3_0p0",
                        "cpq3_0p25",
                        "cpq3_0p5",
                        "cpq3_0p75",
                        "cpq3_1p0",
                        "cpq3_1p25",
                        "cpq3_1p5",
                        "cpq3_1p75",
                        "cpq3_2p0",
                        "cpq3_2p25",
                        "cpq3_2p5",
                        "cpq3_2p75",
                        "cpq3_3p0",
                        "cpq3_3p25",
                        "cpq3_3p5",
                        "cpq3_3p75",
                        "cpq3_4p0",
                        "cpq3_4p25",
                        "cpq3_4p5",
                        "cpq3_4p75",
                        "cpq3_5p0",
                        "cqq13_-3p0",
                        "cqq13_-2p75",
                        "cqq13_-2p5",
                        "cqq13_-2p25",
                        "cqq13_-2p0",
                        "cqq13_-1p75",
                        "cqq13_-1p5",
                        "cqq13_-1p25",
                        "cqq13_-1p0",
                        "cqq13_-0p75",
                        "cqq13_-0p5",
                        "cqq13_-0p25",
                        "cqq13_0p0",
                        "cqq13_0p25",
                        "cqq13_0p5",
                        "cqq13_0p75",
                        "cqq13_1p0",
                        "cqq13_1p25",
                        "cqq13_1p5",
                        "cqq13_1p75",
                        "cqq13_2p0",
                        "cqq13_2p25",
                        "cqq13_2p5",
                        "cqq13_2p75",
                        "cqq13_3p0",
                        "cqq13_3p25",
                        "cqq13_3p5",
                        "cqq13_3p75",
                        "cqq13_4p0",
                        "cqq13_4p25",
                        "cqq13_4p5",
                        "cqq13_4p75",
                        "cqq13_5p0",
                        "cqq13_5p25",
                        "cqq13_5p5",
                        "cqq13_5p75",
                        "cqq13_6p0",
                        "cqq13_6p25",
                        "cqq13_6p5",
                        "cqq13_6p75",
                        "cqq13_7p0",
                        "cqq83_-8p0",
                        "cqq83_-7p6",
                        "cqq83_-7p2",
                        "cqq83_-6p8",
                        "cqq83_-6p4",
                        "cqq83_-6p0",
                        "cqq83_-5p6",
                        "cqq83_-5p2",
                        "cqq83_-4p8",
                        "cqq83_-4p4",
                        "cqq83_-4p0",
                        "cqq83_-3p6",
                        "cqq83_-3p2",
                        "cqq83_-2p8",
                        "cqq83_-2p4",
                        "cqq83_-2p0",
                        "cqq83_-1p6",
                        "cqq83_-1p2",
                        "cqq83_-0p8",
                        "cqq83_-0p4",
                        "cqq83_0p0",
                        "cqq83_0p4",
                        "cqq83_0p8",
                        "cqq83_1p2",
                        "cqq83_1p6",
                        "cqq83_2p0",
                        "cqq83_2p4",
                        "cqq83_2p8",
                        "cqq83_3p2",
                        "cqq83_3p6",
                        "cqq83_4p0",
                        "cqq83_4p4",
                        "cqq83_4p8",
                        "cqq83_5p2",
                        "cqq83_5p6",
                        "cqq83_6p0",
                        "cqq83_6p4",
                        "cqq83_6p8",
                        "cqq83_7p2",
                        "cqq83_7p6",
                        "cqq83_8p0",
                        "ctw_-3p37_ctwi_-0p01_cpq3_4p4_cqq13_6p5_cqq83_-5p76",
                        "ctw_0p32_ctwi_1p2_cpq3_-3p54_cqq13_6p49_cqq83_-4p11",
                        "ctw_-3p36_ctwi_0p06_cpq3_0p05_cqq13_4p73_cqq83_5p65",
                        "ctw_-1p73_ctwi_2p08_cpq3_0p62_cqq13_2p35_cqq83_1p59",
                        "ctw_-0p54_ctwi_-1p28_cpq3_4p3_cqq13_4p63_cqq83_2p81",
                        "ctw_3p86_ctwi_1p04_cpq3_-4p29_cqq13_-0p72_cqq83_4p85",
                        "ctw_-1p54_ctwi_4p4_cpq3_1p89_cqq13_0p8_cqq83_-0p12",
                        "ctw_-4p48_ctwi_-0p52_cpq3_-4p2_cqq13_-0p07_cqq83_6p87",
                        "ctw_-3p52_ctwi_1p77_cpq3_1p51_cqq13_3p93_cqq83_4p71",
                        "ctw_3p64_ctwi_-1p43_cpq3_-1p29_cqq13_-1p26_cqq83_2p23",
                        "ctw_-2p88_ctwi_-3p04_cpq3_4p9_cqq13_3p2_cqq83_3p69",
                        "ctw_0p02_ctwi_-3p67_cpq3_4p1_cqq13_3p51_cqq83_-7p21",
                        "ctw_4p14_ctwi_-0p32_cpq3_0p98_cqq13_0p72_cqq83_-1p65",
                        "ctw_2p78_ctwi_-2p66_cpq3_-4p91_cqq13_0p13_cqq83_7p77",
                        "ctw_-3p8_ctwi_1p0_cpq3_3p65_cqq13_2p34_cqq83_1p5",
                        "ctw_0p37_ctwi_-1p78_cpq3_-1p37_cqq13_6p27_cqq83_3p39",
                        "ctw_1p43_ctwi_-0p43_cpq3_4p11_cqq13_5p77_cqq83_7p47",
                        "ctw_0p44_ctwi_-0p91_cpq3_-3p08_cqq13_-0p27_cqq83_5p88",
                        "ctw_-3p26_ctwi_4p46_cpq3_0p86_cqq13_5p47_cqq83_0p48",
                        "ctw_2p55_ctwi_0p87_cpq3_3p28_cqq13_2p14_cqq83_3p59",
                        "ctw_3p04_ctwi_-1p27_cpq3_-3p97_cqq13_3p75_cqq83_5p6",
                        "ctw_3p34_ctwi_-3p73_cpq3_2p84_cqq13_1p88_cqq83_0p74",
                        "ctw_-3p68_ctwi_2p7_cpq3_-2p61_cqq13_3p99_cqq83_-2p64",
                        "ctw_4p57_ctwi_2p12_cpq3_-1p45_cqq13_2p48_cqq83_-4p99",
                        "ctw_-0p27_ctwi_4p31_cpq3_4p29_cqq13_0p2_cqq83_-0p93",
                        "ctw_0p25_ctwi_0p39_cpq3_4p9_cqq13_-2p14_cqq83_-4p84",
                        "ctw_1p08_ctwi_-1p64_cpq3_-1p22_cqq13_3p95_cqq83_-4p81",
                        "ctw_-4p71_ctwi_4p71_cpq3_-1p84_cqq13_-2p0_cqq83_-7p11",
                        "ctw_1p8_ctwi_-0p76_cpq3_-1p6_cqq13_6p8_cqq83_6p87",
                        "ctw_2p15_ctwi_4p03_cpq3_-3p34_cqq13_6p72_cqq83_-7p86",
                        "ctw_0p97_ctwi_4p21_cpq3_4p1_cqq13_4p04_cqq83_6p82",
                        "ctw_1p09_ctwi_-3p44_cpq3_4p21_cqq13_-1p06_cqq83_-0p75",
                        "ctw_-0p55_ctwi_0p02_cpq3_-1p64_cqq13_-0p31_cqq83_6p22",
                        "ctw_-1p08_ctwi_1p86_cpq3_-0p06_cqq13_-0p52_cqq83_-3p42",
                        "ctw_4p7_ctwi_4p04_cpq3_-3p42_cqq13_-0p67_cqq83_7p34",
                        "ctw_4p21_ctwi_-1p46_cpq3_0p44_cqq13_1p3_cqq83_-0p89",
                        "ctw_0p5_ctwi_2p8_cpq3_2p18_cqq13_4p17_cqq83_-6p87",
                        "ctw_-2p56_ctwi_1p46_cpq3_2p14_cqq13_1p85_cqq83_4p47",
                        "ctw_4p57_ctwi_1p92_cpq3_3p32_cqq13_0p27_cqq83_5p08",
                        "ctw_0p18_ctwi_3p06_cpq3_1p82_cqq13_0p03_cqq83_4p32",
                        "ctw_3p0_ctwi_0p01_cpq3_-3p57_cqq13_0p75_cqq83_-5p17",
                        "ctw_-3p63_ctwi_-1p99_cpq3_-0p12_cqq13_5p01_cqq83_1p62",
                            };

            std::vector<std::string> eft_weight_names_large = {
                "SM",
                "ctW_-5p0",
                "ctW_-4p75",
                "ctW_-4p5",
                "ctW_-4p25",
                "ctW_-4p0",
                "ctW_-3p75",
                "ctW_-3p5",
                "ctW_-3p25",
                "ctW_-3p0",
                "ctW_-2p75",
                "ctW_-2p5",
                "ctW_-2p25",
                "ctW_-2p0",
                "ctW_-1p75",
                "ctW_-1p5",
                "ctW_-1p25",
                "ctW_-1p0",
                "ctW_-0p75",
                "ctW_-0p5",
                "ctW_-0p25",
                "ctW_0p0",
                "ctW_0p25",
                "ctW_0p5",
                "ctW_0p75",
                "ctW_1p0",
                "ctW_1p25",
                "ctW_1p5",
                "ctW_1p75",
                "ctW_2p0",
                "ctW_2p25",
                "ctW_2p5",
                "ctW_2p75",
                "ctW_3p0",
                "ctW_3p25",
                "ctW_3p5",
                "ctW_3p75",
                "ctW_4p0",
                "ctW_4p25",
                "ctW_4p5",
                "ctW_4p75",
                "ctW_5p0",
                "ctWI_-5p0",
                "ctWI_-4p75",
                "ctWI_-4p5",
                "ctWI_-4p25",
                "ctWI_-4p0",
                "ctWI_-3p75",
                "ctWI_-3p5",
                "ctWI_-3p25",
                "ctWI_-3p0",
                "ctWI_-2p75",
                "ctWI_-2p5",
                "ctWI_-2p25",
                "ctWI_-2p0",
                "ctWI_-1p75",
                "ctWI_-1p5",
                "ctWI_-1p25",
                "ctWI_-1p0",
                "ctWI_-0p75",
                "ctWI_-0p5",
                "ctWI_-0p25",
                "ctWI_0p0",
                "ctWI_0p25",
                "ctWI_0p5",
                "ctWI_0p75",
                "ctWI_1p0",
                "ctWI_1p25",
                "ctWI_1p5",
                "ctWI_1p75",
                "ctWI_2p0",
                "ctWI_2p25",
                "ctWI_2p5",
                "ctWI_2p75",
                "ctWI_3p0",
                "ctWI_3p25",
                "ctWI_3p5",
                "ctWI_3p75",
                "ctWI_4p0",
                "ctWI_4p25",
                "ctWI_4p5",
                "ctWI_4p75",
                "ctWI_5p0",
                "cpQ3_-5p0",
                "cpQ3_-4p75",
                "cpQ3_-4p5",
                "cpQ3_-4p25",
                "cpQ3_-4p0",
                "cpQ3_-3p75",
                "cpQ3_-3p5",
                "cpQ3_-3p25",
                "cpQ3_-3p0",
                "cpQ3_-2p75",
                "cpQ3_-2p5",
                "cpQ3_-2p25",
                "cpQ3_-2p0",
                "cpQ3_-1p75",
                "cpQ3_-1p5",
                "cpQ3_-1p25",
                "cpQ3_-1p0",
                "cpQ3_-0p75",
                "cpQ3_-0p5",
                "cpQ3_-0p25",
                "cpQ3_0p0",
                "cpQ3_0p25",
                "cpQ3_0p5",
                "cpQ3_0p75",
                "cpQ3_1p0",
                "cpQ3_1p25",
                "cpQ3_1p5",
                "cpQ3_1p75",
                "cpQ3_2p0",
                "cpQ3_2p25",
                "cpQ3_2p5",
                "cpQ3_2p75",
                "cpQ3_3p0",
                "cpQ3_3p25",
                "cpQ3_3p5",
                "cpQ3_3p75",
                "cpQ3_4p0",
                "cpQ3_4p25",
                "cpQ3_4p5",
                "cpQ3_4p75",
                "cpQ3_5p0",
                "cQq13_-3p0",
                "cQq13_-2p75",
                "cQq13_-2p5",
                "cQq13_-2p25",
                "cQq13_-2p0",
                "cQq13_-1p75",
                "cQq13_-1p5",
                "cQq13_-1p25",
                "cQq13_-1p0",
                "cQq13_-0p75",
                "cQq13_-0p5",
                "cQq13_-0p25",
                "cQq13_0p0",
                "cQq13_0p25",
                "cQq13_0p5",
                "cQq13_0p75",
                "cQq13_1p0",
                "cQq13_1p25",
                "cQq13_1p5",
                "cQq13_1p75",
                "cQq13_2p0",
                "cQq13_2p25",
                "cQq13_2p5",
                "cQq13_2p75",
                "cQq13_3p0",
                "cQq13_3p25",
                "cQq13_3p5",
                "cQq13_3p75",
                "cQq13_4p0",
                "cQq13_4p25",
                "cQq13_4p5",
                "cQq13_4p75",
                "cQq13_5p0",
                "cQq13_5p25",
                "cQq13_5p5",
                "cQq13_5p75",
                "cQq13_6p0",
                "cQq13_6p25",
                "cQq13_6p5",
                "cQq13_6p75",
                "cQq13_7p0",
                "cQq83_-8p0",
                "cQq83_-7p6",
                "cQq83_-7p2",
                "cQq83_-6p8",
                "cQq83_-6p4",
                "cQq83_-6p0",
                "cQq83_-5p6",
                "cQq83_-5p2",
                "cQq83_-4p8",
                "cQq83_-4p4",
                "cQq83_-4p0",
                "cQq83_-3p6",
                "cQq83_-3p2",
                "cQq83_-2p8",
                "cQq83_-2p4",
                "cQq83_-2p0",
                "cQq83_-1p6",
                "cQq83_-1p2",
                "cQq83_-0p8",
                "cQq83_-0p4",
                "cQq83_0p0",
                "cQq83_0p4",
                "cQq83_0p8",
                "cQq83_1p2",
                "cQq83_1p6",
                "cQq83_2p0",
                "cQq83_2p4",
                "cQq83_2p8",
                "cQq83_3p2",
                "cQq83_3p6",
                "cQq83_4p0",
                "cQq83_4p4",
                "cQq83_4p8",
                "cQq83_5p2",
                "cQq83_5p6",
                "cQq83_6p0",
                "cQq83_6p4",
                "cQq83_6p8",
                "cQq83_7p2",
                "cQq83_7p6",
                "cQq83_8p0",
                "ctW_-3p37_ctWI_-0p01_cpQ3_4p4_cQq13_6p5_cQq83_-5p76",
                "ctW_0p32_ctWI_1p2_cpQ3_-3p54_cQq13_6p49_cQq83_-4p11",
                "ctW_-3p36_ctWI_0p06_cpQ3_0p05_cQq13_4p73_cQq83_5p65",
                "ctW_-1p73_ctWI_2p08_cpQ3_0p62_cQq13_2p35_cQq83_1p59",
                "ctW_-0p54_ctWI_-1p28_cpQ3_4p3_cQq13_4p63_cQq83_2p81",
                "ctW_3p86_ctWI_1p04_cpQ3_-4p29_cQq13_-0p72_cQq83_4p85",
                "ctW_-1p54_ctWI_4p4_cpQ3_1p89_cQq13_0p8_cQq83_-0p12",
                "ctW_-4p48_ctWI_-0p52_cpQ3_-4p2_cQq13_-0p07_cQq83_6p87",
                "ctW_-3p52_ctWI_1p77_cpQ3_1p51_cQq13_3p93_cQq83_4p71",
                "ctW_3p64_ctWI_-1p43_cpQ3_-1p29_cQq13_-1p26_cQq83_2p23",
                "ctW_-2p88_ctWI_-3p04_cpQ3_4p9_cQq13_3p2_cQq83_3p69",
                "ctW_0p02_ctWI_-3p67_cpQ3_4p1_cQq13_3p51_cQq83_-7p21",
                "ctW_4p14_ctWI_-0p32_cpQ3_0p98_cQq13_0p72_cQq83_-1p65",
                "ctW_2p78_ctWI_-2p66_cpQ3_-4p91_cQq13_0p13_cQq83_7p77",
                "ctW_-3p8_ctWI_1p0_cpQ3_3p65_cQq13_2p34_cQq83_1p5",
                "ctW_0p37_ctWI_-1p78_cpQ3_-1p37_cQq13_6p27_cQq83_3p39",
                "ctW_1p43_ctWI_-0p43_cpQ3_4p11_cQq13_5p77_cQq83_7p47",
                "ctW_0p44_ctWI_-0p91_cpQ3_-3p08_cQq13_-0p27_cQq83_5p88",
                "ctW_-3p26_ctWI_4p46_cpQ3_0p86_cQq13_5p47_cQq83_0p48",
                "ctW_2p55_ctWI_0p87_cpQ3_3p28_cQq13_2p14_cQq83_3p59",
                "ctW_3p04_ctWI_-1p27_cpQ3_-3p97_cQq13_3p75_cQq83_5p6",
                "ctW_3p34_ctWI_-3p73_cpQ3_2p84_cQq13_1p88_cQq83_0p74",
                "ctW_-3p68_ctWI_2p7_cpQ3_-2p61_cQq13_3p99_cQq83_-2p64",
                "ctW_4p57_ctWI_2p12_cpQ3_-1p45_cQq13_2p48_cQq83_-4p99",
                "ctW_-0p27_ctWI_4p31_cpQ3_4p29_cQq13_0p2_cQq83_-0p93",
                "ctW_0p25_ctWI_0p39_cpQ3_4p9_cQq13_-2p14_cQq83_-4p84",
                "ctW_1p08_ctWI_-1p64_cpQ3_-1p22_cQq13_3p95_cQq83_-4p81",
                "ctW_-4p71_ctWI_4p71_cpQ3_-1p84_cQq13_-2p0_cQq83_-7p11",
                "ctW_1p8_ctWI_-0p76_cpQ3_-1p6_cQq13_6p8_cQq83_6p87",
                "ctW_2p15_ctWI_4p03_cpQ3_-3p34_cQq13_6p72_cQq83_-7p86",
                "ctW_0p97_ctWI_4p21_cpQ3_4p1_cQq13_4p04_cQq83_6p82",
                "ctW_1p09_ctWI_-3p44_cpQ3_4p21_cQq13_-1p06_cQq83_-0p75",
                "ctW_-0p55_ctWI_0p02_cpQ3_-1p64_cQq13_-0p31_cQq83_6p22",
                "ctW_-1p08_ctWI_1p86_cpQ3_-0p06_cQq13_-0p52_cQq83_-3p42",
                "ctW_4p7_ctWI_4p04_cpQ3_-3p42_cQq13_-0p67_cQq83_7p34",
                "ctW_4p21_ctWI_-1p46_cpQ3_0p44_cQq13_1p3_cQq83_-0p89",
                "ctW_0p5_ctWI_2p8_cpQ3_2p18_cQq13_4p17_cQq83_-6p87",
                "ctW_-2p56_ctWI_1p46_cpQ3_2p14_cQq13_1p85_cQq83_4p47",
                "ctW_4p57_ctWI_1p92_cpQ3_3p32_cQq13_0p27_cQq83_5p08",
                "ctW_0p18_ctWI_3p06_cpQ3_1p82_cQq13_0p03_cQq83_4p32",
                "ctW_3p0_ctWI_0p01_cpQ3_-3p57_cQq13_0p75_cQq83_-5p17",
                "ctW_-3p63_ctWI_-1p99_cpQ3_-0p12_cQq13_5p01_cQq83_1p62",
            };

            std::vector<std::string> eft_weight_names = {} ; //eft_weight_names_small;

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
            lheDebug=false;
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

