#ifndef TL1EVENTCLASS_H
#define TL1EVENTCLASS_H

#include <string>
#include <vector>

#include <TMath.h>
#include <TVector2.h>

#ifndef __CINT__
#include "DataFormats/L1Trigger/interface/EtSum.h"
#endif // __CINT__
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "TL1PrimitiveEventClass.h"

class TL1EventClass
{
    public:
        TL1EventClass(std::vector<std::string> inDir);

        bool Next();
        void GetEntry(int i);
        void GetDerivatives();
        void Reset();

        TL1PrimitiveEventClass const * GetPEvent() const;

        // Get L1
        double fL1Met, fL1Mht, fL1Ett, fL1Htt, fL1MetPhi, fL1MhtPhi;
        double fL1Mex, fL1Mey;
        double fL1EmuMet, fL1EmuMht, fL1EmuEtt, fL1EmuHtt, fL1EmuMetPhi, fL1EmuMhtPhi;
        double fL1EmuMex, fL1EmuMey;
        double fL1MetHF, fL1MhtHF, fL1EttHF, fL1HttHF, fL1MetPhiHF, fL1MhtPhiHF;
        double fL1EmuMetHF, fL1EmuMhtHF, fL1EmuEttHF, fL1EmuHttHF, fL1EmuMetPhiHF, fL1EmuMhtPhiHF;
        std::vector<double> fL1JetEt, fL1JetPhi, fL1JetEta;
        std::vector<double> fL1EmuJetEt, fL1EmuJetPhi, fL1EmuJetEta;

        // Filter flags
        bool fMuonFilterPassFlag, fMetFilterPassFlag;
        std::vector<bool> fJetFilterPassFlags;
        
        struct MissingEt{
            TVector2 fMET;
            double phi()const{return fMET.Phi();}
            double met()const{return fMET.Mod();}
            void AddEt(double Et, double phi){
                TVector2 temp(0.0,0.0);
                temp.SetMagPhi(Et,phi);
                fMET-=temp;
            }
            void Reset(){ fMET.Set(0.,0.); }
        };

        // Recalc L1 MET sums
        double fRecalcL1Met, fRecalcL1MetPhi;
        double fRecalcL1MetHF, fRecalcL1MetPhiHF;
        MissingEt fRecalcL1Met28Only, fRecalcL1MetNot28, fRecalcL1MetNot28HF;

        MissingEt fRecalcL1EmuMet;
        MissingEt fRecalcL1EmuMetHF;
        MissingEt fRecalcL1EmuMet28Only;
        MissingEt fRecalcL1EmuMetNot28;
        MissingEt fRecalcL1EmuMetPUS;
        MissingEt fRecalcL1EmuMetPUS28;
        MissingEt fRecalcL1EmuMetPUSHF;
        MissingEt fRecalcL1EmuMetPUSThresh;
        MissingEt fRecalcL1EmuMetPUSThreshHF;

        // Recalc L1 Ht/Et Sums
        double fRecalcL1Mht, fRecalcL1MhtPhi;
        int fNJetL1Mht;
        double fRecalcL1Ett;

        // Recalc Reco Ht/Et Sums
        double fRecalcRecoMht, fRecalcRecoMhtPhi, fRecalcRecoHtt;
        int fNJetRecoMht;
        bool fMhtPassFlag;
        double fRecalcRecoEtt;

        // Jets
        bool fIsLeadingRecoJet, fIsMatchedL1Jet, fIsMatchedL1EmuJet;
        int fLeadingRecoJetIndex, fMatchedL1JetIndex, fMatchedL1EmuJetIndex;

    private:
        TL1PrimitiveEventClass * fPrimitiveEvent;

        // Get L1
        void GetL1Jets();
        void GetL1EmuJets();
        void GetL1Sums();
        void GetL1EmuSums();

        // Filters
        void MuonFilter();
        void JetFilter();
        void SumsFilter();

        // Recalc L1 MET sum
        void GetRecalcL1Met();
        void GetRecalcL1EmuMet();

        // Recalc L1 Ht/Et Sums
        void GetRecalcL1Mht();
        void GetRecalcL1Ett();

        // Recalc Reco Ht/Et Sums
        void GetRecalcRecoHtSums();
        void GetRecalcRecoEtt();

        // Jets
        void GetLeadingRecoJet();
        void GetMatchedL1Jet(const std::string & l1Type);

};

TL1EventClass::TL1EventClass(std::vector<std::string> inDir) :
    fPrimitiveEvent(new TL1PrimitiveEventClass(inDir)),
    fMuonFilterPassFlag(true), fMetFilterPassFlag(true),
    fMhtPassFlag(true),
    fIsLeadingRecoJet(true), fIsMatchedL1Jet(true), fIsMatchedL1EmuJet(true)
{
}

bool TL1EventClass::Next()
{
    bool isNext = fPrimitiveEvent->Next();
    this->GetDerivatives();
    return isNext;
}

void TL1EventClass::GetEntry(int i)
{
    Reset();
    fPrimitiveEvent->GetEntry(i);
    this->GetDerivatives();
}

void TL1EventClass::Reset(){
        fL1Met=fL1Mht=fL1Ett=fL1Htt=fL1MetPhi=fL1MhtPhi=0;
        fL1Mex=fL1Mey=0;
        fL1EmuMet=fL1EmuMht=fL1EmuEtt=fL1EmuHtt=fL1EmuMetPhi=fL1EmuMhtPhi=0;
        fL1EmuMex=fL1EmuMey=0;
        fL1MetHF=fL1MhtHF=fL1EttHF=fL1HttHF=fL1MetPhiHF=fL1MhtPhiHF=0;
        fL1EmuMetHF=fL1EmuMhtHF=fL1EmuEttHF=fL1EmuHttHF=fL1EmuMetPhiHF=fL1EmuMhtPhiHF=0;

        // Filter flags
        fMuonFilterPassFlag=fMetFilterPassFlag=false;
        fJetFilterPassFlags.clear();

        // Recalc L1 MET sums
        fRecalcL1Met=fRecalcL1MetPhi=0;
        fRecalcL1MetHF=fRecalcL1MetPhiHF=0;
        fRecalcL1Met28Only.Reset();fRecalcL1MetNot28.Reset();fRecalcL1MetNot28HF.Reset();

        fRecalcL1EmuMet.Reset();
        fRecalcL1EmuMetHF.Reset();
        fRecalcL1EmuMet28Only.Reset();
        fRecalcL1EmuMetNot28.Reset();
        fRecalcL1EmuMetPUS.Reset();
        fRecalcL1EmuMetPUS28.Reset();
        fRecalcL1EmuMetPUSHF.Reset();
        fRecalcL1EmuMetPUSThresh.Reset();
        fRecalcL1EmuMetPUSThreshHF.Reset();


        // Recalc L1 Ht/Et Sums
        fRecalcL1Mht=fRecalcL1MhtPhi=0;
        fNJetL1Mht=0;
        fRecalcL1Ett=0;

        // Recalc Reco Ht/Et Sums
        fRecalcRecoMht=fRecalcRecoMhtPhi=fRecalcRecoHtt=0;
        fNJetRecoMht=0;
        fMhtPassFlag=false;
        fRecalcRecoEtt=0;

        // Jets
        fIsLeadingRecoJet=fIsMatchedL1Jet=fIsMatchedL1EmuJet=false;
        fLeadingRecoJetIndex=fMatchedL1JetIndex=fMatchedL1EmuJetIndex=0;
}

void TL1EventClass::GetDerivatives()
{
    // L1
    bool isUpgrade = fPrimitiveEvent->fIsUpgrade;
    bool isEmuUpgrade = fPrimitiveEvent->fIsEmuUpgrade;
    if( isUpgrade ) this->GetL1Jets();
    if( isUpgrade ) this->GetL1Sums();
    if( isEmuUpgrade ) this->GetL1EmuJets();
    if( isEmuUpgrade ) this->GetL1EmuSums();

    // Filter
    bool isJetReco = fPrimitiveEvent->fIsJetReco;
    if( fPrimitiveEvent->fIsMuonReco ) this->MuonFilter();
    if( isJetReco ) this->JetFilter();
    if( fPrimitiveEvent->fIsMetFilterReco ) this->SumsFilter();

    // Recalc
    if( fPrimitiveEvent->fIsCaloTower ) this->GetRecalcL1Met();
    if( fPrimitiveEvent->fIsEmuCaloTower ) this->GetRecalcL1EmuMet();
    //this->GetRecalcL1Mht();
    //this->GetRecalcL1Ett();
    
    //this->GetRecalcRecoHtSums();
    //this->GetRecalcRecoEtt();

    // Jets
    if( isJetReco ) this->GetLeadingRecoJet();
    if( isJetReco ) this->GetMatchedL1Jet("HW");
    if( isJetReco ) this->GetMatchedL1Jet("EMU");
}

TL1PrimitiveEventClass const * TL1EventClass::GetPEvent() const
{
    return fPrimitiveEvent;
}

void TL1EventClass::GetL1Jets()
{
    fL1JetEt.clear();
    fL1JetEta.clear();
    fL1JetPhi.clear();
    for(unsigned iJet=0; iJet<fPrimitiveEvent->fUpgrade->nJets; ++iJet)
    {
        if( fPrimitiveEvent->fUpgrade->jetBx[iJet] == 0 )
        {
            fL1JetEt.push_back(fPrimitiveEvent->fUpgrade->jetEt[iJet]);
            fL1JetEta.push_back(fPrimitiveEvent->fUpgrade->jetEta[iJet]);
            fL1JetPhi.push_back(fPrimitiveEvent->fUpgrade->jetPhi[iJet]);
        }
    }
}

void TL1EventClass::GetL1EmuJets()
{
    fL1EmuJetEt.clear();
    fL1EmuJetEta.clear();
    fL1EmuJetPhi.clear();
    for(unsigned iJet=0; iJet<fPrimitiveEvent->fEmuUpgrade->nJets; ++iJet)
    {
        if( fPrimitiveEvent->fEmuUpgrade->jetBx[iJet] == 0 )
        {
            fL1EmuJetEt.push_back(fPrimitiveEvent->fEmuUpgrade->jetEt[iJet]);
            fL1EmuJetEta.push_back(fPrimitiveEvent->fEmuUpgrade->jetEta[iJet]);
            fL1EmuJetPhi.push_back(fPrimitiveEvent->fEmuUpgrade->jetPhi[iJet]);
        }
    }
}

void TL1EventClass::GetL1Sums()
{
    auto upgrades = fPrimitiveEvent->fUpgrade;
    for(int iter=0; iter<upgrades->nSums; ++iter)
    {
        double et = upgrades->sumEt[iter];
        double phi = upgrades->sumPhi[iter];
        if( upgrades->sumBx[iter] == 0 )
        {
            int sumType = upgrades->sumType[iter];
            if(sumType == l1t::EtSum::EtSumType::kTotalEt)   fL1Ett = et;
            else if(sumType == l1t::EtSum::EtSumType::kTotalHt)   fL1Htt = et;
            else if(sumType == l1t::EtSum::EtSumType::kMissingEt)
            {
                fL1Met = et;
                fL1MetPhi = phi;
            }
            else if(sumType == l1t::EtSum::EtSumType::kMissingHt)
            {
                fL1Mht = et;
                fL1MhtPhi = phi;
            }
            else if(sumType == l1t::EtSum::EtSumType::kTotalEtHF) fL1EttHF = et;
            else if(sumType == l1t::EtSum::EtSumType::kTotalHtHF) fL1HttHF = et;
            else if(sumType == l1t::EtSum::EtSumType::kMissingEtHF)
            {
                fL1MetHF = et;
                fL1MetPhiHF = phi;
            }
            else if(sumType == l1t::EtSum::EtSumType::kMissingHtHF)
            {
                fL1MhtHF = et;
                fL1MhtPhiHF = phi;
            }
            else if(sumType == l1t::EtSum::EtSumType::kTotalEtx) fL1Mex = et;
            else if(sumType == l1t::EtSum::EtSumType::kTotalEty) fL1Mey = et;
        }
    }
}

void TL1EventClass::GetL1EmuSums()
{
    auto emuUpgrades = fPrimitiveEvent->fEmuUpgrade;
    for(int iter=0; iter<emuUpgrades->nSums; ++iter)
    {
        double et = emuUpgrades->sumEt[iter];
        double phi = emuUpgrades->sumPhi[iter];

        int sumType = emuUpgrades->sumType[iter];
        if(sumType == l1t::EtSum::EtSumType::kTotalEt)   fL1EmuEtt = et;
        if(sumType == l1t::EtSum::EtSumType::kTotalHt)   fL1EmuHtt = et;
        if(sumType == l1t::EtSum::EtSumType::kMissingEt)
        {
            fL1EmuMet = et;
            fL1EmuMetPhi = phi;
        }
        if(sumType == l1t::EtSum::EtSumType::kMissingHt)
        {
            fL1EmuMht = et;
            fL1EmuMhtPhi = phi;
        }
        if(sumType == l1t::EtSum::EtSumType::kTotalEtHF) fL1EmuEttHF = et;
        if(sumType == l1t::EtSum::EtSumType::kTotalHtHF) fL1EmuHttHF = et;
        if(sumType == l1t::EtSum::EtSumType::kMissingEtHF)
        {
            fL1EmuMetHF = et;
            fL1EmuMetPhiHF = phi;
        }
        if(sumType == l1t::EtSum::EtSumType::kMissingHtHF)
        {
            fL1EmuMhtHF = et;
            fL1EmuMhtPhiHF = phi;
        }
        if(sumType == l1t::EtSum::EtSumType::kTotalEtx) fL1EmuMex = et;
        if(sumType == l1t::EtSum::EtSumType::kTotalEty) fL1EmuMey = et;
    }
}

void TL1EventClass::MuonFilter()
{
    bool pass(false);
    for(unsigned jMuon=0; jMuon<fPrimitiveEvent->fMuons->nMuons; ++jMuon)
    {
        double pt = fPrimitiveEvent->fMuons->pt[jMuon];
        double iso = fPrimitiveEvent->fMuons->iso[jMuon];
        double isLoose = fPrimitiveEvent->fMuons->isLooseMuon[jMuon];
        if( pt>=20.0 && iso<=0.1 && isLoose==1.0 ){
            pass = true;
            break;
        }
    }
    fMuonFilterPassFlag = pass;
}

void TL1EventClass::JetFilter()
{
    fJetFilterPassFlags.clear();
    for(unsigned iJet=0; iJet<fPrimitiveEvent->fJets->nJets; ++iJet)
    {
        auto recoJets = fPrimitiveEvent->fJets;

        // Muon mult filter
        int muMult  = recoJets->muMult[iJet];
        if( muMult != 0 ){ fJetFilterPassFlags.push_back(false); continue; }

        // No hf jets
        double eta  = recoJets->eta[iJet];
        if( TMath::Abs(eta) > 3. ){ fJetFilterPassFlags.push_back(false); continue; }

        // nhef
        double nhef = recoJets->nhef[iJet];
        if( nhef >= 0.9 ){ fJetFilterPassFlags.push_back(false); continue; }

        // pef
        double pef  = recoJets->pef[iJet];
        if( pef >= 0.9 ){ fJetFilterPassFlags.push_back(false); continue; }

        // mef
        double mef  = recoJets->mef[iJet];
        if( mef >= 0.8 ){ fJetFilterPassFlags.push_back(false); continue; }

        // sumMult
        int sumMult = recoJets->chMult[iJet] + recoJets->chMult[iJet] + recoJets->elMult[iJet];
        int nhMult  = recoJets->nhMult[iJet];
        int phMult  = recoJets->phMult[iJet];
        if( sumMult+nhMult+phMult <= 1 ){ fJetFilterPassFlags.push_back(false); continue; }

        if( TMath::Abs(eta)>2.4 ){ fJetFilterPassFlags.push_back(true); continue; }

        // sumMult
        if( sumMult <= 0 ){ fJetFilterPassFlags.push_back(false); continue; }

        // chef
        double chef = recoJets->chef[iJet];
        if( chef <= 0 ){ fJetFilterPassFlags.push_back(false); continue; }

        // eef
        double eef  = recoJets->eef[iJet];
        if( eef >= 0.9 ){ fJetFilterPassFlags.push_back(false); continue; }

        // pass
        fJetFilterPassFlags.push_back(true);
    }
    if( fJetFilterPassFlags.size() != fPrimitiveEvent->fJets->nJets )
        std::cerr << "ERROR: Size of fJetFilterPassFlags does not match the number of reco Jets in fPrimitiveEvent" << std::endl;
}

void TL1EventClass::SumsFilter()
{
    fMetFilterPassFlag = fPrimitiveEvent->fMetFilters->hbheNoiseFilter;
}

void TL1EventClass::GetRecalcL1Mht()
{
    TVector2 * mht(new TVector2(0.0,0.0));
    int jetCount=0;
    for(unsigned iJet=0; iJet<fL1JetEt.size(); ++iJet)
    {
        if( TMath::Abs(fL1JetEta[iJet]) > 3.0 ) continue;
        if( fL1JetEt[iJet] < 30.0 ) continue;
        TVector2 * jet(new TVector2(0.0,0.0));
        jet->SetMagPhi(fL1JetEt[iJet], fL1JetPhi[iJet]);
        ++jetCount;
        *(mht) -= *(jet);
    }
    fRecalcL1Mht = mht->Mod();
    fRecalcL1MhtPhi = mht->Phi();
    fNJetL1Mht = jetCount;
}

void TL1EventClass::GetRecalcL1Ett()
{
    int iEtaMax = 28;
    Double_t sumEt = 0.0;
    int jetCount = 0;
    for(int jTower=0; jTower<fPrimitiveEvent->fCaloTowers->nTower; ++jTower)
    {
        int ieta = fPrimitiveEvent->fCaloTowers->ieta[jTower];
        int iet  = fPrimitiveEvent->fCaloTowers->iet[jTower];
        if( TMath::Abs(ieta) <= iEtaMax )
        {
            sumEt += 0.5 * (double)iet;
            ++jetCount;
        }
    }
    fRecalcL1Ett = sumEt;
    //fNJetsL1Ett = jetCount;
}

void TL1EventClass::GetRecalcL1Met()
{
    MissingEt met, metHF;
    fRecalcL1Met28Only.Reset(); fRecalcL1MetNot28.Reset(); fRecalcL1MetNot28HF.Reset();
    auto caloTowers = fPrimitiveEvent->fCaloTowers;
    for(int jTower=0; jTower<caloTowers->nTower; ++jTower)
    {
        int ieta = caloTowers->ieta[jTower];
        double phi = (TMath::Pi()/36.0) * (double)caloTowers->iphi[jTower];
        double et = 0.5 * (double)caloTowers->iet[jTower];

        if( abs(ieta) < 28 ) {
            met.AddEt(et,phi);
            fRecalcL1MetNot28.AddEt(et,phi);
            fRecalcL1MetNot28HF.AddEt(et,phi);
        } else if( abs(ieta) == 28 ){
            met.AddEt(et,phi);
            fRecalcL1Met28Only.AddEt(et,phi);
        } else{
            fRecalcL1MetNot28HF.AddEt(et,phi);
        }
        metHF.AddEt(et,phi);
    }
    fRecalcL1Met = met.met();
    fRecalcL1MetPhi = met.phi();
    fRecalcL1MetHF = metHF.met();
    fRecalcL1MetPhiHF = metHF.phi();
}

void TL1EventClass::GetRecalcL1EmuMet() {
    MissingEt met,metHF;
    auto caloTowers = fPrimitiveEvent->fEmuCaloTowers;
    const int nTow = caloTowers->nTower;

    int nTT4=0;
    for(int jTower=0; jTower<caloTowers->nTower; ++jTower)
        if(abs(caloTowers->ieta[jTower])<5) ++nTT4;

    static const double etaSizes[] = 
            {1.72, 2.05, 1.72, 4.02, 3.29, 2.01, 2.02, 2.00,  
             2.02, 2.  , 2.03, 1.99, 2.01, 2.03, 1.99, 3.47}; 

    for(int jTower=0; jTower<caloTowers->nTower; ++jTower) {
        const int ieta = caloTowers->ieta[jTower];
        const double phi = (TMath::Pi()/36.0) * (double)caloTowers->iphi[jTower];
        const double et = 0.5 * (double)caloTowers->iet[jTower];
        double et_pus28only=et;
        double et_pusAll=et;

        // Calculate PUS
        if(abs(ieta)==28){
            //et_pus28only = 1 + et - 0.5*round((1./8)*pow((et-0.5),0.5)*(floor(nTow/128)));
            et_pus28only = et - 0.5*round(7*(nTT4/100.)); 
            if(et_pus28only<0) et_pus28only = 0;
        }
        if( abs(ieta) > 24 ){
            et_pusAll = et - 0.5*round((etaSizes[abs(ieta)-25]/4)*7*(nTT4/100.));
            if(et_pusAll < 0.) et_pusAll = 0.;
        }

        // MET sums
        const double threshold_et=0.5*(round(abs(ieta)/10.)+1) ;
        if( abs(ieta) <= 28 ){
            fRecalcL1EmuMet.AddEt(et,phi);
            fRecalcL1EmuMetPUS.AddEt(et_pusAll,phi);
            fRecalcL1EmuMetPUS28.AddEt(et_pus28only,phi);
            if (et > threshold_et){
                fRecalcL1EmuMetPUSThresh.AddEt(et_pusAll,phi);
            }
            if( abs(ieta) < 28 ){
                fRecalcL1EmuMetNot28.AddEt(et,phi);
            } else if( abs(ieta) == 28 ){
                fRecalcL1EmuMet28Only.AddEt(et,phi);
            }
        }
        fRecalcL1EmuMetHF.AddEt(et,phi);
        fRecalcL1EmuMetPUSHF.AddEt(et_pusAll,phi);
        if (et > threshold_et){
            fRecalcL1EmuMetPUSThreshHF.AddEt(et_pusAll,phi);
        }
    }
}

void TL1EventClass::GetRecalcRecoHtSums()
{
    TVector2 * mht(new TVector2(0.,0.));
    double jetEt(0.0);
    unsigned jetCount(0);
    fMhtPassFlag = true;
    for(int iJet=0; iJet<fPrimitiveEvent->fJets->nJets; ++iJet)
    {
        bool passJetFilter = fJetFilterPassFlags[iJet];
        if( fPrimitiveEvent->fJets->etCorr[iJet] < 30.0 ) continue;
        if( !passJetFilter )
        {
            fMhtPassFlag = false;
            break;
        }
        TVector2 * jet(new TVector2(0.,0.));
        jet->SetMagPhi(fPrimitiveEvent->fJets->etCorr[iJet], fPrimitiveEvent->fJets->phi[iJet]);
        jetEt += fPrimitiveEvent->fJets->etCorr[iJet];
        ++jetCount;
        *(mht) -= *(jet);
    }
    fRecalcRecoMht = mht->Mod();
    fRecalcRecoMhtPhi = mht->Phi();
    fRecalcRecoHtt = jetEt;
    fNJetRecoMht = jetCount;
}

void TL1EventClass::GetRecalcRecoEtt()
{
    double jetEt(0.0);
    int jetCount = 0;
    for(int iJet=0; iJet<fPrimitiveEvent->fJets->nJets; ++iJet)
    {
        if( fPrimitiveEvent->fJets->etCorr[iJet] >= 30.0 )
        {
            jetEt += fPrimitiveEvent->fJets->etCorr[iJet];
            ++jetCount;
        }
    }
    fRecalcRecoEtt = jetEt;
    //fNJetsRecoEtt = jetCount;
}

void TL1EventClass::GetLeadingRecoJet()
{
    int iLeadingReco = 0;
    double leadingRecoEt = 10.0;
    fIsLeadingRecoJet = false;
    for(int iReco=0; iReco<fJetFilterPassFlags.size(); ++iReco)
    {
        // If central jet, require tightLepVeto and muMult==0
        auto recoJet = fPrimitiveEvent->fJets;
        if( TMath::Abs(recoJet->eta[iReco])<=3.0 )
            if( !fJetFilterPassFlags[iReco] ) continue;
        if( recoJet->etCorr[iReco] < leadingRecoEt ) continue;
        iLeadingReco = iReco;
        leadingRecoEt = recoJet->etCorr[iReco];
    }
    if( leadingRecoEt > 10.0 )
    {
        fIsLeadingRecoJet = true;
        fLeadingRecoJetIndex = iLeadingReco;
    }
}

void TL1EventClass::GetMatchedL1Jet(const std::string & l1Type)
{
    // No leading jet, no match
    if( l1Type == "HW" )
        fIsMatchedL1Jet = false;
    else if( l1Type == "EMU" )
        fIsMatchedL1EmuJet = false;
    if( !fIsLeadingRecoJet ) return;

    // Keep track of the current matched l1 jet. Set minDeltaR to 0.3
    int iMatchedL1 = 0;
    double minDeltaR = 0.3;

    // Get Leading jet parameters
    auto recoJet = fPrimitiveEvent->fJets;
    double leadingRecoJetEta = recoJet->eta[fLeadingRecoJetIndex]; 
    double leadingRecoJetPhi = recoJet->phi[fLeadingRecoJetIndex];

    std::vector<double> l1JetEts, l1JetEtas, l1JetPhis;
    if(l1Type == "HW")
    {
        l1JetEts = fL1JetEt;
        l1JetEtas = fL1JetEta;
        l1JetPhis = fL1JetPhi;
    }
    else if(l1Type == "EMU")
    {
        l1JetEts = fL1EmuJetEt;
        l1JetEtas = fL1EmuJetEta;
        l1JetPhis = fL1EmuJetPhi;
    }
    for(int iL1=0; iL1<l1JetEts.size(); ++iL1)
    {
        // Get l1 jet parameters
        double l1JetEta = l1JetEtas[iL1];
        double l1JetPhi = l1JetPhis[iL1];

        // Check dEta > dR
        double dEta = leadingRecoJetEta - l1JetEta;
        if( dEta > minDeltaR ) continue;
        
        // Check dPhi > dR
        double dPhi = leadingRecoJetPhi - l1JetPhi;
        if( dPhi > minDeltaR ) continue;

        // Check dR < mindR
        double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        if( dR < minDeltaR )
        {
            // Update new minDeltaR and matched L1 Jet
            minDeltaR = dR;
            iMatchedL1 = iL1;
        }
    }
    // If minDeltaR has changed then we have a match
    if( minDeltaR < 0.3 )
    {
        if( l1Type == "HW" )
        {
            fIsMatchedL1Jet = true;
            fMatchedL1JetIndex = iMatchedL1;
        }
        else if( l1Type == "EMU" )
        {
            fIsMatchedL1EmuJet = true;
            fMatchedL1EmuJetIndex = iMatchedL1;
        }
    }
}

#endif
