#include "mmRawTree.h"
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <iostream>

mmRawTree::mmRawTree(TTree *tree, TString infile) : fChain(0) {
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.

    nsamples = 0; // use BuildMap to get the proper number
    napvs = 0;
    cc2module = NULL;
    cc2strip = NULL;

    if (tree == 0) {
        TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(infile);
        if (!f || !f->IsOpen()) {
            f = new TFile(infile);
        }
        f->GetObject("raw", tree);
    }
    Init(tree);
}

mmRawTree::~mmRawTree() {
    if (!fChain)
        return;
    delete fChain->GetCurrentFile();
    if (smona.size() > 0) { // at list one chamber/module available -> means the map as been built
        delete cc2strip;
        delete cc2module;
    }
}

Int_t mmRawTree::GetEntry(Long64_t entry) {
    // Read contents of entry.
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}
Long64_t mmRawTree::LoadTree(Long64_t entry) {
    // Set the environment to read one entry
    if (!fChain)
        return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0)
        return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void mmRawTree::Init(TTree *tree) {
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    apv_fecNo = 0;
    apv_id = 0;
    apv_ch = 0;
    mm_id = 0;
    mm_readout = 0;
    mm_strip = 0;
    apv_q = 0;
    // Set branch addresses and branch pointers
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("apv_evt", &apv_evt, &b_apv_evt);
    fChain->SetBranchAddress("time_s", &time_s, &b_time_s);
    fChain->SetBranchAddress("time_us", &time_us, &b_time_us);
    fChain->SetBranchAddress("apv_fecNo", &apv_fecNo, &b_apv_fecNo);
    fChain->SetBranchAddress("apv_id", &apv_id, &b_apv_id);
    fChain->SetBranchAddress("apv_ch", &apv_ch, &b_apv_ch);
    fChain->SetBranchAddress("mm_id", &mm_id, &b_mm_id);
    fChain->SetBranchAddress("mm_readout", &mm_readout, &b_mm_readout);
    fChain->SetBranchAddress("mm_strip", &mm_strip, &b_mm_strip);
    fChain->SetBranchAddress("apv_q", &apv_q, &b_apv_q);
    fChain->SetBranchAddress("apv_presamples", &apv_presamples, &b_apv_presamples);
    Notify();
}

Bool_t mmRawTree::Notify() {
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void mmRawTree::Show(Long64_t entry) {
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain)
        return;
    fChain->Show(entry);
}
Int_t mmRawTree::Cut(Long64_t entry) {
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    std::cout << "Don't know what to do with this function yet, entry " << entry << " was given\n";
    return 1;
}

/***
 * build map from raw TTree data; assume NO ZERO suppression!
 * assume APV indeces are consecutive and less then napv
 *
 * input par:
 *   nentries: number of TTree entries to consider to extract the map
 * return number of chamber/modules found
 *
 */

Int_t mmRawTree::BuildMap(long int nentries) {

    int nach = 128; // number of channels on each APV

    nsamples = 0; // initial number of samples

    printf("### Extract Map from TTree raw data from first %ld entries\n", nentries);

    for (Long64_t i = 0; i < nentries; i++) { // loop on entries
        Long64_t j = LoadTree(i);
        if (j < 0) {
            break;
        }
        GetEntry(i);

        Long64_t ecount = apv_id->size();

        if (smona.size() == 0) { // first loop, initialize variables
            // find number of APV cards
            int napv0 = 999999;
            int napv1 = -1;

            for (Long64_t k = 0; k < ecount; k++) {
                int jjj = (int)apv_id->at(k);
                napv0 = (jjj < napv0) ? jjj : napv0;
                napv1 = (jjj > napv1) ? jjj : napv1;
            }
            printf(" Min / Max APV card indices = %d %d\n", napv0, napv1);
            napvs = napv1 - napv0 + 1;
            cc2module = new int[nach * napvs];
            cc2strip = new int[nach * napvs];
            for (int i = 0; i < nach * napvs; i++) {
                cc2module[i] = -1;
                cc2strip[i] = -1;
            }
            printf(" Mapping vector initialized\n");
        }

        for (Long64_t k = 0; k < ecount; k++) { // loop on all strips (apv id and relative channels)

            unsigned int jid = apv_id->at(k);       // APV id
            unsigned int jch = apv_ch->at(k);       // APV ch
            unsigned int jstrip = mm_strip->at(k);  // Strip
            unsigned int jaxis = mm_readout->at(k); // Axis  (assume 2 axis max)
            TString jmod = mm_id->at(k);            // module name

            int nnss = (int)apv_q->at(k).size(); // get number of samples
            if (nsamples == 0) {
                nsamples = nnss;
            } else {
                if (nsamples != nnss) {
                    printf(" WARNING: current number of samples %d does not match previous values %d\n", nnss, nsamples);
                }
            }

            int lin = jid * 128 + jch;

            if (cc2strip[lin] < 0) {
                cc2strip[lin] = jstrip;
            } else {
                if (cc2strip[lin] != (int)jstrip) {
                    printf(" ERROR: strip %d for apv id/ch %d %d does not mach previous assignment %d\n", jstrip, jid, jch, cc2strip[lin]);
                }
            }

            int idx = -1;
            for (int m = 0; m < (int)smona.size(); m++) {
                if (jmod.EqualTo(smona[m]) && (jaxis == raxis[m])) {
                    idx = m;
                    break;
                }
            }
            if (idx < 0) {
                smona.push_back(jmod);
                raxis.push_back(jaxis);
                idx = smona.size() - 1;
            }
            cc2module[lin] = idx;
        }
    }

    printf(" Found %d samples/event\n", nsamples);

    printf(" Found %zu modules (chambers):", smona.size());
    int nmod = smona.size();
    stripl = new int[nmod];
    striph = new int[nmod];
    for (int i = 0; i < (int)nmod; i++) {
        stripl[i] = 999999;
        striph[i] = -1;
    }
    printf("\n");

    // check all channels are mapped and evaluate strips limits on each chamber (module)
    for (int i = 0; i < nach * napvs; i++) {
        int imod = cc2module[i];
        int istr = cc2strip[i];
        if (imod >= 0) {
            stripl[imod] = (stripl[imod] > istr) ? istr : stripl[imod];
            striph[imod] = (striph[imod] < istr) ? istr : striph[imod];
        } else {
            printf(" WARNING: APV %d Channel %d does not map to chamber\n", (int)i / 128, i % 128);
        }
        if (istr < 0) {
            printf(" WARNING: APV %d Channel %d does not map to strip\n", (int)i / 128, i % 128);
        }
    }

    for (int i = 0; i < nmod; i++) {
        printf("  Chamber %s axis %d has strips from %d to %d\n", smona[i].Data(), raxis[i], stripl[i], striph[i]);
    }
    return smona.size();
}

void mmRawTree::Loop() {
    //   In a ROOT session, you can do:
    //      root> .L mmRawTree.C
    //      root> mmRawTree t
    //      root> t.GetEntry(12); // Fill t data members with entry number 12
    //      root> t.Show();       // Show values of entry 12
    //      root> t.Show(16);     // Read and show values of entry 16
    //      root> t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    // by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0)
        return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        // if (Cut(ientry) < 0) continue;
    }
}
