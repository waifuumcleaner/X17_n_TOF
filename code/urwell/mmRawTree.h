//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 31 21:58:54 2022 by ROOT version 6.16/00
// from TTree raw/rawapvdata
// found on file: run13.root
//////////////////////////////////////////////////////////

#ifndef mmRawTree_h
#define mmRawTree_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class mmRawTree {
  private:
    int *cc2strip;                   // APV card, APV ch to strip index
    int *cc2module;                  // APV card, APV ch to chamber/module
    std::vector<TString> smona;      // module (chamber) string name
    std::vector<unsigned int> raxis; // radout axis in module
    int nsamples;                    // number of time slots (samples) per event
    int napvs;                       // number of apvs
    int *stripl, *striph;            // min and max strip index for each module

  public:
    TTree *fChain;  //! pointer to the analyzed TTree or TChain
    Int_t fCurrent; //! current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    UInt_t apv_evt;
    Int_t time_s;
    Int_t time_us;
    std::vector<unsigned int> *apv_fecNo;
    std::vector<unsigned int> *apv_id;
    std::vector<unsigned int> *apv_ch;
    std::vector<std::string> *mm_id;
    std::vector<unsigned int> *mm_readout;
    std::vector<unsigned int> *mm_strip;
    std::vector<std::vector<short>> *apv_q;
    UInt_t apv_presamples;

    // List of branches
    TBranch *b_apv_evt;        //!
    TBranch *b_time_s;         //!
    TBranch *b_time_us;        //!
    TBranch *b_apv_fecNo;      //!
    TBranch *b_apv_id;         //!
    TBranch *b_apv_ch;         //!
    TBranch *b_mm_id;          //!
    TBranch *b_mm_readout;     //!
    TBranch *b_mm_strip;       //!
    TBranch *b_apv_q;          //!
    TBranch *b_apv_presamples; //!

    mmRawTree(TTree *tree = 0, TString infile = "runMMRT.root");
    virtual ~mmRawTree();
    virtual Int_t Cut(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
#ifdef mmRawTest_cxx
    virtual void Loop();
#endif
    virtual Bool_t Notify();
    virtual void Show(Long64_t entry = -1);

    // ----

    TString buildModNameAxis(TString modid, unsigned int iax) { // module id from mmdaq, iax = 0, 1 axis index
        TString saxis[4] = {'x', 'y', 'u', 'v'};                // u,v not actually used
        TString name = modid + saxis[iax];
        return name;
    }

    Int_t BuildMap(long int ne = 10);
    // !!! need to call BuildMap before using the following !!!
    // return module/chamber numeric index from module string name (mod_id)
    Int_t getModId(TString smod) { // smod = module_name + axis name (x/y)
        for (int i = 0; i < (int)smona.size(); i++) {
            if (smod.EqualTo(smona[i])) {
                return i;
            }
        }
        return -1;
    };
    TString getModName(int idx) { return smona[idx]; };
    unsigned int getAxis(int idx) { return raxis[idx]; };

    Int_t numMods() { return (Int_t)smona.size(); };
    Int_t numAPVs() { return napvs; };
    Int_t numSamples() { return nsamples; };

    // return chamber integer index corresponding to apv_id and apv_ch
    Int_t getModule(Int_t apvid, Int_t apvch) { return cc2module[apvid * 128 + apvch]; };

    // return strip corresponding to apv_id and apv_ch
    Int_t getStrip(Int_t apvid, Int_t apvch) { return cc2strip[apvid * 128 + apvch]; };

    Int_t getStripMin(int idx) { return stripl[idx]; } // idx = module index, return lower strip index of module
    Int_t getStripMax(int idx) { return striph[idx]; } // idx = module index, return higher strip index of module
};

#endif

#ifdef mmRawTree_cxx
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
                    printf(" ERROR: strip %d for apv id/ch %d %d does not mach previous assignment %d\n",
                           jstrip, jid, jch, cc2strip[lin]);
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

#endif // #ifdef mmRawTree_cxx
