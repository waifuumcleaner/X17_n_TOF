//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 31 21:58:54 2022 (!!!) by ROOT version 6.16/00
// from TTree raw/rawapvdata
// found on file: run13.root
//////////////////////////////////////////////////////////

#ifndef MMRAWTREE_H
#define MMRAWTREE_H

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TString.h>

#include <vector>

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
    TTree *fChain; //! pointer to the analyzed TTree or TChain
    int fCurrent;  //! current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    uint apv_evt;
    int time_s;
    int time_us;
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
    virtual int Cut(Long64_t entry);
    virtual int GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
    virtual void Loop();
    virtual Bool_t Notify();
    virtual void Show(Long64_t entry = -1);

    // ----

    TString buildModNameAxis(TString modid, unsigned int iax) { // module id from mmdaq, iax = 0, 1 axis index
        TString saxis[4] = {'x', 'y', 'u', 'v'};                // u,v not actually used
        TString name = modid + saxis[iax];
        return name;
    }

    int BuildMap(long int ne = 10);
    // !!! need to call BuildMap before using the following !!!
    // return module/chamber numeric index from module string name (mod_id)
    int getModId(TString smod) { // smod = module_name + axis name (x/y)
        for (int i = 0; i < (int)smona.size(); i++) {
            if (smod.EqualTo(smona[i])) {
                return i;
            }
        }
        return -1;
    };
    TString getModName(int idx) { return smona[idx]; };
    unsigned int getAxis(int idx) { return raxis[idx]; };

    int numMods() { return (int)smona.size(); };
    int numAPVs() { return napvs; };
    int numSamples() { return nsamples; };

    // return chamber integer index corresponding to apv_id and apv_ch
    int getModule(int apvid, int apvch) { return cc2module[apvid * 128 + apvch]; };

    // return strip corresponding to apv_id and apv_ch
    int getStrip(int apvid, int apvch) { return cc2strip[apvid * 128 + apvch]; };

    int getStripMin(int idx) { return stripl[idx]; } // idx = module index, return lower strip index of module
    int getStripMax(int idx) { return striph[idx]; } // idx = module index, return higher strip index of module
};

#endif
