/*
 * FERS + Scintillators bars utility functions
 *
 * Version: 0.1
 * Date: Oct/2023
 */

#ifndef SCINTIFERS_H
#define SCINTIFERS_H

#include "TH2F.h"
#include "TParameter.h"
#include <TObjString.h>
#include <ctime>
#include <set>
#include <unordered_map>

#include "./ntofDAQ.h"

// bars are 500 x 17 mm2
// SiPM have 16 pixels 3 mm each
// 1 SiPM (= 4 channels) couples to 3 bars
constexpr float half_bar_length = 250.;                                     // [mm] half bar length
constexpr float si_channel_length = 12.;                                    // [mm] aggregated SiPM pixels length (4 pixelx3mm)
constexpr float bar_width = 17.;                                            // [mm] scintillation bar width
constexpr int max_sig = 512;                                                // max number of channels/signals in a single trigger/time reference (for a single board the maximum channels in a single event is 64, both in timing and spectroscopy mode!)
constexpr float light_speed = 299.792458 / 1.5;                             // [mm/ns] light speed in the medium (refractive index is actually 1.58, but for code purposes it is better to approximate)
constexpr float bar_time_width = (half_bar_length * 2.) / light_speed;      // [ns] time equivalent of a scintillator bar length
constexpr int true_n_bar_portions = static_cast<int>(bar_time_width / 0.5); // number of portions in which the scintillator bar can be divided, given 0.5 ns timing resolution
/*
 * MAPPING
 */

/*
Draws a red square around the crossing bars region in a canvas with xy hit map
*/
void draw_square(TCanvas *canvas, const std::string &plot_type);

std::vector<std::vector<float>> get_ToT_pedestal(const std::vector<int> act_ch_v, // vector of active channels
                                                 const float def_rms = 5.,        // [ns] assumed rms where there is no channel data
                                                 const TString basepath = "../../test_231020/data", const std::vector<int> ped_runs = std::vector<int>{13, 22, 33});

/*
 * Return the active channels either from the static "db" or from the provided TTree
 */
std::vector<int> getActiveChannels(TTree *tl = NULL) {
    std::vector<int> vret;
    std::set<int> unique_chs;

    if (tl == NULL) {
        for (int i = 0; i < 32; i++) { // assume active ch 0,1  4,5,  8,9 ...
            int cc = (i / 2) * 4 + (i % 2);
            vret.push_back(cc);
        }
    } else {
        tl->Draw("chan", "Entry$<1000", "goff");
        int nn = tl->GetSelectedRows();
        for (int i = 0; i < nn; ++i) {
            int ch = (int)tl->GetV1()[i];
            unique_chs.insert(ch);
            /*(std::find(vret.begin(), vret.end(), ch) == vret.end())
            {
              vret.push_back(cur);
            }*/
        }
    }
    vret.assign(unique_chs.begin(), unique_chs.end());
    return vret;
}

/*
 * Return axis: 0=x, 1=y
 */

int c2axis(int ch) {
    int axis = (ch / 32); // 0=x (vertical bars), 1=y (horizontal bars)
    return axis;
};

// determine the side of the sipm channel
// return -1 if channel is masked, or 0: top, 1: right, 2: bottom, 3: left
int c2side(int ch) {
    int maskmap[2][4] = {{0, 0, 1, 1}, {0, 0, 1, 1}}; // [axis][0..3] 1 masked, 0 running ;
    int lrudmap[2][2] = {{2, 0}, {1, 3}};             // [axis][even,odd] <- TO BE VERIFIED

    int axis = (ch / 32); // 0 x, 1 y
    int lrud = ch % 2;    // odd: top or left, even: bottom, right
    int masked = maskmap[axis][ch % 4];

    int rval = lrudmap[axis][lrud];
    if (masked)
        rval = -1;

    return rval;
};

// Return the index position (not coordinate) along the side (similar to x,y coordinate)
// do not check active/inactive/masked channels
int c2idx(int ch) {
    int axis = (ch / 32);
    // int lrud = ch % 2;

    int rval = (ch - axis * 32) / 4;
    //  if (lrud) rval = 7-rval;
    return rval;
}

// mapping channel to x / y sipm position
float c2x(int ch) {
    float x = 0;
    float axis = (float)(ch / 32);
    float ipos = 0;
    if (axis > 0) {
        x = half_bar_length - 2 * half_bar_length * (ch % 2);
    } else {
        ipos = (float)(ch / 4);
        x = ipos * si_channel_length - si_channel_length * 3.5;
    }
    return x;
};

float c2y(int ch) {
    float y = 0;
    float axis = (float)(ch / 32);
    float ipos = 0;
    if (axis == 0) {
        y = -half_bar_length + 2 * half_bar_length * (ch % 2);
    } else {
        ipos = (float)((ch - 32) / 4); // cern test, previous iss test 32 -> 34
        y = -ipos * si_channel_length + si_channel_length * 3.5;
    }
    return y;
};

void fill_xy_scintibar(const int bar_index, TH2F &xy_map);

void timing_analysis(const std::string basepath = "../../test_231020/data", const std::vector<int> sig_runs = {29}, const std::vector<int> ped_runs = {13, 22, 33}, float tot_default_rms = 12., float tot_rms_n = 3., float time_window = bar_time_width);

// get the companion paired channel on the other scinti bar end
// assume paired channels are consecutive
int pairedChannel(int ch) {

    int rch = ch + 1;
    if ((c2side(ch - 1) > 0) && (c2idx(ch - 1) == c2idx(ch))) {
        rch = ch - 1;
    }

    return rch;
};

// just a test of the previous functions
int testa() {

    // int a[16]={0,1,4,5,8,9,12,13,16,17,20,21,24,25,28,29};
    int a[16] = {0, 1, 2, 3, 4, 5, 6, 7, 32, 33, 34, 35, 36, 37, 38, 39};
    for (int i = 0; i < 16; i++) {
        int idx = a[i];
        float x = c2x(idx);
        float y = c2y(idx);
        printf("xx %d %f %f  %d \n", idx, x, y, c2side(idx));
        idx = idx + 34;
        x = c2x(idx);
        y = c2y(idx);
        printf("yy %d %f %f  %d\n", idx, x, y, c2side(idx));
    }

    return 0;
};

/*
 * compute coincidence of each side of scinti-bars
 * bsig: each bit corresponds to a channel, if firing is set to 1
 * assume sipm channels on the two side of the scintillation bars are consecutive channels in FERS
 * and only 2 out of 4 channels are connected (e.g. 0,1,4,5,8,9 ...)
 * return position (see ch2idx) of firing channel coincidences ( one bit per positional coincidence)
 */

uint32_t evalCoincSingle(uint32_t axx) {
    uint32_t even = axx & 0x22222222; // 0xAAAAAAAA if all channels are used
    uint32_t odd = axx & 0x11111111;  // 0x55555555 if all channels used
    uint32_t coinc = (even >> 1) & odd;
    uint32_t cpos = 0;                           // bit 0 = coinc first pairs, bit 1 = second pairs ...
    for (int i = 0; i < 8; i++) {                // if all channels used 8 -> 16
        cpos |= ((coinc >> (4 * i) & 0x1) << i); // if all channels used 4 -> 2
    }
    return cpos;
};

/*
 * as previous eval
 * first 32 channels for one axis, second 32 channels on the other axis
 * return the two coincidence words (each bit is a coincidence of the two opposite channels in the bar) of each axis
 */

std::vector<uint32_t> evalCoincidence(uint64_t bsig) {

    std::vector<uint32_t> vret;
    for (int i = 0; i < 2; i++) { // axis
        uint32_t axx = (uint32_t)((bsig >> 32 * i) & 0xffffffff);
        uint32_t coinc = evalCoincSingle(axx);
        //    if (bsig > (((uint64_t) 1) << 32)) { printf("%d %lx -> %x %x %x -> %x\n",i,bsig, axx,even,odd,coinc); }
        vret.push_back(coinc);
    }

    return vret;
}

inline std::unordered_map<int, int> runDelays = {{15, 0}, {16, 500}, {17, 1000}, {18, 1500}, {19, 2000}, {20, 2500}, {21, 3000}, {28, 3500}, {29, 4000}, {26, 3000}, {27, 1000}, {30, 10000}, {31, 0}, {44, 500}, {45, 2000}, {46, 10000}, {48, 10000}, {49, 2000}, {53, 2000}, {54, 2000}, {55, 0}, {56, 2000}, {57, 2000}, {58, 2000}};

struct coinc_counter {
    int run_n;
    int ns_delay;
    int coinc_n;
    int double_coinc_n;
    int cross_coinc_n;
    coinc_counter() : run_n(0), ns_delay(0), coinc_n(0), double_coinc_n(0), cross_coinc_n(0) {}
    coinc_counter(int run_number, int gflash_delay_in_ns, int coincidences_n, int double_coincidences_n, int cross_coincidences_n) : run_n(run_number), ns_delay(gflash_delay_in_ns), coinc_n(coincidences_n), double_coinc_n(double_coincidences_n), cross_coinc_n(cross_coincidences_n) {}
};
/*
 * Collect summary data, generally from a run
 */

class summaryList {
  private:
    TList *iList; // list of integers
    TList *fList; // list of floats ... to be improved merging them
    TParameter<float> *fPar;
    TParameter<int> *iPar;
    std::vector<coinc_counter> counters;

  public:
    summaryList(TString name = "list") {
        iList = new TList();
        fList = new TList();
        iList->SetName("int_" + name);
        fList->SetName("float_" + name);
    };
    ~summaryList() {
        if (iList != NULL) {
            delete iList;
            delete fList;
        }
    };
    void resetLists() {
        iList->Clear();
        fList->Clear();
    };
    void addFloat(float val, TString name) {
        fPar = new TParameter<float>(name, val);
        fList->Add(fPar);
    };
    void addFloat(TParameter<float> *p) { fList->Add(p); };

    void addInt(int val, TString name) {
        iPar = new TParameter<int>(name, val);
        iList->Add(iPar);
    };
    void addInt(TParameter<int> *p) { iList->Add(p); };

    void add_counter(const coinc_counter &counter) { counters.push_back(counter); }
    const std::vector<coinc_counter> &get_counter() const { return counters; }

    std::string getNames(int format = 0) { // 0=just space between names, 1=use TTree format
        std::string rstring = "# ";
        std::string sform = " %s";
        if (format == 1)
            sform = " %s/I";

        for (int i = 0; i < iList->GetSize(); i++) {
            iPar = (TParameter<int> *)iList->At(i);
            rstring += Form(sform.data(), iPar->GetName());
            if (format == 1) {
                sform = ":%s/I";
            }
        }
        if (format == 1)
            sform = ":%s/F";
        for (int i = 0; i < fList->GetSize(); i++) {
            fPar = (TParameter<float> *)fList->At(i);
            rstring += Form(sform.data(), fPar->GetName());
        }
        return rstring;
    };

    std::string getValues() {
        std::string rstring = " ";
        for (int i = 0; i < iList->GetSize(); i++) {
            iPar = (TParameter<int> *)iList->At(i);
            rstring += Form(" %d", iPar->GetVal());
        }
        for (int i = 0; i < fList->GetSize(); i++) {
            fPar = (TParameter<float> *)fList->At(i);
            rstring += Form(" %f", fPar->GetVal());
        }
        return rstring;
    };

    TParameter<float> *getFloat(int idx) {
        if (idx < fList->GetSize()) {
            return (TParameter<float> *)fList->At(idx);
        } else {
            return NULL;
        }
    };

    TParameter<int> *getInt(int idx) {
        if (idx < iList->GetSize()) {
            return (TParameter<int> *)iList->At(idx);
        } else {
            return NULL;
        }
    }

    void merge(summaryList *sl) {
        int i = 0;
        while ((fPar = sl->getFloat(i)) != NULL) {
            addFloat(fPar);
            i++;
        }
        i = 0;
        while ((iPar = sl->getInt(i)) != NULL) {
            addInt(iPar);
            i++;
        }
    };
};

/*
 * This class do not work yet ...
 */

class DAQpar : public TNamed { // public TObject {

  private:
    int run_number;
    int acqmode;
    int histoch;
    float timebin; // [ns]
    TString start_time;

  public:
    DAQpar() { // const char *name, const char *title) {
        //    this->SetName(name);
        //    this->SetTitme(title);
        run_number = -1;
        acqmode = -1;
        histoch = -1;
        timebin = -1;
        start_time = "";
    };

    ~DAQpar(){};

    int runNumber(int run = -1) {
        if (run >= 0) {
            run_number = run;
        }
        return run_number;
    };

    int acqMode(int acm = -1) {
        if (acm >= 0) {
            acqmode = acm;
        }
        return acqmode;
    };

    int histoChannels(int hc = -1) {
        if (hc > 0) {
            histoch = hc;
        }
        return histoch;
    };

    float timeBinWidth(float tb = -1) {
        if (tb > 0) {
            timebin = tb;
        }
        return timebin;
    };

    TString startTime(TString st = "") {
        if (st.Length() > 0) {
            start_time = st;
        }
        return start_time;
    };

    // overidd virtual
    void Print(Option_t *opt = "") const override {
        if (strcmp(opt, "") != 0) {
            std::cout << "Option " << opt << " will be ignored\n";
        }
        std::cout << "DAQ paramereters" << std::endl;
        std::cout << "| Run Number : " << run_number << std::endl;
        std::cout << "| Start Time : " << start_time << std::endl;
        std::cout << "| Acquisition Mode : " << acqmode << " (0:spect, 1:timing, 2:spect_timing)" << std::endl;
        std::cout << "| Energy histo bins: " << histoch << std::endl;
        std::cout << "| Time bin width   : " << timebin << std::endl;
    };

    ClassDef(DAQpar, 1)
};

/*
class manageChain() {
 private:
  TChain *cate;
  int nfiles;

 public:
  manageChain(TChain* chain) {
    cate = chain;

    TObjArray * tob = cate->GetListofFiles();
    nfiles = tob->GetEntries();
    for (int i=0;i<nfiles;i++) {
      TTree

  }

};
*/

/*
 *
 */
class manageTree {

  private:
    // Tree related variables
    TFile *fout;
    TTree *tlist;

    double xb, yb;

    int countEvents;

    int amode; // acquisition mode

    // ttree global variable, to be improved

    float timeus;
    int trgid;
    int novt;  // number of over threshold signals or with charge
    int *chan; // board*64 + channel
    int *lgain;
    int *hgain;
    int *tot;
    int *toa;
    double *xc, *yc; // centroid on each side/axis
    double *qt;      // total charge on each side/axis

    float beamDelta; // beam delta time of arrival relative tu timeus
    float beamPulse; // beam pulse from ntof data

    TList *userlist;

  public:
    manageTree(TString sofile) {
        fout = new TFile(sofile, "recreate");
        tlist = new TTree("tlist", "Event Lists");
        countEvents = 0;

        chan = new int[max_sig];
        lgain = new int[max_sig];
        hgain = new int[max_sig];
        tot = new int[max_sig];
        toa = new int[max_sig];

        xc = new double[4];
        yc = new double[4];
        qt = new double[4];

        novt = 0;
    }

    ~manageTree() {
        delete[] chan;
        delete[] lgain;
        delete[] hgain;
        delete[] tot;
        delete[] toa;
        delete[] xc;
        delete[] yc;
        delete[] qt;
        //    delete tlist;
    }

    void allocBranches(int acqmode) { // 0: spectroscopy, 1: timing 2:spect_timing

        amode = acqmode;

        tlist->Branch("timeus", &timeus, "timeus/F");
        tlist->Branch("trigger", &trgid, "trgid/I");
        tlist->Branch("novt", &novt, "novt/I"); // number of signals
        tlist->Branch("chan", chan, "chan[novt]/I");

        if (amode != 1) {
            tlist->Branch("lgain", lgain, "lgain[novt]/I");
            tlist->Branch("hgain", hgain, "hgain[novt]/I");
            tlist->Branch("xcen", xc, "xc[4]/D");
            tlist->Branch("ycen", yc, "yc[4]/D");
            tlist->Branch("charge", qt, "qt[4]/D");
        }
        if (amode != 0) {
            tlist->Branch("toa", toa, "toa[novt]/I");
            tlist->Branch("tot", tot, "tot[novt]/I");
        }
        tlist->Branch("bPulse", &beamPulse, "bPulse/F");
        tlist->Branch("bDTime", &beamDelta, "bDTime/F");

        //    branch[6] = tlist->Branch("ybal", &yb, "yb/I");

        printf("  tree branches allocated\n");
    };

    // return true if channel is masked
    bool isChannelMasked(int ch) { return (c2side(ch) == -1) ? true : false; };

    void evalCentroids() {

        for (int i = 0; i < 64; i++) {
            double charge = (double)hgain[i]; // single board @@@
            // if (charge>5000) printf(" ... still %d %d %f\n",i,hgain[i],charge);
            int cyc = c2side(i);
            if (cyc < 0)
                continue;
            double x = c2x(i);
            double y = c2y(i);
            xc[cyc] += x * charge;
            yc[cyc] += y * charge;
            qt[cyc] += charge;
        }
        for (int j = 0; j < 4; j++) {
            if (qt[j] > 0) {
                xc[j] = xc[j] / qt[j];
                yc[j] = yc[j] / qt[j];
            }
        }

        return;
    };

    void resetVars() {
        novt = 0;
        for (int j = 0; j < 4; j++) {
            xc[j] = 0;
            yc[j] = 0;
            qt[j] = 0;
        }
    };

    void setTimeTrg(float time, float trg) {
        timeus = time;
        trgid = (int)trg;
    };

    void setChannel(float board, float channel) {
        if (novt >= max_sig) {
            printf("WARNING: number of signals/channels exceed current limit %d, no further signals considered\n", novt);
        } else {
            chan[novt] = ((int)board) * 32 + ((int)channel);
            novt += 1;
        }
    };

    void setADCs(float lgadc, float hgadc) {
        lgain[novt - 1] = (int)lgadc;
        hgain[novt - 1] = (int)hgadc;
    };

    void setTime(float toav, float totv) {
        toa[novt - 1] = (int)toav;
        tot[novt - 1] = (int)totv;
    };

    void setVars(float *arr) {
        setChannel(arr[0], arr[1]);
        switch (amode) {
        case 0: // spectroscopy
            setADCs(arr[2], arr[3]);
            break;
        case 1: // timing
            setTime(arr[2], arr[3]);
            break;
        case 2: // spect_timing
            setADCs(arr[2], arr[3]);
            setTime(arr[4], arr[5]);
            break;
        }
    };

    void setBeamInfo(float pulse, float dtime) {
        beamPulse = pulse;
        beamDelta = dtime;
    };

    int proFill() {
        //    printf(" fill %d %d %f\n",novt,chan[0], timeus);
        if (novt <= 0) {
            return 0;
        }
        if ((amode % 2) == 0)
            evalCentroids(); // hgain[0]); // only board 0 @@@
        tlist->Fill();
        resetVars();
        countEvents += 1;
        return 0;
    };

    int run_number;
    int acqmode;
    int histoch;
    float timebin; // [ns]
    TString start_time;

    void setUserData(int runnum, float timebin, int histoch, TString startime) {
        userlist = (TList *)tlist->GetUserInfo();
        TParameter<int> *tp0 = new TParameter<int>("RunNumber", runnum);
        userlist->Add(tp0);
        TParameter<std::time_t> *tp0a = new TParameter<std::time_t>("StartTime", convertTime(startime));
        userlist->Add(tp0a);
        TParameter<int> *tp0b = new TParameter<int>("AcqMode", amode);
        userlist->Add(tp0b);
        TParameter<float> *tp1 = new TParameter<float>("TimeBinWidth", timebin);
        userlist->Add(tp1);
        TParameter<int> *tp2 = new TParameter<int>("HistoBins", histoch);
        userlist->Add(tp2);
        tlist->GetUserInfo()->Print();
        tlist->GetCurrentFile()->Write();
    };

    TTree *getTree() { return tlist; };

    void Print() { printf(" Trg: %d at time [ms] %f\n", trgid, timeus / 1000); };

    void Save() {
        fout->Write();
        fout->Close();
    };
};

/*
 * Parse FERS output run list text file and save data into root file
 * ntofRFile: root file with ntof pulse data extracted by ntofData.cpp macro
 * return number of trigger loaded/parsed (<0 on error)
 *
 // spectroscopy gain data should go from 0 to up to "Energy N Channels" parameter, maximum is 8192 - 13 bits) but in output ascii file there are values larger than this limit and also approaching 2^15; probbly a bug!
  *
 */

int parseListFile(int run, TString basepath, TString prefix = "/fers/Run", TString ntofRoot = "__NONE__") {

    TString fname = basepath + prefix + Form("%d_list.txt", run);
    TString ntofRFile = basepath + ntofRoot;

    printf("Parse Text file %s\n", fname.Data());

    std::ifstream infile(fname.Data());

    if (infile.fail()) {
        printf("ERROR: file %s not found\n", fname.Data());
        return 0;
    }

    // get ntof root file if available
    manageNTOF *ntof = new manageNTOF(ntofRFile);

    // prepare ttree
    TString ofile = basepath + prefix + Form("%d.root", run);
    manageTree *mT = new manageTree(ofile);

    TString line;
    TObjString *obs;

    TString ss;
    // int headflag = 0;
    TObjArray *oa;

    TString sitem[4] = {"Acquisition Mode:", "Run start time:", "ToA/ToT LSB:", "Energy Histogram Channels:"}; // first element shall be the Acquisition Mode
    TString sival[4] = {"", "", "", ""};                                                                       // values of items in header
    TString spname;                                                                                            // names of the data variables (columns)

    TString acqmode[3] = {"Spectroscopy", "Timing", "Spect_Timing"};
    int acqcol[3] = {6, 6, 8}; // number of columns with Tstamp, depending on acq mode
    int acqidx = -1;           // acq mode index

    int counttrg = 0;

    float adummy[8]; // temporary storage of data

    int ncols = -1; // number of data columns of first line of new events (include Tstam_us and TrgID)

    //  int novt=0;
    // parse text file and fill data in root file

    Long64_t sdati; // start data and time as unix time

    while (line.ReadLine(infile)) {
        oa = NULL;
        if (line.Index("//") >= 0) {      // parse Header
            for (int j = 0; j < 4; j++) { // search relevant items
                int idx = line.Index(sitem[j]);
                if (idx > 0) {
                    sival[j] = line(idx + sitem[j].Length() + 1, 9999);
                    printf("  %s %s\n", sitem[j].Data(), sival[j].Data());
                    if (j == 0) { // load acquisition mode
                        for (int k = 0; k < 3; k++) {
                            if (sival[j] == acqmode[k]) {
                                acqidx = k;
                                ncols = acqcol[k];
                                printf("  %s acquisition mode selected (%d) -> cols %d\n", sival[j].Data(), acqidx, ncols);
                                mT->allocBranches(acqidx);
                                break;
                            }
                        }
                    }
                    if (j == 1) { // start time
                        sdati = convertTime(sival[1]);
                        ntof->setStartDaTi(sdati);
                    }
                    break;
                }
            }
            continue;
        } // end of header parsing

        // data (or titles of columns)

        oa = line.Tokenize(" ");
        int ntos = oa->GetEntries();

        obs = (TObjString *)oa->First();
        ss = obs->GetString();

        if (ss.IsFloat()) {                         // data
            for (int i = (ntos - 1); i >= 0; i--) { // loop on token
                TObjString *obsl = (TObjString *)oa->At(i);
                TString ssl = obsl->GetString();
                // ssl.ReplaceAll("\n","").ReplaceAll("\t","").ReplaceAll(" ","");
                sscanf(ssl.Data(), "%f\n", &adummy[i]);
            }
            int nskip = 0;
            if (ntos == ncols) { // new event with time and trigger id * depend on modality *
                mT->proFill();
                mT->setTimeTrg(adummy[0], adummy[1]);

                std::vector<float> pulse = ntof->getBeamPulse(adummy[0]); // return beam pulse related to the time "adummy[0]"
                mT->setBeamInfo(pulse[0], pulse[1]);
                nskip = 2;
                if ((counttrg % 1000) == 0) {
                    printf("  parsed events so far: %d\n", counttrg);
                }
                counttrg += 1;
            } // end of new event
            //      int board= (int) adummy[nskip]; // not currently used
            //      int channel= (int) adummy[nskip+1];

            mT->setVars(&adummy[nskip]);
        } else { // parameters names (before data)
            spname = line;
            printf("  Par data Names: %s\n", spname.Data()); // first row of data contains the name of the variables
            if (acqidx < 0) {
                printf("ERROR: acquisition mode %s is not supported\n", sival[0].Data());
                return -1;
            }
        }

    } // loop on input file lines

    mT->proFill(); // last event

    printf(" Total Events parsed: %d\n", counttrg);

    // extract and add daq parameters to ttre user info
    int idum;
    sscanf(sival[3].Data(), "%d", &idum);
    float fdum;
    sscanf(sival[2].Data(), "%f", &fdum);

    mT->setUserData(run, fdum, idum, sival[1]);
    mT->Save();

    delete ntof;
    delete mT;

    return counttrg;
}

/*
 * Check if run already in root format,
 * if not try to parse text file (FERS output) and convert to root
 * to force parsing and conversion, use a non existing prefix
 *
 * the optional ntofRFile can be used to synch with the Cube/nTOF data (time and beam pulse to distinguish main and parassitic pulses, TO BE DONE: add the info on the Cube signal itself)
 *
 */

int checkAndConvertText2Root(TString basepath, std::vector<int> vrun, TString prefix = "/fers/Run", TString ntofRFile = "__NONE__") {

    std::cout << "Check and parse text data if needed\n";

    int count = 0;
    // convert each ASCII file to root format if not yet done
    for (uint i = 0; i != vrun.size(); ++i) {
        if (vrun[i] >= 0) {
            TString ifile = basepath + prefix + Form("%d.root", vrun[i]);
            TFile *fin = new TFile(ifile, "read");
            if (fin->IsZombie()) {
                // printf(" %d %s %s\n", vrun[i], basepath.Data(), ntofRFile.Data());
                int nevt = parseListFile(vrun[i], basepath, prefix, ntofRFile);
                if (nevt > 0)
                    ++count;
            } else {
                fin->Close();
            }
        }
    }

    printf(" Parsed %d (out of %zu) text file and converted to root\n", count, vrun.size());

    return count;
};

#endif
