/*
 * methods for SiPM CAEN readout output
 *
 * First Version Jul/2023
 *
 */

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <span>
#include <stdexcept>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TEntryList.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>
#include <TTree.h>

#include "../utils.h"
#include "sipm_analysis.h"

// FROM SHELL:  g++ anasipm_timing.cpp  -Wall -Wextra -g -O3 `root-config --cflags --libs`
//              ./a.out

void set_local_style() {
    // gStyle->SetTitleW(0.4);
    // gStyle->SetTitleH(0.07);
    // gStyle->SetTitleSize(0.06, "xyzt");
    // gStyle->SetTitleXOffset(.8);
    // gStyle->SetTitleYOffset(.8);
    // gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadTopMargin(.1);
    gStyle->SetPadLeftMargin(.1);
    gStyle->SetPadRightMargin(.05);
    gStyle->SetPadBottomMargin(.1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
}

coinc_counter::coinc_counter() : run_n(0), ns_delay(0), trig_n(0), coinc_n(0), double_coinc_n(0), cross_coinc_n(0) {}

coinc_counter::coinc_counter(int run_number, int gflash_delay_in_ns, int trigger_n, int coincidences_n, int double_coincidences_n, int cross_coincidences_n) : run_n(run_number), ns_delay(gflash_delay_in_ns), trig_n(trigger_n), coinc_n(coincidences_n), double_coinc_n(double_coincidences_n), cross_coinc_n(cross_coincidences_n) {}

coincidence_Info::coincidence_Info() : counters(std::vector<coinc_counter>()) {}

void coincidence_Info::add_counter(const coinc_counter &counter) { counters.push_back(counter); }

const std::vector<coinc_counter> &coincidence_Info::get_counter() const { return counters; }

summaryList::summaryList(TString name) {
    iList = new TList();
    fList = new TList();
    iList->SetName("int_" + name);
    fList->SetName("float_" + name);
}
summaryList::~summaryList() {
    if (iList != NULL) {
        delete iList;
        delete fList;
    }
}
void summaryList::resetLists() {
    iList->Clear();
    fList->Clear();
}
void summaryList::addFloat(float val, TString name) {
    fPar = new TParameter<float>(name, val);
    fList->Add(fPar);
}
void summaryList::addFloat(TParameter<float> *p) { fList->Add(p); }
void summaryList::addInt(int val, TString name) {
    iPar = new TParameter<int>(name, val);
    iList->Add(iPar);
}
void summaryList::addInt(TParameter<int> *p) { iList->Add(p); }
std::string summaryList::getNames(int format) { // 0=just space between names, 1=use TTree format
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
}
std::string summaryList::getValues() {
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
}
TParameter<float> *summaryList::getFloat(int idx) {
    if (idx < fList->GetSize()) {
        return (TParameter<float> *)fList->At(idx);
    } else {
        return NULL;
    }
}
TParameter<int> *summaryList::getInt(int idx) {
    if (idx < iList->GetSize()) {
        return (TParameter<int> *)iList->At(idx);
    } else {
        return NULL;
    }
}
void summaryList::merge(summaryList *sl) {
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
}

std::vector<int> getActiveChannels(TTree *tl) {
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
}

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
}

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
}

manageTree::manageTree(TString sofile) {
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

manageTree::~manageTree() {
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

void manageTree::allocBranches(int acqmode) { // 0: spectroscopy, 1: timing 2:spect_timing

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
}

bool manageTree::isChannelMasked(int ch) { return (c2side(ch) == -1) ? true : false; }

void manageTree::evalCentroids() {

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
}

void manageTree::resetVars() {
    novt = 0;
    for (int j = 0; j < 4; j++) {
        xc[j] = 0;
        yc[j] = 0;
        qt[j] = 0;
    }
}

void manageTree::setTimeTrg(float time, float trg) {
    timeus = time;
    trgid = (int)trg;
}

void manageTree::setChannel(float board, float channel) {
    if (novt >= max_sig) {
        printf("WARNING: number of signals/channels exceed current limit %d, no further signals considered\n", novt);
    } else {
        chan[novt] = ((int)board) * 32 + ((int)channel);
        novt += 1;
    }
}

void manageTree::setADCs(float lgadc, float hgadc) {
    lgain[novt - 1] = (int)lgadc;
    hgain[novt - 1] = (int)hgadc;
}

void manageTree::setTime(float toav, float totv) {
    toa[novt - 1] = (int)toav;
    tot[novt - 1] = (int)totv;
}

void manageTree::setVars(float *arr) {
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
}

void manageTree::setBeamInfo(float pulse, float dtime) {
    beamPulse = pulse;
    beamDelta = dtime;
}

int manageTree::proFill() {
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
}

void manageTree::setUserData(int runnum, float timebin, int histoch, TString startime) {
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
}

TTree *manageTree::getTree() { return tlist; }

void manageTree::Print() { printf(" Trg: %d at time [ms] %f\n", trgid, timeus / 1000); }

void manageTree::Save() {
    fout->Write();
    fout->Close();
}

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

int checkAndConvertText2Root(TString basepath, std::vector<int> vrun, TString prefix, TString ntofRFile) {

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
}

int find_first_empty_bin_after_max(TH1F *histogram) {
    int max_bin = histogram->GetMaximumBin();
    int nbins = histogram->GetNbinsX();

    for (int i = max_bin; i != nbins; ++i) {
        if (histogram->GetBinContent(i + 1) == 0) {
            return i;
        }
    }
    // if it gets to here it means that no drop to zero was found!
    return nbins;
}

/*
void adjust_y_scale(TH1F *h, TPad *pad) {
    double h_max = h->GetMaximum();
    double pad_max = pad->GetUymax(); // Get upper y-coordinate of the pad
    if (h_max > pad_max) {
        h->GetYaxis()->SetRangeUser(0, h_max * 1.2); // Set y-axis range
    }
}
*/

std::vector<std::vector<float>> get_ToT_pedestal(const std::vector<int> act_ch_v, const float def_rms, const std::string basepath, const std::vector<int> ped_runs) {
    int ch = 0;
    float mean = 0.;
    float rms = 0.;
    int entries = 0;
    const int nchs = act_ch_v.size();
    std::vector<std::vector<float>> ped_v(4, std::vector<float>(nchs)); // channel, pedestal mean, rms and total events

    const std::string timing_ped_path = basepath + "/fers/timing_pedestal.root";
    std::ifstream file_check(timing_ped_path);
    if (file_check.good()) {
        std::cout << "\nReading pedestal from root file . . .\n";
        std::unique_ptr<TFile> ped_root_file(TFile::Open(timing_ped_path.c_str(), "READ"));
        if (!ped_root_file || ped_root_file->IsZombie()) {
            std::cerr << "Error opening file" << '\n';
            exit(-1);
        }
        auto ped_tree = ped_root_file->Get<TTree>("ped_tree");
        ped_tree->SetBranchAddress("ch", &ch);
        ped_tree->SetBranchAddress("mean", &mean);
        ped_tree->SetBranchAddress("rms", &rms);
        ped_tree->SetBranchAddress("entries", &entries);
        for (int i = 0; ped_tree->LoadTree(i) >= 0; ++i) {
            ped_tree->GetEntry(i);
            ped_v[0][i] = ch;
            ped_v[1][i] = mean;
            ped_v[2][i] = rms;
            ped_v[3][i] = entries;
        }
    } else {
        if (!ped_runs.empty()) {
            std::cout << "\nMaking pedestal and writing in root file . . .\n";
            const int mincount = 20; // at least this number of counts in the channel for pedestal estimation, otherwise use default
            std::unique_ptr<TFile> ped_root_file(TFile::Open(timing_ped_path.c_str(), "RECREATE"));
            auto ped_tree = std::make_unique<TTree>("ped_tree", "Pedestal ToT TTree");
            ped_tree->Branch("ch", &ch);
            ped_tree->Branch("mean", &mean);
            ped_tree->Branch("rms", &rms);
            ped_tree->Branch("entries", &entries);

            TChain *pedestal_runs = new TChain("tlist"); // class that inherits from TTree to access separate files’ tree parts as one large tree
            const std::string prefix = "/fers/Run";
            for (int prun : ped_runs) {
                const std::string pedrun_root = basepath + prefix + std::to_string(prun) + ".root";
                std::cout << pedrun_root << '\n';
                if (std::filesystem::exists(pedrun_root)) {
                    pedestal_runs->Add(pedrun_root.c_str());
                } else {
                    throw std::invalid_argument("Invalid pedestal run root file name provided: " + pedrun_root + " does not exist\n");
                }
            }

            TCanvas *c_ped = new TCanvas("c_ped", "Single Channel ToT - Pedestal", 0, 0, 1920, 1000);
            c_ped->Divide(2, 8, 0.0001, 0.002); // [x,y][pos]
            c_ped->Update();
            TCanvas *c_temp = new TCanvas("c_temp", "Temp canvas for tot pedestal extraction", 0, 0, 1920, 1000);
            c_temp->Divide(2, 8, 0.0001, 0.002); // [x,y][pos]
            c_temp->Update();
            TText *text = new TText(0, 0, "ch_mean_rms_entries");
            text->SetTextColor(1);
            text->SetTextSize(0.15);

            for (int i = 0; i != nchs; ++i) {
                ch = act_ch_v[i];
                int couple_number = ch / 4; // same number for couples at the end of scintillator bars
                                            // (0 for 0 and 1, 1 for 4 and 5, 2 for 8 and 9, ...)
                c_temp->cd(couple_number + 1);
                gPad->SetLogy();

                pedestal_runs->Draw("tot>>htot", Form("chan==%d", ch), "histo");
                entries = pedestal_runs->GetSelectedRows();
                TH1F *hped = (TH1F *)((TH1F *)gPad->GetPrimitive("htot"))->Clone("hped");

                int ped_edge_bin = 0;
                if (entries > mincount) {
                    ped_edge_bin = find_first_empty_bin_after_max(hped);
                }
                TH1F *hcutped = (TH1F *)hped->Clone("hcutped");
                hcutped->Reset();
                for (int i = 0; i != ped_edge_bin; ++i) {
                    hcutped->SetBinContent(i, hped->GetBinContent(i));
                }
                c_ped->cd(couple_number + 1); // plot couples at the end of scintillator
                                              // bars in the same graph
                gPad->SetLogy();

                std::string title = "ToT Distribution, channel " + std::to_string(couple_number * 4) + ((ch != 28) ? (" and " + std::to_string(couple_number * 4 + 1)) : (" (" + std::to_string(couple_number * 4 + 1) + " is dead)")) + "\n";
                hcutped->SetTitle(title.c_str());
                hcutped->SetLineColor((ch % 2) + 1);
                hcutped->Draw("histosame");
                hped->SetTitle(title.c_str());
                hped->SetLineColor((ch % 2) + 1);
                hped->SetLineStyle(2);
                hped->Draw("histosame");

                mean = 0;
                rms = def_rms;
                if (entries > mincount) {
                    mean = hcutped->GetMean();
                    rms = hcutped->GetRMS();
                }

                text->SetTextColor((ch % 2) + 1);
                float yn = 0.5 + 0.2 * static_cast<float>(1 - (ch % 2));
                text->DrawTextNDC(0.7, yn, Form("%d : %.1f +/- %.1f (%d)", ch, mean, rms, entries));

                ped_v[0][i] = ch;
                ped_v[1][i] = mean;
                ped_v[2][i] = rms;
                ped_v[3][i] = entries;
                ped_tree->Fill();
            }

            ped_tree->Write();
            c_ped->Update();
            c_ped->Write();
            c_ped->SaveAs("../code/sipm/ToT_cut.png");
            // ped_root_file->Close();
        } else {
            std::cout << "\nWarning: no pedestal runs were provided, setting "
                         "pedestals to defalut values\n";
            for (int i = 0; i != nchs; ++i) {
                ped_v[0][i] = act_ch_v[i];
                ped_v[1][i] = 0;
                ped_v[2][i] = def_rms;
                ped_v[3][i] = 0;
            }
        }
    };

    return ped_v;
};

void plot_timing_data(TH1F *h_time_diff, TH1F *h_coinc_per_trigger, std::vector<TH1F *> dead_t_ch_even, std::vector<TH1F *> dead_t_ch_odd, TH2F *h_toa_corr, TH1F *h_toa, TGraph *tot_global_corr, TH1F *h_tot_prod, std::vector<TGraph *> &tot_corr_ch, TH2F *hxy_channels, TH2F *hxy_bars, TH2F *hxy_bars_toa, TH2F *hxy_bars_tot, TH2F *hxy_cross_bars, const float thr_toa, const int run_number) {
    TCanvas *cv_time_diff = new TCanvas("cv_time_diff", "Time differences canvas", 0, 0, 1920, 1000);
    h_time_diff->Draw("HIST");
    h_time_diff->SetFillColor(kBlue);
    h_time_diff->SetLineColor(kBlack);
    h_time_diff->SetLineWidth(2);
    h_time_diff->SetStats(kTRUE);
    gStyle->SetOptStat(1110);
    TLine *coinc_thr_line = new TLine(thr_toa, gPad->GetUymin(), thr_toa, gPad->GetUymax());
    coinc_thr_line->SetLineColor(kRed);
    coinc_thr_line->SetLineWidth(4);
    coinc_thr_line->Draw();
    // gPad->SetLogx();
    // gPad->SetLogy();
    cv_time_diff->Update();
    cv_time_diff->SaveAs(Form("../code/sipm/Run_%d/Run_%d_timediff.png", run_number, run_number));
    delete h_time_diff;
    delete cv_time_diff;

    TCanvas *cv_coinc_per_trigger = new TCanvas("cv_coinc_per_trigger", "Coincidences per trigger", 0, 0, 1920, 1000);
    h_coinc_per_trigger->Draw("HIST");
    h_coinc_per_trigger->SetFillColor(kBlue);
    h_coinc_per_trigger->SetLineColor(kBlack);
    h_coinc_per_trigger->SetLineWidth(2);
    h_coinc_per_trigger->SetStats(kTRUE);
    cv_coinc_per_trigger->Update();
    cv_coinc_per_trigger->SaveAs(Form("../code/sipm/Run_%d/Run_%d_coinc.png", run_number, run_number));
    delete h_coinc_per_trigger;
    delete cv_coinc_per_trigger;

    TCanvas *cv_dead_t_even = new TCanvas("cv_dead_t_even", "Time between consecutives hit in a channel - right and bottom channels", 0, 0, 1920, 1000);
    cv_dead_t_even->Divide(4, 4, 0.0001, 0.002);
    for (int i = 0; i != 16; ++i) {
        cv_dead_t_even->cd(i + 1);
        dead_t_ch_even[i]->SetFillColor(kBlue);
        dead_t_ch_even[i]->SetLineColor(kBlack);
        dead_t_ch_even[i]->SetLineWidth(2);
        gPad->SetBottomMargin(0.15);
        // gPad->SetLeftMargin(0.15);
        dead_t_ch_even[i]->Draw("HIST");
    }
    cv_dead_t_even->Update();
    cv_dead_t_even->SaveAs(Form("../code/sipm/Run_%d/Run_%d_dead_time_even.png", run_number, run_number));
    for (auto h : dead_t_ch_even)
        delete h;
    delete cv_dead_t_even;

    TCanvas *cv_dead_t_odd = new TCanvas("cv_dead_t_odd", "Time between consecutives hit in a channel - left and top channels", 0, 0, 1920, 1000);
    cv_dead_t_odd->Divide(4, 4, 0.0001, 0.002);
    for (int i = 0; i != 16; ++i) {
        cv_dead_t_odd->cd(i + 1);
        dead_t_ch_odd[i]->SetFillColor(kBlue);
        dead_t_ch_odd[i]->SetLineColor(kBlack);
        dead_t_ch_odd[i]->SetLineWidth(2);
        gPad->SetBottomMargin(0.15);
        // gPad->SetLeftMargin(0.15);
        dead_t_ch_odd[i]->Draw("HIST");
    }
    cv_dead_t_odd->Update();
    cv_dead_t_odd->SaveAs(Form("../code/sipm/Run_%d/Run_%d_dead_time_odd.png", run_number, run_number));
    for (auto h : dead_t_ch_odd)
        delete h;
    delete cv_dead_t_odd;

    TCanvas *cv_hit_toa = new TCanvas("cv_hit_toa", "Information from Time of Arrival", 0, 0, 2400, 1200);
    cv_hit_toa->SetMargin(0.12, 0.05, 0.15, 0.05);
    cv_hit_toa->SetFrameBorderMode(0);
    cv_hit_toa->Divide(2, 1);

    cv_hit_toa->cd(1);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.19);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.13);
    gPad->SetBorderSize(2);
    gPad->SetFrameBorderMode(0);
    h_toa_corr->SetTitle("");
    h_toa_corr->GetXaxis()->SetLabelFont(42);
    h_toa_corr->GetXaxis()->SetLabelSize(0.04);
    h_toa_corr->GetXaxis()->SetTitleSize(0.04);
    h_toa_corr->GetXaxis()->SetTitleOffset(1.5);
    h_toa_corr->GetYaxis()->SetLabelFont(42);
    h_toa_corr->GetYaxis()->SetLabelSize(0.04);
    h_toa_corr->GetYaxis()->SetTitleSize(0.04);
    h_toa_corr->GetYaxis()->SetTitleOffset(2.5);
    h_toa_corr->Draw("colz text");
    gPad->Modified();
    /*TLine *bisector = new TLine(h_toa_corr->GetXaxis()->GetXmin(),
    h_toa_corr->GetXaxis()->GetXmin(),
    h_toa_corr->GetXaxis()->GetXmax(),
    h_toa_corr->GetXaxis()->GetXmax()); bisector->SetLineColor(kRed);
    bisector->SetLineWidth(3);
    bisector->Draw();*/

    cv_hit_toa->cd(2);
    gPad->Range(-840, -4.9, 8465, 34.4);
    gPad->SetBorderSize(2);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    h_toa->SetTitle("");
    h_toa->GetXaxis()->SetLabelFont(42);
    h_toa->GetXaxis()->SetLabelSize(0.04);
    h_toa->GetXaxis()->SetTitleSize(0.04);
    h_toa->GetXaxis()->SetTitleOffset(1.5);
    h_toa->GetYaxis()->SetLabelFont(42);
    h_toa->GetYaxis()->SetLabelSize(0.04);
    h_toa->GetYaxis()->SetTitleSize(0.04);
    h_toa->GetYaxis()->SetTitleOffset(1.3);
    h_toa->SetFillColor(kBlue);
    h_toa->SetLineColor(kBlack);
    h_toa->SetLineWidth(2);
    h_toa->SetStats(kTRUE);
    h_toa->Draw("HIST");
    gPad->Modified();

    cv_hit_toa->cd();
    cv_hit_toa->Modified();
    cv_hit_toa->SetSelected(cv_hit_toa);
    cv_hit_toa->SaveAs(Form("../code/sipm/Run_%d/Run_%d_hit_toa.png", run_number, run_number));
    delete h_toa_corr;
    delete h_toa;
    delete cv_hit_toa;

    TCanvas *cv_tot_info = new TCanvas("cv_tot_info", "Extracted information from Time over Threshold", 0, 0, 2400, 1200);
    cv_tot_info->SetMargin(0.12, 0.05, 0.15, 0.05);
    cv_tot_info->SetFrameBorderMode(0);
    cv_tot_info->Divide(2, 1);

    cv_tot_info->cd(1);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.19);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.13);
    gPad->SetBorderSize(2);
    gPad->SetFrameBorderMode(0);
    /*TPaveText *pt = new TPaveText(0.32, 0.92, 0.8, 0.993, "blNDC");
    pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(23);
    pt->SetTextFont(42);
    gStyle->SetTitleH(0.1);
    (void)pt->AddText(tot_global_corr->GetTitle());
    pt->Draw();*/
    tot_global_corr->SetTitle("");
    tot_global_corr->GetXaxis()->SetLabelFont(42);
    tot_global_corr->GetXaxis()->SetLabelSize(0.04);
    tot_global_corr->GetXaxis()->SetTitleSize(0.04);
    tot_global_corr->GetXaxis()->SetTitleOffset(1.5);
    tot_global_corr->GetYaxis()->SetLabelFont(42);
    tot_global_corr->GetYaxis()->SetLabelSize(0.04);
    tot_global_corr->GetYaxis()->SetTitleSize(0.04);
    tot_global_corr->GetYaxis()->SetTitleOffset(2);
    // tot_global_corr->GetYaxis()->SetRangeUser(0., 80);
    tot_global_corr->SetMarkerStyle(8);
    tot_global_corr->SetMarkerSize(0.5);
    tot_global_corr->Draw("AP");
    gPad->Modified();

    cv_tot_info->cd(2);
    gPad->Range(-840, -4.9, 8465, 34.4);
    gPad->SetBorderSize(2);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    /*pt = new TPaveText(0.2, 0.925, 0.8, 0.995, "blNDC");
    pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextFont(42);
    (void)pt->AddText(h_tot_prod->GetTitle());
    pt->Draw();*/
    h_tot_prod->SetTitle("");
    h_tot_prod->GetXaxis()->SetLabelFont(42);
    h_tot_prod->GetXaxis()->SetLabelSize(0.04);
    h_tot_prod->GetXaxis()->SetTitleSize(0.04);
    h_tot_prod->GetXaxis()->SetTitleOffset(1.5);
    h_tot_prod->GetYaxis()->SetLabelFont(42);
    h_tot_prod->GetYaxis()->SetLabelSize(0.04);
    h_tot_prod->GetYaxis()->SetTitleSize(0.04);
    h_tot_prod->GetYaxis()->SetTitleOffset(1.2);
    h_tot_prod->SetFillColor(kBlue);
    h_tot_prod->SetLineColor(kBlack);
    h_tot_prod->SetLineWidth(2);
    h_tot_prod->SetStats(kTRUE);
    h_tot_prod->Draw("HIST");
    gPad->Modified();

    cv_tot_info->cd();
    cv_tot_info->Modified();
    cv_tot_info->SetSelected(cv_tot_info);
    cv_tot_info->SaveAs(Form("../code/sipm/Run_%d/Run_%d_tot_info.png", run_number, run_number));
    delete tot_global_corr;
    delete h_tot_prod;
    delete cv_tot_info;

    TCanvas *cv_tot_corr = new TCanvas("cv_tot_corr", "Time over Threshold correlation between channels", 0, 0, 1920, 1000);
    cv_tot_corr->Divide(3, 4, 0.0001, 0.002);
    for (int i = 0; i != 12; ++i) {
        cv_tot_corr->cd(i + 1);
        tot_corr_ch[i]->SetMarkerStyle(8);
        tot_corr_ch[i]->SetMarkerSize(0.5);
        gPad->SetBottomMargin(0.15);
        // gPad->SetLeftMargin(0.15);
        tot_corr_ch[i]->Draw("AP");
    }
    cv_tot_corr->Update();
    cv_tot_corr->SaveAs(Form("../code/sipm/Run_%d/Run_%d_tot_correlation.png", run_number, run_number));
    for (auto h : tot_corr_ch)
        delete h;
    delete cv_tot_corr;

    TCanvas *cv_chan_coinc = new TCanvas("cv_chan_coinc", "Channel coincidences xy hit map", 0, 0, 1000, 1000);
    hxy_channels->Draw("colz text");
    // gPad->SetLogz();
    draw_square(cv_chan_coinc, "Channels");
    cv_chan_coinc->Update();
    cv_chan_coinc->SaveAs(Form("../code/sipm/Run_%d/Run_%d_hxy_channels.png", run_number, run_number));
    delete hxy_channels;
    delete cv_chan_coinc;

    TCanvas *cv_bar_coinc = new TCanvas("cv_bar_coinc", "Scintillator bar coincidences xy hit map", 0, 0, 1000, 1000);
    hxy_bars->Draw("colz text");
    // gPad->SetLogz();
    draw_square(cv_bar_coinc, "Bars");
    cv_bar_coinc->Update();
    cv_bar_coinc->SaveAs(Form("../code/sipm/Run_%d/Run_%d_hxy_bars.png", run_number, run_number));
    delete hxy_bars;
    delete cv_bar_coinc;

    TCanvas *cv_bar_coinc_toa = new TCanvas("cv_bar_coinc_toa", "Scintillator bar coincidences xy hit map using ToA", 0, 0, 1000, 1000);
    hxy_bars_toa->Draw("colz text");
    // gPad->SetLogz();
    draw_square(cv_bar_coinc_toa, "Bars");
    cv_bar_coinc_toa->Update();
    cv_bar_coinc_toa->SaveAs(Form("../code/sipm/Run_%d/Run_%d_hxy_bars_toa.png", run_number, run_number));
    delete hxy_bars_toa;
    delete cv_bar_coinc_toa;

    TCanvas *cv_bar_coinc_tot = new TCanvas("cv_bar_coinc_tot", "Scintillator bar coincidences xy hit map using ToT", 0, 0, 1000, 1000);
    hxy_bars_tot->Draw("colz text");
    // gPad->SetLogz();
    draw_square(cv_bar_coinc_tot, "Bars");
    cv_bar_coinc_tot->Update();
    // cv_bar_coinc_tot->SaveAs(Form("Run_%d/Run_%d_hxy_bars_tot.png", run_number, run_number));
    delete hxy_bars_tot;
    delete cv_bar_coinc_tot;

    TCanvas *cv_bar_cross_coinc = new TCanvas("cv_bar_cross_coinc", "Scintillator bar cross coincidences xy hit map", 0, 0, 1000, 1000);
    hxy_cross_bars->Draw("colz text");
    // gPad->SetLogz();
    draw_square(cv_bar_cross_coinc, "Bars");
    cv_bar_cross_coinc->Update();
    cv_bar_cross_coinc->SaveAs(Form("../code/sipm/Run_%d/Run_%d_hxy_cross_bars.png", run_number, run_number));
    delete hxy_cross_bars;
    delete cv_bar_cross_coinc;
}

void draw_square(TCanvas *canvas, const std::string &plot_type) {
    float half_side_length;
    int multiplier;
    if (plot_type == "Bars") {
        half_side_length = bar_width;
        multiplier = 3;
    } else if (plot_type == "Channels") {
        half_side_length = si_channel_length;
        multiplier = 4;
    } else {
        std::cerr << "Invalid plot type: " << plot_type << " - No square will be drawn\n";
        return;
    }
    half_side_length *= multiplier;

    TLine *line1 = new TLine(-half_side_length, -half_side_length, -half_side_length, +half_side_length);
    line1->SetLineColor(kRed);
    line1->SetLineWidth(4);
    line1->Draw();

    TLine *line2 = new TLine(-half_side_length, +half_side_length, +half_side_length, +half_side_length);
    line2->SetLineColor(kRed);
    line2->SetLineWidth(4);
    line2->Draw();

    TLine *line3 = new TLine(+half_side_length, +half_side_length, +half_side_length, -half_side_length);
    line3->SetLineColor(kRed);
    line3->SetLineWidth(4);
    line3->Draw();

    TLine *line4 = new TLine(+half_side_length, -half_side_length, -half_side_length, -half_side_length);
    line4->SetLineColor(kRed);
    line4->SetLineWidth(4);
    line4->Draw();

    canvas->Update();
}

int check_couple_or_adjacent(const int chj, const int ch) {
    if ((chj != ch) && (chj / 4 == ch / 4)) {
        // check if the channels are exactly one opposite to the other end of a
        // scintibar, e.g. 0 with 1
        return 0;
    } else if ((chj / 16 == ch / 16) && ((std::abs(chj - ch) == 3) || (std::abs(chj - ch) == 5))) {
        // check if the channels see, at least in part, the same scintibar at
        // opposite ends, e.g. 0 with 5 or 1 with 4 BUT NOT 13 with 16 or 12 with 17
        return 1;
    } else {
        return -1;
    }
}

int chan_to_bar_index(const int chj, const int ch) {
    const int share_scintibar = check_couple_or_adjacent(chj, ch);
    int bar_index = -1;       // default value, channels don't see any common scintillator bar
    if (share_scintibar == 0) // channels are one opposite to the other
    {
        bar_index = (((ch / 4) + 1) / 2) + (ch / 16);
        // std::cout << "Channel " << chj << " and " << ch << " are coupled to bar
        // no. " << bar_index << " and " << bar_index << "\n";
    } else if (share_scintibar == 1) // channels partially see the same scintillator bar
    {
        bar_index = ((ch + chj - 5) / 8) - (ch / 16);
        // std::cout << "Channel " << chj << " and " << ch << " are coupled to bar
        // no. " << bar_index << " and " << bar_index << "\n";
    }
    return bar_index;
}

void fill_xy_opposite_channels_coinc(const int ch, TH2F &xy_map) {
    if ((ch > 61) || (ch / 4 >= 16)) {
        throw std::out_of_range("Channel is invalid (valid range: 0-61)");
    }
    if constexpr (!std::is_same_v<std::remove_reference_t<decltype(xy_map)>, TH2F>) {
        throw std::invalid_argument("Invalid histogram type. Expected TH2F.");
    }
    if (((ch / 4) >= 0) && ((ch / 4) < 8)) { // Fill vertical channels, assuming they have indices 0-7
        float left_bin_center = -4.0 * si_channel_length + 0.5 * si_channel_length;
        float x = left_bin_center + (ch / 4) * si_channel_length;
        float vertical_edge_bin_center = xy_map.GetYaxis()->GetBinCenter(xy_map.GetNbinsY());
        for (float y = -vertical_edge_bin_center; y <= +vertical_edge_bin_center; y += si_channel_length) {
            xy_map.Fill(x, y);
        }
    } else if (((ch / 4) >= 8) && ((ch / 4) < 16)) { // Fill horizontal channels, assuming they have
                                                     // indices 8-15
        float top_bin_center = 4.0 * si_channel_length - 0.5 * si_channel_length;
        float y = top_bin_center - ((ch / 4) - 8) * si_channel_length;
        float horizontal_edge_bin_center = xy_map.GetXaxis()->GetBinCenter(xy_map.GetNbinsX());
        for (float x = -horizontal_edge_bin_center; x <= +horizontal_edge_bin_center; x += si_channel_length) {
            xy_map.Fill(x, y);
        }
    }
}

void fill_xy_scintibar(const int bar_index, TH2F &xy_map) {
    if (bar_index < 0 || bar_index > 11) {
        throw std::out_of_range("Scintillator bar index is out of range (valid range: 0-11)");
    }
    if constexpr (!std::is_same_v<std::remove_reference_t<decltype(xy_map)>, TH2F>) {
        throw std::invalid_argument("Invalid histogram type. Expected TH2F.");
    }
    if ((bar_index >= 0 && bar_index < 6)) { // Fill vertical bars, assuming they have indices 0-5
        float left_bin_center = -3.0 * bar_width + 0.5 * bar_width;
        float x = left_bin_center + bar_index * bar_width;
        float vertical_edge_bin_center = xy_map.GetYaxis()->GetBinCenter(xy_map.GetNbinsY());
        for (float y = -vertical_edge_bin_center; y <= +vertical_edge_bin_center; y += bar_width) {
            xy_map.Fill(x, y);
        }
    } else if ((bar_index >= 6 && bar_index < 12)) { // Fill horizontal bars, assuming they have
                                                     // indices 6-11
        float top_bin_center = 3.0 * bar_width - 0.5 * bar_width;
        float y = top_bin_center - (bar_index - 6) * bar_width;
        float horizontal_edge_bin_center = xy_map.GetXaxis()->GetBinCenter(xy_map.GetNbinsX());
        for (float x = -horizontal_edge_bin_center; x <= +horizontal_edge_bin_center; x += bar_width) {
            xy_map.Fill(x, y);
        }
    }
}

int toa_coordinate(const int ch_i, const float t_i, const int n_bar_portions) {
    if ((ch_i > 61) || (ch_i / 4 >= 16)) {
        throw std::out_of_range("Channel is invalid (valid range: 0-61)");
    }
    if (t_i < 0 || t_i > (n_bar_portions * 0.5)) {
        throw std::out_of_range("Time given to function is invalid (valid range: 0 - selected time window), " + std::to_string(t_i) + " was given\n");
    }
    const float true_t = t_i / ((float)(n_bar_portions) / (float)(true_n_bar_portions));
    int coord = -1;
    const int furthest_coord = (int)(bar_time_width / 0.5) - 1;
    const float time_offset = (true_t <= (bar_time_width / 2.)) ? +0.1 : -0.1;
    // channel is either at the top or to the left of a bar if ch_i is odd, at the bottom or to the right if even
    coord = (ch_i % 2 == 1) ? (int)((true_t + time_offset) * 2) : furthest_coord - (int)((true_t + time_offset) * 2);
    if (coord < 0 || coord >= (int)(bar_time_width / 0.5)) {
        throw std::runtime_error("Something went wrong with coordinate assignment using time of arrival\n");
    }
    return coord;
}

void fill_xy_scintibar_toa(const int bar_index, const int ch_i, const float t_i, const int n_bar_portions, TH2F &xy_map) {
    const int toa_coord = toa_coordinate(ch_i, t_i, n_bar_portions);
    if (bar_index < 0 || bar_index > 11) {
        throw std::out_of_range("Scintillator bar index is out of range (valid range: 0-11)");
    }
    if constexpr (!std::is_same_v<std::remove_reference_t<decltype(xy_map)>, TH2F>) {
        throw std::invalid_argument("Invalid histogram type. Expected TH2F.");
    }
    const int nbinsx = xy_map.GetNbinsX();
    const int nbinsy = xy_map.GetNbinsY();
    if ((bar_index >= 0 && bar_index < 6)) { // Fill vertical bars, assuming they have indices 0-5
        float left_bin_center = -3.0 * bar_width + 0.5 * bar_width;
        float x = left_bin_center + bar_index * bar_width;
        float vertical_edge_bin_center = xy_map.GetYaxis()->GetBinCenter(nbinsy) - toa_coord * bar_width * (nbinsy / true_n_bar_portions);
        for (float y = vertical_edge_bin_center - bar_width * ((nbinsy / true_n_bar_portions) - 1); y <= +vertical_edge_bin_center; y += bar_width) {
            xy_map.Fill(x, y);
        }
    } else if ((bar_index >= 6 && bar_index < 12)) { // Fill horizontal bars, assuming they have indices 6-11
        float top_bin_center = 3.0 * bar_width - 0.5 * bar_width;
        float y = top_bin_center - (bar_index - 6) * bar_width;
        float horizontal_edge_bin_center = xy_map.GetXaxis()->GetBinCenter(nbinsx) - toa_coord * bar_width * (nbinsx / true_n_bar_portions);
        for (float x = horizontal_edge_bin_center - bar_width * ((nbinsx / true_n_bar_portions) - 1); x <= +horizontal_edge_bin_center; x += bar_width) {
            xy_map.Fill(x, y);
        }
    }
}

void fill_xy_cross_scintibar(const int bar_index1, const int bar_index2, TH2F &xy_map) {
    if (bar_index1 < 0 || bar_index1 > 11 || bar_index2 < 0 || bar_index2 > 11) {
        throw std::out_of_range("Scintillator bar index is out of range (valid range: 0-11)");
    }
    if constexpr (!std::is_same_v<std::remove_reference_t<decltype(xy_map)>, TH2F>) {
        throw std::invalid_argument("Invalid histogram type. Expected TH2F.");
    }
    float left_bin_center = -3.0 * bar_width + 0.5 * bar_width;
    float top_bin_center = 3.0 * bar_width - 0.5 * bar_width;
    if ((bar_index1 >= 0 && bar_index1 < 6) && (bar_index2 >= 6 && bar_index2 < 12)) {
        float x = left_bin_center + bar_index1 * bar_width;
        float y = top_bin_center - (bar_index2 - 6) * bar_width;
        xy_map.Fill(x, y);
    } else if ((bar_index2 >= 0 && bar_index2 < 6) && (bar_index1 >= 6 && bar_index1 < 12)) {
        float x = left_bin_center + bar_index2 * bar_width;
        float y = top_bin_center - (bar_index1 - 6) * bar_width;
        xy_map.Fill(x, y);
    } else {
        std::cerr << "Something went wrong: maybe bars are oriented in the same way?\n";
    }
}

summaryList *channelTimeProcess(const int run_number, TTree *tl, coincidence_Info &coinc_info, const float thr_toa, const float nsigma, const std::vector<std::vector<float>> &vped) {

    const int n_bar_portions = (thr_toa) / 0.5; // number of portions in which the scintillator bar can be divided, according to the given time window
    const float toa_space_resolution = (half_bar_length * 2.) / n_bar_portions;
    std::cout << "\nWith 0.5 ns resolution on timing measurements, a "
                 "scintillator bar can be divided in "
              << n_bar_portions << " sections -> x \"resolution\" using time of arrival information is " << toa_space_resolution << " mm\n";
    summaryList *rpl = new summaryList("channelTimeProcess");

    tl->Draw("timeus/1e6", "", "goff"); // s
    double acqtime = TMath::MaxElement(tl->GetSelectedRows(), tl->GetV1()) - TMath::MinElement(tl->GetSelectedRows(), tl->GetV1());
    std::cout << "Total Acquisition time [s]: " << acqtime << '\n';

    // find unique indeces of acquired channels
    std::vector<int> act_ch_v = getActiveChannels(tl);
    const int nchs = act_ch_v.size();
    const int n_chs_ped = vped[0].size();
    if (nchs != n_chs_ped) {
        std::cerr << "Error: getActiveChannels detects " << nchs << " active channels while pedestal vector has " << n_chs_ped << " !\n";
    }
    std::cout << "Number of active SiPM channels from data: " << nchs << '\n';

    //  getToTPedestal(tl, act_ch_v, thr_tot);

    int ntref = tl->GetEntries();
    printf("Total Reference Times (triggers) : %d\n", ntref);

    // find time coincidences
    TCanvas *cv = new TCanvas("cv");
    cv->Divide(2, 2);
    cv->Update();

    TH2F *hxy = new TH2F("hxy", "x/y coinc Hit map", 20, -60, 60, 20, -60, 60);
    hxy->SetXTitle("x centroid [mm]");
    hxy->SetYTitle("y centroid [mm]");

    TH1F *hx = new TH1F("hx", "x Hit Map", 20, -60, 60);
    hx->SetXTitle("x centroid [mm]");
    hx->SetLineColor(1);
    TH1F *hix = new TH1F("hix", "Hit Map", 8, -0.5, 7.5);
    hix->SetXTitle("axis position index");
    hix->SetLineColor(1);

    TH1F *hy = new TH1F("hy", "y Hit Map", 20, -60, 60);
    hy->SetXTitle("y centroid [mm]");
    hy->SetLineColor(2);
    TH1F *hiy = new TH1F("hiy", "Hit map", 8, -0.5, 7.5);
    hiy->SetXTitle("axis position index");
    hiy->SetLineColor(2);

    TString scut = "";
    TString sand = "";
    std::cout << "\nPedestal: active channels, mean ToT, RMS ToT and entries: \n";
    std::cout << std::setw(10) << std::left << "Channel" << std::setw(10) << std::left << "Mean ToT" << std::setw(10) << std::left << "RMS ToT" << std::setw(10) << std::left << "Entries" << "\n";
    for (uint i = 0; i < vped[0].size(); ++i) {       // can be improved including the "and" between opposite scinti-bar
                                                      // ends ToT channels
        int ch = static_cast<int>(vped[0][i]);        // one of the 31 active SiPM channels
        float mean_ToT = vped[1][i];                  // pedestal mean
        float rms_ToT = vped[2][i];                   // pedestal rms
        int entries = static_cast<int>(vped[3][i]);   // pedestal entries
        float thr = vped[1][i] + nsigma * vped[2][i]; // threshold is pedestal mean + nsigma
                                                      // (around 3) * pedestal rms
        scut += Form("%s(tot>=%.2f && chan == %d)", sand.Data(), thr,
                     ch); // creating the string for TCut: channel i must have tot >
                          // threshold
        sand = "||";
        std::printf("%-10d%-10.2f%-10.2f%-10d\n", ch, mean_ToT, rms_ToT, entries);
    }
    std::cout << "\nOut of " << n_chs_ped << " channels, " << std::count(vped[1].begin(), vped[1].end(), 0.) << " have not enough entries to compute pedestals\n";

    std::cout << '\n';
    TCut pedcut = scut.Data();
    pedcut.Print(); // print TCut string

    TCut cbepu = "bPulse<6e12"; // low beam pulse cut
    tl->Draw("tot", "", "goff");
    int nsig_tot = tl->GetSelectedRows();
    tl->Draw("tot", pedcut, "goff");
    int nsig_over_thr = tl->GetSelectedRows();
    std::cout << "N. of ToT signals over threshold: " << nsig_over_thr << " out of " << nsig_tot << " total signals ( " << 100 * nsig_over_thr / nsig_tot << "% ) \n";

    TH1F *h_time_diff = new TH1F("h_time_diff", "Time Differences [ns]", 60, 0, 30);
    h_time_diff->GetXaxis()->SetTitle("Time differences [ns]");
    int tottimediff = 0;
    int negtimediff = 0;

    const int trg_rbin = run_delays.at(run_number) != 0 ? 5 : 25;
    TH1F *h_coinc_per_trigger = new TH1F("h_coinc_per_trigger", "Coincidences per trigger", trg_rbin, 0, trg_rbin);
    h_coinc_per_trigger->GetXaxis()->SetTitle("# of coincidences per trigger");

    std::vector<TH1F *> dead_t_ch_even(16); // For estimation of dead time in right and bottom channels
    for (int i = 0; i != 16; ++i) {
        dead_t_ch_even[i] = new TH1F(Form("dead_t_even_ch_%d", i * 4), Form("Channel no. %d", i * 4), 100, 0, 500);
        dead_t_ch_even[i]->GetXaxis()->SetTitle("Time diff [ns]");
    }

    std::vector<TH1F *> dead_t_ch_odd(16); // For estimation of dead time in right and bottom channels
    for (int i = 0; i != 16; ++i) {
        dead_t_ch_odd[i] = new TH1F(Form("dead_t_odd_ch_%d", i * 4 + 1), Form("Channel no. %d", i * 4 + 1), 100, 0, 500);
        dead_t_ch_odd[i]->GetXaxis()->SetTitle("Time diff [ns]");
    }

    TH2F *h_toa_corr = new TH2F("h_toa_corr", "Time of Arrival correlation between coincidences", n_bar_portions * 1.5, 0., thr_toa / 1.9, n_bar_portions * 1.5, thr_toa / 2.2, thr_toa * 1.1);
    h_toa_corr->SetXTitle("t_{i} [ns]");
    h_toa_corr->SetYTitle("t_{j} [ns]");

    TH1F *h_toa = new TH1F("h_toa", "Coincidences hit time distribution ", 36, 0, 36000);
    h_toa->GetXaxis()->SetTitle("Hit time [ns]");
    h_toa->GetYaxis()->SetTitle("Occurrences");

    TGraph *tot_global_corr = new TGraph();
    tot_global_corr->SetName("tot_global_corr");
    tot_global_corr->SetTitle("Correlation of Time over Threshold between coincidences channels");
    tot_global_corr->GetXaxis()->SetTitle("tot_{i}");
    tot_global_corr->GetYaxis()->SetTitle("tot_{j}");

    TH1F *h_tot_prod = new TH1F("h_tot_prod", "Distribution of tot_i * tot_j of a bar coincidence", 50, 0, 4000);
    h_tot_prod->GetXaxis()->SetTitle("tot_{i} * tot_{j} [ns^{2}]");
    h_tot_prod->GetYaxis()->SetTitle("Occurrences");

    std::vector<TGraph *> tot_corr_ch(12); // Correlation between ToT of coincidences in a bar
    for (int i = 0; i != 12; ++i) {
        tot_corr_ch[i] = new TGraph();
        tot_corr_ch[i]->SetName(Form("tot_corr_couple_%d", i));
        tot_corr_ch[i]->SetTitle(Form("Bar no. %d", i));
        tot_corr_ch[i]->GetXaxis()->SetTitle(Form("tot_{%s}", i < 6 ? "bottom" : "right"));
        tot_corr_ch[i]->GetYaxis()->SetTitle(Form("tot_{%s}", i < 6 ? "up" : "left"));
    }

    TH2F *hxy_channels = new TH2F("hxy_channels", "x/y coinc hit map for channels (only same bar ends)", 14, -si_channel_length * 7, si_channel_length * 7, 14, -si_channel_length * 7, si_channel_length * 7);
    hxy_channels->SetXTitle("x [mm]");
    hxy_channels->SetYTitle("y [mm]");

    TH2F *hxy_bars = new TH2F("hxy_bars", "x/y Coincidences hit map for scintillation bars", 10, -bar_width * 5, bar_width * 5, 10, -bar_width * 5, bar_width * 5);
    hxy_bars->SetXTitle("x [mm]");
    hxy_bars->SetYTitle("y [mm]");
    int bar_coinc_counter = 0;

    TH2F *hxy_bars_toa = new TH2F("hxy_bars_toa", "x/y Coincidences hit map for scintillation bars, using ToA", 30, -bar_width * 15, bar_width * 15, 30, -bar_width * 15, bar_width * 15);
    hxy_bars_toa->SetXTitle("x [mm]");
    hxy_bars_toa->SetYTitle("y [mm]");

    TH2F *hxy_bars_tot = new TH2F("hxy_bars_tot", "x/y Coincidences hit map for scintillation bars, using ToT", 10, -bar_width * 5, bar_width * 5, 10, -bar_width * 5, bar_width * 5);
    hxy_bars_tot->SetXTitle("x [mm]");
    hxy_bars_tot->SetYTitle("y [mm]");

    TH2F *hxy_cross_bars = new TH2F("hxy_cross_bars", "x/y Coincidences hit map for crossed scintillation bars", 10, -bar_width * 5, bar_width * 5, 10, -bar_width * 5, bar_width * 5);
    hxy_cross_bars->SetXTitle("x [mm]");
    hxy_cross_bars->SetYTitle("y [mm]");
    int bar_2coinc_counter = 0;
    int bar_cross_coinc_counter = 0;

    for (int trigger_idx = 0; trigger_idx < ntref; ++trigger_idx) { // loop on reference times (triggers)
        float timestamp;
        tl->SetBranchAddress("timeus", &timestamp);
        tl->GetEntry(trigger_idx);
        TCut ecut = Form("(Entry$==%d)", trigger_idx);

        tl->Draw("toa:tot:chan", ecut && pedcut,
                 "goff");                 // goff is the option to NOT draw anything
        int nsig = tl->GetSelectedRows(); // 0 when no event meets the required constraints
        if (nsig <= 0) {
            continue;
        }

        std::span<const Double_t> toa_v(tl->GetV1(), nsig);      // Vector to store time of arrivals for current trigger event
        std::span<const Double_t> tot_v(tl->GetV2(), nsig);      // Vector to store time over threshold for current trigger event
        std::span<const Double_t> channels_v(tl->GetV3(), nsig); // Vector to store channels for current trigger event

        if (toa_v.size() != channels_v.size() || tot_v.size() != channels_v.size()) {
            std::cerr << "\n\n ATTENTION: Size mismatch between ToA/TOT and channels "
                         "vectors \n\n";
            break;
        }

        int coinc_per_trigger = 0;
        // std::vector<float> time_coincidences;                 // Vector to store time coincidences (i.e. < thr_toa) for each trigger event
        // std::vector<std::vector<float>> ch_coincidences(16);  // Vector to store time coincidences (i.e. < thr_toa) for each trigger event for different channels
        // std::vector<std::vector<float>> bar_coincidences(12); // Vector to store time coincidences (i.e. < thr_toa) for each trigger event for different bars
        std::vector<float> coinc_hit_toa; // Vector to store the average time of arrival of a bar coincidence for each trigger event
        std::vector<int> coinc_bar_index; // Vector to store the scintillator bar index of when a coincidence happens

        for (size_t i = 0; i != toa_v.size(); ++i) {
            int ch_i = channels_v[i];
            float toa_i = toa_v[i] * 0.5; // *0.5 to have it in [ns]
            float tot_i = tot_v[i] * 0.5; // *0.5 to have it in [ns]

            // Iterate over other channels to find adjacent channels and select events within the time window
            for (size_t j = i + 1; j != toa_v.size(); ++j) {
                // OPEN QUESTION: appena trovo una coincidenza escludo la coppia o cerco
                // se nella stessa finestra un canale vicino è stato colpito? Come
                // gestisco poi la cosa?
                int ch_j = channels_v[j];
                float toa_j = toa_v[j] * 0.5; // *0.5 to have it in [ns]
                if ((ch_i == ch_j) && (toa_i != toa_j)) {

                    float time_difference = (toa_j - toa_i);
                    ch_i % 2 == 0 ? dead_t_ch_even[ch_i / 4]->Fill(time_difference) : dead_t_ch_odd[ch_i / 4]->Fill(time_difference);
                }
                const int share_scintibar = check_couple_or_adjacent(ch_j, ch_i);
                if (share_scintibar != -1) // Check if channels are a couple
                {
                    float time_difference = (toa_j - toa_i);
                    ++tottimediff;
                    if (time_difference >= 0) {
                        h_time_diff->Fill(time_difference);
                        if (time_difference < thr_toa) {
                            // coincidence!
                            ++bar_coinc_counter;
                            ++coinc_per_trigger;
                            // time_coincidences.push_back(time_difference);
                            const int bar_index = chan_to_bar_index(ch_j, ch_i);
                            fill_xy_scintibar(bar_index, *hxy_bars);
                            // bar_coincidences[bar_index].push_back(time_difference);
                            if (share_scintibar == 0) { // channels are at exactly opposite ends of a scintillator bar -> ch / 4 == chj / 4
                                fill_xy_opposite_channels_coinc(ch_j, *hxy_channels);
                                // ch_coincidences[ch_j / 4].push_back(time_difference);
                            }
                            float tot_j = tot_v[j] * 0.5; // *0.5 to have it in [ns]

                            // compute coincidence between bars

                            // OPEN QUESTION: richiedo che le scintibar siano perpendicolari?
                            // Oppure così perdo eventi nella giunzione delle scintibar? In
                            // generale le scintibar si parlano? Cioè la luce da una barra
                            // all'altra si propaga? OPEN QUESTION: se proprio succede che c'è
                            // più di una coincidenza tra barre, prendo quella minima o le
                            // prendo tutte?

                            /*
                            Estimate of true hit time t_hit (called hit_toa in the code):
                            toa_i = t_hit + t_i
                            toa_j = t_hit + t_j
                            t_i + t_j = bar_length / c
                            ===> toa_i + toa_j = bar_length / c + 2 * t_hit
                            only thing is that bar_length / c could be higher than real value, because of the chosen time window
                            (which in principle should be bar_time_width but could be higher to account for smearings and delays)
                            */
                            float hit_toa = ((toa_j + toa_i) - (thr_toa)) / 2.; // [ns] estimate of true hit time

                            h_toa_corr->Fill((toa_i - hit_toa), (toa_j - hit_toa));
                            h_toa->Fill(hit_toa);
                            tot_global_corr->AddPoint(tot_i, tot_j);
                            h_tot_prod->Fill(tot_i * tot_j);
                            tot_corr_ch[bar_index]->AddPoint((ch_i % 2 == 0) ? tot_i : tot_j, (ch_i % 2 == 0) ? tot_j : tot_i);

                            for (size_t k = 0; k < coinc_hit_toa.size(); ++k) {
                                const float hit_toa_k = coinc_hit_toa[k];
                                if (std::fabs(hit_toa_k - hit_toa) < (thr_toa)) {
                                    // 2 bars coincidence!
                                    ++bar_2coinc_counter;
                                    const int bar_index_k = coinc_bar_index[k];
                                    if ((bar_index_k / 6) != (bar_index / 6)) {
                                        // cross coincidence!
                                        ++bar_cross_coinc_counter;
                                        if (run_delays.at(run_number) != 0) {
                                            std::cout << "\nCross coincidence found at trigger " << trigger_idx << " (timestamp " << timestamp << ") : scintibars " << bar_index << " and " << bar_index_k << " fired at " << hit_toa << " and " << hit_toa_k << " respectively\n\n";
                                        }
                                        fill_xy_cross_scintibar(bar_index_k, bar_index, *hxy_cross_bars);
                                    }
                                }
                            }
                            coinc_hit_toa.push_back(hit_toa);
                            coinc_bar_index.push_back(bar_index);

                            // trying to do xy mapping with tot (epic fail)

                            // trying to do xy mapping with toa
                            fill_xy_scintibar_toa(bar_index, ch_i, (toa_i - hit_toa), n_bar_portions, *hxy_bars_toa);
                        } else {
                            ++negtimediff;
                        }
                    }
                }
            }
        }

        h_coinc_per_trigger->Fill(coinc_per_trigger);

        // cisbani_processing(tl, thr_toa, ntref, trigger_idx, nsig, cv, hx, hy, hix, hiy, hxy);
    }

    plot_timing_data(h_time_diff, h_coinc_per_trigger, dead_t_ch_even, dead_t_ch_odd, h_toa_corr, h_toa, tot_global_corr, h_tot_prod, tot_corr_ch, hxy_channels, hxy_bars, hxy_bars_toa, hxy_bars_tot, hxy_cross_bars, thr_toa, run_number);

    std::cout << "\nTotal anomalous negative time differences of this run: " << negtimediff << " out of " << tottimediff << " ( " << 100. * (float)negtimediff / (float)tottimediff << "% )\n";
    std::cout << std::setw(30) << std::left << "Total coincidences" << std::setw(30) << std::left << "Total double coincidences" << std::setw(30) << std::left << "Total cross coincidences" << "\n";
    std::cout << std::setw(30) << std::left << bar_coinc_counter << std::setw(30) << std::left << bar_2coinc_counter << std::setw(30) << std::left << bar_cross_coinc_counter << "\n";
    //
    //
    //
    coinc_counter counter(run_number, run_delays.at(run_number), ntref, bar_coinc_counter, bar_2coinc_counter, bar_cross_coinc_counter);
    coinc_info.add_counter(counter);
    // estimate rates
    float scale = 1. / ((float)ntref);

    float totxy = hxy->Integral();
    float totx = hx->Integral();
    float toty = hy->Integral();
    printf("Statistics:     xy      x      y\n");
    printf("Sum of signals: %.4f %.4f %.4f\n", totxy, totx, toty);
    printf("Signal/Trigger: %.4f %.4f %.4f\n", totxy * scale, totx * scale, toty * scale);
    printf("Signal/Time[s]: %.4f %.4f %.4f\n", totxy / acqtime, totx / acqtime, toty / acqtime);

    std::vector<float> rval;
    rval.push_back(acqtime);
    rval.push_back((float)ntref);
    rpl->addFloat(acqtime, "AcqTime_s");
    rpl->addInt(ntref, "NEntries");

    for (int i = 0; i < 8; i++) {
        float count = hix->GetBinContent(i + 1);
        rval.push_back(count);
        rpl->addFloat(count, Form("xCoinc[%d]", i));
    }
    for (int i = 0; i < 8; i++) {
        float count = hiy->GetBinContent(i + 1);
        rval.push_back(count);
        rpl->addFloat(count, Form("yCoinc[%d]", i));
    }

    // hxy->Scale(scale);
    hy->Scale(scale);
    hy->SetYTitle("Signals/Triggers");
    hx->Scale(scale);
    hx->SetYTitle("Signals/Triggers");
    cv->cd(1);
    hxy->Draw("colz");
    gPad->SetLogz(1);
    cv->Update();
    cv->cd(2);
    hy->Draw();
    cv->Update();
    cv->cd(3);
    hx->Draw();
    cv->cd(4);
    hiy->Draw();
    hix->Draw("same");
    cv->Update();

    return rpl;
}

summaryList *plotTimeData(TChain *chainT) {

    // process and display data stored in TTree

    TCut cbepu = "bPulse<6e12"; // low beam pulse cut

    // list of "summary" parameters

    summaryList *spl = new summaryList("plotTimeData");

    //  TParameter<float> *pl[20];
    // for (int i=0;i<20;i++) { pl[i] = NULL; }

    chainT->LoadTree(0);
    chainT->GetTree()->GetUserInfo()->Print();
    //  int run = daqp->runNumber();
    TList *userlist = (TList *)chainT->GetTree()->GetUserInfo();
    TParameter<int> *par = (TParameter<int> *)userlist->At(0); // first par is run number
    int run = par->GetVal();
    printf("Plot data from run %d\n", run);
    float timebw = ((TParameter<float> *)userlist->FindObject("TimeBinWidth"))->GetVal();
    std::time_t start_time = ((TParameter<std::time_t> *)userlist->FindObject("StartTime"))->GetVal();

    std::cout << "START TIME : " << start_time << std::endl;

    TCanvas *cc1 = new TCanvas("cc1", Form("Timing mode correlation (run %d)", run));
    cc1->Divide(1, 1);
    cc1->Update();
    cc1->cd(1);
    chainT->Draw(Form("tot*%f:toa*%f", timebw, timebw), "", "colz");
    TH2F *hdummy2 = (TH2F *)gPad->GetPrimitive("htemp");
    hdummy2->SetTitle(Form("Timing / Amplitude correlation (run %d)", run));
    hdummy2->SetXTitle("Time of Arrival [ns]");
    hdummy2->SetYTitle("Time over Threshold [ns]");
    cc1->Update();

    TCanvas *cc2 = new TCanvas("cc2", "Ref Trigger Time Distribution");
    //  cc2->Divide(2,2);
    //  cc2->Update();
    TPad *pad1 = new TPad("pad1", "Pad1", 0.05, 0.51, 0.98, 0.97);
    TPad *pad2 = new TPad("pad2", "Pad2", 0.05, 0.02, 0.98, 0.49);
    pad1->Draw();
    pad2->Draw();
    pad2->cd();
    TPad *pad21 = new TPad("pad21", "Pad21", 0.02, 0.05, 0.48, 0.98);
    TPad *pad22 = new TPad("pad22", "Pad22", 0.52, 0.05, 0.98, 0.98);
    pad21->Draw();
    pad22->Draw();

    // cc2->cd(1);
    pad1->cd();
    chainT->Draw("novt:timeus/1000.", "", "prof");
    int nn = chainT->GetSelectedRows();

    float smean = -1;
    float srms = -1;
    if (nn > 0) {
        //    pl[0] = new
        //    TParameter<float>("Signal_mean",TMath::Mean(nn,chainT->GetV1())); //
        //    signal over thr time mean ~ charge mean pl[1] = new
        //    TParameter<float>("Signal_rms", TMath::StdDev(nn,chainT->GetV1())); //
        //    signal over thr. time rms ~ charge rms
        smean = TMath::Mean(nn, chainT->GetV1());
        srms = TMath::StdDev(nn, chainT->GetV1());
    }
    spl->addFloat(smean, "Signal_mean");
    spl->addFloat(srms, "Signal_rms");

    float tstart = TMath::MinElement(nn, chainT->GetV2());
    float tstop = TMath::MaxElement(nn, chainT->GetV2());
    float acqtime = (tstop - tstart) / 1000.; // [s]
    float scale = 1. / acqtime;

    printf(" number of loaded events: %d  acqtime %f [ms] -> event rate %f [Hz]\n", nn, acqtime, ((float)nn) / acqtime);

    hdummy2 = (TH2F *)gPad->GetPrimitive("htemp");
    hdummy2->SetTitle(Form("OverThreshold Signals vs Trigger Time since start of run %d", run));
    hdummy2->SetXTitle("Trigger Time [ms]");
    hdummy2->SetYTitle("Number OverThr");

    float min = 1e99;
    float max = -1e99;

    for (int i = 1; i < nn; i++) {
        float delta = chainT->GetV2()[i] - chainT->GetV2()[i - 1];
        min = (delta < min) ? delta : min;
        max = (delta > max) ? delta : max;
    }
    float dmm = (max - min) / 99.;

    printf(" delta time min/max %f %f [ms]\n", min, max);

    TH1F *hdt = new TH1F("hdt", "Time Distance between triggers", 200, min - dmm, max + dmm);
    hdt->SetXTitle("Delta Time [ms]");

    for (int i = 1; i < nn; i++) {
        float delta = chainT->GetV2()[i] - chainT->GetV2()[i - 1];
        hdt->Fill(delta);
    }
    //  cc2->cd(2)->SetLogy();
    pad21->cd()->SetLogy();
    hdt->Draw();

    pad22->cd();
    chainT->Draw("novt"); // number of overtreshold channels / event
    TH1F *hdummy = (TH1F *)gPad->GetPrimitive("htemp");
    hdummy->SetTitle(Form("Acquired number of channels/event (run %d)", run));
    hdummy->SetXTitle("NumChannels/Event");
    TF1 *fpe = new TF1("fpe", "pol2(0)+expo(3)", 2, 62); // fit to find minimum
    fpe->SetRange(2, 62);
    hdummy->Fit(fpe, "R");
    float ovt_thr = fpe->GetMinimumX(2, 60);
    // float ovt_low = hdummy->Integral(0, ovt_thr);
    spl->addFloat(ovt_thr, "Count_thr");

    //    pl[2] = new TParameter<float>("Count_Thr",ovt_thr); // per event

    //  pl[3] = new TParameter<float>("Ovt_count_Low", ovt_low);
    //  pl[4] = new TParameter<float>("Ovt_count_High",
    //  hdummy->Integral(0,64)-ovt_low);
    cc2->Update();

    TCanvas *cc3 = new TCanvas("cc3", "Compare to beam pulses");
    TPad *cc3p1 = new TPad("cc3p1", "Pad31", 0.05, 0.51, 0.98, 0.97);
    TPad *cc3p2 = new TPad("cc3p2", "Pad32", 0.05, 0.02, 0.98, 0.49);
    cc3p1->Draw();
    cc3p2->Draw();
    cc3p1->cd();
    TPad *cc3p11 = new TPad("cc3p11", "Pad11", 0.02, 0.05, 0.48, 0.98);
    TPad *cc3p12 = new TPad("cc3p12", "Pad312", 0.52, 0.05, 0.98, 0.98);
    cc3p11->Draw();
    cc3p12->Draw();

    cc3p11->cd();
    chainT->SetMarkerStyle(1);
    chainT->Draw("novt:bPulse");
    chainT->SetMarkerStyle(7);
    chainT->SetMarkerColor(2);
    chainT->Draw("novt:bPulse", cbepu, "same");
    float mean_novt0 = TMath::Mean(chainT->GetSelectedRows(), chainT->GetV1());
    chainT->SetMarkerColor(1);
    chainT->Draw("novt:bPulse", !cbepu, "same");
    float mean_novt1 = TMath::Mean(chainT->GetSelectedRows(), chainT->GetV1());

    spl->addFloat(mean_novt0, "Mean_novt_low_beam");
    spl->addFloat(mean_novt1, "Mean_novt_high_beam");

    cc3p12->cd();
    chainT->Draw("bDTime");
    // chainT->Draw(Form("tot*%f:bPulse",timebw),Form("%f",scale),"colz");
    chainT->Draw(Form("toa*%f:bPulse", timebw), Form("%f", scale), "colz");

    cc3p2->cd();
    chainT->Draw("novt:timeus/1000000");
    chainT->SetMarkerColor(2);
    chainT->Draw("novt:timeus/1000000", "bPulse<6e12", "same");
    chainT->SetMarkerColor(1);
    cc3->Update();

    chainT->Draw("novt", Form("novt<=%f", ovt_thr), "goff");
    int ndu = chainT->GetSelectedRows();
    //  pl[3] = new TParameter<float>("NSignal_Low", ((float) ndu)); // total
    //  signals below threshold
    spl->addInt(ndu, "NSignal_Low");

    smean = srms = -1;
    if (ndu > 0) {
        //    pl[4] = new TParameter<float>("Signal_mean_Low",
        //    TMath::Mean(ndu,chainT->GetV1())); // signal/event pl[5] = new
        //    TParameter<float>("Signal_rms_Low",
        //    TMath::StdDev(ndu,chainT->GetV1()));
        smean = TMath::Mean(ndu, chainT->GetV1());
        srms = TMath::StdDev(ndu, chainT->GetV1());
    }
    spl->addFloat(smean, "Signal_mean_Low");
    spl->addFloat(srms, "Signal_rms_Low");

    chainT->Draw("novt", Form("novt>%f", ovt_thr), "goff");
    ndu = chainT->GetSelectedRows();
    // pl[6] = new TParameter<float>("NSignal_high", ((float) ndu));
    spl->addInt(ndu, "NSignal_high");

    smean = srms = -1;
    if (ndu > 0) {
        //    pl[7] = new TParameter<float>("Signal_mean_High",
        //    TMath::Mean(ndu,chainT->GetV1())); pl[8] = new
        //    TParameter<float>("Signal_rms_High",
        //    TMath::StdDev(ndu,chainT->GetV1()));
        smean = TMath::Mean(ndu, chainT->GetV1());
        srms = TMath::StdDev(ndu, chainT->GetV1());
    }
    spl->addFloat(smean, "Signal_mean_High");
    spl->addFloat(srms, "Signal_rms_High");

    printf(" scale %f\n", scale);
    TCanvas *cc = new TCanvas("cc", Form("Timing mode (run %d)", run));
    cc->Divide(1, 2);
    cc->Update();

    cc->cd(1);
    chainT->Draw(Form("tot*%f", timebw), Form("%f", scale));

    hdummy = (TH1F *)gPad->GetPrimitive("htemp");
    hdummy->SetTitle(Form("Signal Time Lenght - Overthreshold (run %d)", run));
    hdummy->SetXTitle("Time Over Threshold [ns]");

    //  pl[9] = new TParameter<float>("ToT_mean",hdummy->GetMean());
    //  pl[10] = new TParameter<float>("ToT_rms",hdummy->GetRMS());
    spl->addFloat(hdummy->GetMean(), "ToT_mean");
    spl->addFloat(hdummy->GetRMS(), "ToT_rms");

    cc->cd(2)->SetLogy();
    chainT->Draw(Form("toa*%f", timebw), Form("%f", scale));
    float reftwin = TMath::MaxElement(chainT->GetSelectedRows(), chainT->GetV1());
    //  float reft0 = TMath::MinElement(chainT->GetSelectedRows(),
    //  chainT->GetV1());
    printf(" Ref Time Window : %f\n", reftwin);
    hdummy = (TH1F *)gPad->GetPrimitive("htemp");

    // pl[11] = new TParameter<float>("ToA_mean",hdummy->GetMean());
    // pl[12] = new TParameter<float>("ToA_rms",hdummy->GetRMS());
    spl->addFloat(hdummy->GetMean(), "ToA_mean");
    spl->addFloat(hdummy->GetRMS(), "ToA_rms");

    hdummy->SetTitle(Form("Signal Arrival Time run %d", run));
    hdummy->SetXTitle("Time of Arrival [ns]");
    cc->Update();

    //  pl[13] = new TParameter<float>("AcqTime_s",acqtime);
    //  pl[14] = new TParameter<float>("Events",((float) chainT->GetEntries()));
    spl->addFloat(acqtime, "pAcqTime_s");
    spl->addInt(chainT->GetEntries(), "pEvents");
    /*
        TList *rpl = new TList();
        for (int i=0;i<20;i++) { if (pl[i]) { rpl->Add(pl[i]); } }
        rpl->Print();
    */

    return spl;
}

summaryList *processTEvents(const std::string basepath, std::vector<int> srun, std::vector<int> prun, float totsig, float nsigma, float time_window, TString ntoFile) {

    const std::string prefix = "/fers/Run";
    set_local_style();

    summaryList *srList = NULL; // list with run summary parameters

    // printf("##### Check and parse text data if needed\n");
    srun.insert(srun.end(), prun.begin(),
                prun.end()); // for simplicity parse all runs (signal/pedestal) in
                             // the same way

    checkAndConvertText2Root(basepath, srun, prefix, ntoFile);

    std::cout << "Plotting signal data and extracting active channels\n";
    srun.erase(srun.end() - prun.size(),
               srun.end()); // remove pedestal runs from list of signal runs

    for (uint i = 0; i != srun.size(); ++i) {
        std::cout << "Analysing run " << srun[i] << "\n";
    }

    TChain *ctlist;

    if (srun[0] != -1) {
        ctlist = new TChain("tlist");

        for (uint i = 0; i != srun.size(); ++i) {
            const std::string ifile = basepath + prefix + std::to_string(srun[i]) + ".root";
            ctlist->Add(ifile.c_str());
        }

        srList = plotTimeData(ctlist);
    } else {
        ctlist = NULL;
    } // look at the pedestal only

    std::vector<int> ach = getActiveChannels(ctlist); // this requires the signal ttree and not the pedestal ttree

    ctlist->Delete();

    // get pedestals

    std::cout << "\nGetting pedestal . . .\n";
    std::vector<std::vector<float>> ped_v = get_ToT_pedestal(ach, totsig, basepath, prun);

    // now reopen the signal chain
    printf("Processing signal data . . . \n");

    std::vector<float> agd;
    summaryList *spList = NULL;
    if (srun[0] != -1) {
        ctlist = new TChain("tlist");
        for (uint i = 0; i < srun.size(); ++i) {
            const std::string ifile = basepath + prefix + std::to_string(srun[i]) + ".root";
            std::cout << ifile << '\n';
            ctlist->Add(ifile.c_str());
        }
        coincidence_Info coinc_info;
        if (ctlist != NULL) {
            spList = channelTimeProcess(srun[0], ctlist, coinc_info, time_window, nsigma, ped_v);
            std::cout << "\n Run " << srun[0] << " analysis completed\n";
            for (uint i = 0; i < agd.size(); i++) {
                std::cout << agd[i] << '\n';
            }
        }
    }

    // append the two list objects
    if (srList != NULL) {
        spList->merge(srList);
        delete srList;
        /*
          for (int i=0;i<srList->GetSize();i++) {
          TParameter<float> *pp;
          pp = (TParameter<float> *) srList->At(i);
          agd.push_back(pp->GetVal());
          }
        */
    }
    return spList;
}

void timing_analysis(const std::string basepath, const std::vector<int> sig_runs, const std::vector<int> ped_runs, float tot_default_rms, float tot_rms_n, float time_window) {
    set_local_style();
    const std::string prefix = "/fers/Run";
    const std::string ntoFile = "/cube/run_pkup_sall.root";
    checkAndConvertText2Root(basepath, sig_runs, prefix, ntoFile);
    checkAndConvertText2Root(basepath, ped_runs, prefix, ntoFile);

    TChain *runs_tree = new TChain("tlist");
    for (const int run_n : sig_runs) {
        const std::string run_filename = basepath + prefix + std::to_string(run_n) + ".root";
        runs_tree->Add(run_filename.c_str());
    }
    std::vector<int> active_chs = getActiveChannels(runs_tree);
    std::vector<std::vector<float>> ped_v = get_ToT_pedestal(active_chs, tot_default_rms, basepath, ped_runs);

    summaryList *summ_list = NULL;
    coincidence_Info coinc_info;
    for (const int run_n : sig_runs) {
        const std::string run_filename = basepath + prefix + std::to_string(run_n) + ".root";
        std::unique_ptr<TFile> root_file(TFile::Open(run_filename.c_str(), "READ"));
        if (!root_file || root_file->IsZombie()) {
            std::cerr << "Error opening file" << '\n';
            exit(-1);
        }
        auto run_tree = root_file->Get<TTree>("tlist");
        summ_list = channelTimeProcess(run_n, run_tree, coinc_info, time_window, tot_rms_n, ped_v);
        std::cout << "\nRun " << run_n << " analysis completed\n";
        std::cout << summ_list->getNames() << '\n';
        std::cout << summ_list->getValues() << '\n';
    }
    const std::vector<coinc_counter> coincidences = coinc_info.get_counter();

    TH1F *h_coinc = new TH1F("h_coinc", "Coincidences distribution normalized to total triggers", 60, 0, 12000);
    h_coinc->GetXaxis()->SetTitle("Delay from #gamma flash [ns]");
    h_coinc->GetYaxis()->SetTitle("100 * coincidences / triggers");
    h_coinc->SetFillColor(kBlue);
    h_coinc->SetLineColor(kBlack);
    h_coinc->SetLineWidth(2);
    gPad->SetLogy();

    TH1F *h_2coinc = new TH1F("h_2coinc", "Double coincidences distribution normalized to total triggers", 60, 0, 12000);
    h_2coinc->GetXaxis()->SetTitle("Delay from #gamma flash [ns]");
    h_2coinc->GetYaxis()->SetTitle("100 * coincidences / triggers");
    h_2coinc->SetFillColor(kBlue);
    h_2coinc->SetLineColor(kBlack);
    h_2coinc->SetLineWidth(2);
    gPad->SetLogy();

    TH1F *h_xcoinc = new TH1F("h_xcoinc", "Cross coincidences distribution normalized to total triggers", 60, 0, 12000);
    h_xcoinc->GetXaxis()->SetTitle("Delay from #gamma flash [ns]");
    h_xcoinc->GetYaxis()->SetTitle("100 * coincidences / triggers");
    h_xcoinc->SetFillColor(kBlue);
    h_xcoinc->SetLineColor(kBlack);
    h_xcoinc->SetLineWidth(2);
    gPad->SetLogy();

    for (const auto &counter : coincidences) {
        std::cout << "Run " << counter.run_n << ": " << counter.double_coinc_n << " double coincidences and " << counter.cross_coinc_n << " cross coincidences among " << counter.trig_n << " triggers\n";
        h_coinc->Fill(counter.ns_delay, 100 * (float)counter.coinc_n / (float)counter.trig_n);
        h_2coinc->Fill(counter.ns_delay, 100 * (float)counter.double_coinc_n / (float)counter.trig_n);
        h_xcoinc->Fill(counter.ns_delay, 100 * (float)counter.cross_coinc_n / (float)counter.trig_n);
    }

    TCanvas *cv_coinc = new TCanvas("cv_coinc", "Coincidences distribution", 800, 600);
    gPad->SetLogy();
    h_coinc->Draw("histo");
    cv_coinc->SaveAs("../code/sipm/Coincidences_distribution.png");

    TCanvas *cv_2coinc = new TCanvas("cv_2coinc", "Double coincidences distribution", 800, 600);
    gPad->SetLogy();
    h_2coinc->Draw("histo");
    cv_2coinc->SaveAs("../code/sipm/Double_coincidences_distribution.png");

    TCanvas *cv_xcoinc = new TCanvas("cv_xcoinc", "Cross coincidences distribution", 800, 600);
    gPad->SetLogy();
    h_xcoinc->Draw("histo");
    cv_xcoinc->SaveAs("../code/sipm/Cross_coincidences_distribution.png");
}

/* ==================================================
 * analyze multiple runs in timing mode
 * (adapted to 2310 CERN test)
 *
 */
void runTimingMode(const std::string ipath = "../test_231020/data") {

    std::vector<std::vector<int>> run = {{15, 31},
                                         {16},
                                         {17, 27}, // at least one for each delay
                                                   // and other run variables
                                         {18},
                                         {19},
                                         {20},
                                         {21, 26},
                                         {28},
                                         {29},
                                         {30},
                                         {44},
                                         {45},
                                         {46},
                                         {48},
                                         {49, 53},
                                         {55},
                                         {54},
                                         {56},
                                         {57},
                                         {58}};

    std::vector<int> delay = {0,   500,  1000,  1500,  2000, 2500, 3000, 3500, 4000, 10000, // [ns] delay from cube trigger signal
                              500, 2000, 10000, 10000, 2000, 500,  0,    2000, 2000, 2000};

    std::vector<float> tdgap = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 10, 10, 10, 30, 30, 30, 30, 20, 20, 20}; // [cm] target detector distance

    std::vector<int> target = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, // target type 0: light 0.3mm, 1:
                                  // light-vertical, 2:medium-vertical 0.5mm,
                                  // 3:heavy-vertical 0.5mm+reinforced dome
                               0, 0, 0, 0, 0, 0, 0, 1, 2, 3};

    std::vector<int> prun = {13, 22, 33}; // pedestal runs

    std::vector<float> ctime = {50, 100, 200, 350, 500}; // coincidence gate window time

    //  std::vector<std::vector<float>> vdb;
    std::vector<std::string> vdb;
    TString stitle = "";

    uint nruns = run.size();

    summaryList *sld;

    for (uint i = 0; i < nruns; i++) {
        printf("- - - - first signal run %d\n", run[i][0]);
        for (uint j = 0; j < ctime.size(); j++) {
            sld = processTEvents(ipath, run[i], prun, 6.6, 3., ctime[j]);
            sld->addFloat(ctime[j], "CoincGate_par");
            sld->addFloat(tdgap[i], "TargDetGap_par");
            sld->addInt(delay[i], "gammaDelay_par");
            sld->addInt(target[i], "Target_par");
            sld->addInt(run[i][0], "Run_par");
            vdb.push_back((sld->getValues()).data());
            if (stitle == "") {
                stitle = sld->getNames(1);
            }
            delete sld;
            /*
        v.insert(v.begin(),ctime[j]);
        v.insert(v.begin(),tdgap[i]);
        v.insert(v.begin(),(float) delay[i]);
        vdb.push_back(v);
            */
        }
    }

    printf("----- Summary\n%s\n", stitle.Data());

    for (uint i = 0; i < vdb.size(); i++) {
        printf("%s\n", vdb[i].data());
    }

    /*
    for (uint i=0;i<vdb.size();i++) {
      for (uint j=0;j<vdb[i].size(); j++) {
        printf("%.2f ",vdb[i][j]);
      }
      printf("\n");
    }
    */
};

int Check_Coinc() // std::vector<int> runs)
{

    /*
    This is the way in which time info is managed in the run_pkup_sall.root file
    ("gtime" and "gdate" are from the data of pkup and cube) ctime = gtime; cdate
    = gdate; utime = getUTime(gdate, gtime) - 2*3600; // convert to UTC, ntof
    seems to be CET, FERS use UTC; in Oct/23 CET = UTC+2h ## IT MAY CHANGE FOR
    OTHER PERIOD!
    */

    /*
    //This is the function getUTime:

    //convert ntof time to unix-time (https://cplusplus.com/reference/ctime/tm/)

    Long64_t getUTime(int date, int time) { // date and time in ntof format
      time_t timer;
      struct tm y2k; // = {0,0,0,0,0,0,0,0,0,0};
      y2k.tm_sec = time % 100; // 0-59
      y2k.tm_min = ((int) time / 100 ) % 100; // 0-59
      y2k.tm_hour = ((int) time / 10000); // 0-23
      date = date - 1000000;
      y2k.tm_mday = date % 100; // day of the month starting from 1
      y2k.tm_mon = ((int) (date/100))%100 - 1; // months since January, starting
    from 0! y2k.tm_year = ((int) (date/10000)) + 100; // year from 1900

      timer = mktime(&y2k); // unix time with Daylight Saving Time ??

      return (Long64_t) timer;

    }
    */

    // Get old file, old tree and set top branch address
    TString pkup_dir = "../test_231020/data/cube";
    gSystem->ExpandPathName(pkup_dir);
    const auto pkup_filename = gSystem->AccessPathName(pkup_dir) ? "./run_pkup_sall.root" : "../test_231020/data/cube/run_pkup_sall.root";

    TString fers_dir = "../test_231020/data/fers";
    gSystem->ExpandPathName(fers_dir);
    const auto fers_filename = gSystem->AccessPathName(fers_dir) ? "./Run29.root" : "../test_231020/data/fers/Run29.root";

    TFile pkup(pkup_filename);
    TTree *pkup_tree;
    pkup.GetObject("tsel", pkup_tree);

    TFile fers(fers_filename);
    TTree *fers_tree;
    fers.GetObject("tlist", fers_tree);

    const auto nentries = pkup_tree->GetEntries();

    Long64_t pkup_utime;
    pkup_tree->SetBranchAddress("utime", &pkup_utime);

    Float_t fers_utime;
    fers_tree->SetBranchAddress("timeus", &fers_utime);

    // Create a new file + a clone of pkup tree in new file
    TFile coinc_file("../test_231020/Coincidences/coinc.root", "recreate");
    auto coinc_tree = pkup_tree->CloneTree(0);

    for (int i = 0; i < nentries; ++i) {
        pkup_tree->GetEntry(i);

        // if (pkup_utime == fers_utime)
        coinc_tree->Fill();
    }

    coinc_tree->Print();
    coinc_file.Write();

    return 0;
}

int main() {

    /*std::vector<int> run_numbers;
    for (const auto &pair : run_delays)
    {
        run_numbers.push_back(pair.first);
    }
    std::sort(run_numbers.begin(), run_numbers.end());*/
    const std::vector<int> run_numbers = {15, 16, 17, 18, 19, 20, 21, 28, 29, 30}; // 30 has 10000ns delay but very heavy file
    const std::vector<int> ped_run_numbers = {13, 22, 33};

    try {
        std::ofstream output_file("../code/sipm/output.txt");
        (void)!freopen("../code/sipm/output.txt", "w", stdout);
        timing_analysis("../test_231020/data", run_numbers, ped_run_numbers, 12., 3., bar_time_width * 4);
        // processTEvents("../test_231020/data", run_numbers, ped_run_numbers, 5., 5., bar_time_width * 4., "cube/run_pkup_sall.root");
        output_file.close();
    } catch (const std::exception &e) {
        handle_exception(e);
        return 1;
    }
    return 0;
}