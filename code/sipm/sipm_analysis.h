/*
 * FERS + Scintillators bars utility functions
 *
 * Version: 0.1
 * Date: Oct/2023
 */

#ifndef SIPM_ANALYSIS_H
#define SIPM_ANALYSIS_H

#include "TChain.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TParameter.h"
#include <TObjString.h>
#include <ctime>
#include <set>
#include <unordered_map>

#include "../cube/ntofDAQ.h"

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

inline std::unordered_map<int, int> run_delays = {{15, 0}, {16, 500}, {17, 1000}, {18, 1500}, {19, 2000}, {20, 2500}, {21, 3000}, {28, 3500}, {29, 4000}, {26, 3000}, {27, 1000}, {30, 10000}, {31, 0}, {44, 500}, {45, 2000}, {46, 10000}, {48, 10000}, {49, 2000}, {53, 2000}, {54, 0}, {55, 500}, {56, 2000}, {57, 2000}, {58, 2000}};

struct coinc_counter {
    int run_n;
    int ns_delay;
    int trig_n;
    int coinc_n;
    int double_coinc_n;
    int cross_coinc_n;

    coinc_counter();
    coinc_counter(int run_number, int gflash_delay_in_ns, int trigger_n, int coincidences_n, int double_coincidences_n, int cross_coincidences_n);
};

class coincidence_Info {
  private:
    std::vector<coinc_counter> counters;

  public:
    coincidence_Info();
    void add_counter(const coinc_counter &counter);
    const std::vector<coinc_counter> &get_counter() const;
};

/**
 * @class summaryList
 * @brief Collects summary data, generally from a run
 * @see coinc_counter
 */
class summaryList {
  private:
    TList *iList; // list of integers
    TList *fList; // list of floats ... to be improved merging them
    TParameter<float> *fPar;
    TParameter<int> *iPar;
    std::vector<coinc_counter> counters;

  public:
    summaryList(TString name = "list");
    ~summaryList();
    void resetLists();
    void addFloat(float val, TString name);
    void addFloat(TParameter<float> *p);
    void addInt(int val, TString name);
    void addInt(TParameter<int> *p);
    std::string getNames(int format = 0);
    std::string getValues();
    TParameter<float> *getFloat(int idx);
    TParameter<int> *getInt(int idx);
    void merge(summaryList *sl);
};

/**
 * @brief
 * Returns the active channels either from the static "db" or from the provided TTree
 */
std::vector<int> getActiveChannels(TTree *tl = NULL);

/**
 * @brief
 * Determines the side of the sipm channel
 * @return -1 if channel is masked, or 0: top, 1: right, 2: bottom, 3: left
 *
 */
int c2side(int ch);

/**
 * @brief
 * Maps channel to x / y sipm position
 */
float c2x(int ch);

float c2y(int ch);

/**
 * @class manageTree
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
    manageTree(TString sofile);
    ~manageTree();
    void allocBranches(int acqmode);

    // returns true if channel is masked
    bool isChannelMasked(int ch);

    void evalCentroids();
    void resetVars();
    void setTimeTrg(float time, float trg);
    void setChannel(float board, float channel);
    void setADCs(float lgadc, float hgadc);
    void setTime(float toav, float totv);
    void setVars(float *arr);
    void setBeamInfo(float pulse, float dtime);

    int proFill();
    int run_number;
    int acqmode;
    int histoch;
    float timebin; // [ns]

    TString start_time;
    void setUserData(int runnum, float timebin, int histoch, TString startime);
    TTree *getTree();
    void Print();
    void Save();
};

/**
 * @brief
 * Parses FERS output run list text file and saves data into root file.
 *
 * @param ntofRoot: root file with ntof pulse data extracted by ntofData.cpp macro
 *
 * @return Number of trigger loaded/parsed (< 0 on error)
 *
 * spectroscopy gain data should go from 0 to up to "Energy N Channels" parameter, maximum is 8192 - 13 bits) but in output ascii file there are values larger than this limit and also approaching 2^15; probably a bug!
 */
int parseListFile(int run, TString basepath, TString prefix, TString ntofRoot);

/**
 * @brief
 * Checks if run is already in root format,
 * if not tries to parse text file (FERS output) and convert to root.
 * To force parsing and conversion, use a non-existing prefix.
 * The optional ntofRFile can be used to sync with the Cube/nTOF data (time and beam pulse to distinguish main and parassitic pulses, TO BE DONE: add the info on the Cube signal itself)
 */
int checkAndConvertText2Root(TString basepath, std::vector<int> vrun, TString prefix, TString ntofRFile);

/**
 * @brief
 * Setting ROOT plot style
 */
void set_local_style();

/**
 * @brief
 * Used to cut pedestal time over threshold distribution
 *
 * @return
 * Bin index of pedestal time over threshold histogram
 */
int find_first_empty_bin_after_max(TH1F *histogram);

/**
 * @brief
 * Gets pedestal of time over threshold.
 * Can be used to plot time over threshold distributions of each channels;
 * useful to evaluate equalization of channels responses
 */
std::vector<std::vector<float>> get_ToT_pedestal(const std::vector<int> act_ch_v, // vector of active channels
                                                 const float def_rms = 5.,        // [ns] assumed rms where there is no channel data
                                                 const TString basepath = "../../test_231020/data", const std::vector<int> ped_runs = std::vector<int>{13, 22, 33});

/**
 * @brief
 * Plots all relevant information gathered in channelTimeProcess
 */
void plot_timing_data(TH1F *h_time_diff, TH1F *h_coinc_per_trigger, std::vector<TH1F *> dead_t_ch_even, std::vector<TH1F *> dead_t_ch_odd, TH2F *h_toa_corr, TH1F *h_toa, TGraph *tot_global_corr, TH1F *h_tot_prod, std::vector<TGraph *> &tot_corr_ch, TH2F *hxy_channels, TH2F *hxy_bars, TH2F *hxy_bars_toa, TH2F *hxy_bars_tot, TH2F *hxy_cross_bars, const float thr_toa, const int run_number);

/**
 * @brief
 * Draws a red square around the crossing bars region in a canvas with xy hit map
 */
void draw_square(TCanvas *canvas, const std::string &plot_type);

/**
 * @return
 * 0 if the channels are at exactly opposite ends of the bar(s), 1 if they see the same bar(s) only partially, -1 otherwise
 */
int check_couple_or_adjacent(const int chj, const int ch);

/**
 * @return
 * Index of the scintillator bar seen by the channels, assuming a coincidence happened:
 * 0-5 for vertical bars (left to right), 6-11 for horizontal bars (top to bottom).
 * If the channels don't see the same scintillator bar returns -1.
 */
int chan_to_bar_index(const int chj, const int ch);

/**
 * @brief
 * Fills a TH2F xy map with opposite channels coincidences
 */
void fill_xy_opposite_channels_coinc(const int ch, TH2F &xy_map);

/**
 * @brief
 * Fills a TH2F xy map with scintillator bar hits.
 * OPEN QUESTION: a che barra associare la coincidenza di due canali sulle
 * giunzioni tra due barre? Per ora a quella che vedono per 2/3 a priori
 */
void fill_xy_scintibar(const int bar_index, TH2F &xy_map);

/**
 * @return
 * An integer between 0 and (max-1), where max is the number of sections the scintibar can be
 * divided in. The returned integer represents the portion of the bar that has been hit.
 * 0 represents top if the bar is vertical, left if the bar is horizontal.
 */
int toa_coordinate(const int ch_i, const float t_i, const int n_bar_portions);

/**
 * @brief
 * Fills a TH2F xy map with scintillator bar hits, using toa info -> not whole bar fires.
 * OPEN QUESTION: a che barra associare la coincidenza di due canali sulle
 * giunzioni tra due barre? Per ora a quella che vedono per 2/3 a priori
 */
void fill_xy_scintibar_toa(const int bar_index, const int ch_i, const float t_i, const int n_bar_portions, TH2F &xy_map);

/**
 * @brief
 * TO BE IMPROVED (CODE EFFICIENCY AND CONCISENESS).
 * Fills a TH2F xy map with scintillator bar hits when 2 bars have been hit.
 */
void fill_xy_cross_scintibar(const int bar_index1, const int bar_index2, TH2F &xy_map);

/**
 * @brief
 * Main function of sipm analysis. Extracts timing information and looks for coincidences for a specific run.
 *
 * @param run_number Number of FERS run
 * @param tl TTree of the chosen run
 * @param coincidence_Info Coincidence information class to be filled in the analysis
 * @param thr_toa Chosen time window (in ns). Default: time traveled by light in the scintillator (~2.5 ns)
 * @param nsigma Number of pedestal time over threshold std deviations after which accept signals
 * @param vped pedestal time over threshold vector
 */
summaryList *channelTimeProcess(const int run_number, TTree *tl, coincidence_Info &coinc_info, const float thr_toa = bar_time_width, const float nsigma = 5., const std::vector<std::vector<float>> &vped = std::vector<std::vector<float>>());

/**
 * @brief
 * TO BE ADJUSTED AND IMPROVED.
 * Plots relevant (hopefully) data and produce x/y hit map
 */
summaryList *plotTimeData(TChain *chainT);

/**
 * @brief
 * Reads and processes timing events. Optionally correlates each event to the cube root data.
 *
 * @param basepath Path to data files
 * @param srun Signal runs, requires "{}" brackets, if <0 looks at pedestal only
 * @param prun Pedestal runs, can be cumulated, requires "{}" brackets, if <0 not considered
 * @param totsig Default time over threshold std deviation [ns]
 * @param nsigma Number of pedestal time over threshold std deviations after which accept signals
 * @param time_window Chosen time window (in ns). Default: time traveled by light in the scintillator (~2.5 ns)
 * @param ntoFile Cube data root file
 */
summaryList *processTEvents(const std::string basepath = "../../test_231020/data", std::vector<int> srun = {-1}, std::vector<int> prun = {13, 22, 33}, float totsig = 5., float nsigma = 5., float time_window = bar_time_width, TString ntoFile = "/cube/run_pkup_sall.root");

/**
 * @brief
 * Analyses multiple runs in timing mode (adapted to 2310 CERN test)
 */
void timing_analysis(const std::string basepath = "../../test_231020/data", const std::vector<int> sig_runs = {29}, const std::vector<int> ped_runs = {13, 22, 33}, float tot_default_rms = 12., float tot_rms_n = 3., float time_window = bar_time_width);

int sipm_analysis_main();

#endif
