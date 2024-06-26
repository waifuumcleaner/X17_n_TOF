#ifndef SIPM_CLASSES_H
#define SIPM_CLASSES_H

#include "TChain.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TParameter.h"
#include <TObjString.h>
#include <ctime>
#include <optional>
#include <set>
#include <unordered_map>

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

// TO BE REMOVED!!
inline std::unordered_map<int, int> run_delays = {{15, 0}, {16, 500}, {17, 1000}, {18, 1500}, {19, 2000}, {20, 2500}, {21, 3000}, {28, 3500}, {29, 4000}, {26, 3000}, {27, 1000}, {30, 10000}, {31, 0}, {44, 500}, {45, 2000}, {46, 10000}, {48, 10000}, {49, 2000}, {53, 2000}, {54, 0}, {55, 500}, {56, 2000}, {57, 2000}, {58, 2000}};

enum class TargetType { Light, Medium, Heavy };
enum class Mode { Pedestal, Physical };

class RunsInfo {
  private:
    std::vector<int> m_numbers;
    std::vector<int> m_delays;             // ns delay from g-flash
    std::vector<TargetType> m_targetTypes; // light, medium or heavy
    std::vector<Mode> m_modes;             // pedestal or physical run

  public:
    RunsInfo(std::vector<int> numbers = {}, std::vector<int> delays = {}, std::vector<TargetType> targetTypes = {}, std::vector<Mode> modes = {});

    void addRun(int number, int delay, TargetType target, Mode mode);

    std::vector<int> getNumbers(const std::optional<std::vector<int>> &subset = std::nullopt) const;

    std::vector<int> getDelays(const std::optional<std::vector<int>> &subset = std::nullopt) const;

    std::vector<TargetType> getTargetTypes(const std::optional<std::vector<int>> &subset = std::nullopt) const;

    std::vector<Mode> getModes(const std::optional<std::vector<int>> &subset = std::nullopt) const;

    std::vector<int> getPedestalRuns() const;
    std::vector<int> getPhysicalRuns() const;
    void setNumbers(const std::vector<int> &numbers);
    void setDelays(const std::vector<int> &delays);
    void setTargetTypes(const std::vector<TargetType> &targetTypes);
    void setModes(const std::vector<Mode> &modes);
    void setMode(int run_number, Mode mode);
    void setTargetType(int run_number, TargetType targetType);
};

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

#endif // SIPM_CLASSES_H
