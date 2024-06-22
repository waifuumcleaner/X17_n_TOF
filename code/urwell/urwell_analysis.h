/*
 * urwell_analysis.cpp header file, containing basic struct and class definitions
 */
#ifndef URWELL_ANALYSIS_H
#define URWELL_ANALYSIS_H

#include <vector>

#include <TFile.h>
#include <TProfile.h>
#include <TString.h>

constexpr int n_apv = 10;    // number of APV cards
constexpr int n_sample = 30; // max number of time samples per channel

// single strip event data
struct STRIPS {
    int idc, ids; // chamber and strip indeces

    int nsam; // number of samples over threshold (multiplicity)

    int asam;   // max adjacent samples
    int adjs;   // first sample in max adjacent slot
    float adjq; // total charge of max adjacent samples
    float adjc; // centroid of max adjacent samples
    // fit parameters of adjacent samples

    float totq;   // total charge (over threshold)
    float intq;   // integrated charge (over and below threshold)
    float centro; // sample centroid (all over threshold)
    float first;  // first sample over threshold
};

// event info
struct EVENTINFO {
    int evt;    // event index
    float time; // event time in us from first event (first event = 0)
};

// time cluster: consecutive samples over threshold
struct TCLUSTER {
    int ch;   // absolute channel index (from 0 = APV0,CH0 to 1024 = APV7,CH127)
    int nt;   // number of consecutive samples
    int t0;   // first sample
    float qt; // total charge (pedestal subtracted)
    float ct; // sample centroid
    // fit parameters of adjacent samples
    int lsc; // link to spatial cluster index
};

// spatial cluster: adjacent strips above threshold in time coincidence
struct SCLUSTER {
    int ur;   // chamber index
    int ns;   // number of adjacent strips forming the cluster
    float qs; // total charge (pedestal subtracted)
    float cs; // spatial centroid of cluster
    float ms; // propagation of space (x = m * t + x0) where x is the strip - store the atan(m)
    float mq; // proagation of charge (q = m * t + q0) q is the charge / store the atan(m)
              // fit parameters of adjacent samples
};

// time-space cluster
struct TSCLUSTER {
    int ur;   // chamber
    int nn;   // number of time-space points
    int s0;   // first strip in chamber
    int s1;   // last strip
    int t0;   // first sample
    int t1;   // last sample
    float q;  // total charge (pedestal subtracted)
    float ct; // sample centroid
    float cs; // spatial (strip) centroid
    float rs; // time-space cross correlation coefficient (Pearson)
    float ms; // propagation of space (x = m * t + x0) where x is the strip - store the atan(m)
    float mq; // propagation of charge (q = m * t + q0) q is the charge / store the atan(m)
};

struct PEDES { // pedestal structure
    std::vector<float> mean[n_apv];
    std::vector<float> sigma[n_apv]; // pedestal sigma (if negative -> hot channel, to be masked)
};

// from apv and ch (electronics) indeces to absolute channel index (or strip index)
int e2a(int iapv, int irch);

// from absolute channel to apv index
int a2apv(int iach);

// from absolute channel to apv channel index
int a2ch(int iach);

// from absolute channel to chamber index
int a2c(int iach);

// from absolute channel to strip index in chamber
int a2s(int iach);

// from chamber and strip to absolute channel
int c2a(int icha, int istrip);

class mmRawTree;

/**
 * @brief
 * Plots pedestals and detects hot channels (std > hothr).
 * The sigma of hot channels is set to a negative value.
 *
 * @param rawT Pointer to the raw data tree
 * @param infile Input file name (default is "run0.root")
 * @param fIn Pointer to the input file (default is null)
 * @param hothr Threshold for detecting hot channels (default is 100)
 * @return PEDES The resulting pedestal data
 */
PEDES getPedestal(mmRawTree *rawT, // used to get the module names
                  TString infile = "run0.root", TFile *fIn = 0, float hothr = 100);

/**
 *
 * @brief
 * Subtracts noise on sample series of given channel.
 * Gets raw ADC data, subtracts pedestal and common noise,
 * sets to 0 values below nsig*sigma_pedestal
 *
 * @param adc Vector of raw ADC data
 * @param iapv apv index
 * @param ich Channel index
 * @param ped Pedestal data vector
 * @param cnoise Common noise for the selected channel
 * @param nsig Number of pedestal std deviations (sigma) defining the noise level
 * @param hpo TProfile to be filled, to see charge correlation with channels
 *
 * @return Sample signal vector, set to 0 if below nsigma
 */
std::vector<float> pedestalSubtract(std::vector<short> adc, int iapv, int ich, PEDES ped, std::vector<float> cnoise, float nsig = 3, TProfile *hpo = NULL);

/**
 * @brief
 * Single Strip sample analysis.
 * Searches for temporal clustering on single strip.
 *
 * @param adc Vector of ADC signal data, with pedestal and common noise subtracted
 * @param iapv apv index
 * @param ich Channel index
 *
 * @return Temporal cluster vector
 */
std::vector<TCLUSTER> samplesAnalysis(std::vector<float> adc, // signal, noise suppressed
                                      int iapv, int ich);

class SpaceCluster {
  private:
    std::vector<float> last_strip;
    std::vector<float> last_sample;

    std::vector<float> sum_xi; // time
    std::vector<float> sum_yi; // space
    std::vector<float> sum_xiyi;
    std::vector<float> sum_xi2;
    std::vector<float> sum_zi; // charge
    std::vector<float> sum_xizi;

    std::vector<SCLUSTER> cluster;

    float sigma;      // sigma threshold for simultaneous track [sample units]
    int thr_nsample;  // number of consecutive samples forming a signal time cluster
    int max_adj_mask; // max number of adjacent masked strip separating a single cluster; if two strips are separated by more than this number, they cannot form a single cluster

  public:
    SpaceCluster(float sigmat = 2, int nsample = 4, int maxmask = 3);

    ~SpaceCluster();

    /**
     * @brief
     * Search for spatial clusters in time coincidence.
     *
     * @param tc Sample cluster structure vector
     * @param deltas Next strip step, 0 means new chamber, typically 1; 2 or more in case of previously masked channels
     *
     * @return Just 0 for some reason
     */
    int searchSpatialCluster(std::vector<TCLUSTER> &tc, int deltas);

    /**
     * @brief
     * Computes centroid and other parameters, to be called at the end of loop on strips (before filling tree)
     */
    void finalizeClustering();

    uint size();

    SCLUSTER getCluster(uint i);

    void clear();
};

/**
 * @brief
 * Estimates common noise for all time samples and APV.
 * Skips hot channels (impact is negligible)
 *
 * @param tt Pointer to the raw data tree
 * @param ped Pedestal data vector
 * @param method Selects the method (0: median (default), ...)
 *
 * @return Vector[napv][nsample]
 */
std::vector<std::vector<float>> commonNoise(mmRawTree *tt, PEDES ped, int method = 0);

/**
 * @brief
 * BRIEF TO BE DONE
 *
 * @param phys_run Physics run number
 * @param ped_run Pedestal run number. If <0 assumes phys_run is already pedestal subtracted
 * @param nsigma Number of sigma for signal threshold
 * @param min_ncsample Min number of consecutive samples forming a signal time cluster
 * @param max_ncsample Max number of consecutive sample that cannot be exceeded to form a cluster (TBC)
 * @param verbose Both for verbosity and event display: if>0 show event every "verbose" delta; if 0, no event display and minimal verbosity, if <0 affect only verbosity as absolute value
 * @param outType Data written on output file. 0: cluster information, 1: track and cluster information ...
 * @param dpath Data path of raw data
 *
 * @return TO BE DONE
 */
int decodePhysRun(int phys_run = 1, int ped_run = -1, float nsigma = 3, int min_ncsample = 6, int max_ncsample = 20, int verbose = 0, int outType = 0, TString dpath = "../../test_231020/data/srs");

/**
 * @brief
 * BRIEF TO BE DONE
 *
 * @param nrun Physics run number
 * @param thr_q Minimum charge of cluster
 * @param thr_ns Number of minimum strips forming a cluster
 * @param thr_nt Number of consecutive sample above threshold forming a signal (number of samples)
 * @param dpath Data path of raw data
 *
 * @return TO BE DONE
 */
int readAna(int nrun = 1, Float_t thr_q = 1000, Int_t thr_ns = 1, Int_t thr_nt = 5, TString dpath = "../../test_231020/data/srs");

/**
 * @brief
 * Sort of template of post processing (of decoded SRS data)
 *
 * @param run Physics run number
 * @param dpath Data path of raw data
 */
void readClusterOut(int run = 1, TString dpath = "../../test_231020/data/srs");

/**
 * @brief
 * Sort of template of post processing (of decoded SRS data)
 *
 * @param run Physics run number
 * @param dpath Data path of raw data
 */
void readTrackOut(int run = 1, TString dpath = "../../test_231020/data/srs");

int urwell_analysis_main();
#endif