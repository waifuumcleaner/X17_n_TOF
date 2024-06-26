/*
 * Utility functions, classes and structures
 *
 */
#ifndef UTILF_H
#define UTILF_H

#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <iostream>

// signal on single strip with time samples: consecutive samples over threshold
struct uTSIGNAL {
    int imod;             // chamber (module) index
    int istrip;           // chamber strip
    int nt;               // number of consecutive samples
    int t0;               // first sample
    float qt;             // total charge (pedestal subtracted)
    std::vector<float> q; // charge on each sample
    float centroid;       // sample centroid (time centroid)
    float rms;            // sample rms (time rms)
    int clidx;            // cluster index the signal belong to
    int fitidx;           // index of the fitting function
    float fitpar[10];     // fit parameters, up to 10 (?)
    float tstart;         // signals starting time (from fit)

    uTSIGNAL() : imod(0), istrip(0), nt(0), t0(0), qt(0.0f), centroid(0.0f), rms(0.0f), clidx(0), fitidx(0), tstart(0.0f) { std::fill(std::begin(fitpar), std::end(fitpar), 0.0f); }
};

// cluster collecting spatially adjacent strips with temporal compatible signals (having temporal overlap)
typedef struct {
    std::vector<int> its; // vector of indices of uTSIGNAL strip with time-signals belonging to the cluster
    float tstart;         // time (sample) start of cluster
    float tstop;          // time (sample) end of cluster
    float tdstart;        // max time difference between the signal time starts (this should be "constant" for a given chamber gap)
    float centroid;       // along strips (space centroid)
    float rms;            // space rms
    float totcharge;      // total charge of all signals in cluster
    float sfirst, slast;  // first and last strips
    float m, q;           // s = m * t + q  # fit of strip vs start time of signals on strips (at lest 2 strips in cluster)
} uCLUSTER;

/* = = = = = = = = = = = = = = = = = = =
 *
 * fvts - vector of signals
 * fvcl[cluster].its[signal] = index of fvts
 *
 */

class tsCluster {
  private:
    double *qts; // charge versus channel versus sample

    std::vector<uTSIGNAL> fvts; // vector of strips with time signal
    //  std::vector<std::vector<int>> fvcl; // cluster vector, contain vector of indices of fvts
    std::vector<uCLUSTER> fvcl;

    std::vector<std::vector<int>> vtrack; // vector of tracks, each element contains the indices of the participating clusters

    int nurw;    // number of chambers
    int ncha;    // channels (strips) / chamber
    int nsam;    // samples / channel
    int *maskCh; // array of the masked channel (if 1 channel is masked)
    int countmask;
    float fchargemax;

    TF1 *fRT_pulse; // APV25 fitting shapes

    int fverbose; // verbose printing

    int getIdx(int ur, int ch, int sa) {
        int idx = ur * ncha * nsam + ch * nsam + sa;
        return idx;
    };

  public:
    tsCluster(int nurwell, int nchannel, int nsample, int verbose) { // module/chamber, strips/module, number of time samples/channel
        nurw = nurwell;
        ncha = nchannel;
        nsam = nsample;

        fverbose = TMath::Abs(verbose);

        qts = new double[nurw * ncha * nsam]; // charge versus chambers x channels x samples
        maskCh = new int[nurw * ncha];
        memset(maskCh, 0, nurw * ncha * sizeof(int));
        countmask = 0;

        //    fRT_pulse = new TF1("fRT_pulse","[2]*(1.-exp(-(x-[3])/[0]))*exp(-(x-[3])/[1])*(TMath::Sign(0.5,x-[3])+0.5)",-0.5,(nsample-0.5));
        //    fRT_pulse[0] = new TF1("fRT_pulse","[2]*TMath::Landau(x,[0],[1])",-0.5,nsam-0.5);
        // fRT_pulse[1] = new TF1("fRT_pulse","[2]*TMath::Landau(x,[0],[1])+[5]*TMath::Landau(x,[3],[4])",-0.5,nsam-0.5);
        fRT_pulse = new TF1("fRT_pulse", "[2]*TMath::Landau(x,[0],[1])+[5]*TMath::Landau(x,[3],[4])+[8]*TMath::Landau(x,[6],[7])", -0.5, nsam - 0.5);
    };

    ~tsCluster() { delete qts; };

    void maskChannel(int iurw, int ich) { // chamber, strip
        if (fverbose > 1)
            printf(" mask channel %d %d\n", iurw, ich);
        maskCh[iurw * ncha + ich] = 1;
        countmask++;
    };

    // set charge-array to 0
    void reset() {
        memset(qts, 0, nurw * ncha * nsam);
        fvcl.clear();
        fvts.clear();
        vtrack.clear();
    };

    // fill charge array with vector of sample of given strip/chamber
    void fillQ(int iurw, int ich, std::vector<float> tsample) {
        uint nn = tsample.size();

        for (uint i = 0; i < nn; i++) { // loop on sample
            qts[getIdx(iurw, ich, i)] = tsample[i];
        }
    };

    void checkMask() {
        int count = 0;
        for (int i = 0; i < nurw * ncha; i++) {
            if (maskCh[i] > 0)
                count++;
        }
        if (count > countmask) {
            if (fverbose > 1)
                printf("Masked channel in db: %d\n", count);
        }
    };

    // get functions

    float chargeMax() { return fchargemax; };
    int numClusters() { return (int)fvcl.size(); };
    int numSignals(int cl) { return (int)fvcl[cl].its.size(); }; // num signals (strips) in given cluster
    int numTracks() { return (int)vtrack.size(); };

    int clusterModule(int idx) { // return chamber associated at the given cluster index
        int isig = fvcl[idx].its[0];
        return fvts[isig].imod;
    }

    float clusterSlope(int idx) { // return cluster linear-fitted strip-time slope (m)
        return fvcl[idx].m;
    };

    float clusterConst(int idx) { // return cluster linear-fitted strip-time constant term (q)
        return fvcl[idx].q;
    };

    uTSIGNAL signal(int cl, int is) { // cluster and signal index
        int idx = fvcl[cl].its[is];
        return fvts[idx];
    };

    // return std::vector of uTSIGNAL associated to the given module (chamber)
    std::vector<uTSIGNAL> signalInMod(int mod) {
        std::vector<uTSIGNAL> vret;
        int n = fvts.size();
        for (int i = 0; i < n; i++) {
            if (fvts[i].imod == mod) {
                vret.push_back(fvts[i]);
            }
        }
        return vret;
    };

    TF1 *getSignalFitFun(int idx = 0) {
        std::cout << "Don't know what to do with given idx " << idx << " yet\n";
        return fRT_pulse;
    };

    /*
     */

    /***
     * Estimate cluster main parameters:
     *    total charge (sum of each signal charge)
     *    space centroid and rms (weighted for signal charge)
     *    time_start, time_stop (by fitting Landau functions and interpolating fitted function)
     *    (the 3 possible fitting functions are Landau, Landau+Landau, Landau+Landau+Landau), chi2r is pretty high!
     *
     * Find track:
     *    search clusters on adjacent chambers with time starts within given delta_time
     *    track shall have signal (cluster) on each chamber
     *    procedure is currently limited and biased by cluster orders ... need improvements!!!
     *
     * fvcl : vector of clusters with indices to fvts
     * fvts : vector of time signals (for each channel)
     *
     * sample_period : [ns] the time between samples
     * chi2r_thr : max chi2r a "good" fit can get
     * dtstart_max : max difference between cluster starts difference on the same tracks on adjacent chambers
     *
     */

    int clusterAnalysis(float chi2r_thr = 5000., float dtstart_max = 4) {

        // float sample_period = 25.; // [ns]
        //  float funTau0 = 200. / sample_period; // raising time constant (sample unit of 25 ns) // it is 8
        //  float funTau1 = 900. / sample_period; // it is 36

        TGraph *gr = new TGraph(nsam);

        std::vector<int> mod2cl[nurw]; // list of clusters beloging to each module (chamber/urwell)

        if (fverbose > 1)
            printf(" Cluster (and track) analysis start\n");
        for (int i = 0; i < (int)fvcl.size(); i++) { // loop on identified clusters
            if (fverbose > 2)
                printf("  Cluster %d: ", i);

            int imod = -1;

            float xcentro = 0;
            float xrms = 0;
            float allcharge = 0;

            float thmin = 9e99;    // cluster first start time
            float tstart_max = -1; // cluster latest start time
            float thmax = -1;      // cluster latest "ending" time

            int numsig = (int)fvcl[i].its.size();

            //      TGraph *gtimestrip = new TGraph(numsig); // start time vs strip position (index)
            std::vector<float> st0; // start time (from fit) of each signal
            std::vector<float> sx0; // strip position of signal
            std::vector<float> sqt; // charge weight
            int strip0 = 999999;    // first and last strips
            int strip1 = -1;

            for (int j = 0; j < numsig; j++) { // loop on time signals (sample series) of the cluster

                int isig = fvcl[i].its[j]; // signal index

                fvts[isig].clidx = i; // cluster index the signal will belong to

                int strip = fvts[isig].istrip;
                strip0 = strip < strip0 ? strip : strip0;
                strip1 = strip > strip1 ? strip : strip1;

                imod = fvts[isig].imod;
                if (j == 0) {
                    if (fverbose > 2)
                        printf(" module %d\n", imod);
                }

                float qtot = fvts[isig].qt; // use total charge, max charge in single sample should be better

                xcentro += qtot * ((float)strip);
                xrms += qtot * ((float)(strip * strip));
                allcharge += qtot;

                int nst = fvts[isig].nt;
                int ns0 = fvts[isig].t0;
                std::vector<float> charge = fvts[isig].q;

                float funT0 = (float)ns0;
                float funT1 = (float)(ns0 + nst);

                float centroid = fvts[isig].centroid; // time centroid
                float rms = fvts[isig].rms;

                int countSig = 0;
                for (int k = 0; k < nsam; k++) {
                    float y = 0;
                    int h = k - ns0;
                    if ((h >= 0) && (h < nst)) {
                        y = charge[h];
                    }
                    if (y > 0) {
                        countSig += 1;
                    }
                    gr->SetPoint(k, (float)k, y);
                }

                if (fverbose > 3) {
                    printf(" sample evolution:\n");
                    gr->Print();
                }

                float chi2r;
                int kf;
                for (int kp = 3; kp < 9; kp++) {
                    fRT_pulse->FixParameter(kp, 0.);
                }

                float chi2rmin = 9e99;
                float t0half = -1; // begin and end time at "half maximum"
                float t1half = -2;

                int fparMax = countSig / 10 + 1;       // a fitting function cannot exceed 3 parameters for 10 points to be fitted
                fparMax = (fparMax < 3) ? fparMax : 3; // max up to 9 pars for the fitting function
                for (kf = 0; kf < fparMax; kf++) {     // try 3 fitting functions (sums of Landau), choose the one with minimum chi2r
                    float p1 = centroid + funT0 / (kf + 1) * kf;
                    p1 = (p1 > funT1) ? funT1 : p1;
                    fRT_pulse->SetParameter(0 + 3 * kf, p1);
                    fRT_pulse->SetParLimits(0 + 3 * kf, (float)ns0, (float)(ns0 + nst));
                    fRT_pulse->SetParameter(1 + 3 * kf, rms);
                    fRT_pulse->SetParLimits(1 + 3 * kf, rms * .1, rms * 5.);
                    fRT_pulse->SetParameter(2 + 3 * kf, qtot / (kf + 1)); // normalization
                    fRT_pulse->SetParLimits(2 + 3 * kf, 0., qtot * 5);

                    gr->Fit(fRT_pulse, "Q");
                    chi2r = fRT_pulse->GetChisquare() / ((float)fRT_pulse->GetNDF());
                    if (chi2r < chi2rmin) {
                        if (fverbose > 2)
                            printf("  fit %d, signal: %d (%d %d): ", kf, isig, imod, strip); //%.2f %.2f %.1f (%.1f)\n",isig, imod, istrip, fRT_pulse->GetParameter(0), fRT_pulse->GetParameter(1), fRT_pulse->GetParameter(2), chi2r); // corresponding signal index
                        for (int k = 0; k < 9; k++) {
                            fvts[isig].fitpar[k] = fRT_pulse->GetParameter(k);
                            if (k < (kf + 1) * 3) {
                                if (fverbose > 2)
                                    printf(" %.2f", fRT_pulse->GetParameter(k));
                            }
                        }
                        float fmax = fRT_pulse->GetMaximum(funT0, funT1);
                        float tmax = fRT_pulse->GetMaximumX(funT0, funT1);
                        float fmin0 = fRT_pulse->GetMinimum(funT0, tmax);
                        float fmin1 = fRT_pulse->GetMinimum(tmax, funT1);
                        float fthr = fmin0 + (fmax - fmin0) / 2.;
                        if (tmax < 1) {          // signal started before first sample and only tail acquired, better remove from analysis
                            t0half = funT1 + 1;  // skip this signal
                            t1half = t0half - 1; // NEED IMPROVEMENT THERE ARE STILL INF VALUES @@@@@
                            break;
                        } else {
                            t0half = fRT_pulse->GetX(fthr, funT0, tmax); // half maximum - assumed as starting time of the signal
                        }

                        fthr = fmin1 + (fmax - fmin1) / 2;
                        t1half = fRT_pulse->GetX(fthr, tmax, funT1); // half maximum - assumed as ending time of the signal

                        fvts[isig].tstart = t0half;

                        thmin = (thmin > t0half) ? t0half : thmin;
                        thmax = (thmax < t1half) ? t1half : thmax;

                        tstart_max = (tstart_max < t0half) ? t0half : tstart_max;

                        if (fverbose > 2)
                            printf(" (chi2r %.1f) %f %f\n", chi2r, thmin, thmax); // t0half, tmax, t1half, funT0, funT1);
                        fvts[isig].fitidx = kf;
                        chi2rmin = chi2r;
                    }

                    if (chi2r < chi2r_thr) {
                        break;
                    } // good fit if chi2r below given threshold
                } // end fit loop

                if (t0half > t1half) { // fit did not succeded
                    if (fverbose > 0)
                        printf("WARNING: bad sample fit for signal %d, not saved\n", isig);
                } else { // add data for TPC fit @@ should pass some check for cleaner data @@
                    st0.push_back(t0half);
                    sx0.push_back((float)strip);
                    sqt.push_back(qtot); // not used ... should represent the weight of the fit
                }
            } // end loop on time signals in cluster

            // try linear fit time vs strip for single signal
            fvcl[i].m = 0;                                        // invalid value
            fvcl[i].q = -1e5;                                     // invalid value
            if ((st0.size() > 1) and (tstart_max - thmin) > 0.) { // at least two strips with "good" sample fit and deltas between signals start times is larger than 0
                TGraph *gdummy = new TGraph(st0.size(), &st0[0], &sx0[0]);

                TFitResultPtr r = gdummy->Fit("pol1", "SQ");
                fvcl[i].m = r->Parameter(1); // NOTE unit is [strip index]/[25 ns]!!!!
                fvcl[i].q = r->Parameter(0); // unit is [strip index]
                delete gdummy;

                if (fverbose > 1)
                    printf("     Linear Fit (%d): %.2f %.2f\n", (int)st0.size(), fvcl[i].m, fvcl[i].q);
                st0.clear();
                sx0.clear();
                sqt.clear();
            }

            mod2cl[imod].push_back(i); // add cluster to module (chamber)

            fvcl[i].sfirst = strip0;
            fvcl[i].slast = strip1;

            fvcl[i].tstart = thmin;
            fvcl[i].tstop = thmax;
            fvcl[i].tdstart = tstart_max - thmin;
            fvcl[i].centroid = xcentro / allcharge;
            fvcl[i].totcharge = allcharge;
            if (numsig > 1) {
                fvcl[i].rms = TMath::Sqrt(xrms / allcharge - xcentro * xcentro / allcharge / allcharge);
            } else {
                fvcl[i].rms = 0.;
            }

            if (fverbose > 1)
                printf("   --> centroid %.2f +/- %.2f (tstart: %.2f)\n", fvcl[i].centroid, fvcl[i].rms, fvcl[i].tstart);

        } // end cluster loop

        // find tracks (clusters on two or more chambers, with consistent timing) - need improvement
        // NOTE: not applicable if chambers are not stackered for tracking
        // for (vector<uCLUSTER>::iterator vc = fvcl.begin(); vc != fvcl.end(); vc++) {
        std::vector<int> cktrack; // vector of cluster indices, one for each chamber
        // int ntrack = 0;
        // int nchamb = 0;
        for (int ic = 0; ic < (int)mod2cl[0].size(); ic++) { // loop on clusters of the first chamber
            int clm0 = mod2cl[0][ic];
            cktrack.push_back(clm0);

            float tstart_prev = fvcl[clm0].tstart;
            for (int im = 1; im < nurw; im++) {                          // loop on other chambers
                for (int icm = 0; icm < (int)mod2cl[im].size(); icm++) { // loop on clusters of the "im" chamber
                    int clm1 = mod2cl[im][icm];
                    float t0d = TMath::Abs(fvcl[clm1].tstart - tstart_prev);
                    if (t0d < dtstart_max) { // assume clusters on two chambers belong to the same track if their starting times differ for less then dtstart_max
                        cktrack.push_back(clm1);
                        tstart_prev = fvcl[clm1].tstart;
                        break; // go to next chamber
                    }
                }
            }
            if ((int)cktrack.size() == nurw) { // all chambers have a signal cluster
                vtrack.push_back(cktrack);
                if (fverbose > 2) {
                    printf(" Track Clusters: ");
                    for (std::vector<int>::iterator cki = cktrack.begin(); cki != cktrack.end(); cki++) {
                        printf(" %d", *cki);
                    }
                    printf("\n");
                }
                cktrack.clear();
            }
        }

        delete gr;

        return 0;
    };

    /*
     * Output the track data on provided out stream file
     *
     * if evt<0 return the title string
     * time_from_start : time in second from start of run
     *
     */

    TString trackOutput(int evt, float time_from_start) {

        TString sret = "";

        // event id
        if (evt < 0) { // print title
            sret = "# EventId Time TrackId ChamberId TimeStart TimeDeltaStart TimeEnd Charge StripFirst StripLast StripCentroid StripRMS mlFit qlFit\n";
            return sret;
        }

        for (int i = 0; i < (int)vtrack.size(); i++) { // loop on tracks

            int track_id = i;

            std::vector<int> cc = vtrack[i]; // cluster lists

            for (int j = 0; j < (int)cc.size(); j++) {

                int k = cc[j]; // actual cluster index
                int chamber_id = clusterModule(k);
                float time_start = fvcl[k].tstart;
                float time_delta = fvcl[k].tdstart;
                float time_end = fvcl[k].tstop;
                // float time_centroid;
                float charge_total = fvcl[k].totcharge;
                int strip_first = fvcl[k].sfirst;
                int strip_last = fvcl[k].slast;
                float strip_centroid = fvcl[k].centroid;
                float strip_rms = fvcl[k].rms;
                float xmt = fvcl[k].m;
                float xqt = fvcl[k].q;

                sret = sret + Form("%d %f %d %d %.3f %.3f %.3f %.2f %d %d %.3f %.3f %.3f %.2f\n", evt, time_from_start, track_id, chamber_id, time_start, time_delta, time_end, charge_total, strip_first, strip_last, strip_centroid, strip_rms, xmt, xqt);
            }

            // printf("%s\n",sret.Data());
        }

        return sret; // (int) vtrack.size();
    };

    /*
     * Output the processed event cluster-data on provided out stream file
     *
     * if evt<0 return the title string
     */

    TString eventOutput(int evt, float time_from_start) {

        TString sret = "";

        // event id
        if (evt < 0) { // print title
            sret = "# EventId Time ChamberId NumStrips TimeStart TimeDeltaStart TimeEnd Charge StripCentroid StripRMS mlFit qlFit\n";
            return sret;
        }

        for (uint k = 0; k < fvcl.size(); k++) { // loop on clusters

            int chamber_id = clusterModule(k);
            float time_start = fvcl[k].tstart;
            float time_delta = fvcl[k].tdstart;
            float time_end = fvcl[k].tstop;
            // float time_centroid;
            float charge_total = fvcl[k].totcharge;
            float strip_centroid = fvcl[k].centroid;
            float strip_rms = fvcl[k].rms;
            float xmt = fvcl[k].m;
            float xqt = fvcl[k].q;
            int nstrip = numSignals(k);

            if (isinf(time_start) || isinf(time_end))
                continue; // do not print clusters with inf **TO BE DETECTED during process**
            sret = sret + Form("%d %f %d %d %.3f %.3f %.3f %.2f %.3f %.3f %.3f %.3f\n", evt, time_from_start, chamber_id, nstrip, time_start, time_delta, time_end, charge_total, strip_centroid, strip_rms, xmt, xqt);
        }

        return sret;
    };

    /*
     * Search for time and space clusters (consecutive samples and adjacent strips)
     * at least "thrsample" and no more than "maxsample" consecutive number of samples - TO BE IMPROVED
     *
     */

    int searchClustersX(int thrsample = 8, int maxsample = 20) { // search space-time clusters, shall have thrsample consecutive samples with charge>0 and cannot exceed maxSample contiguous samples

        // int idxClu = 0; // running cluster index
        // int count = 0;

        uTSIGNAL tsd;
        if (fverbose > 1) {
            printf(" Search clusters (time and spatial) signals\n");
        }

        std::vector<float> vcharge;

        // search space-time cluster
        fchargemax = 0;
        for (int iu = 0; iu < nurw; iu++) { // loop on chambers

            for (int ic = 0; ic < ncha; ic++) { // loop on channels (strips)
                if (maskCh[iu * ncha + ic] > 0)
                    continue;        // masked channel
                int countsample = 0; // consecutive samples
                float totcharge = 0;
                float centroid = 0;
                float rms = 0;
                // search time developing signals
                for (int is = 0; is <= nsam; is++) { // loop on samples (time); extended by 1 unit beyond last sample to consider truncated signals
                    double charge = 0;
                    if (is < nsam) {
                        int ij = getIdx(iu, ic, is);
                        charge = qts[ij];
                    }

                    if (charge > 0) {
                        fchargemax = (charge > fchargemax) ? charge : fchargemax;
                        vcharge.push_back(charge);
                        centroid += charge * is;
                        rms += charge * is * is;
                        totcharge += charge;
                        countsample += 1;
                    } else {
                        if ((countsample > thrsample) && (countsample < maxsample)) { // new time signal
                            tsd.imod = iu;
                            tsd.istrip = ic;
                            tsd.t0 = is - countsample; //+1;
                            tsd.nt = countsample;
                            tsd.qt = totcharge;
                            tsd.q = vcharge;
                            tsd.centroid = centroid / totcharge;
                            tsd.rms = TMath::Sqrt(rms / totcharge - tsd.centroid * tsd.centroid);
                            fvts.push_back(tsd);
                            if (fverbose > 2)
                                printf(" New time signal: %d %d : %d %d : %.1f [ %.2f +/- %.2f ]\n", tsd.imod, tsd.istrip, tsd.t0, tsd.nt, tsd.qt, tsd.centroid, tsd.rms);
                        }
                        vcharge.clear();
                        centroid = 0;
                        rms = 0;
                        totcharge = 0;
                        countsample = 0; // reset sample count
                    }
                } // end sample loop
            } // end strip loop
        } // end chamber loop

        // now combine time signals from adjacent strips and congruent timing, taking into account masked channels
        // assume ordered strips
        int im = -1;
        int is = -1;
        int it = -1;
        int dt = 0;

        int notadj = 1; // not adjacent
        uCLUSTER ucl;
        std::vector<int> vdd;
        for (int ii = 0; ii < (int)fvts.size(); ii++) { // loop on time signals

            int im1 = fvts[ii].imod;
            int is1 = fvts[ii].istrip;
            int it1 = fvts[ii].t0;
            int dt1 = fvts[ii].nt;

            if (im1 == im) {                                                                     // same chamber
                if (TMath::Abs(is1 - is) <= 2) {                                                 // signal in adjacent strips (+/2 strips) @@ NEED TO ADD CHECK ON MASKED CHANNELS !!! + maskCh[im*ncha+is+1]
                    if (((it1 >= it) && (it1 < it + dt)) || ((it >= it1) && (it < it1 + dt1))) { // signals overlap in time
                        notadj = 0;
                    }
                }
            }

            if ((notadj) && (ii > 0)) { // not adjacent means new cluster (or first index in loop)
                ucl.its = vdd;
                fvcl.push_back(ucl);
                vdd.clear();
            }

            im = im1;
            is = is1;
            it = it1;
            dt = dt1;
            vdd.push_back(ii);
            notadj = 1;

        } // end loop on time signals

        if (vdd.size() > 0) {
            ucl.its = vdd;
            fvcl.push_back(ucl);
        } // last cluster

        clusterAnalysis();
        /*
        for (int i=0;i<fvcl.size();i++) { // loop on clusters
          printf(" Cluster %d:",i);
          for (int j=0;j<fvcl[i].size();j++) { // loop on signals of the cluster
            int isig = fvcl[i][j];
            fvts[isig].clidx = i; // cluster index the signal will belong to
            printf(" %d",isig); // corresponding signal index
          }
          printf("\n");
        }
        */

        /*
        time - centroid and rms of each strip
        m,q of linear fit of centroid/rms vs strip
        */

        return 0;
    };
};

/****
 *
 *
 */

int alignmentOptimization(tsCluster *tc) {

    // static float ns1 = 0;
    // static float ds1 = 0;

    int nt = tc->numTracks(); // current event tracks

    if (nt != 1) {
        return 0;
    } // consider "clean events" only

    // float s1 = 0;

    return 0;
};

/*****
 * Manage output of data into a file (for the moment is a text file)
 */

class dataOut {

  private:
    FILE *fostream;
    TString foname;

  public:
    dataOut(TString outfname, TString info) {
        foname = outfname;
        fostream = fopen(foname.Data(), "w+");
        if (fostream == NULL) {
            printf("ERROR: Cannot open %s file\n", foname.Data());
            exit(-1);
        } else {
            fprintf(fostream, "%s", info.Data());
            printf(" File %s opened for writing\n", foname.Data());
        }
    };

    ~dataOut() {
        if (fostream) {
            fclose(fostream);
            printf(" Output text file %s closed\n", foname.Data());
        }
    };

    void comment(TString info) {
        fprintf(fostream, "%s", info.Data());
        fflush(fostream);
    };

    /*
     * type: the provided output data
     */
    void event(int evt, float time, tsCluster *ts, int type = 0) {
        TString sdum;
        switch (type) {
        case 0:
            sdum = ts->eventOutput(evt, time);
            break;
        case 1:
            sdum = ts->trackOutput(evt, time);
            break;
        default:
            return;
        }
        fprintf(fostream, "%s", sdum.Data());
        fflush(fostream);
    };
};

#endif
