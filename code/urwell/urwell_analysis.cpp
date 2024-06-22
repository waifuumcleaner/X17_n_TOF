/*
 * uRwell SRS APV root file analysis
 *
 * first version: May/2022
 *
 */

#include <limits.h>
#include <string>
#include <time.h>

#include <Riostream.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TDatime.h>
#include <TEntryList.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h>
#include <TLinearFitter.h>
#include <TMarker.h>
#include <TMath.h>
#include <TPolyMarker.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TText.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TVirtualFFT.h>

#include "mmRawTree.h"

#include "../cube/ntofDAQ.h"
#include "../utils.h"
#include "eventDisplay.h"
#include "urwell_analysis.h"
#include "utilf.h"

// FROM SHELL: g++ urwana.cpp  -Wall -Wextra -g -O3 `root-config --cflags --libs`
//              ./a.out
// FROM ROOT : .L urwana.cpp+
//             decodePhysRun(29, 27, 5, 6, 18)

// from apv and ch (electronics) indeces to absolute channel index (or strip index)
int e2a(int iapv, int irch) { return iapv * 128 + irch; }

// from absolute channel to apv index
int a2apv(int iach) { return (int)iach / 128; }

// from absolute channel to apv channel index
int a2ch(int iach) { return (iach % 128); }

// from absolute channel to chamber index
int a2c(int iach) { return (int)iach / 512; }

// from absolute channel to strip index in chamber
int a2s(int iach) { return (iach % 512); }

// from chamber and strip to absolute channel
int c2a(int icha, int istrip) { return (icha * 512 + istrip); }

PEDES getPedestal(mmRawTree *rawT, // used to get the module names
                  TString infile, TFile *fIn, float hothr) {

    // TText *testo = new TText();

    //  const char cname[2][10]={"BOT","TOP"};
    // const char cname[2][10]={"uRwA","uRwB"};
    //  const char cname[5][10]={"B0x","B1x","C0x","D3x","D3y"};

    if (fIn == 0) {
        fIn = new TFile(infile.Data(), "READ");
        if (fIn->IsZombie() != 0) {
            printf("Root pedestal run file %s not found or broken\n", infile.Data());
            exit(0);
        }
    }
    fIn->ls();

    printf("###### Read Pedestal file %s\n", infile.Data());

    TTree *fTp = (TTree *)fIn->Get("pedestals");

    int nchamb = rawT->numMods();

    TCanvas *cc0 = new TCanvas("cc0", "Chamber - Pedestal Mean and RMS");
    cc0->Divide(nchamb, 2);
    cc0->Update();

    TCanvas *cc1 = new TCanvas("cc1", "APV Card - Pedestal Mean");
    cc1->Divide(n_apv / 2, 2);
    cc1->Update();

    TCanvas *cc2 = new TCanvas("cc2", "APV Card - Pedestal RMS");
    cc2->Divide(n_apv / 2, 2);
    cc2->Update();

    for (int i = 0; i < nchamb; i++) { // loop on chambers
        TString modid = rawT->getModName(i);
        unsigned int axis = rawT->getAxis(i);
        printf(" pedestal process %s %d\n", modid.Data(), axis);

        cc0->cd(i + 1);
        fTp->SetMarkerStyle(20 + i);
        fTp->SetMarkerColor(1 + i);
        fTp->Draw("apv_pedmean:mm_strip", Form("(mm_id==\"%s\")&&(mm_readout==%d)", modid.Data(), axis), "p");

        // cc1->cd(i+1);
        cc0->cd(i + 1 + nchamb);
        fTp->SetMarkerStyle(20 + i);
        fTp->SetMarkerColor(1 + i);
        fTp->Draw("apv_pedstd:mm_strip", Form("(mm_id==\"%s\")&&(mm_readout==%d)", modid.Data(), axis), "p");
        //    fTp->SetMarkerStyle(24);
        //    fTp->Draw("apv_pedstd:mm_strip", Form("mm_id==\"%s\"", cname[i]),"p,same");
        //    cc1->Update();
        cc0->Update();
    }

    PEDES vped; // apv_id, apv_ch
    printf(" Hot channels from pedestal data (APV, ch):\n");
    for (int i = 0; i < n_apv; i++) { // APV cards
        cc1->cd(i + 1);
        fTp->SetMarkerStyle(1);
        fTp->SetMarkerColor(0);
        fTp->Draw("apv_pedmean:apv_ch", "", "p");
        fTp->SetMarkerStyle(20 + i / 2);
        fTp->SetMarkerColor(1 + i / 2);
        fTp->Draw("apv_pedmean:apv_ch", Form("apv_id==%d", i), "p,same");
        for (int k = 0; k < fTp->GetSelectedRows(); k++) {
            vped.mean[i].push_back(fTp->GetV1()[k]);
        }
        cc1->Update();

        cc2->cd(i + 1);
        fTp->SetMarkerStyle(1);
        fTp->Draw("apv_pedstd:apv_ch", "", "p");

        fTp->SetMarkerStyle(20 + i / 2);
        fTp->SetMarkerColor(1 + i / 2);
        fTp->Draw("apv_pedstd:apv_ch", Form("apv_id==%d", i), "p,same");
        for (int k = 0; k < fTp->GetSelectedRows(); k++) {
            vped.sigma[i].push_back(fTp->GetV1()[k]);
        }
        fTp->Draw("apv_pedstd:apv_ch", Form("(apv_id==%d)&&(apv_pedstd>%f)", i, hothr), "p,same"); // hot channels
        for (int k = 0; k < fTp->GetSelectedRows(); k++) {
            int idxch = (int)fTp->GetV2()[k];
            vped.sigma[i][idxch] = -vped.sigma[i][idxch]; // masked channels set to their negative values
            printf("   %d %d (ped: %f)\n", i, idxch, vped.sigma[i][idxch]);
        }
        cc2->Update();
    }

    cc2->cd(8);
    /*
    testo->SetTextColor(1);
    testo->DrawTextNDC(.5,.7,Form("%s : black", cname[0]));
    testo->SetTextColor(2);
    testo->DrawTextNDC(.5,.8,Form("%s : red", cname[1]));
    */
    cc2->Update();

    //  fTp->Draw("mm_strip","apv_id==7");

    return vped;
}

std::vector<float> pedestalSubtract(std::vector<short> adc, int iapv, int ich, PEDES ped, std::vector<float> cnoise, float nsig, TProfile *hpo) {

    int nsamples = adc.size();

    std::vector<float> vsig;

    float mped = ped.mean[iapv][ich];
    float sped = ped.sigma[iapv][ich];

    if (sped <= 0) {
        printf("WARNING: %d %d sigma pedestal wrong : %f\n", iapv, ich, sped);
    }

    /*
      int saold=-1; // previous over thr sample
      int adjsize=1;
      int adj0=-1;
      float adjqp=0;
      float adjcc=0;
    */

    // float mcharge=0;
    // float scharge=0;

    for (int is = 0; is < nsamples; is++) {

        float charge = mped - (float)adc[is] - cnoise[is]; // raw signals are negative (relative to baseline) -> now changed to positive

        if (hpo)
            hpo->Fill(ich + is * 128, charge);

        //    mcharge += charge;
        //    scharge += charge*charge;

        charge = (charge > nsig * sped) ? charge : 0; // suppress noise

        vsig.push_back(charge);
    }

    return vsig;
};

std::vector<TCLUSTER> samplesAnalysis(std::vector<float> adc, // signal, noise suppressed
                                      int iapv, int ich) {

    int nsamples = adc.size();

    std::vector<TCLUSTER> vtc; // returned vector of next tc elements
    TCLUSTER tc;               // working cluster structure

    tc.ch = e2a(iapv, ich);

    int saold = -1; // previous over thr sample
    int adjsize = 1;
    int adj0 = -1;   // first sample
    float adjqp = 0; // total charge
    float adjcc = 0; // sample centroid
    // int adjmax = -1; // sample with max

    for (int is = 0; is < nsamples; is++) {

        float charge = adc[is];

        if (charge > 0) { // signal
            if (adj0 < 0)
                adj0 = is;
            if ((saold + 1) == is) { // next time sample, append to previous cluster
                adjsize += 1;
                adjqp += charge;
                adjcc += charge * ((float)is);
            } else { // sample is not consecutive, save previous time-cluster if any and "open" new
                if (adjqp > 0) {
                    tc.qt = adjqp;
                    tc.ct = adjcc / adjqp;
                    tc.t0 = adj0;
                    tc.nt = adjsize;
                    vtc.push_back(tc);
                    adjqp = 0;
                }
                adjsize = 1; // new time-cluster
                adjqp = charge;
                adjcc = charge * ((float)is);
                adj0 = is;
            }
            saold = is;
        }

        if (is == (nsamples - 1)) { // last sample (probably can be moved out of the loop)
            if (adjqp > 0) {        // need to save last group
                tc.qt = adjqp;
                tc.ct = adjcc / adjqp;
                tc.t0 = adj0;
                tc.nt = adjsize;
                vtc.push_back(tc);
            }
        }
    }

    /*
      printf(" VTC size %d : ", vtc.size());
      for (int i=0;i<vtc.size();i++) {
      printf(" %f",vtc[i].qt);
      }
      printf("\n");
    */

    return vtc;
};

SpaceCluster::SpaceCluster(float sigmat, int nsample, int maxmask) : sigma(sigmat), thr_nsample(nsample), max_adj_mask(maxmask){};

SpaceCluster::~SpaceCluster(){};

int SpaceCluster::searchSpatialCluster(std::vector<TCLUSTER> &tc, int deltas) {

    uint ntc = tc.size();

    uint ncl = cluster.size(); // cluster size may change in the next lines, but we are interested in previous clusters only, not on current strip (... confusing)

    if ((deltas == 0) || (deltas > max_adj_mask))
        ncl = 0; // force new chamber, previous clusters are "closed"

    for (uint i = 0; i < ntc; i++) { // loop on time cluster in strip

        if (tc[i].nt < thr_nsample) {
            continue;
        } // skip time clusters with number of consecutive samples < thr_nsample
        float ctime = (tc[i].ct + tc[i].t0) / 2; // assumed as "signal start time"

        int aflag = 0;

        printf("time sample - %.1f\n", ctime);
        for (uint j = 0; j < ncl; j++) { // loop on already defined spatial clusters

            printf("    %d : %f\n", j, last_strip[j]);
            if (                                              // (tc[i].nt >= thr_nsample) &&   // current time cluster has thr_nsample at least
                ((last_strip[j] + deltas) == tc[i].ch) &&     // belongs to adjacent strip
                (TMath::Abs(ctime - last_sample[j]) < sigma)) // is time correlated (in TPC should be a loose constraint)
            {

                printf(" ... add\n");
                cluster[j].qs += tc[i].qt;
                cluster[j].cs += ((float)a2s(tc[i].ch)) * tc[i].qt;
                cluster[j].ns += 1;
                last_sample[j] = ctime;
                last_strip[j] = tc[i].ch;

                sum_xi[j] += ctime; // used for linear regression estimation of projection speed
                sum_xi2[j] += (ctime * ctime);
                sum_yi[j] += (float)a2s(tc[i].ch);
                sum_xiyi[j] += ctime * ((float)a2s(tc[i].ch));
                sum_zi[j] += tc[i].qt;
                sum_xizi[j] += ctime * tc[i].qt;

                tc[i].lsc = (Int_t)j; // link to spatial cluster

                aflag = 1; // tc[i] linked to existing cluster

                break; // TBC - in principle more than 1 time-cluster could associate to the next strip time-cluster
            }
        } // end loop on already defined clusters

        if (aflag == 0) { // new spatial cluster
            SCLUSTER dummy;

            dummy.ur = a2c(tc[i].ch);
            dummy.qs = tc[i].qt;
            dummy.cs = ((float)a2s(tc[i].ch)) * tc[i].qt;
            dummy.ns = 1;
            dummy.ms = 0;
            dummy.mq = 0;
            //	printf("%d :charge : %e %d\n",i,dummy.qs, cluster.size());

            cluster.push_back(dummy);

            last_sample.push_back(ctime);
            last_strip.push_back(tc[i].ch);

            sum_xi.push_back(ctime); // used for linear regression estimation of projection speed
            sum_xi2.push_back(ctime * ctime);
            sum_yi.push_back((float)a2s(tc[i].ch));
            sum_xiyi.push_back(ctime * ((float)a2s(tc[i].ct)));
            sum_zi.push_back((tc[i].qt));
            sum_xizi.push_back(ctime * tc[i].qt);

            tc[i].lsc = (Int_t)(cluster.size() - 1); // link to spatial cluster
        }
    }

    return 0;
};

void SpaceCluster::finalizeClustering() {

    uint n = cluster.size();

    int inq = 0;
    int ims = 0;
    for (uint i = 0; i < n; i++) { // loop on clusters

        float q = cluster[i].qs;
        if (q <= 0) {
            inq++;
            continue;
        }
        cluster[i].cs = cluster[i].cs / q;

        float nn = (float)cluster[i].ns;
        printf(" cluster %d : charge %.1f, strip %f space %.1f\n", i, q, nn, cluster[i].cs);
        float den = (nn * sum_xi2[i] - sum_xi[i] * sum_xi[i]);
        if (den == 0) {
            if (nn > 1)
                ims++;
        }
        cluster[i].ms = TMath::ATan2((nn * sum_xiyi[i] - sum_xi[i] * sum_yi[i]), den);
        cluster[i].mq = TMath::ATan2((nn * sum_xizi[i] - sum_xi[i] * sum_zi[i]), den);
    }

    printf("finalizeClustering: charge0 %d, den0 %d (%d) of %d\n", inq, ims, inq + ims, n);
};

uint SpaceCluster::size() { return cluster.size(); };

SCLUSTER SpaceCluster::getCluster(uint i) { return cluster[i]; };

void SpaceCluster::clear() {
    cluster.clear();
    last_strip.clear();
    last_sample.clear();

    sum_xi.clear();
    sum_xi2.clear();
    sum_xiyi.clear();
    sum_yi.clear();
    sum_xizi.clear();
    sum_zi.clear();
};

std::vector<std::vector<float>> commonNoise(mmRawTree *tt, PEDES ped, int method) {

    int napv = tt->numAPVs();
    int nsample = tt->numSamples();

    Long64_t nn = tt->apv_id->size();

    std::vector<float> vdum;

    // assume max n_sample samples per channel
    float **vcharge = new float *[napv * nsample]; // apv, sample, channel
    int *ncount = new int[napv];                   // count unmasked channels on each APV

    std::vector<std::vector<float>> cono; // returned common noise vector [apv idx][sample idx]

    for (int i = 0; i < napv; i++) {
        ncount[i] = 0;
    }

    for (int i = 0; i < napv * nsample; i++) {
        vcharge[i] = new float[128]; // single APV charges at fixed sample
        for (int j = 0; j < 128; j++) {
            vcharge[i][j] = 0;
        }
    }

    int maxsample = 0;                  // retrieve charge from each channel and sample
    for (Long64_t k = 0; k < nn; k++) { // loop on cards and channels

        std::vector<unsigned int> *vid = tt->apv_id;
        std::vector<unsigned int> *vch = tt->apv_ch;

        unsigned int jid = vid->at(k); // apv
        unsigned int jch = vch->at(k); // channel

        if (ped.sigma[jid][jch] < 0)
            continue; // skip hot channels

        int icnt = ncount[jid];

        std::vector<short> vsq = tt->apv_q->at(k); // charge of given channels, all samples

        int nsamp = (int)vsq.size();
        maxsample = (maxsample > nsamp) ? maxsample : nsamp;

        for (int h = 0; h < nsamp; h++) { // loop on samples
            vcharge[(jid * nsample) + h][icnt] = ped.mean[jid][jch] - (float)vsq[h];
        }
        ncount[jid] += 1;
    }

    // computer common noise for each APV card
    for (int k = 0; k < napv; k++) {
        std::vector<float> dum;
        for (int h = 0; h < maxsample; h++) {

            float val = 0;
            switch (method) {
            case 0:
                val = TMath::Median(ncount[k], vcharge[k * nsample + h]);
                break;
            }
            dum.push_back(val);
        }

        cono.push_back(dum);
        dum.clear();
    }

    delete[] vcharge;
    delete[] ncount;

    return cono;
};

/*
 * ############################################
 * MAIN METHODS
 * ############################################
 */

int decodePhysRun(int phys_run, int ped_run, float nsigma, int min_ncsample, int max_ncsample, int verbose, int outType, TString dpath) {

    TString infile = Form("%s/run%d.root", dpath.Data(), phys_run); // mmdaq root file in physics mode
    TString ofile = Form("%s/run%d_ana.root", dpath.Data(), phys_run);

    TString otxtfile;
    switch (outType) {
    case 0:
        otxtfile = Form("%s/run%d_cluster.txt", dpath.Data(), phys_run);
        break;
    case 1:
        otxtfile = Form("%s/run%d_track.txt", dpath.Data(), phys_run);
        break;
    default:
        printf(" no output data type selected, use cluster data\n");
        otxtfile = Form("%s/run%d_cluster.txt", dpath.Data(), phys_run);
    }

    // TText *testo = new TText();
    //  ISS  const char cname[2][10]={"BOT","TOP"};
    //  const char cname[2][10]={"uRwA","uRwB"};

    // open raw file and load ttree
    TFile *fIn = new TFile(infile.Data(), "READ");
    if (fIn->IsZombie() != 0) {
        printf("Root physics run file %s not found or broken\n", infile.Data());
        exit(0);
    }

    fIn->ls();

    printf("Pedestal processed\n");

    TTree *tree;
    fIn->GetObject("raw", tree);

    tree->SetEstimate(-1); // all entries will be considered to estimate variable limits

    // link tree to proper object
    mmRawTree *rT = new mmRawTree(tree);

    long int nentries = rT->fChain->GetEntries();

    printf("Total Entries (Events) : %ld\n", nentries);
    // rT->Show(0);

    rT->BuildMap(10); // map apv_id,apv_ch to chamber,strip (and retrieve vector of module names)

    // process pedestals either from raw data file or as separate file
    PEDES pdata; // apv_id, apv_ch -> pedestal
    if (ped_run >= 0) {
        TString pedfile = Form("%s/run%d.root", dpath.Data(), ped_run);
        pdata = getPedestal(rT, pedfile, 0);
    } else {
        pdata = getPedestal(rT, "null", fIn);
    }

    // prepare output trees
    // STRIPS sso; // why is this not used??
    TCLUSTER tcln;
    SCLUSTER scln;
    EVENTINFO evto;

    //  TString ofile = Form("%s_ana.root",infile.Data());
    TFile *fOut = new TFile(ofile.Data(), "RECREATE");
    printf("Output data on file: %s\n", ofile.Data());

    TTree *tout = new TTree("tana", Form("Time Samples from run %d (ped: %d)", phys_run, ped_run));
    tout->Branch("event", &evto, "evt/I:time/F");
    tout->Branch("tclu", &tcln, "ch/I:nt/I:t0/I:qt/F:ct/F:lsc/I");

    TTree *tous = new TTree("sana", Form("Space Clusters from run %d (ped: %d)", phys_run, ped_run));
    tous->Branch("event", &evto, "evt/I:time/F");
    tous->Branch("sclu", &scln, "ur/I:ns/I:qs/F:cs/F:ms/F:mq/F");

    // new ? (June/2023) output format To be improved
    dataOut *dOut = new dataOut(otxtfile, Form("# Command = decodePhysRun(%d %d %f %d %d %d %d %s)\n", phys_run, ped_run, nsigma, min_ncsample, max_ncsample, verbose, outType, dpath.Data()));
    dOut->comment(Form("# Output data type = %d\n# Signal run = %d\n# Pedestal run = %d\n", outType, phys_run, ped_run));

    printf(" Allocate event display\n");
    eventDisplay *eD = new eventDisplay(phys_run, rT, verbose);

    int lNAPV = rT->numAPVs();
    int lNSAMPLE = rT->numSamples();

    printf("Allocate histograms for %d APVs\n", lNAPV);

    // allocate histograms for processing
    TProfile *hpcono[lNAPV]; // common noise profiles
    // TProfile *hposc[lNAPV];   // single event "oscilloscope"
    for (int i = 0; i < lNAPV; i++) {
        printf("new tprofile %d\n", i);
        hpcono[i] = new TProfile(Form("hpc%d", i), Form("Common Noise APV %d", i), lNSAMPLE, -0.5, 29.5);
        // hposc[i] = new TProfile(Form("hpox%d",i),Form("Charge vs Strip+TimeSample APV %d",i),128*n_sample,-0.5,128*n_sample-0.5);
    }

    printf("Histo et other variables allocated ...\n");

    // preliminary plots and params
    /*
      TCanvas *cc0 = new TCanvas("cc0","Raw Spectra");
      cc0->Divide(4,2);
      cc0->Update();

      TCanvas *cc0a = new TCanvas("cc0a","ADC vs Time");
      cc0a->Divide(4,2);
      cc0a->Update();

      TH2F *hcorr = new TH2F("hcorr","BOT vs TOP",512,-0.5,511.5, 512, -0.5,511.5);

      double srmean[8];   // mean charge
      double srrms[8];    // charge RMS
      for (int i=0;i<8;i++) {
      cc0->cd(i+1)->SetLogy();
      tree->Draw("apv_q",Form("apv_id==%d",i));
      printf(" elements : %lld\n",tree->GetSelectedRows());
      //    srmean[i] = TMath::Mean(tree->GetSelectedRows(),tree->GetV1());
      //    srrms[i] = TMath::RMS(tree->GetSelectedRows(),tree->GetV1());
      cc0a->cd(i+1);
      tree->Draw("apv_q:apv_evt",Form("apv_id==%d",i),"colz");
      break; // to-be-removed, for speed up
      }
    */

    // int deltastrip = 1; // used to manage masked channel for spatial cluster identification //not actually used??

    long int start_time_s = -1;
    TString string_stime;

    tsCluster *tsC = new tsCluster(rT->numMods(), 512, rT->numSamples(), verbose); // chambers, channels/chamber, samples/channel, verbosity

    // set masked channels
    for (int i = 0; i < rT->numAPVs(); i++) {
        for (int j = 0; j < 128; j++) {
            int imod = rT->getModule(i, j); // chamber index
            int istr = rT->getStrip(i, j);  // strip inde in chamber
            if (pdata.sigma[i][j] < 0) {    // masked channel
                tsC->maskChannel(imod, istr);
            }
        }
    }
    //  nentries = 10;
    // start processing, event by event
    for (Long64_t i = 0; i < nentries; i++) { // loop on entries
        Long64_t j = rT->LoadTree(i);
        if (j < 0) {
            break;
        }
        rT->GetEntry(i);

        Long64_t ecount = rT->apv_id->size();

        if (start_time_s < 0) {
            start_time_s = rT->time_s; // time in second since epoch
            time_t sss = (time_t)rT->time_s;
            string_stime = ctime(&sss);
            dOut->comment(Form("# Run Start Time (s) = %ld\n", start_time_s));
            dOut->event(-1, 0, tsC, outType); // write comment line with variable names that will be stored
        }

        evto.evt = rT->apv_evt;
        evto.time = ((float)(rT->time_s - start_time_s)) + ((float)rT->time_us) / 1e6; // time from run start in sec with us precision
        //    printf(" %lld : evt %d @ %ld s %d us: ecount %lld\n",i,rT->apv_evt,rT->time_s-start_time_s, rT->time_us, ecount);

        if ((verbose != 0) || ((evto.evt % 100) == 0))
            printf(" - - - - - - - - - - - - EVENT n. %d out of %ld - - - - - - - - - -\n", evto.evt - 1, nentries);
        std::vector<std::vector<float>> vcn = commonNoise(rT, pdata, 0);

        for (int k = 0; k < lNAPV; k++) {
            for (int h = 0; h < (int)vcn[k].size(); h++) {
                hpcono[k]->Fill(h, vcn[k][h]);
            }
        }

        tsC->reset();

        for (Long64_t k = 0; k < ecount; k++) { // loop on all strips (apv id and relative channels)

            std::vector<unsigned int> *vid = rT->apv_id; // @@@@@ put BEFORE start of loop ??? @@@@@
            std::vector<unsigned int> *vch = rT->apv_ch;

            unsigned int jid = vid->at(k);             // APV id
            unsigned int jch = vch->at(k);             // APV ch
            unsigned int jstrip = rT->mm_strip->at(k); // Strip

            std::vector<short> vsq = rT->apv_q->at(k); // sampled charges of single strip

            int jmod = rT->getModule(jid, jch); // chamber (module) index

            if (pdata.sigma[jid][jch] <= 0) { // channel is masked
                continue;
            }
            // subtract pedestal and common noise
            std::vector<float> vsig = pedestalSubtract(vsq, jid, jch, pdata, vcn[jid], nsigma, eD->getPOsc(jid));

            tsC->fillQ(jmod, jstrip, vsig);

        } // loop on all id,ch data

        tsC->searchClustersX(min_ncsample, max_ncsample);
        dOut->event(evto.evt, evto.time, tsC, outType);

        eD->plotx(i, tsC);
        if (eD->parserKB() < 0) { // quit
            delete eD;
            delete dOut;
            return -1;
            // ...
        }

        /*
        tsC->searchClusters();
        TCanvas *ce1=new TCanvas("ce1","Test");
        TH2F *hgf = tsC->plotCluster(1);
        ce1->Update();
        */

        // tous->Fill();

    } // end loop on entries

    printf("Loop on events completed\n");

    TCanvas *cc3 = new TCanvas("cc3", "Common noise");
    cc3->Divide(2, n_apv / 2);
    cc3->Update();
    for (int i = 0; i < lNAPV; i++) {
        cc3->cd(i + 1);
        hpcono[i]->SetMarkerStyle(20);
        hpcono[i]->DrawCopy("P");
    }

    cc3->Update();

    tout->Write();
    tous->Write();

    tout->Print();
    tous->Print();

    fOut->Write();
    fOut->Close();

    printf("Output data file %s closed\n", ofile.Data());

    delete dOut;

    return 0;
}

int readAna(int nrun, Float_t thr_q, Int_t thr_ns, Int_t thr_nt, TString dpath) {

    TString infile = Form("%s/run%d_ana.root", dpath.Data(), nrun);

    TFile *fIn = new TFile(infile.Data(), "READ");
    if (fIn->IsZombie() != 0) {
        printf("Root physics run file %s not found or broken\n", infile.Data());
        exit(0);
    }

    fIn->ls();

    TTree *tree;
    TTree *tres;
    fIn->GetObject("tana", tree);
    fIn->GetObject("sana", tres);

    TCanvas *ccp = new TCanvas("ccp", "Common Noise (events averaged)");
    ccp->Divide(2, n_apv);
    ccp->Update();

    TProfile *hpc[n_apv];
    for (int i = 0; i < n_apv; i++) {
        fIn->GetObject(Form("hpc%d", i), hpc[i]);
        ccp->cd(i + 1);
        hpc[i]->Draw();
    }
    ccp->Update();

    tree->Print();

    // sample cluster data
    TCanvas *cc0 = new TCanvas("cc0", "Sample/Time Clusters");
    cc0->Divide(2, 3);
    cc0->Update();

    TCut cnsam = Form("(nt>=%d)", thr_nt);
    TCut cnq = Form("(qt>=%f)", thr_q);

    tree->SetMarkerStyle(7);
    cc0->cd(1);
    tree->Draw("ct", cnsam && cnq);
    cc0->cd(2);
    tree->Draw("t0", cnsam && cnq);
    cc0->cd(3);
    //  tree->Draw("nt");
    tree->Draw("ch");
    cc0->cd(4);
    tree->Draw("ct-t0", cnsam && cnq);
    cc0->cd(5)->SetLogy();
    tree->Draw("qt", cnsam);
    cc0->cd(6);
    tree->Draw("qt:nt", cnsam && cnq, "colz");

    cc0->Update();

    tres->Print();

    TCut cnspa = Form("(ns>=%d)", thr_ns);
    TCut cnqs = Form("(qs>=%f)", thr_q);
    TCut curw; // chamber index cut

    TCanvas *ccs[2];
    for (int i = 0; i < 2; i++) { // loop on chambers

        curw = Form("ur==%d", i);
        ccs[i] = new TCanvas(Form("ccs%d", i), Form("Spatial Cluster, chamber %d", i));
        ccs[i]->Divide(2, 3);
        ccs[i]->Update();

        ccs[i]->cd(1)->SetLogy();
        tres->Draw("ns", curw);
        ccs[i]->cd(2)->SetLogy();
        tres->Draw("qs", cnspa && curw);
        ccs[i]->cd(3);
        tres->Draw("cs", cnspa && curw && cnqs);
        ccs[i]->cd(4)->SetLogy();
        tres->Draw("ms*TMath::RadToDeg()", cnspa && curw && cnqs);
        ccs[i]->Update();
        ccs[i]->cd(5);
        tres->Draw("ms*TMath::RadToDeg():cs", cnspa && curw && cnqs, "colz");
        ccs[i]->cd(6)->SetLogy();
        //    tres->Draw("ms*TMath::RadToDeg():mq*TMath::RadToDeg()", cnspa&&curw,"colz");
        tres->Draw("qs:ns", curw, "colz");
        ccs[i]->Update();
    }

    return 0;
};

void readClusterOut(int run, TString dpath) {

    int nchamber = 3;

    TString infile = Form("%s/run%d_cluster.txt", dpath.Data(), run);

    // get ntof root-ttree file if available
    TString ntofRFile = Form("%s/../cube/run_pkup_sall.root", dpath.Data()); // This is the summary file of the nTOF/CUBE scintillator data and related beam puk (optional)
    manageNTOF *ntof = new manageNTOF(ntofRFile);                            // internal tree pointer is null if file cannot be retrived (managed within the class)

    Long64_t sdati = getLPar("Run Start Time", 0, infile); // read start date-time in unix format
    ntof->setStartDaTi(sdati);

    TFile *tempfile = TFile::Open("tempout.root", "recreate");

    TTree *tclu = new TTree("tclu", "Cluster Post Analysis");
    tclu->ReadFile(infile, "Event/I:Time/F:Chamber/I:nStrips/I:tStart/F:tDeltaStart/F:tEnd/F:Charge/F:sCentroid/F:sRMS/F:mlFit/F:qlFit/F");

    tclu->Write(); // need to put the tree on a file

    TTreeReader reader("tclu", tempfile);

    TTreeReaderValue<int> bevt(reader, "Event");
    TTreeReaderValue<int> bchamb(reader, "Chamber");
    TTreeReaderValue<float> btime(reader, "Time");
    TTreeReaderValue<float> bcharge(reader, "Charge");

    int nentries = reader.GetEntries();
    printf(" Entries : %d\n", nentries);

    int count = 0;
    int ppe = nentries / 20; // progress period

    reader.SetEntry(count - 1);

    TH1F *h1 = new TH1F("h1", "Beam - SRS Time difference [s]", 100, -10., 10.);
    TGraph *gqp[nchamber]; // total charge per event in single chamber side
    float totcharge[nchamber];

    for (int i = 0; i < nchamber; i++) {
        totcharge[i] = 0;
        gqp[i] = new TGraph();
        gqp[i]->SetMarkerStyle(20 + i);
        gqp[i]->SetTitle(Form("Chamber %d", i));
    }

    int old_event = -1;
    int ientry = 0;
    float qbeam = 0;

    while (reader.Next()) {

        float time_us = (*btime) * 1e6; // convert to us
        int ievt = *bevt;

        if (ievt != old_event) { // new event

            std::vector<float> pulse = ntof->getBeamPulse(time_us); // get beam pulse charge
            qbeam = pulse[0] / 1e9;
            h1->Fill(pulse[1]);

            if (ientry > 0) { // at least one event already processed
                for (int i = 0; i < nchamber; i++) {
                    gqp[i]->SetPoint(ientry - 1, qbeam, totcharge[i]);
                    totcharge[i] = 0;
                }
            }
            ientry++;
            old_event = ievt;
        }

        float ucharge = *bcharge;
        int ich = *bchamb; // chamber index

        totcharge[ich] += ucharge;

        if ((count % ppe) == 0) {
            printf(" %6d (of %d): %f %f\n", count, nentries, time_us, ucharge);
        }

        count = count + 1;
        //    if (count>1000) break;
    }

    TCanvas *ccc = new TCanvas("ccc", "beam - srs");
    ccc->Divide(2, 2);
    ccc->Update();

    ccc->cd(1);
    h1->DrawCopy();

    for (int i = 0; i < nchamber; i++) {
        ccc->cd(i % nchamber + 2);
        gqp[i]->Draw("PAW");
    }
    ccc->Update();

    tempfile->Close();

    return;

    // ---- original stuff

    TCanvas *cb = new TCanvas("cb", "Charge et al");
    cb->Divide(3, 2);
    cb->Update();

    TCanvas *cc = new TCanvas("cc", "Timing");
    cc->Divide(3, 2);
    cc->Update();

    TCanvas *ce = new TCanvas("ce", "Fit");
    ce->Divide(3, 2);
    ce->Update();

    for (int i = 0; i < 3; i++) { // loop on chambers
        TCut ccut = Form("(Chamber==%d)", i);
        TCut cns = "(nStrips==1)";
        cb->cd(1 + i)->SetLogy();
        tclu->Draw("Charge", ccut);
        cb->cd(4 + i);
        tclu->Draw("sCentroid", ccut);

        cc->cd(1 + i);
        tclu->SetLineColor(1);
        tclu->Draw("tDeltaStart", ccut);
        cc->cd(4 + i);
        tclu->Draw("tEnd-tStart", ccut);
        tclu->SetLineColor(2);
        tclu->Draw("tEnd-tStart", ccut && cns, "same");

        ce->cd(1 + i);
        tclu->Draw("mlFit", ccut && (!cns));
        ce->cd(4 + i);
        tclu->Draw("qlFit", ccut && (!cns));
    }

    cb->Update();
    cc->Update();
    ce->Update();

    TCanvas *cf = new TCanvas("cf", "Comp");
    cf->Divide(2, 2);
    cf->Update();

    tclu->SetLineColor(1);
    cf->cd(1)->SetLogy();
    tclu->Draw("Charge");
    cf->cd(2)->SetLogy();
    tclu->Draw("tDeltaStart", "nStrips>1");
    cf->cd(3)->SetLogy();
    tclu->Draw("tEnd-tStart");
    cf->cd(4)->SetLogy();
    tclu->Draw("mlFit", "nStrips>1");

    tclu->SetLineWidth(2);
    for (int i = 0; i < 3; i++) {
        TCut ccut = Form("(Chamber==%d)", i);
        TCut cns = "(nStrips>1)";
        tclu->SetLineColor(i + 2);
        //    tclu->SetFillColor(i+2);
        //    tclu->SetFillStyle(3305+i*20);

        cf->cd(1)->SetLogy();
        tclu->Draw("Charge", ccut, "same");
        cf->cd(2)->SetLogy();
        tclu->Draw("tDeltaStart", ccut && cns, "same");
        cf->cd(3)->SetLogy();
        tclu->Draw("tEnd-tStart", ccut, "same");
        cf->cd(4)->SetLogy();
        tclu->Draw("mlFit", ccut && cns, "same");
    }
    cf->Update();
};

void readTrackOut(int run, TString dpath) {

    TString infile = Form("%s/run%d_track.txt", dpath.Data(), run);

    TTree *tclu = new TTree("tclu", "Track data Post Analysis");
    tclu->ReadFile(infile, "Event/I:Track/I:Chamber/I:tStart/F:tDeltaStart/F:tEnd/F:Charge/F:StripFirst/I:StripLast/I:StripCentroid/F:StripRMS/F:mlFit/F:qlFit/F");

    TCut cns = Form("StripLast>StripFirst");

    /*
    tclu->SetLineColor(18); //light gray
    cb->cd(1);
    tclu->Draw("tStart");
    cb->cd(2);
    tclu->Draw("tDeltaStart",cns);
    cb->cd(3);
    tclu->Draw("tEnd-tStart");

    tclu->SetLineWidth(2);
    for (int i=0;i<3;i++) { // loop on chamber
      tclu->SetLineColor(1+i);

      TCut cchamber = Form("Chamber==%d",i);
      cb->cd(1);
      tclu->Draw("tStart",cchamber,"same");
      cb->cd(2);
      tclu->Draw("tDeltaStart",cchamber&&cns,"same");
      cb->cd(3);
      tclu->Draw("tEnd-tStart",cchamber,"same");
    }
    */

    TCanvas *cc = new TCanvas("cc", "Space - Charge");
    cc->Divide(2, 2);
    cc->Update();

    TCanvas *cf = new TCanvas("cf", "Timing");
    cf->Divide(3, 2);
    cf->Update();

    tclu->SetLineColor(18);
    tclu->SetMarkerStyle(7);
    tclu->SetMarkerColor(18);

    cc->cd(1)->SetLogy();
    tclu->Draw("Charge");
    cc->cd(3);
    tclu->Draw("StripCentroid");
    cc->cd(2)->SetLogy();
    tclu->Draw("StripRMS");
    cc->cd(4);
    tclu->Draw("(StripLast-StripFirst):StripRMS", "", "p");

    TCut cmfit = "TMath::Abs(mlFit)<10000"; // ARBITRARY TO BE CHECKED
    TCut cqfit = "TMath::Abs(qlFit)<10000"; // ARBITRARY TO BE CHECKED

    cf->cd(1);
    tclu->Draw("Charge:(tEnd-tStart)", "", "p");
    cf->cd(2)->SetLogy();
    tclu->Draw("tDeltaStart", cns);
    cf->cd(5);
    tclu->Draw("(tEnd-tStart):tDeltaStart", cns, "p");
    cf->cd(4)->SetLogy();
    tclu->Draw("tEnd-tStart");
    cf->cd(3)->SetLogy();
    tclu->Draw("mlFit", cns && cmfit);
    cf->cd(6)->SetLogy();
    tclu->Draw("qlFit", cns && cqfit);

    tclu->SetLineWidth(2);
    for (int i = 0; i < 3; i++) {
        TCut ccut = Form("(Chamber==%d)", i);
        tclu->SetLineColor(i + 1);
        tclu->SetMarkerColor(i + 1);
        tclu->SetMarkerStyle(i + 20);

        //    tclu->SetFillColor(i+2);
        //    tclu->SetFillStyle(3305+i*20);

        cc->cd(1)->SetLogy();
        tclu->Draw("Charge", ccut, "same");
        cc->cd(3);
        tclu->Draw("StripCentroid", ccut, "same");
        cc->cd(2)->SetLogy();
        tclu->Draw("StripRMS", ccut, "same");
        cc->cd(4);
        tclu->Draw("(StripLast-StripFirst):StripRMS", ccut, "p,same");

        cf->cd(1);
        tclu->Draw("Charge:(tEnd-tStart)", ccut, "p,same");
        cf->cd(2)->SetLogy();
        tclu->Draw("tDeltaStart", ccut && cns, "same");
        cf->cd(5);
        tclu->Draw("(tEnd-tStart):tDeltaStart", cns && ccut, "p,same");
        cf->cd(3)->SetLogy();
        tclu->Draw("mlFit", ccut && cns, "same");

        cf->cd(4)->SetLogy();
        tclu->Draw("tEnd-tStart", ccut, "same");
        cf->cd(6)->SetLogy();
        tclu->Draw("qlFit", ccut && cns, "same");
    }
    cf->Update();
    cc->Update();
};

int main() {

    try {
        std::ofstream output_file("output.txt");
        (void)!freopen("output.txt", "w", stdout);
        decodePhysRun(29, 27, 5, 6, 18, 1, 0);
        output_file.close();
    } catch (const std::exception &e) {
        handle_exception(e);
        return 1;
    }
    return 0;
}
