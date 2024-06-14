/*
 * Event display class
 *
 */
#include "utilf.h"
#include <TCanvas.h>
#include <TCanvasImp.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TString.h>
#include <stdlib.h>

// #include "mmRawTree.h"

class eventDisplay {

  private:
    mmRawTree *frT; // used for mapping

    Int_t fNDH; // number of detector head
    Int_t fRun; // run number

    Int_t nfe;      // number of front end cards (APV)
    Int_t nchamber; // number of chamber (URWELL)
    Int_t nsample;  // number of samples

    Int_t flastEvt; // last event plotted

    std::vector<TProfile *> fPOsc; // single event "oscilloscope"

    std::vector<TH2F *> fHSvs;     // strip vs time-sample for each chamber
    std::vector<TH1F *> fHSample;  // dummy for time-samples of each chamber
    std::vector<TCanvas *> fCEdis; // canvases

    Int_t fEvt; // current event

    Int_t fDisable;  // disable event display
    Int_t fDeltaEvt; // event increment step (or absolute event if >1)
    Int_t fSavePlot; // save plot
    Int_t fFitMode;  // fit mode

    Int_t fFlag_Verbose;

    TString foutFilePref;

    Int_t fMarker_Size;
    Float_t fHisto_Max;
    Int_t fEvent_Direction;
    Int_t fFlag_Fill;
    Int_t fSample0;
    Int_t fDSample;
    Int_t fSolver_Method;
    Float_t fProj_Thr;  // charge projection threshold for printing channel info
    Float_t fCrossTalk; // set cross talk correction factor (if <0 do not apply croo talk correction)

    Float_t fnSigma; // number of sigma for noise suppression

    Float_t fxmin, fxmax; // if centroid in x range, event is shown
    Float_t fymin, fymax; // if centroid in y range, event is shown

    Int_t fNColor;    // number of color in custom palette
    Int_t *fCPalette; // custom color palette

    // from event data
    Float_t *fCharge_max;
    Float_t *fCharge_tot;
    Float_t *fNhit;
    Float_t *fChargeMaxThr; // charge max threshold (below this value event is not shown)

    Int_t crossIdx[128]; // cross talk channels
    Int_t c2s[128];      // apv channel to strip

  public:
    // napv: number of front-end cards (APV, assume 128 channel each)
    // nurw: number of chambers (URWELL...)
    // nts : number of time-samples
    // devt : event delta (if <= 0 no event display)

    eventDisplay(Int_t irun, mmRawTree *rT, Int_t devt = 1) {

        Int_t napv = rT->numAPVs();
        Int_t nurw = rT->numMods();
        Int_t nts = rT->numSamples();

        frT = rT;

        fEvent_Direction = 1;
        fSample0 = 0;
        fDSample = 12;
        fSolver_Method = 0;

        fFlag_Verbose = (devt > 0) ? 1 : 0;
        fFlag_Fill = 0;

        fProj_Thr = 1000.; // charge thr for verbose printing
        fCrossTalk = 0.0;  // default cross talk correction

        foutFilePref = "plot_event"; // outfile plot file prefix

        fDeltaEvt = (devt > 0) ? devt : 0;
        nfe = napv;
        nchamber = nurw;
        nsample = nts;

        int nstrip = 512; // TBC assume same number of strips in each chamber

        for (int i = 0; i < nfe; i++) {
            printf("eventdisplay: allocate hpo tprofile %d\n", i);
            fPOsc.push_back(new TProfile(Form("hpo%d", i), Form("Charge vs Strip+TimeSample APV %d", i), 128 * 30, -0.5, 128 * 30 - 0.5));
            fPOsc[i]->SetXTitle("sample+strip");
            fPOsc[i]->SetYTitle("charge");
        }

        for (int i = 0; i < nchamber; i++) {
            printf("eventdisplay: allocate hsvs tprofile %d\n", i);
            fHSvs.push_back(new TH2F(Form("hsvs%d", i), Form("Charge vs Strip and Time Sample Chamber %s", rT->getModName(i).Data()),
                                     nstrip, -0.5, nstrip - 0.5, nsample, -0.5, nsample - 0.5));
            fHSvs[i]->SetXTitle("strip");
            fHSvs[i]->SetYTitle("time-sample");
            fHSample.push_back(new TH1F(Form("hsam%d", i), Form("Time Samples Chamber %s", rT->getModName(i).Data()),
                                        10, -0.5, nsample - 0.5)); //, 10,-0.5,1000));
            fHSample[i]->SetXTitle("time sample [25ns]");
            fHSample[i]->SetYTitle("Charge (a.u.)");
        }

        printf("eventdisplay: allocate TCanvases for %d %d\n", nfe, nchamber);

        TCanvas *cC; // delete existing canvases
        for (int i = 0; i < 3; i++) {
            cC = (TCanvas *)gROOT->FindObject(Form("cedis%d", i));
            if (cC)
                delete cC;
        }

        fCEdis.push_back(new TCanvas("cedis0", Form("Oscilloscope Frames (Run %d)", irun)));
        fCEdis.push_back(new TCanvas("cedis1", Form("Charge vs Strips and Time (Run %d)", irun))); // 2d representation of cedis0
        fCEdis.push_back(new TCanvas("cedis2", Form("Charge vs Time (Run %d)", irun)));
        //    fCEdis.push_back( new TCanvas("cedis3",Form("Channels Samples Frames (Run %d)",irun))); // same info in cedis0, splitted on different pads

        fCEdis[0]->Divide(2, TMath::Ceil(((float)nfe) / 2.), .0001, .00005);
        fCEdis[0]->Update();

        fCEdis[1]->Divide(1, nchamber);
        fCEdis[1]->Update();

        fCEdis[2]->Divide(2, nchamber);
        fCEdis[2]->Update();

        /*
          fCEdis[3]->Divide(6, 5, -1, -1); // single chamber
          fCEdis[3]->Update();
        */

        // mapping!
        // cross talk mapping (see crosstalk.map) ... to be improved reading the map file

        for (int i = 0; i < 96; i++) {
            crossIdx[i] = 32 + i;
        }
        for (int i = 0; i < 24; i++) {
            crossIdx[i + 96] = i + 8;
        }
        for (int i = 0; i < 7; i++) {
            crossIdx[i + 120] = i + 1;
        }
        crossIdx[127] = -1; // no cross talk

        /*
          printf("Cross Talk Map (%f)\n",fCrossTalk);
          for (int ch=0;ch<128;ch++) {
          printf("%d : %d\n",ch,crossIdx[ch]);
          }
        */

        for (int ch = 0; ch < 8; ch++) { // channel to strip, APV mapping (To Be Doubled Checked)
            c2s[ch] = 2 * ch;
            c2s[120 + ch] = ch * 2 + 1;
            c2s[8 + ch] = 16 + ch * 2;
            c2s[96 + ch] = 16 + ch * 2 + 1;
            c2s[16 + ch] = 32 + ch * 2;
            c2s[104 + ch] = 32 + ch * 2 + 1;
            c2s[24 + ch] = 48 + ch * 2;
            c2s[112 + ch] = 48 + ch * 2 + 1;
            c2s[32 + ch] = 64 + ch * 2;
            c2s[64 + ch] = 64 + ch * 2 + 1;
            c2s[40 + ch] = 80 + ch * 2;
            c2s[72 + ch] = 80 + ch * 2 + 1;
            c2s[48 + ch] = 96 + ch * 2;
            c2s[80 + ch] = 96 + ch * 2 + 1;
            c2s[56 + ch] = 112 + ch * 2;
            c2s[88 + ch] = 112 + ch * 2 + 1;
        }
    };

    ~eventDisplay() {
        for (uint i = 0; i < fPOsc.size(); i++) {
            delete fPOsc[i];
        }

        for (uint i = 0; i < fHSvs.size(); i++) {
            delete fHSvs[i];
        }

        for (uint i = 0; i < fCEdis.size(); i++) {
            delete fCEdis[i];
        }
        //    if (fCPalette) delete fCPalette;
    };

    Int_t next() { // ??
        Int_t ret;
        if ((fDeltaEvt > 1) || (fDeltaEvt < -1) || (fDeltaEvt == 0)) {
            ret = TMath::Abs(fDeltaEvt);
            fDeltaEvt = 1;
        } else {
            ret = (fEvt + fDeltaEvt);
        }
        return ret;
    }; // get next event

    Int_t clearHisto() {
        for (uint i = 0; i < fPOsc.size(); i++) {
            fPOsc[i]->Reset();
            if (i < fHSvs.size()) {
                fHSvs[i]->Reset();
            }
        }

        return 0;
    };

    //
    void crossCorr(Float_t cfactor = 0.1) { // cross talk correction

        if (cfactor <= 0) {
            return;
        }

        for (uint i = 0; i < fPOsc.size(); i++) { // loop on single APV time-ch profile

            Int_t nnx = fPOsc[i]->GetXaxis()->GetNbins();
            Int_t nsam = nnx / 128; // number of samples

            TProfile *pdummy = (TProfile *)fPOsc[i]->Clone("hpdum"); // keep the original signals

            for (Int_t jch = 0; jch < 128; jch++) { // channel on APV

                Int_t ghch = crossIdx[jch]; // cross-talk channel
                //	Float_t ghcorr = crossCor[jch]/100.; // cross-talk correction
                if (ghch < 0) {
                    continue;
                }

                for (Int_t isa = 0; isa < nsam; isa += 1) { // sample

                    Float_t charge = pdummy->GetBinContent(1 + isa * 128 + jch);       // bin start from 1!!
                    Float_t ghostCharge = pdummy->GetBinContent(1 + isa * 128 + ghch); // bin start from 1!!
                    if ((charge <= 0) || (ghostCharge <= 0))
                        continue;
                    Float_t ghostCorr = ghostCharge - charge * cfactor;

                    fPOsc[i]->SetBinContent(1 + ghch + isa * 128, ghostCorr);
                    printf("%d %d : %f %f = %f %d\n", jch, isa, charge, ghostCharge, ghostCorr, ghch);
                }
            }
            delete pdummy;
        }

        return;
    };

    //
    Int_t plotx(Int_t evt, tsCluster *tc) {

        flastEvt = evt;

        if (fFlag_Verbose == 0)
            return 0;

        gStyle->SetOptStat(0);

        crossCorr(fCrossTalk);

        TH2F *hdum;

        fCEdis[0]->SetTitle(Form("Oscilloscope Frames (Event = %d)", evt));
        for (uint i = 0; i < fPOsc.size(); i++) { // loop on single APV
            fCEdis[0]->cd(i + 1);
            fPOsc[i]->Draw();

            /* */
            Int_t nnx = fPOsc[i]->GetXaxis()->GetNbins();
            Int_t nsam = nnx / 128;                                                // number of samples
            for (Int_t isa = 0; isa < nsam; isa += 1) {                            // sample
                for (Int_t jch = 0; jch < 128; jch++) {                            // channel on APV
                    Float_t charge = fPOsc[i]->GetBinContent(1 + isa * 128 + jch); // bin start from 1!!

                    Int_t kcha = frT->getModule(i, jch); // chamber

                    Int_t strip = frT->getStrip(i, jch); // strip

                    fHSvs[kcha]->Fill(strip, isa, charge);
                }
            }
        }
        fCEdis[0]->Update();

        float uppercharge = tc->chargeMax();

        TGraph *ghd = new TGraph(); // charge vs sample for each strip
        /*
          TGraph *gts = new TGraph();  // strips vs start time
          gts->SetMarkerStyle(34);
          gts->SetMarkerSize(2);
          gts->SetMarkerColor(1);
        */
        TF1 *fpol1 = new TF1("p1", "[0]+[1]*x", -100., 100);

        TF1 *f1 = tc->getSignalFitFun(0);

        int mcolor[22] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 17, 20, 28, 29, 32, 38, 40, 42, 45, 46, 47};
        f1->SetLineWidth(1);
        for (uint i = 0; i < fHSvs.size(); i++) { // loop on chambers (one plot per chamber)
            fCEdis[1]->cd(i + 1);
            fHSvs[i]->Draw("colz");

            // fill time-signal plots
            fCEdis[2]->cd(2 * i + 1);
            fHSample[i]->SetMaximum(uppercharge * 1.05);
            //      fHSample[i]->SetMinimum(0.);
            fHSample[i]->Draw();
            std::vector<uTSIGNAL> vs = tc->signalInMod(i);
            int ns = (int)vs.size();
            float ystrip[ns];
            float xstime[ns];
            int cluidx[ns];
            std::vector<int> icl;          // unique cluster index
            for (int k = 0; k < ns; k++) { // loop on signals in given module
                ystrip[k] = vs[k].istrip;
                xstime[k] = vs[k].tstart;
                int nsam = vs[k].nt;
                int t0 = (float)vs[k].t0;
                float vx[nsam];
                for (int is = 0; is < nsam; is++) {
                    vx[is] = (float)(t0 + is);
                }
                std::vector<float> vq = vs[k].q;
                float *vy = &vq[0];
                int clidx = vs[k].clidx; // cluster signal belong to
                cluidx[k] = clidx;
                int ff = 1;                                 // temp flag
                for (int h = 0; h < (int)icl.size(); h++) { // check cluster already stored
                    if (icl[h] == clidx) {
                        ff = 0;
                        break;
                    }
                }
                if (ff)
                    icl.push_back(clidx);

                ghd->SetMarkerStyle(k % 30 + 20);
                ghd->SetMarkerColor(mcolor[clidx % 22]);
                ghd->SetLineColor(mcolor[clidx % 22]);
                ghd->DrawGraph(nsam, vx, vy, "p");
                for (int jp = 0; jp < f1->GetNpar(); jp++) {
                    f1->SetParameter(jp, vs[k].fitpar[jp]);
                }
                f1->SetLineColor(mcolor[clidx % 22]);
                f1->DrawCopy("same");
            } // loop on signals
            fCEdis[2]->cd(2 * i + 2);
            float dy0 = TMath::MinElement(ns, ystrip);
            float dy1 = TMath::MaxElement(ns, ystrip);
            float dx0 = TMath::MinElement(ns, xstime);
            float dx1 = TMath::MaxElement(ns, xstime);
            float ddx = (dx1 - dx0) / 10.;
            float ddy = (dy1 - dy0) / 10.;
            hdum = new TH2F(Form("hdum%d", i), "Strip vs Start Time (> 1 strip)", 5, dx0 - ddx, dx1 + ddx, 5, dy0 - ddy, dy1 + ddy);
            hdum->SetXTitle("Start Time [x25 ns]");
            hdum->SetYTitle("Strip index");

            hdum->DrawCopy();
            delete hdum;
            //      gts->DrawGraph(ns,xstime,ystrip,"p");
            TMarker *tmark = new TMarker();
            tmark->SetMarkerSize(1.5);
            for (int k = 0; k < ns; k++) { // loop on single strip with sample fit data
                if (tc->numSignals(cluidx[k]) < 2) {
                    continue;
                } // plot only cluster with at least 2 strips
                tmark->SetMarkerStyle(k + 20);
                tmark->SetMarkerColor(mcolor[cluidx[k] % 22]);
                tmark->DrawMarker(xstime[k], ystrip[k]);
            }
            for (int h = 0; h < (int)icl.size(); h++) {
                float m = tc->clusterSlope(icl[h]);
                float q = tc->clusterConst(icl[h]);
                fpol1->SetParameter(0, q);
                fpol1->SetParameter(1, m);
                fpol1->SetLineColor(mcolor[icl[h] % 22]);
                fpol1->DrawCopy("same");
            }
            icl.clear();
        } // end loop on chambers
        fCEdis[1]->Update();
        fCEdis[2]->Update();

        //    delete gts;
        delete ghd;
        delete fpol1;

        return 0;
    };

    //
    //
    Int_t parserKB() { // return number of keys pressed (return itself is excluded), -1 if quit pressed

        UInt_t cm_w, cm_h;
        Int_t flag_next;
        TString str, str1;
        TCanvasImp *cimp;

        TString outFile = Form("%s_%d.pdf", foutFilePref.Data(), flastEvt);
        int nnn;
        std::vector<int> v1, v2, v3;
        if (fFlag_Verbose == 0) {
            return 0;
        }

        do {

            flag_next = 1;

            printf("Press [key] RETURN (key = h for help) : ");
            str.ReadLine(std::cin, kFALSE);
            printf("Pressed: %s (%d)\n", str.Data(), str.Length());

            Int_t cc = 0;
            Float_t mfac = 1;

            str.Append("_"); // append a dummy character to avoid index overflow, this character cannot be used as command
            while (cc < str.Length()) {
                switch (str[cc]) {
                case 'h': // h - help
                    printf(" Available commands ([sequence of keys] and RETURN):\n");
                    printf(" General:\n");
                    printf("  h  : this help\n");
                    printf("  q  : quit\n");
                    printf("  v  : verbosity on/off (%d)\n", fFlag_Verbose);
                    printf("  RET: next event\n");

                    printf(" Cuts (character followd by + or -)\n");
                    printf("  s : number of sigma for noise threshold (%4.1f)\n", fnSigma);

                    printf(" Plots and graphics:\n");
                    printf("  c : clear (reset) histogram with statistics, and move to next event\n");
                    printf("  p : histogram points (markers) larger or smaller by 25%% / 20%% respectively\n");
                    printf("  x : strip profiles (print channel for charge larger than %f)\n", fProj_Thr);
                    printf("  y : histogram maximum increase/decrease (%5.1f)\n", fHisto_Max);
                    printf("  w : write/save plts on file prefix %s\n", outFile.Data());
                    printf("  z : zoom canvas\n");
                    printf("  o : histogram zoom around max and charge threshold for projection\n");

                    printf(" Analysis:\n");
                    printf("  e : invert event direction (%d)\n", fEvent_Direction);
                    printf("  f : toogle fill-tree flag (%d)\n", fFlag_Fill);
                    printf("  i : first sample in analysis (inc/dec by 1) (%d)\n", fSample0);
                    printf("  l : number of samples used in analysis (inc/dec by 1) (%d)\n", fDSample);
                    printf("  m : solver method cycle (%d)\n", fSolver_Method);
                    printf("  t : cross talk correction factor (%f)\n", fCrossTalk);

                    printf(" NOTE: Many commands increse the respective parameter by 25%% or decrease it by 20%% if followed by the \"-\" sign\n");
                    printf("       for example: s-RET will decrease the number of sigma by 20%%\n");

                    flag_next = 0;
                    break;
                case 'e':
                    fEvent_Direction = -fEvent_Direction;
                    printf("Event Direction is now: %d\n", fEvent_Direction);
                    flag_next = 0;
                    break;
                case 'c':
                    clearHisto();
                    break;
                case 'f':
                    fFlag_Fill = 1 - fFlag_Fill;
                    printf("Fill-tree flga is now: %d (0=no fill)\n", fFlag_Fill);
                    flag_next = 0;
                    break;
                case 'i':
                    cc++;
                    fSample0 += (str[cc] == '-') ? -1 : 1;
                    fSample0 = fSample0 % nsample;
                    break;
                case 'l':
                    cc++;
                    fDSample += (str[cc] == '-') ? -1 : 1;
                    if (fDSample < 2) {
                        fDSample = 2;
                    };
                    fDSample = fDSample % nsample;
                    break;
                case 'm':
                    fSolver_Method = (fSolver_Method + 1) % 4;
                    printf("Solver method is now: %d\n", fSolver_Method);
                    break;
                case 'o':
                    printf("Provide first and last bins and intensity threshold: ");
                    str1.ReadLine(std::cin, kFALSE);
                    Int_t first, last;
                    sscanf(str1.Data(), "%d %d %f", &first, &last, &fProj_Thr);
                    for (uint ic = 0; ic < fHSvs.size(); ic++) {
                        fCEdis[1]->cd(1 + ic);
                        fHSvs[ic]->GetXaxis()->SetRange(first, last);
                        fHSvs[ic]->Draw("colz");
                    }
                    fCEdis[1]->Update();
                    flag_next = 0;
                    break;
                case 'p':
                    cc++;
                    fMarker_Size *= (str[cc] == '-') ? 0.80 : 1.25;
                    break;
                case 'q': // q - quit
                    return -1;
                    break;
                case 's':
                    cc++;
                    fnSigma += (str[cc] == '-') ? -0.5 : 0.5;
                    break;
                case 't':
                    printf("Provide Cross Talk Correction Factor (if <=0, no correction applied on next events): ");
                    str1.ReadLine(std::cin, kFALSE);
                    sscanf(str1.Data(), "%f", &fCrossTalk);
                    break;
                case 'x':
                    v1.clear();
                    v2.clear();
                    v3.clear();
                    for (uint ic = 0; ic < fHSvs.size(); ic++) { // chamber
                        fCEdis[1]->cd(1 + ic);
                        TH1D *h1d = (TH1D *)fHSvs[ic]->ProjectionX();
                        h1d->Draw();
                        h1d->SetXTitle("strip");
                        h1d->SetYTitle("time-cumulated charge");
                        for (int iv = 0; iv < h1d->GetNbinsX(); iv++) {
                            float y = h1d->GetBinContent(iv);
                            if (y > fProj_Thr) {
                                float x = h1d->GetBinCenter(iv);
                                v1.push_back((int)x);  // strip
                                v2.push_back((int)y);  // charge
                                v3.push_back((int)ic); // chamber
                            }
                        }
                    }
                    nnn = v1.size();
                    if (nnn > 0) {
                        printf("Chamber Strip Tot-Charge\n");
                        for (int ic = 0; ic < nnn; ic++) {
                            printf(" %d %d %d\n", v3[ic], v1[ic], v2[ic]);
                        }
                    }
                    fCEdis[1]->Update();
                    flag_next = 0;
                    break;
                case 'v':
                    fFlag_Verbose = 1 - fFlag_Verbose;
                    printf("Verbosity is now: %d (0 = off)\n", fFlag_Verbose);
                    flag_next = 0;
                    break;
                case 'y':
                    cc++;
                    fHisto_Max *= (str[cc] == '-') ? 0.80 : 1.25;
                    break;
                case 'w':
                    printf("Write plots of event %d on file %s\n", flastEvt, outFile.Data());
                    fCEdis[0]->Print(outFile + "(");
                    fCEdis[1]->Print(outFile);
                    fCEdis[2]->Print(outFile + ")");
                    flag_next = 0;
                    break;
                case 'z': // z - zoom canvas
                    cc++;
                    mfac = (str[cc] == '-') ? 0.80 : 1.25;
                    for (uint i = 0; i < fCEdis.size(); i++) {
                        cm_w = fCEdis[i]->GetWindowWidth();
                        cm_h = fCEdis[i]->GetWindowHeight();
                        printf("Zoom window w / h : %d %d by %f\n", cm_w, cm_h, mfac);
                        cm_h = cm_h * mfac;
                        cm_w = cm_w * mfac;
                        cimp = fCEdis[i]->GetCanvasImp();
                        cimp->SetWindowSize(cm_w, cm_h);
                        cimp->ForceUpdate();
                    }
                    flag_next = 0;
                    break;
                default:
                    break;
                }
                cc++;
            } // while on input string characters

        } while (flag_next == 0);

        return (str.Length() - 1); // return number of characters (return excluded)
    };

    Int_t getColor(Int_t n = 0) {
        Int_t i = 0;
        i = (n < fNColor) ? n : (fNColor - 1);
        if (i < 0)
            i = 0;
        return fCPalette[i];
    }

    Int_t getColor(Float_t z, Float_t zmax, Float_t z0 = 0.) {
        Int_t i;
        if ((zmax - z0) != 0) {
            i = (z - z0) / (zmax - z0) * fNColor;
        } else {
            i = fNColor;
        }
        return getColor(i);
    }

    Int_t isDisabled() { return fDisable; };

    TProfile *getPOsc(Int_t idx) { return fPOsc[idx]; }; // used to fill fPOsc

    TH2F *getHSvs(Int_t idx) { return fHSvs[idx]; };
};
