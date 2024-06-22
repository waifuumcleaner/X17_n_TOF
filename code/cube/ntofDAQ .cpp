/*
 * utility classes and methods to manage the ntof DAQ root files and correlate ntof data with fers and srs DAQ
 */

#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <iomanip>
#include <iostream>

#include "ntofDAQ.h"

TString getSPar(TString parname, Int_t idx, TString inFile, TString defval) {

    TString cmd;
    TString line;
    TString tok;
    Ssiz_t from = 0;
    Int_t match = 0;

    Int_t count = 0;

    TString retval = defval;

    cmd = Form("grep \"%s\" %s", parname.Data(), inFile.Data());

    // cmd = cmd + " | grep -v \"#\""; // line with #
    cmd = cmd + " | grep \"#\"";
    FILE *fp = gSystem->OpenPipe(cmd.Data(), "r");

    while (line.Gets(fp)) {
        count = 0;
        if (!line.BeginsWith("#"))
            continue;
        from = line.First("=");
        if (from < 0)
            continue;
        while (line.Tokenize(tok, from, "[ \t]+")) { // look at values after the "="
            //    printf("%d out: %s\n",count,tok.Data());
            if (count > 0) {
                retval = tok;
            }
            if (count == (idx + 1)) {
                match = 1;
                break;
            }
            count++;
        }
        if (match)
            break;
    }

    gSystem->ClosePipe(fp);

    if (match) {
        retval = tok;
    } else {
        printf("WARNING (get*Par): parameter %s ", parname.Data());
        if (count == 0) { // parameter is absent
            printf("does not exist, return default value %s\n", defval.Data());
            retval = defval;
        }
        if (count == 1) { // parameter is present but no value is associated
            printf("exists but no value is assigned, return default value %s\n", defval.Data());
            retval = defval;
        }
        if (count > 1) {
            printf("exists but assigned values do not extend to selected index %d, return last useful value %s\n", idx, retval.Data());
        }
    }
    printf("%s %d : %s\n", parname.Data(), idx, retval.Data());
    return retval;
}

Float_t getFPar(TString parname, Int_t idx, TString inFile, Float_t defval) {

    Float_t retval = getSPar(parname, idx, inFile, Form("%f", defval)).Atof();
    return retval;
}

Int_t getIPar(TString parname, Int_t idx, TString inFile, Int_t defval) {

    Int_t retval = getSPar(parname, idx, inFile, Form("%d", defval)).Atoi();
    return retval;
}

Long64_t getLPar(TString parname, Int_t idx, TString inFile, Long64_t defval) {

    Int_t retval = getSPar(parname, idx, inFile, Form("%lld", defval)).Atoll();
    return retval;
}

std::time_t convertTime(TString stime) { // GMT=UTC time

    // https://stackoverflow.com/questions/17681439/convert-string-time-to-unix-timestamp

    std::string s{stime.Data()};
    std::tm t{};
    std::istringstream ss(s);
    ss >> std::get_time(&t, "%a %b %d %H:%M:%S %Y"); // do not consider "UTC"
    if (ss.fail()) {
        throw std::runtime_error{"failed to parse time string"};
    }

    std::time_t time_stamp = mktime(&t); // convert to time since epoch
    std::cout << time_stamp << " current seconds since epoch" << std::endl;

    //  cout << std::put_time(&t,"%Y-%m-%d %s") << endl;
    //  std::time_t current = std::time(&time_stamp);

    return time_stamp;
};

int *convert2DateTime(TString stime) {

    int *vret = new int[2];

    std::time_t ts = convertTime(stime);
    tm *tmg = gmtime(&ts);
    vret[0] = 1000000 + (tmg->tm_year - 100) * 10000 + (tmg->tm_mon + 1) * 100 + tmg->tm_mday; // date
    vret[1] = tmg->tm_hour * 10000 + tmg->tm_min * 100 + tmg->tm_sec;

    return vret;
};

manageNTOF::manageNTOF(TString filepath) {
    centry = 0;
    nentries = 0;
    reader = NULL;
    fin = new TFile(filepath, "read");
    if (!fin->IsZombie()) {
        reader = (TTree *)fin->Get("tsel");
        reader->SetBranchAddress("time", &btime);
        reader->SetBranchAddress("date", &bdate);
        reader->SetBranchAddress("utime", &butime);
        reader->SetBranchAddress("pulse", &bpulse);
        nentries = reader->GetEntries();
        printf(" NTOF Tree entries :%lld\n", nentries);
    }
    start_time = 0;
}

manageNTOF::~manageNTOF() {
    if (fin->IsZombie()) {
        return;
    }
    fin->Close();
}

void manageNTOF::setStartDaTi(Long64_t udatetime) { // set start time (reference run time) of fers/srs run under analysis
    start_time = udatetime;
    printf(" Start Time : %lld\n", start_time);
    centry = 0;
}

std::vector<float> manageNTOF::getBeamPulse(float current_us, float dtc_s) {

    std::vector<float> rbp = {-1, -1};

    if (reader == NULL) {
        return rbp;
    }

    Long64_t actual_time = start_time + ((Long64_t)round(current_us / 1e6)); // second resolution

    Long64_t mini = 99999999; // minimum time difference between pulse timing and fers/srs timing in seconds

    int count = centry; //-1; // start from last search entry (expected date are read in increasing order)

    while (count < nentries) { // find min time
        reader->GetEntry(count);
        Long64_t delta = abs(butime - actual_time);
        // printf("%d  %d %d %d\n",count, delta, butime, actual_time);
        if (mini > delta) {
            mini = delta;
            centry = count;
            if (mini < delta) {
                break;
            } // minimum is before
            if (mini <= dtc_s) {
                break;
            } // this should be the matching condition
        }
        count = count + 1;
    }
    reader->GetEntry(centry);

    rbp[0] = bpulse;
    rbp[1] = (float)mini;

    // printf("  beam pulse %d : %lld %d %d %.3f %.1f\n", centry, actual_time, btime, bdate, rbp[0], rbp[1]);
    return rbp;
}

int manageNTOF::matchedEntry() { return centry; } // return last matched entry
