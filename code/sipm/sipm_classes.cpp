#include "TFile.h"

#include "sipm_classes.h"
#include "../cube/ntofDAQ.h"
RunsInfo::RunsInfo(std::vector<int> numbers, std::vector<int> delays, std::vector<TargetType> targetTypes, std::vector<Mode> modes) : m_numbers{numbers}, m_delays{delays}, m_targetTypes{targetTypes}, m_modes{modes} {
    if (m_numbers.size() != m_delays.size() || m_numbers.size() != m_targetTypes.size() || m_numbers.size() != m_modes.size()) {
        throw std::invalid_argument("All vectors must have the same size.\n");
    }
}
void RunsInfo::addRun(int number, int delay, TargetType target, Mode mode) {
    m_numbers.push_back(number);
    m_delays.push_back(delay);
    m_targetTypes.push_back(target);
    m_modes.push_back(mode);
}
std::vector<int> RunsInfo::getNumbers(const std::optional<std::vector<int>> &subset) const {
    if (!subset)
        return m_numbers;
    std::vector<int> result;
    for (int num : *subset) {
        if (std::find(m_numbers.begin(), m_numbers.end(), num) != m_numbers.end()) {
            result.push_back(num);
        }
    }
    return result;
}

std::vector<int> RunsInfo::getDelays(const std::optional<std::vector<int>> &subset) const {
    if (!subset)
        return m_delays;
    std::vector<int> result;
    for (int num : *subset) {
        auto it = std::find(m_numbers.begin(), m_numbers.end(), num);
        if (it != m_numbers.end()) {
            result.push_back(m_delays[std::distance(m_numbers.begin(), it)]);
        }
    }
    return result;
}

std::vector<TargetType> RunsInfo::getTargetTypes(const std::optional<std::vector<int>> &subset) const {
    if (!subset)
        return m_targetTypes;
    std::vector<TargetType> result;
    for (int num : *subset) {
        auto it = std::find(m_numbers.begin(), m_numbers.end(), num);
        if (it != m_numbers.end()) {
            result.push_back(m_targetTypes[std::distance(m_numbers.begin(), it)]);
        }
    }
    return result;
}

std::vector<Mode> RunsInfo::getModes(const std::optional<std::vector<int>> &subset) const {
    if (!subset)
        return m_modes;
    std::vector<Mode> result;
    for (int num : *subset) {
        auto it = std::find(m_numbers.begin(), m_numbers.end(), num);
        if (it != m_numbers.end()) {
            result.push_back(m_modes[std::distance(m_numbers.begin(), it)]);
        }
    }
    return result;
}

std::vector<int> RunsInfo::getPedestalRuns() const {
    std::vector<int> pedestalRuns;
    std::copy_if(m_numbers.begin(), m_numbers.end(), std::back_inserter(pedestalRuns), [this, i = 0](int) mutable { return m_modes[i++] == Mode::Pedestal; });
    return pedestalRuns;
}
std::vector<int> RunsInfo::getPhysicalRuns() const {
    std::vector<int> physicalRuns;
    std::copy_if(m_numbers.begin(), m_numbers.end(), std::back_inserter(physicalRuns), [this, i = 0](int) mutable { return m_modes[i++] == Mode::Physical; });
    return physicalRuns;
}

void RunsInfo::setNumbers(const std::vector<int> &numbers) { m_numbers = numbers; }
void RunsInfo::setDelays(const std::vector<int> &delays) { m_delays = delays; }
void RunsInfo::setTargetTypes(const std::vector<TargetType> &targetTypes) { m_targetTypes = targetTypes; }
void RunsInfo::setModes(const std::vector<Mode> &modes) { m_modes = modes; }
void RunsInfo::setMode(int run_number, Mode mode) {
    auto it = std::find(m_numbers.begin(), m_numbers.end(), run_number);
    if (it != m_numbers.end()) {
        m_modes[std::distance(m_numbers.begin(), it)] = mode;
    }
}
void RunsInfo::setTargetType(int run_number, TargetType targetType) {
    auto it = std::find(m_numbers.begin(), m_numbers.end(), run_number);
    if (it != m_numbers.end()) {
        m_targetTypes[std::distance(m_numbers.begin(), it)] = targetType;
    }
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
