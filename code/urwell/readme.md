# uRwell X17 analysis tools

Based on the SRS-APV25 electronics and the mmdaq ATLAS acquistion kindly provided by Mauro Iodice

WORK IN PROGRESS

# Software Requirements

root 6
( python3 with numpy, matplotlib, lmfit ) 

# Codes

- urwana.cpp : main methods (to be clean!), at the moment produce the runXX_track.txt file with tracks information
- mmRawTree.C .h : the interface to read and do the first processing of the TTrees of the mmdaq output root files
- utilf.h : utility functions (e.g. cluster search class ...)
- eventDisplay.h : event display functions
- multisetfit.py : plots the tracks data produced by urwana.cpp and tries to compute the alignment parameters
- combine.cpp : Combine nTOF, FERS and SRS data, based on trigger time
- getSRStime.awk : awk script to extract timing information (et al.) from the processed text data file

# Usage

enter root:

> .L urwana.cpp+
>  decodePhysRun(9, 12, 3, 5, 0)  (if the last par is 1, each event is displayed = event display)
(more details in the urwana.cpp comment before the method definition)

NOTE: There are some hardcoded parameters:
  NAPV: number of APVs used - this shall be set to the proper values used in the acquisition
  DATAPATH: path of the acquired data (this is actually the last default parameter in the macro call)

The output is a text file that can be easily read by the TTree::ReadFile() method

# How it work (concisely)

in *urwana* and *utilf* for each event:
 1. get charge hits (vs time samples and chamber/strips) 
 2. subtract pedestal
 3. if larger than `nsigma * pedestal-rms` assume as signal
 4. for each strip with signal, determine the "time" clusters, that is a minimum number (selectable) of adjacent samples with some charge
 5. for each chamber determine space clusters; each cluster has adjacent strips with signal in overlapping time slots
 6. for each time-space cluster determine:
   1. total charge,
   2. strips involved, strip centroide and standard deviation, 
   3. start and end time (by Landau fit function at half maximum) on each clusters (combining sample serie and adjacent strips)
   4. for clusters with multiple strips: evaluate the linear fit between strip position and related signal start_time (`strip = s0 + m * start_time`);
      `m` is the `vdrift/tan(theta)` where `theta` is the track entrance angle into the chamber (on the plane perpendicular to the direction of the strips);
	  `vdrift` is the electron drift velocity in the chamber, depending on the electric field and gas mixture (~4 cm/us for Ar/CO2 70/30 and ~1.5 kV/cm); 
	  one can estimate vdrift from the relation: `tend-tstart * vdrift = dgap` where: tstart and tend are the first and last `start_time` 
	  (not trivial on a single strip), and dgap is the cathode-anode(top) distance (which is known by construction).
 7. in case of multiple chambers, search for tracks, that is at least 2 clusters on two different chambers, 
    whose start times differ for less than a give delta_time (selectable); 
	the search starts from the first chamber and currently all chambers shall have a cluster belonging to the track
 8. the tracks information are saved into a text file, each line has the following parameters:
     EventId TrackId ChamberId TimeStart TimeEnd Charge StripFirst StripLast StripCentroid StripRMS *mlFit* qlFit

then in the *multisetfit* (not used for CERN test):
 1. select events with single track
 2. reconstruct the slope of the track from the stripCentroid and distance between chambers; 
    for 3 chambers 2 slopes can be reconstructed for a single track.
	The mean of these two slopes is assumed as the slope of the track which shall be related to *mlFit* from the "TPC" approach.
 3. for straight track, the estimated slopes shall be very similar (to be improved by estimating the chamber misalignments); 
    then select events (tracks) which have "similar" estimated slopes
 4. different plots are produced; the most relevant are expected to be the scatter plot
    between *mean reconstructed track slope by strip centroid* vs the *mlFit* slope evaluated by TPC method; 
	both are converted to angles (relative to the vertical axis). T

*multisetfit* should also be able to estimate the alignment parameters trying to minimize the distance between the reconstructed slopes on different chamber couples. 
Currently the leastsquare method does not converge properly ... probably some free parameters are cross-related 

# Examples of processing flow

Single Event display mode:

enter root then:
>.L urwana.cpp+
> decodePhysRun(29, 27, 5, 6, 18, 1, 0)   // for noisy events

simple post processing of cluster data (optionally combined to CUBE data):

> readClusterOut(29)



# To DO
Move from txt output file to a root file in decodePhysRun (using the *_ana.root file already created for histograms)

Clean up the code from methods no longer used and/or replaced

-----

---
Some post processing with:

>  python3 multisetfit.py 

(originally developed for cosmic post processing analysis; no longer used ...)



TO BE REMOVED:
DA TRASFERIRE SU CODICE per scrivere parametri nel root file di output
 //TParameter<char> *SH_Atypea = new TParameter<char>("type_SH",SH_Atype);

  TParameter<float> *poffx_SH_type = new TParameter<float>("OffsetX_SH_type",offsetx_SH_type);
  TParameter<float> *poffy_SH_type = new TParameter<float>("OffsetY_SH_type",offsety_SH_type);
  TParameter<int> *nanodex_SH_type = new TParameter<int>("NumAnodeX_SH_type",paramx_SH_type);
  TParameter<int> *nanodey_SH_type = new TParameter<int>("NumAnodeY_SH_type",paramy_SH_type);

  TParameter<float> *poffx_LH_parh = new TParameter<float>("OffsetX_LH_parh",offsetx_LH_parh);
  TParameter<float> *poffy_LH_parh = new TParameter<float>("OffsetY_LH_parh",offsety_LH_parh);
  TParameter<int> *nanodex_LH_parh = new TParameter<int>("NumAnodeX_LH_parh",paramx_LH_parh);
  TParameter<int> *nanodey_LH_parh = new TParameter<int>("NumAnodeY_LH_parh",paramy_LH_parh);

  TParameter<float> *pixsize_SH_type = new TParameter<float>("PixelSize_SH_type",pmt_pixel_size_SH_type);
  TParameter<float> *pixsize_LH_parh = new TParameter<float>("PixelSize_LH_parh",pmt_pixel_size_LH_parh);
  TParameter<int> *pnevt = new TParameter<int>("SimuEvents",evtID);
   
   TTree *tma = new TTree("tma","Anode Charge Matrix");
   
  TList *ui = (TList*) tma->GetUserInfo();

  if (!ui) { printf("cannot get user info from digitized root tree ... this is very strange\n"); }

  ui->SetName(Form("Digitization from simulated data file %s\n",inputfile));

  //ui->Add(SH_Atypea);
  
  ui->Add(poffx_SH_type);
  ui->Add(poffy_SH_type);
  ui->Add(nanodex_SH_type);
  ui->Add(nanodey_SH_type);

  ui->Add(poffx_LH_parh);
  ui->Add(poffy_LH_parh);
  ui->Add(nanodex_LH_parh);
  ui->Add(nanodey_LH_parh);

  ui->Add(pixsize_SH_type);
  ui->Add(pixsize_LH_parh);
  ui->Add(pnevt);
  
  tma->Write();
  fTdata->Write();

  fout->Flush();
  fout->Close();
