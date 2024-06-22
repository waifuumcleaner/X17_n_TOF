#
# Extract times from SRS _cluster processed data
#
#
# awk -v run=29 -v dtime=3.0 -f getSRStime.awk ../../test_231020/data/srs/run29_cluster.txt >out.txt
#

BEGIN {
    evtold=-1;
    printf("# SRS_Run_Number= %s\n",run);
    printf("# Delay_from_trigger_us= %s\n",dtime)
    qtot=0 # total charger for single event
    crec=0  # number of records in single event
    etime="-1";  # event time to be printed
}

(NF==6) {
    if ($2 == "Run") {
	printf("# RunStartTime= %s\n",$6);
    }
}

(NF==13) {
    printf("# TStamp_s nClusters TotCharge\n");
}

(NF==12) {  # records with data
    if ($1!=evtold) { # new event
	if (evtold>0) {
	    printf("%s %d %f\n",etime, crec, qtot);
	}
	evtold=$1;
	etime=$2 # sec
	crec=0;
	qtot=0;
    }
    crec=crec+1;  # count records in single event
    qtot=qtot+$8; # sum charge
}

END { # last event
    printf("%s %d %f\n",etime, crec, qtot);
}
