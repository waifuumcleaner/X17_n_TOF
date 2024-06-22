#
# Try to estimate raising and trailing time constant of the GEM-APV25 signals
# using "golden signals" and fitting simultaneously, keeping the two parameter identical to all signals
#
# inspired by: https://stackoverflow.com/questions/20339234/python-and-lmfit-how-to-fit-multiple-datasets-with-shared-parameters
#
# version 2.0 Jun/2023
# derived from jlab/tracker/bprog/analysis

import sys
import random
import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit

#### APV25 typical function
# amp : amplitude
# sta : starting point
# rup : raising constant
# tdo : trailing constant
# npow: power index

def fapv3(x, amp, sta, tau):
    "three pars APV function (Integral from x-sta=0 to infty is amp"
    xx = (x - sta)  
    return amp*0.5*(np.sign(xx)+1)*xx/(tau*tau)*np.exp(-xx/tau)

def fapv4(x, amp, sta, tau, npow):
    "4 pars APV function generalized by Roberto Perrino, integral should be: npow! * amp * tau ... from wikipedia"
    xx = (x - sta)
    yy = xx/tau
    return amp*0.5*(np.sign(xx)+1)*np.power(yy,npow)*np.exp(-yy)

def fapv5(x, amp, sta, rup, tdo):
    "basic 4 pars  APV function, integral is amp"
    xx = (x - sta)
    fac = (rup+tdo)/tdo/tdo
    return amp*fac*0.5*(np.sign(xx)+1)*(1-np.exp(-xx/rup))*np.exp(-xx/tdo)
    
###
#

def fapv_dataset(params, i, x, selfunc):
    """calc fapv from params for data set i
    using simple, hardwired naming convention"""
    amp = params['a_%i' % i].value
    sta = params['s_%i' % i].value
    rup = params['lea'].value
    tdo = params['tra'].value
    ff=0
    if (selfunc==3):
        ff=fapv3(x, amp, sta, rup)
    if (selfunc==4):
        ff=fapv4(x, amp, sta, rup, tdo) # tdo is used for "n" power
    if (selfunc==5):
        ff=fapv5(x, amp, sta, rup, tdo)

    return ff

def objective(params, x, data, sdata, selfunc):
    """ calculate total residual for fits to several data sets held
    in a 2-D array, and modeled by fapv functions"""
    ndata, nx = data.shape
#    print "objective: %d %d" % (ndata,nx)
    resid = 0.0*data[:]
    # compute residuals per data set, with weighted error
    for i in range(ndata):
        #  print " %d %f " % (i,params['r_%i' % (i+1)].value)
        resid[i, :] = (data[i, :] - fapv_dataset(params, i, x, selfunc)) / sdata[i, :]

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

##
# From chamber reference to lab reference
#  lab reference is defined by the top chamber 
#     x axis is vertical, y axis is horizontal
#     origin of x axis is the first strip of the top chamber
#
#  ystrip: [mm] the 1D strip coordinate in the chamber (data)
#  xdist: [mm] vertical distance of the origin of the chamber from the origin of the lab; from config file or alignment tuning
#  angle: [deg] the angle of the chamber relative to vertical downward axis; 90 deg is horizontal chamber
#
# return the lab coordinates

def chamber2lab(ystrip, xdist, angle):
    arad = angle * 3.1415/180.
    xx = xdist + ystrip * np.cos(arad)  # cos(90+theta) = - cos(theta)
    yy = ystrip * np.sin(arad)

    return xx,yy

####
def histo2Stat(data):

    distr,xbin = np.histogram(data, bins=100)
    imax = np.argmax(distr)
    mids = (xbin[1:]+xbin[:-1])/2.
    vmax = mids[imax]
    mean = np.average(mids, weights=distr)
    std = np.sqrt(np.average((mids-mean)**2, weights=distr))

    return (vmax,mean,std)

###
#
# sel: selected single tracks events
# data: each row is a chamber hits
# cdist: chamber to chamber distances (first element should be 0), evaluated along the vertical axix
# cangle: angle of each chamber relative to vertical axis

def plotAll(run, sel, data, strip_pitch, cdist, cangle, gap):

    ncs=len(cdist)  # number of chambers

    nevent0=len(sel[:,0])
    # select tracks which likely are straight  
    yy=[]
    mm=[]
    for i in range(ncs): # loop on chambers
        sc=data[sel[:,i],9] # hit strip centroid
        yy.append(sc)
        
    figc,axc = plt.subplots(2,2)

    xlab0,ylab0 = chamber2lab(yy[0],cdist[0]*strip_pitch, cangle[0])  # = 0,yy[0]
    xlab1,ylab1 = chamber2lab(yy[1],cdist[1]*strip_pitch, cangle[1])
    xlab2,ylab2 = chamber2lab(yy[2],cdist[2]*strip_pitch, cangle[2])

    dx0 = (ylab1-ylab0)/(xlab1-xlab0)  # angular coefficient of the track
    dx1 = (ylab2-ylab0)/(xlab2-xlab0)

    dangle=np.arctan((dx0+dx1)/2)/3.1415*180. # slope angle (relative to vertical direction) in deg

    angle_distr,xbin = np.histogram(dangle, bins=100)
    
    angmin=np.min(dangle)
    angmax=np.max(dangle)
    angmax = np.argmax(angle_distr)

    print((xbin[angmax]+xbin[angmax+1])/2.)
    
    print("Track slope angular range: ",angmin, angmax)
    
    dd=dx1-dx0
    dd_distr,xbin = np.histogram(dd, bins=100)
    ddamax = np.argmax(dd_distr)
    ddmax = (xbin[ddamax+1]+xbin[ddamax])/2.
    
    print("dd max : ",ddmax)
    
    gg = np.where((dd>(ddmax-1)) & (dd<(ddmax+1)))[0]  # "good" tracks
    nevent1=len(gg)

    # axc[0,0].hist2d(dx1,dx0)
    axc[0,0].scatter(dx1,dx0)
    axc[0,0].scatter(dx1[gg],dx0[gg])
    axc[0,0].set_title("Track slope from chamber 0-1 and chamber 0-2")
    axc[0,0].set_xlabel("Angular Coeff. chamber 0-2")
    axc[0,0].set_ylabel("Angular Coeff. chamber 0-1")
    
    axc[0,1].hist(dangle, bins=50, histtype='bar')
    axc[0,1].hist(dangle[gg], bins=50, histtype='bar')
    axc[0,1].set_title("Track mean angle")
    axc[0,1].set_xlabel("Angle (deg)")

    axc[1,1].hist(dd, bins=50, histtype='bar')
    axc[1,1].hist(dd[gg], bins=50, histtype='bar')
    axc[1,1].set_title("Track slopes difference")
    axc[1,1].set_xlabel("Angular Coeffs Difference")
    
    figa,axa = plt.subplots(ncs,4, sharex='col')
    figb,axb = plt.subplots(ncs,4, sharex='col')
    figd,axd = plt.subplots(ncs,3, sharex='col')

    figa.suptitle("Run {:d} Single Track Events {:d} filtered {:d}".format(run,nevent0,nevent1))
    figb.suptitle("Run {:d} Single Track Events {:d} filtered {:d}".format(run,nevent0,nevent1))
    figc.suptitle("Run {:d} Single Track Events {:d} filtered {:d}".format(run,nevent0,nevent1))
    figd.suptitle("Run {:d} Single Track Events {:d} filtered {:d}".format(run,nevent0,nevent1))

    
    xcmean=[]
    for i in range(ncs): # loop on chambers
        print("Chamber ",i) 
        # TimeStart
        ts = data[sel[:,i],3]
        axa[i,0].hist(ts,bins=50, histtype='bar')
        axa[i,0].hist(ts[gg],bins=50, histtype='bar')
        axa[i,0].set_title("Time Start")

        # deltaTime
        dtstart = data[sel[:,i],4]
        axa[i,1].hist(dtstart,bins=50, histtype='bar')
        axa[i,1].hist(dtstart[gg],bins=50, histtype='bar')
        axa[i,1].set_title("Delta Start Times")
        
        # deltaTime (stop - start)
        dt = data[sel[:,i],5]-data[sel[:,i],3]
        axa[i,2].hist(dt,bins=50, histtype='bar')
        axa[i,2].hist(dt[gg],bins=50, histtype='bar')
        axa[i,2].set_title("Delta Time Start-Stop")
        
        # Charge
        q=data[sel[:,i],6]
        axa[i,3].hist(q,bins=50, histtype='bar')
        axa[i,3].hist(q[gg],bins=50, histtype='bar')
        axa[i,3].set_title("Charge")

        # StripCentroid
        sc=data[sel[:,i],9]
        axb[i,0].hist(sc,bins=50, histtype='bar')
        axb[i,0].hist(sc[gg],bins=50, histtype='bar')
        axb[i,0].set_title("Strip Center")
        xcmean.append(np.mean(sc[gg]))
        
        # StripDelta vs RMS
        sd=data[sel[:,i],8]-data[sel[:,i],7]
        # axb[i,1].hist(sd,bins=60, histtype='bar')
        #        axb[i,1].set_title("Strip Delta")

        # StripRMS
        ss=data[sel[:,i],10]
        # axb[i,1].hist(ss,bins=60, histtype='bar')
        axb[i,1].scatter(ss,sd)
        axb[i,1].scatter(ss[gg],sd[gg])
        axb[i,1].set_title("Strip RMS vs Delta")

        # mlFit = track_slope * v_drift  
        ml=data[sel[:,i],11]*(strip_pitch/10.)*(1000./25.)  # from strip_index/25_ns to cm/us; NOTE ml = m * v_drift where m is tan(theta); v_drift ~ few cm/us; gap is the drift gap in chamber
        axb[i,2].hist(ml,bins=50, histtype='bar')
        axb[i,2].hist(ml[gg],bins=50, histtype='bar')
        axb[i,2].set_title("m Linear Fit")
        print("  mlFit average {:.3f} +/- {:.3f}".format(np.mean(ml),np.std(ml)))
        
        # qlFit
        ql=data[sel[:,i],12]
        axb[i,3].hist(ql,bins=50, histtype='bar')
        axb[i,3].hist(ql[gg],bins=50, histtype='bar')
        axb[i,3].set_title("q Linear Fit")

        # estimte vdrift ... very preliminarely
        gn = np.where(dtstart>1)[0]
        vdrift = (gap[i]/10.)/(dtstart[gn]*0.025) # vdrfit cm/us

        axc[1,0].hist(vdrift, bins=20, histtype='bar')
        axc[1,0].set_title("Estimated Drift Speed (gap/dt)")
        axc[1,0].set_xlabel("Speed cm/us")
        vdrift_mean = np.mean(vdrift)/2./3.
        print("   vdrift {:.3f} cm/us".format(vdrift_mean))

        # track slope scatter...
        eangle = np.arctan(-ml/vdrift_mean)/3.1415*180. + (cangle[i]-90.) # angle estimated from "TPC" method 
        #gm = np.where((eangle[gg]>angmin*2) & (eangle[gg]<angmax*2) & (ql[gg] > -1e4))[0]  # ql=-1e5 is invalid value
        gm = np.where(ql[gg] > -1e4)[0]  # ql=-1e5 is invalid value
        g2 = gg[gm]

        eamax, eamean, eavar = histo2Stat(eangle[g2])
        
        gm = np.where((np.abs(eangle[gg]-eamean)<2.*eavar) & (ql[gg]>-1e4))[0]
        g2 = gg[gm]
        
        print(" Estimated angle max mean std: ",eamax, eamean,eavar)
        #        axd[i,2].hist(eangle[g2],bins=50, range=(eamean-2*eavar,eamean+2*eavar), histtype='bar')
        #        axd[i,2].hist(dangle[g2],bins=50, range=(eamean-2*eavar,eamean+2*eavar), histtype='bar')

        #axd[i,0].hist2d(dangle[g2],eangle[g2]) #,bins=(np.arange(-30,30,1),np.arange(-30,30,1)))
        axd[i,0].scatter(dangle[g2],eangle[g2])
        if (i==0): axd[i,0].set_title("Track Slope: theta TPC vs theta Centroid")
        if (i==(ncs-1)): axd[i,0].set_xlabel("Theta Strip Centroid (deg)")
        axd[i,0].set_ylabel("Theta TPC (deg)")
        par = np.polyfit(dangle[g2],eangle[g2],1)
        print("  Angle_TPC = {:.2f} Angle_Strip + {:.2f} (points {:d})".format(par[0],par[1],len(gm)))

        ares = eangle[g2]-(par[0]*dangle[g2]+par[1])
        armax, armean, arvar = histo2Stat(ares)
        
        axd[i,2].hist(ares,bins=50, histtype='bar')
        print(" reasidual angle = ",armax, armean, arvar)
        
        axd[i,1].scatter(dangle,sd)
        #axd[i,1].hist2d(dangle[gg],sd[gg])
        axd[i,1].set_title("Strip rms vs Track Slope")

        # axd[i,2].scatter(dangle,dtstart)
        # axd[i,2].scatter(dangle[gg],dtstart[gg])
        # axd[i,2].set_title("Delta Start Time vs Track Slope")

        
    print("Average shifts: {:.3f} {:.3f}".format(xcmean[1]-xcmean[0],xcmean[2]-xcmean[0]))

    return gg

###
#

def objTrack(params, x, y):
    """ calculate total residual for fits to several data sets held
    in a 2-D array, and modeled by fapv functions"""

    s1 = params['s_1'].value
    s2 = params['s_2'].value
    theta1 = params['theta_1'].value
    theta2 = params['theta_2'].value
    h1 = params['h_1'].value
    h2 = params['h_2'].value

    sin1=np.sin(theta1)
    sin2=np.sin(theta2)
    cos1=np.cos(theta1)
    cos2=np.cos(theta2)

    nn = len(x) # number of events

    #    resid = np.zeros(nn)
    w1 = x
    w2 = y
    y1 = (w1 + s1)*sin1
    y2 = (w2 + s2)*sin2
    x1 = h1 + (w1 + s1)*cos1
    x2 = h2 + (w2 + s2)*cos2

    resid = (y1/x1 - y2/x2)

    # r2=resid*resid
    # imx = np.argmax(r2)
    # mx=np.max(r2)
    # mean=np.mean(r2)
    # std=np.std(r2)
    # print(mean,std,mx,imx)

    return resid*resid

#
#
#

def ingestConfig(rnum, infconfig):

    print("Ingest config file ", infconfig)
    
    cdata = np.loadtxt(infconfig)
    print(" Config array shape ",cdata.shape)
    
    iln = np.where(cdata[:,0]==rnum)[0] # extract line corresponding to the current run

    if (len(iln)==0):
        print("ERROR: cannot find the current run numer %d in config file" % rnum)
        exit(0)

    cline = cdata[iln[0],1:]  # remove run number
    nchamber = int(len(cline)/5) # 5 parms for each chamber

    idx = np.arange(3)*5
    field=cline[idx]    # electric field in chamber drift 
    xorig=cline[idx+1]  # position of chamber along x-axis (vertical) relative to laboratory frame
    yorig=cline[idx+2]  # y-origin of chamber relative to lab frame
    angle=cline[idx+3]  # angle of chamber relative to vertical x-axis
    gap=cline[idx+4]    # drift gap
    
    print(" returned data: ",field,xorig,yorig,angle)

    return field,xorig,yorig,angle,gap

# -----------------------------------------------------
# MAIN START HERE
# use example:
# python3 multisetfit.py


g_pitch = 0.6 # [mm] strip pitch in chamber

runnum=9
inpath="/home/cisbani/proj/ntof/x17-gustavino/test/2304/"
infilew="run{:d}_track.txt"
configfile="run-config.txt"

if (len(sys.argv)<2):
    print("Syntax: %s cosmic_run_number [input_path run_file_name config_file_name]" % sys.argv[0])
    print(" cosmic_run_number: is inserted in the run_file_name string and is used to extract data from config_file_name [%d]" % runnum)
    print(" input_path: where bith run_file and config_file sits [%s]" % inpath)
    print(" run_file_name: with wildcard for the run number position [%s]" % infilew)
    print(" config_file_name: in input_path, contains information on the run, e.g. chamber configurations [%s] " % configfile)

if (len(sys.argv)>1):
    runnum = int(sys.argv[1])

if (len(sys.argv)>2):
    inpath=sys.argv[2]

if (len(sys.argv)>3):
    infilew = sys.argv[3]

if (len(sys.argv)>4):
    configfile=sys.argv[4]

infile = inpath + "/" + infilew.format(runnum)
confpath = inpath + "/" + configfile

print("input pars:")
print("  Run Number : ",runnum)
print("  Input File : ", infile)
print("  Config File: ", confpath)

g_field, g_hdist, g_yshift, g_angle, g_gap = ingestConfig(runnum,confpath)

g_nchamber=len(g_angle)

print(g_field, g_hdist, g_yshift, g_angle, g_gap)

# ingest data
data = np.loadtxt(infile)

print(" Data shape ",data.shape)

# prepare data for tracking fit
# assume each chamber has cluster (not checked)

goodevt=[]

vhit=[]

# not optimized but more readable code (hopefully)

evtidx=np.unique(data[:,0]) # list of events

for evt in evtidx:
    pts = np.where(data[:,0]==evt)[0]
    sevt=data[pts] # single event data
    ltracks = np.unique(sevt[:,1]) # list of track indices
    if (len(ltracks)>1): # consider events with single tracks only
        continue
    lchambers = sevt[:,2]
    if (len(lchambers)!=g_nchamber): # all chambers shall have one hit
        continue

    w0 = sevt[0,9] # StripCentroid of first chamber is the reference

    wi=[]
    for hit in sevt[1:]: # loop on the other chambers 
        wi.append(hit[9]-w0)

    vhit.append(wi)  # for each event contains the hits coordinate along each chamber relative to the first chamber

    goodevt.append(list(pts)) # list of good, filtered hits
    
ahit = np.asarray(vhit)
x = np.asarray(g_hdist)

print(" ahit shape ",ahit.shape)

goodevt = np.asarray(goodevt)
print(" goodevt shape ", goodevt.shape)

gsel = plotAll(runnum, goodevt, data, g_pitch, g_hdist, g_angle, g_gap)

ahit = ahit[gsel]

xvar=ahit[:,0]
yvar=ahit[:,1]


fit_params = Parameters()
ang=np.asarray(g_angle)/180.*3.1415
flag=False

for i in range(1,g_nchamber):
    fit_params.add( 's_%i' % i, value=0, min=-5., max=5., vary=True) # [strip] shift of i chamber relative to origin of first
    fit_params.add( 'theta_%i' % i, value=ang[i], min=ang[i]-.1, max=ang[i]+.1, vary=True) # [deg] tilt (angle) of i chamber relative to "x" axis (vertical axis)
    fit_params.add( 'h_%i' % i, value=g_hdist[i], min=g_hdist[i]-5, max=g_hdist[i]+5, vary=flag) # [cm] distance along x of i chamber x-intercept respect to origin of first chamber
    flag=True

result = minimize(objTrack, fit_params, args=(xvar,yvar), method='dual_annealing') #method='leastsq')   # x is dummy!!!
report_fit(result.params, min_correl=0.9999)

res = objTrack(result.params, xvar, yvar)
print("Minimization residual: ",np.sum(res*res))

plt.show()

print("Try random minimization: ")

sres_min=9e9
sres_old=9e9
for k in range(100000):
    if (((k%1000)==0)and(sres_min!=sres_old)):
        sdum = " {:10d} {:.3f} ->".format(k,sres_min) 
        for par in fit_min:
            sdum = sdum + " {:7.3f}".format(fit_min[par].value)
        print(sdum)
        sres_old = sres_min

    for pp in fit_params:
        rnd = np.random.rand()
        vmin = fit_params[pp].min
        vdd = fit_params[pp].max - vmin
        fit_params[pp].value = rnd*vdd + vmin
        
    res = objTrack(fit_params, xvar, yvar)
    sres = np.sum(res*res)
    if (sres<sres_min):
        sres_min=sres
        fit_min = fit_params

print(" final Res2min: {:.3f}".format(sres_min))
for par in fit_min:
    print(" {:8s}: {:7.3f}".format(par, fit_min[par].value))

exit(0)

