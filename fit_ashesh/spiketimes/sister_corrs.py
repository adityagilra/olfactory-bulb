# -*- coding: utf-8 -*-

## USAGE: python2.6 sister_corrs.py

import scipy.io
responses = scipy.io.loadmat('sispiketime.mat',struct_as_record=True)
### these are numpy structured arrays
responses = responses['sispiketime']

phases = scipy.io.loadmat('phases.mat',struct_as_record=True)
### these are numpy structured arrays
phases = phases['phases']

## responses[pairnum(size=20),sisnum(size=2)] each of which is 
## an ndarray[odornum(size=42)][trialnum(size=3 or 5)]
## each of which is itself an ndarray of shape numspikes x 1
## Thus responses[19,1][41,2][10,0] will give a spikeindex,
## but responses[19,1,41,2,10,0] will give an error
## Spike time is an index at 10kHz.

from pylab import *
import sys

sys.path.extend(["../../"])

from data_utils import * # for crosscorrgram()

corrbinwidth = 5e-3 # 5ms bin width
halfwindow = 300e-3 # s ## width of the sliding window
samplingf = 10000.0 # Hz float/double, sampling frequency: spike indices are sampled at 10kHz
tcorrlist = arange(-halfwindow,halfwindow+1e-6,corrbinwidth)
min_spikes = 50 # only consider those in peakshift histogram that have at least these many spikes
shift_cutoff = 100e-3 # plot those spike trains with corr peak shifts > this cutoff

def latespike(spikeslists, endtime):
    for spikeslist in spikeslists:
        if len(spikeslist)>0 and spikeslist[-1]>endtime: return True
    return False

def sis_corrs(norm_str):

    air_xcorrgram_list = []
    odor_xcorrgram_list = []
    for pair_response in responses:
        air_xcorrgram_list.append([])
        odor_xcorrgram_list.append([])
        for odornum in range(42):
            ## convert to spiketimes from spike indices at sampling frequency
            mit0_responses = pair_response[0][odornum]/samplingf
            mit1_responses = pair_response[1][odornum]/samplingf
            mit0_spiketimes = [ mit0_responses[trialnum].flatten() \
                for trialnum in range(len(mit0_responses)) ]
            mit1_spiketimes = [ mit1_responses[trialnum].flatten() \
                for trialnum in range(len(mit0_responses)) ]
            ## A set could be 10s air, 5s odor, 5s air OR 10s air, 10s odor, 5s air
            airstarttime = 0.0 # s
            airendtime = 10.0 # s
            odorstarttime = 10.0 # s
            ## if any of the spiketimes > 20s, consider it of the latter variety
            ## all trials will be either the former or the latter, so no complications.
            if latespike(mit0_spiketimes, 20.0) or latespike(mit1_spiketimes, 20.0):
                odorendtime = 20.0 # s
            else: odorendtime = 15.0 # s
                
            air_xcorrgram = crosscorrgram( mit0_spiketimes, mit1_spiketimes,\
                corrbinwidth, halfwindow, airstarttime, airendtime, norm_str )
            odor_xcorrgram = crosscorrgram( mit0_spiketimes, mit1_spiketimes,\
                corrbinwidth, halfwindow, odorstarttime, odorendtime, norm_str )
            air_xcorrgram_list[-1].append(array(air_xcorrgram))
            odor_xcorrgram_list[-1].append(array(odor_xcorrgram))

    ## replace all NaN's by 0's for plotting and averaging
    air_xcorrgram_list = array(air_xcorrgram_list)
    air_xcorrgram_list[where(isnan(air_xcorrgram_list))[0]] = 0
    odor_xcorrgram_list = array(odor_xcorrgram_list)
    odor_xcorrgram_list[where(isnan(odor_xcorrgram_list))[0]] = 0
    return air_xcorrgram_list, odor_xcorrgram_list

def plot_ind_xcorrgrams(air_xcorrgram_list, odor_xcorrgram_list):
    for pairnum in range(len(responses)):
        fig1 = figure(facecolor='none') # 'none' is transparent
        ## A super axes to set common x and y axes labels
        bigAxes1 = fig1.add_axes([0.1,0.1,0.8,0.8],frameon=False) # hide frame
        bigAxes1.set_xticks([])
        bigAxes1.set_yticks([])
        bigAxes1.text(-0.1,0.3,'spiking probability', fontsize=24, rotation='vertical')
        bigAxes1.text(-0.1,-0.11,'tau (s)',\
            fontsize=24, rotation='horizontal')
        bigAxes1.set_title('Cross-correlogram sispair='+str(pairnum),fontsize=30)
        for odornum in range(42):
            ax1 = fig1.add_subplot(6,7,odornum+1) ## 42 odor panels in a fig
            ax1.set_xticks([])
            ax1.plot(tcorrlist,air_xcorrgram_list[pairnum][odornum],'k-,')
            ax1.plot(tcorrlist,odor_xcorrgram_list[pairnum][odornum],'r-,')

def plot_average_xcorrgrams(air_xcorrgram_list, odor_xcorrgram_list):
    ## mean over two axes in two steps
    mean_intermed = mean(air_xcorrgram_list,axis=0)
    mean_airxcorrgram = mean(mean_intermed,axis=0)
    ## integral normalization for the average
    mean_airxcorrgram = mean_airxcorrgram/sum(mean_airxcorrgram)
    ## mean over two axes in two steps
    mean_intermed = mean(odor_xcorrgram_list,axis=0)
    mean_odorxcorrgram = mean(mean_intermed,axis=0)
    ## integral normalization for the average
    mean_odorxcorrgram = mean_odorxcorrgram/sum(mean_odorxcorrgram)
    fig1 = figure(facecolor='none') # 'none' is transparent
    ax1 = fig1.add_subplot(111)
    ax1.set_title('Average cross-correlogram',fontsize=36)
    ax1.plot(tcorrlist,mean_airxcorrgram,'k-,')
    ax1.plot(tcorrlist,mean_odorxcorrgram,'r-,')
    axes_labels(ax1,"time (s)","spike probability")

def plot_peakshift_histogram(air_xcorrgram_list, odor_xcorrgram_list):
    hist_width = 30e-3 #s ## should divide halfwindow integrally
    ## ensure odd number by "div2 +1"
    bins = int(4*halfwindow/hist_width)/2 + 1
    centralbinnum = bins/2 ## integer division
    timeshifthist = array([0.0]*bins)
    thistlist = arange(-halfwindow,halfwindow+1e-6,hist_width)
    multipeaks = 0
    relevant_xcorrs = 0
    shifted_responses = []
    for pairnum in range(len(responses)):
        shifted_responses.append((pairnum,[]))
        for odornum in range(42):
            xcorrgram = odor_xcorrgram_list[pairnum][odornum]
            #xcorrgram = air_xcorrgram_list[pairnum][odornum]
            ## at least min_spikes number of spikes
            if sum(xcorrgram) > min_spikes:
                maxval = max(xcorrgram)
                indices = where(xcorrgram==maxval)[0]
                peaktimes = (indices-1)*corrbinwidth - halfwindow
                if len(indices)>1:
                    print "pairnum,odornum,maxtimes =",\
                        pairnum,odornum,peaktimes
                    multipeaks += 1
                relevant_xcorrs += 1
                peaktime = mean(peaktimes)
                #for peaktime in peaktimes:
                timeshifthist[(peaktime/hist_width)+centralbinnum] += 1
                if peaktime > shift_cutoff:
                    shifted_responses[-1][1].append(odornum)
    print "The number of (pair,odor) combos with spikes >",\
        min_spikes,"is",relevant_xcorrs,"out of",pairnum*odornum
    print "The number of (pair,odor) combos with multiple peaks is",\
        multipeaks,"out of",relevant_xcorrs
    fig1 = figure(facecolor='none') # 'none' is transparent
    ax1 = fig1.add_subplot(111)
    ax1.set_title('Differential peak shift histogram',fontsize=36)
    ax1.plot(thistlist,timeshifthist)
    axes_labels(ax1,"time (s)","number")
    return shifted_responses

def convert_respimage(trialslist, spiketimeslist):
    phaseslist = []
    for trialnum,spikephases in enumerate(trialslist):
        ## to allow responses to be a numpy array,
        ## spikephases is [[spiketime1],[spiketime2],...]
        spikephases = spikephases.flatten()
        oldphase = 0.0
        air_phases = [[]]
        odor_phases = [[]]
        spiketimes = spiketimeslist[trialnum]
        for spikenum,phase in enumerate(spikephases):
            spiketime = spiketimes[spikenum]
            ## air delivery period
            if spiketime < 10.0: # s
                if phase < oldphase: air_phases.append([])
                air_phases[-1].append(phase)
            ## odor delivery period
            elif spiketime < 15.0: # s
                if phase < oldphase: odor_phases.append([])
                odor_phases[-1].append(phase)
            oldphase = phase
            ## discard the later air (or sometimes odor) period:
            ## some trials have odor from 15 to 25s, instead of just till 20s,
            ## but we discard that. air period after odor is unreliable,
            ## so discard that too.
        phaseslist.append((air_phases,odor_phases))
    return phaseslist

def plot_respimages(listof_phaseslists, phasemax, dt):
    """ listof_phaseslists[mit_odor_combonum][trialnum][0=air|1=odor][respnum][phasenum]
    """
    fig = figure(facecolor='w') ## 'none' for transparent
    numplots = len(listof_phaseslists)
    for phaseslistnum,phaseslist in enumerate(listof_phaseslists):
        #ax = fig.add_subplot(2,round(numplots/2.0),phaseslistnum+1)
        ax = fig.add_subplot(numplots,1,phaseslistnum+1)
        ## trialnum is not important, what matters
        ## is how many resp cycles we have plotted
        image = None
        for trialnum, phaselist in enumerate(phaseslist):
            air_phaselist = phaselist[0]
            odor_phaselist = phaselist[1]
            airimage = get_phaseimage(air_phaselist,\
                phasemax, dt, overlay=True, rasterwidth=100)
            if image is None:
                image = airimage
            else:
                image = append(image, airimage, axis=1)
            odorimage = get_phaseimage(odor_phaselist,\
                phasemax, dt, overlay=True, rasterwidth=100)
            image = append(image, odorimage, axis=1)
        imy,imx = image.shape
        ax.imshow(image, cmap=cm.binary, aspect=1.0)
        ax.set_ylim(0,imy)
        axes_off(ax)

def plot_shifted_responses(response_ids):
    for pairnum,odornums in response_ids:
        if len(odornums)==0: continue
        odor_responses = []
        odor_phases = []
        for odornum in odornums:
            ## convert to spiketimes from spike indices at sampling frequency
            mit0_responses = responses[pairnum,0][odornum]/samplingf
            mit1_responses = responses[pairnum,1][odornum]/samplingf
            mit0_spiketimes = [ mit0_responses[trialnum].flatten() \
                for trialnum in range(len(mit0_responses)) ]
            mit1_spiketimes = [ mit1_responses[trialnum].flatten() \
                for trialnum in range(len(mit0_responses)) ]
            odor_responses.append(mit0_spiketimes)
            odor_responses.append(mit1_spiketimes)
            
            ## convert firing phases (0 to 1 along respiration) into 2D array
            mit0_phaseslist = convert_respimage(phases[pairnum,0][odornum], mit0_spiketimes)
            mit1_phaseslist = convert_respimage(phases[pairnum,1][odornum], mit1_spiketimes)
            odor_phases.append( mit0_phaseslist )
            odor_phases.append( mit1_phaseslist )
        ## 20 s runtime (odor from 10 to 15 s, sometimes 20s: see above)
        plot_rasters(odor_responses, 20.0,\
            colorlist = ['r','m','g','y','b','c'], labels=False)
        plot_respimages(odor_phases, 1.0, 10e-3)

if __name__ == "__main__":
    #air_xcorrgram_list, odor_xcorrgram_list = sis_corrs("overall")
    #plot_ind_xcorrgrams(air_xcorrgram_list, odor_xcorrgram_list)

    air_xcorrgram_list, odor_xcorrgram_list = sis_corrs("none")
    
    #plot_average_xcorrgrams(air_xcorrgram_list, odor_xcorrgram_list)
    
    shifted_responses = \
        plot_peakshift_histogram(air_xcorrgram_list, odor_xcorrgram_list)
    plot_shifted_responses(shifted_responses)
    show()
