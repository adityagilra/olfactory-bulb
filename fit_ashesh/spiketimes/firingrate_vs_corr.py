# -*- coding: utf-8 -*-

## USAGE: python2.6 firingrate_vs_corrs.py

import scipy.io, scipy.stats
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

## phases is just like responses,
## just that phase is a respiration phase from 0 to 1.
## To distinguish which phase entry belongs to air / odor period,
## get the indices from the responses and use these on phases.

from pylab import *
import sys

## set path below to where you can find data_utils.py from my OB model files.
sys.path.extend(["../../"])
from data_utils import * # for plotting utility functions

corrbinwidth = 5e-3 # 5ms bin width
halfwindow = 300e-3 # s ## width of the sliding window
samplingf = 10000.0 # Hz float/double, sampling frequency: spike indices are sampled at 10kHz
tcorrlist = arange(-halfwindow,halfwindow+1e-6,corrbinwidth)
min_spikes = 50 # only consider those in peakshift histogram that have at least these many spikes
shift_cutoff = 100e-3 # plot those spike trains with corr peak shifts > this cutoff

numbins = 5
mouse_respirationT = 0.5 # half a second

def latespike(spikeslists, endtime):
    for spikeslist in spikeslists:
        if len(spikeslist)>0 and spikeslist[-1]>endtime: return True
    return False

def flatten2(list_of_lists):
    """ flatten an irregular 2-level list """
    return [ item for list1 in list_of_lists for item in list1 ]

def binPhases(phaselist,numbins):
    """ non-overlapping phase bins, phase is assumed between >=0.0 and <1.0 """
    binlist = [0]*numbins
    for phase in phaselist:
        binlist[int(phase*numbins)] += 1
    ## approximate value taken for mouse_respirationT
    return array(binlist) / (float(mouse_respirationT)/numbins)

def bin_phases(mit_phases_trials,mit_spiketimes, \
        airstarttime,airendtime,odorstarttime,odorendtime):
    mit_phases = [ mit_phases_trials[trialnum].flatten() \
        for trialnum in range(len(mit_spiketimes)) ]
    ## mit_phases and mit_spiketimes are nested lists [trialnum,phasenum]
    ## flatten them before a where, else empty results,
    ## array(...).flatten() won't work as spiketimes/phases are irregular list of lists,
    ## so can't convert to array.
    ## IMPORTANT to convert to array() else where() returns empty ()
    mit_phases = array(flatten2(mit_phases))
    mit_spiketimes = array(flatten2(mit_spiketimes))
    ## Note the () & (), the brackets are a must and so is the bitwise &
    air_indices = where( (mit_spiketimes>airstarttime) & (mit_spiketimes<airendtime) )
    if len(air_indices)>0:
        mit_air = mit_phases[air_indices[0]]
        ## sorting is needed by the plotOverlappingBins() below to bin the phases
        ## all resp cycles over different trials clubbed together
        mit_air = sort(mit_air)
    else: mit_air = ()
    ## Note the () & (), the brackets are a must and so is the bitwise &
    odor_indices = where( (mit_spiketimes>odorstarttime) & (mit_spiketimes<odorendtime) )
    if len(odor_indices)>0:
        mit_odor = mit_phases[odor_indices[0]]
        ## sorting is needed by the plotOverlappingBins() below to bin the phases
        ## all resp cycles over different trials clubbed together
        mit_odor = sort(mit_odor)
    else: mit_odor = ()

    ## many cycles clubbed together! trialnum * (respnums cycles per trial),
    ## (respnums cycles per trial) not known, different for each trial
    ## But doesn't matter for calculating pearson correlation as it is normalized.
    mit_air_binned = binPhases( mit_air, numbins )
    mit_odor_binned = binPhases( mit_odor, numbins )

    return mit_air, mit_odor, mit_air_binned, mit_odor_binned

def bin_phases_withSD(mit_phases_trials,mit_spiketimes_trials,starttime,endtime):
    ## mit_phases_trials and mit_spiketimes_trials are nested lists [trialnum,phasenum]
    num_trials = len(mit_phases_trials)
    num_resps = 1
    binlist = [0]*numbins
    binlist_resps = [binlist]
    toolongISIs = []
    for tnum,spiketimes in enumerate(mit_spiketimes_trials):
        oldphase = 0
        oldspiketime = starttime
        for snum,spiketime in enumerate(spiketimes):
            ## for each trial, only consider spikes between starttime and endtime for each trial
            if spiketime<starttime: continue
            elif spiketime>endtime: break
            phase = mit_phases_trials[tnum][snum]
            ## new resp cycle starts if phase goes back OR
            ## mouse resp period passes between two spikes -- bad method -- not using
            longtimenospike = (spiketime-oldspiketime)>mouse_respirationT
            if phase<oldphase:
                num_resps += 1
                binlist = [0]*numbins
                binlist_resps.append(binlist)
                ## But no spike for a resp period is a bad method
                ## esp if firing rate is close to zero
                ## diagnoze how many such cases are there
                ## fair number of cases with inter-spike-interval > 0.5s, hence use with caution:
                ## low firing air/odor response will be overestimated.
                if longtimenospike: toolongISIs.append(spiketime-oldspiketime)
            oldspiketime = spiketime
            oldphase = phase
            binlist_resps[-1][int(phase*numbins)] += 1
    ## This rate is only approximate as mouse_respirationT is variable over resps and mice
    ## Calculate firing rate by bindt -- approximate as bin is a phasebin!
    binlist_resps = array(binlist_resps)/(mouse_respirationT/float(numbins))
    ## num_resps are the number of respiratory cycle caught.
    ## It is way off since there are no spikes in some cycles
    ## Use approximate 0.5s period and calculate mean
    num_resps_estimate = float((endtime-starttime)/mouse_respirationT)*float(num_trials)
    binnedfratemean = sum(binlist_resps,axis=0) / num_resps_estimate
    ## SD is also a lower estimate as there are many missed cycles!!
    ## Add zero filled arrays for all those cycles that are missed.
    for i in range(int(num_resps_estimate-num_resps)):
        binlist_resps = append(binlist_resps,zeros(5))
    binnedfrateSD = std(binlist_resps,axis=0)
    #print 'ISIs greater than',mouse_respirationT,'s are',\
    #    len(toolongISIs),'with mean',mean(toolongISIs),\
    #    's, num of resp cycles =',num_resps
    return binnedfratemean,binnedfrateSD,num_resps

## binned correlations
def sis_corrs_frates((splpairnum,splodornum,splax)):

    air_corr_list = []
    odor_corr_list = []
    air_frate_list = []
    odor_frate_list = []
    air_var_list = []
    odor_var_list = []
    neg_corr_nums, neg_corr_high_frate_nums = 0,0
    frate_cutoff = 1 # 1 Hz
    ## plot the deviant responses!
    numsubfigs_rows = 12
    numsubfigs_cols = 12
    fig = figure(facecolor='w')
    figtext(0.1,0.94,'Responses with neg corr and both mean frates > '+\
        str(frate_cutoff)+'Hz',fontsize=20)
    delaxes() # delete the main axes to accomodate subpanels

    ## doing for phases,not spiketime -- looking at binned phase decorr
    for pairnum,pair_response in enumerate(responses):
        air_corr_list.append([])
        odor_corr_list.append([])
        air_frate_list.append([])
        odor_frate_list.append([])
        air_var_list.append([])
        odor_var_list.append([])
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

            ## bin air and odor phases for the two sisters
            ## not divided by num resps, needs to be normalized
            mit0_phases_trials = phases[pairnum,0][odornum]
            mit1_phases_trials = phases[pairnum,1][odornum]
            mit0_air, mit0_odor, mit0_air_binned, mit0_odor_binned = \
                bin_phases(mit0_phases_trials,mit0_spiketimes,\
                    airstarttime,airendtime,odorstarttime,odorendtime)
            mit1_air, mit1_odor, mit1_air_binned, mit1_odor_binned = \
                bin_phases(mit1_phases_trials,mit1_spiketimes,\
                    airstarttime,airendtime,odorstarttime,odorendtime)

            ## diagnostic plots to check if resp info can be obtained -- No, can't.
            if pairnum==splpairnum and odornum==splodornum:
                print 'Binning for our special example'
                numtrials = float(len(mit0_spiketimes))
                for tnum,spiketimes in enumerate(mit0_spiketimes):
                    fig2 = figure()
                    ax2 = fig2.add_subplot(111)
                    ax2.plot(spiketimes,mit0_phases_trials[tnum].flatten(),color='r',marker='x')
                    ax2.plot(mit1_spiketimes[tnum],mit1_phases_trials[tnum].flatten(),color='b',marker='x')
                    axes_labels(ax2,'spike time (s)','spike phase',adjustpos=False,fontsize=24)
                    #for snum,spiketime in enumerate(spiketimes):
                    #    phase = mit0_phases_trials[tnum][snum]
                    #    print tnum,snum,spiketime,phase

            ## above binning does not return mean & SD (no num resps info)
            ## do not use the post-odor air periods, as they are tainted with odor!
            mit0_air_binned, mit0_air_binned_SD, numresps_mit0_air = \
                bin_phases_withSD(mit0_phases_trials,mit0_responses,airstarttime,airendtime)
            mit0_odor_binned, mit0_odor_binned_SD, numresps_mit0_odor = \
                bin_phases_withSD(mit0_phases_trials,mit0_responses,odorstarttime,odorendtime)
            mit1_air_binned, mit1_air_binned_SD, numresps_mit1_air = \
                bin_phases_withSD(mit1_phases_trials,mit1_responses,airstarttime,airendtime)
            mit1_odor_binned, mit1_odor_binned_SD, numresps_mit1_odor = \
                bin_phases_withSD(mit1_phases_trials,mit1_responses,odorstarttime,odorendtime)

            ## mean air and odor frate for each pairnum: avg over the pair
            air_frate = float( len(mit0_air)+len(mit1_air) ) / (airendtime-airstarttime) / 2.0
            odor_frate = float( len(mit0_odor)+len(mit1_odor) ) / (odorendtime-odorstarttime) / 2.0
            ## Be careful with this binned mean frate, it is 5x smaller due to mean over 5 bins.
            #air_frate = (mean(mit0_air_binned)+mean(mit1_air_binned))/2.0
            #odor_frate = (mean(mit0_odor_binned)+mean(mit1_odor_binned))/2.0
            
            ## variance in air / odor frate
            air_var = std( flatten2( [mit0_air_binned, mit1_air_binned] ) )
            odor_var = std( flatten2( [mit0_odor_binned, mit1_odor_binned] ) )

            ## BINNED cross correlation
            air_corr = scipy.stats.pearsonr(mit0_air_binned,mit1_air_binned)[0]
            odor_corr = scipy.stats.pearsonr(mit0_odor_binned,mit1_odor_binned)[0]

            air_corr_list[-1].append(air_corr)
            odor_corr_list[-1].append(odor_corr)
            air_frate_list[-1].append(air_frate)
            odor_frate_list[-1].append(odor_frate)
            air_var_list[-1].append(air_var)
            odor_var_list[-1].append(odor_var)
            
            if odor_corr<0.0:
                neg_corr_nums += 1
                if mean(mit0_odor_binned)>frate_cutoff and mean(mit1_odor_binned)>frate_cutoff:
                    print 'sis'+str(pairnum)+'odor'+str(odornum)+' corr='+str(odor_corr)+\
                        ' frate='+str(odor_frate)+' var='+str(odor_var)
                    neg_corr_high_frate_nums += 1
                    ax = fig.add_subplot(numsubfigs_rows,numsubfigs_cols,neg_corr_high_frate_nums)
                    ax.set_title(str(pairnum)+'-'+str(odornum)+'; '+'%1.2f'%(odor_corr))
                    ax.plot(mit0_odor_binned,'r-x')
                    ax.plot(mit0_air_binned,color=(1,0.6,0.6))
                    ax.plot(mit1_odor_binned,'g-o')
                    ax.plot(mit1_air_binned,color=(0.6,1.0,0.6))
                    axes_off(ax,x=True,y=False)
                    
            if pairnum==splpairnum and odornum==splodornum:
                print 'Plotting for our special example'
                xbinmids = arange(1.0/numbins/2.0, 1.0, 1.0/numbins)
                ## do not plot standard error, as # of resp cycles caught
                ## is very different between the two mitrals (see above).
                ## Even the standard deviation is calculated not so cleanly,
                ## so don't plot the error.
                splax.errorbar(x=xbinmids,y=mit1_odor_binned,\
                    color='r',marker='o',ms=marker_size,linewidth=linewidth)
                splax.errorbar(x=xbinmids,y=mit0_odor_binned,\
                    color='b',marker='s',ms=marker_size,linewidth=linewidth,linestyle='dashed')
                splax.set_xlim(0,1)
                splax.set_xticks([0,1])
                splax.set_xticklabels(['0','2$\pi$'])
                xmin,xmax,ymin,ymax = beautify_plot(splax,xticksposn='bottom',yticksposn='left')
                splax.set_yticks([0,ymax])
                axes_labels(splax,'phase','firing rate (Hz)',adjustpos=False,fontsize=label_fontsize)

    print "Number of sis-odor combos with neg corr =",neg_corr_nums
    print "Number of sis-odor combos with neg corr and "\
        "both mits having frates>"+str(frate_cutoff)+"Hz mean =",neg_corr_high_frate_nums
    print "Mean firing rate for odor is",mean(odor_frate_list)
    print "Mean firing rate for air is",mean(air_frate_list)
    return air_corr_list,odor_corr_list,air_frate_list,odor_frate_list,air_var_list,odor_var_list

if __name__ == "__main__":

    fig = figure(figsize=(columnwidth,linfig_height/2.0),dpi=300,facecolor='w') # 'none' is transparent
    ax1 = fig.add_subplot(1,2,1)
    #ax1.text(-0.3, 0.9, 'a', fontweight='bold', fontsize=label_fontsize, transform=ax1.transAxes)
    ax2 = fig.add_subplot(1,2,2)

    ## old example of fig8 feb 2013 was 2-2, May 2013 it is 13-9.
    air_corr_list,odor_corr_list,air_frate_list,odor_frate_list,air_var_list,odor_var_list \
        = sis_corrs_frates((13,9,ax1))

    ax2.hist(flatten2(air_corr_list),10,range=(-1.0,1.0),normed=True,histtype='step',\
        linewidth=linewidth,color='b',ls='dotted')
    ax2.hist(flatten2(odor_corr_list),10,range=(-1.0,1.0),normed=True,histtype='step',\
        linewidth=linewidth,color='r',ls='solid')
    ## ax.transAxes ensures relative to axes size, rather than to data units.
    #ax2.text(-0.3, 0.9, 'b', fontweight='bold', fontsize=label_fontsize, transform=ax2.transAxes)
    beautify_plot(ax2,x0min=False,y0min=True,xticksposn='bottom',yticksposn='left')
    axes_labels(ax2,'correlation','prob. density',adjustpos=False,fontsize=label_fontsize)
    ax2.set_xticks([-1,0,1])
    fig.tight_layout()
    fig_clip_off(fig)
    fig.savefig('../../figures/decorr/ashesh_example.svg', dpi=fig.dpi)
    fig.savefig('../../figures/decorr/ashesh_example.png', dpi=fig.dpi)

    ## mean firing rate vs corr (is there negative corr for low firing rates only - looks not!)
    figure(facecolor='w')
    title('Air mean frate (Hz) vs corr')
    plot(array(air_corr_list).flatten(),array(air_frate_list).flatten(),'x')
    figure(facecolor='w')
    title('Odor mean frate (Hz) vs corr')
    plot(array(odor_corr_list).flatten(),array(odor_frate_list).flatten(),'x')

    ## Approx. SD vs corr (is there negative corr for flat (low variance) firing rates only - looks not!)
    figure(facecolor='w')
    title('Air approx.SD (Hz) vs corr')
    plot(array(air_corr_list).flatten(),array(air_var_list).flatten(),'x')
    figure(facecolor='w')
    title('Odor approx.SD (Hz) vs corr')
    plot(array(odor_corr_list).flatten(),array(odor_var_list).flatten(),'x')

    show()
