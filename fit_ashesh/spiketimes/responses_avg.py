# -*- coding: utf-8 -*-

## USAGE: python2.6 responses_avg.py

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
from data_utils import * # for plotBins()

samplingf = 10000.0 # Hz float/double, sampling frequency: spike indices are sampled at 10kHz
mouse_resp_cycle = 0.5 # s, approx.

## copied recursive_flatten() from
## http://en.wikibooks.org/wiki/Python_Programming/Tips_and_Tricks
def recursive_flatten(seq, list = None):
    """flatten(seq, list = None) -> list 
    Return a flat version of the iterator `seq` appended to `list`
    """
    if list is None:
        list = []
    try:                                    # Can `seq` be iterated over?
        for item in seq:                    # If so then iterate over `seq`
            recursive_flatten(item, list)   # and make the same check on each item.
    except TypeError:                       # If seq isn't iterable
        list.append(seq)                    # append it to the new list.
    return list

def plot_avg_response(odor_responses, titlestr):
    ## total_resps only has those respiration cycles that have at least one spike
    ## see convert_air_odor_responses()
    total_resps_calc = 0
    for pairnum in range(len(odor_responses)):
        for odornum in range(42):
            for mitnum in range(2):
                total_resps_calc += len(odor_responses[pairnum][odornum][mitnum])
    ## Assuming mouse_resp_cycle for mouse respiration, the number of resp cycles is different:
    ## Assuming mouse_resp_cycle for mouse respiration, the number of resp cycles is different:
    if 'air' in titlestr:
        total_resps_est = len(odor_responses)*42*2*(10.0-0.0)/mouse_resp_cycle
    else:
        total_resps_est = len(odor_responses)*42*2*(15.0-10.0)/mouse_resp_cycle
    print "Number of resp cycles with next phase phase < old spike phase is", total_resps_calc
    print "Number of resp cycles dividing time by mouse resp period is", total_resps_est
    merged_responses = recursive_flatten(odor_responses)
    numbins = 20
    ## phase response
    binned_response = plotBins(merged_responses, numbins, runtime=1.0, settletime=0.0)
    avg_response = array(binned_response)/float(total_resps_est)/mouse_resp_cycle
    plot(avg_response,label=titlestr)

def plot_late_responses(odor_responses):
    total_resps = 0
    total_combos = 0
    numbins = 10
    inh_first_responses = zeros(numbins)
    binned_responses = []
    for pairnum in range(len(odor_responses)):
        for odornum in range(42):
            for mitnum in range(2):
                mit_responses = recursive_flatten(odor_responses[pairnum][odornum][mitnum])
                binned_response = plotBins(mit_responses, numbins, runtime=1.0, settletime=0.0)
                maxresponse = max(binned_response)
                if binned_response.index(maxresponse)/float(numbins) >= 0.3:
                    inh_first_responses += binned_response
                    numresps = len(odor_responses[pairnum][odornum][mitnum])
                    total_resps += numresps
                    total_combos += 1
                    #binned_responses.append(array(binned_response)/float(numresps)/mouse_resp_cycle)
                    for resptrain in odor_responses[pairnum][odornum][mitnum]:
                        respflat = recursive_flatten(resptrain)
                        binned_response_oneresp = \
                            plotBins(respflat, numbins, runtime=1.0, settletime=0.0)
                        binned_responses.append(array(binned_response_oneresp)/mouse_resp_cycle)
    ## phase response
    print 'Total number of mit-odor combos is',total_combos
    print 'Total number of respirations is',total_resps
    avg_response = array(inh_first_responses)/float(total_resps)/mouse_resp_cycle
    fig = figure()
    ax = fig.add_subplot(111)
    title('inhibition first responses',fontsize=24)
    plot(avg_response)
    axes_labels(ax,'resp phase bin #','firing rate (Hz)')
    fig = figure()
    ax = fig.add_subplot(111)
    title('Avg response (w/ std error)',fontsize=24)
    avg_response = mean(binned_responses,axis=0)
    ## standard error of the mean, not standard deviation of trials
    response_std = std(binned_responses,axis=0)/sqrt(len(binned_responses))
    errorbar(x=range(numbins),y=avg_response,yerr=response_std)
    axes_labels(ax,'resp phase bin #','firing rate (Hz)')

def convert_air_odor_responses(trialslist, spiketimeslist):
    phaseslist = []
    air_phases = [[]]
    odor_phases = [[]]
    for trialnum,spikephases in enumerate(trialslist):
        ## to allow responses to be a numpy array,
        ## spikephases is [[spiketime1],[spiketime2],...]
        spikephases = spikephases.flatten()
        oldphase = 0.0
        spiketimes = spiketimeslist[trialnum]
        for spikenum,phase in enumerate(spikephases):
            spiketime = spiketimes[spikenum]
            ## air delivery period
            if spiketime < 10.0: # s
                ######## Only those cycle whose startphase < oldphase get selected:
                ######## empty cycles discarded
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
    return (air_phases,odor_phases)

def plot_responses():
    odor_responses = []
    air_phases = []
    odor_phases = []
    for pairnum in range(len(responses)):
        odor_responses.append([])
        air_phases.append([])
        odor_phases.append([])
        for odornum in range(42):
            ## convert to spiketimes from spike indices at sampling frequency
            mit0_responses = responses[pairnum,0][odornum]/samplingf
            mit1_responses = responses[pairnum,1][odornum]/samplingf
            mit0_spiketimes = [ mit0_responses[trialnum].flatten() \
                for trialnum in range(len(mit0_responses)) ]
            mit1_spiketimes = [ mit1_responses[trialnum].flatten() \
                for trialnum in range(len(mit0_responses)) ]
            odor_responses[-1].append( (mit0_spiketimes,mit1_spiketimes) )
            
            ## convert firing phases (0 to 1 along respiration) into 2D array
            mit0_phaseslist = convert_air_odor_responses(phases[pairnum,0][odornum], mit0_spiketimes)
            mit1_phaseslist = convert_air_odor_responses(phases[pairnum,1][odornum], mit1_spiketimes)
            air_phases[-1].append( (mit0_phaseslist[0],mit1_phaseslist[0]) )
            odor_phases[-1].append( (mit0_phaseslist[1],mit1_phaseslist[1]) )
    figure()
    title('Avg phase response')
    ## air/odor_phases[pairnum][odornum][respnum][phasenum]
    plot_avg_response(air_phases,'air')
    plot_avg_response(odor_phases,'odor')
    legend()
    plot_late_responses(odor_phases)

if __name__ == "__main__":
    plot_responses()
    show()
