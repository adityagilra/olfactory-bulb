import scipy.io
air_responses = scipy.io.loadmat('aptc.mat',struct_as_record=True)
odor_responses = scipy.io.loadmat('optc.mat',struct_as_record=True)
### these are numpy structured arrays, like dicts
### obtain the actual data
air_responses = air_responses['aptc'][0]
odor_responses = odor_responses['optc'][0]

### Now air_response and odor_response are mitralnum x odornum x bin = 133 x 42 x 12
### For air, 42 are just repeats: each corresponding to the odor given
### For some mitrals not all 42 odors / air-repeats are present.

from pylab import *

def mysum(response):
    val = 0
    for bin in response:
        val += bin
    return val

figure()
title('air responses')
# mean across all trials of the sum across all bins!
# somehow sum() cribs that axis=1 and axis=2 are out of bounds!
# hence sum the hard way
air_responses_permitral = []
air_responses_permitral_full = []
for mitralnum,mitral_responses in enumerate(air_responses):
    mitral_response_sum = 0
    mitral_response = zeros((len(mitral_responses[0]),))
    for trial_response in mitral_responses:
        mitral_response_sum += mysum(trial_response)
        mitral_response += array(trial_response)
    # append mean of sum_of_spikes in all trials
    air_responses_permitral.append(mitral_response_sum/len(mitral_responses))
    plot(mitral_response/len(mitral_responses))
    air_responses_permitral_full.append(mitral_response/len(mitral_responses))
# sum over all bins for each mitral for each odor
odor_responses_sum = []
for mitral_responses in odor_responses:
    mitral_response_sum = []
    for odor_response in mitral_responses:
        mitral_response_sum.append( mysum(odor_response) )
    odor_responses_sum.append(mitral_response_sum)

##### plot air and odor responses for first 10 mitrals
#for mitnum in range(10):
    #figure()
    #print mitnum,air_responses_permitral[mitnum]
    #for odornum,odorresponse in enumerate(odor_responses[mitnum]):
        #plot(odorresponse)
        #print " |--",odor_responses_sum[mitnum][odornum],odor_responses
    #plot(air_responses_permitral_full[mitnum], 'ko-')

excitatory_list = []
inhibitory_list = []
noresp_list = []

excfig = figure()
inhfig = figure()
norespfig = figure()
excax = excfig.add_subplot(111)
inhax = inhfig.add_subplot(111)
norespax = norespfig.add_subplot(111)

cutoff = 1
for mitralnum,mitral_responses in enumerate(odor_responses_sum):
    airspikesnum = air_responses_permitral[mitralnum]
    for odornum,odorspikesnum in enumerate(mitral_responses):
        if odorspikesnum > airspikesnum + cutoff:
            response = odor_responses[mitralnum][odornum]
            excitatory_list.append(response)
            excax.plot(response)
        elif odorspikesnum < airspikesnum - cutoff:
            response = odor_responses[mitralnum][odornum]
            inhibitory_list.append(response)
            inhax.plot(response)
        else:
            response = odor_responses[mitralnum][odornum]
            noresp_list.append(response)
            norespax.plot(response)

excax.set_title(str(len(excitatory_list))+' Excitatory responses (# odor spikes > # resp spikes + '+str(cutoff)+')')
inhax.set_title(str(len(inhibitory_list))+' Inhibitory responses (# odor spikes < # resp spikes) - '+str(cutoff)+'')
norespax.set_title(str(len(noresp_list))+' No responses (# odor spikes ~= # resp spikes +- '+str(cutoff)+')')

show()
