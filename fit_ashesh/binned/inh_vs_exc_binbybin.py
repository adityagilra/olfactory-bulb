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

cutoff = 2
print "Cutoff / threshold =",cutoff,"spikes per bin."
only_air_responsive = True
print "Do the analysis for only air responsive cells i.e. all air bins>cutoff?", only_air_responsive

exc_air = 0
no_air = 0
mixed_air = 0

figure()
title('Mean air responses of mitrals')
air_responses_permitral_full = []
air_responsive_cells_list = []
for mitralnum,mitral_responses in enumerate(air_responses):
    mitral_response = zeros((len(mitral_responses[0]),))
    for trial_response in mitral_responses:
        mitral_response += array(trial_response)
    mean_air_response = mitral_response/len(mitral_responses)
    if not only_air_responsive: plot(mean_air_response)
    # append mean air response in all trials
    air_responses_permitral_full.append(mean_air_response)
    exc = False
    noresp = False
    for binnum, airresponse_inbin in enumerate(mean_air_response):
        if airresponse_inbin > cutoff:
            exc = True
        else:
            noresp = True
    if exc and noresp:
        mixed_air += 1
    elif exc:
        exc_air+=1
        air_responsive_cells_list.append(mitralnum)
        if only_air_responsive: plot(mean_air_response)
    else:
        no_air += 1
print "Air responses: full-exc =",exc_air,", no response =",no_air,\
    ", mixed response =",mixed_air,"out of",mixed_air+exc_air+no_air

excitatory_list = []
inhibitory_list = []
mixed_list = []
noresp_list = []

excfig = figure()
inhfig = figure()
mixedfig = figure()
excax = excfig.add_subplot(111)
inhax = inhfig.add_subplot(111)
mixedax = mixedfig.add_subplot(111)

if only_air_responsive: cells_list = air_responsive_cells_list
else: cells_list = range(len(odor_responses))

for mitralnum in air_responsive_cells_list:
    mitral_responses = odor_responses[mitralnum]
    airspikes = air_responses_permitral_full[mitralnum]
    for odornum,odorspikes in enumerate(mitral_responses):
        exc = False
        inh = False
        for binnum, odorspikesinbin in enumerate(odorspikes):
            if odorspikesinbin > (airspikes[binnum]+cutoff):
                exc = True
            elif odorspikesinbin < (airspikes[binnum]-cutoff):
                inh = True
        ### Plot only differential responses
        response = odor_responses[mitralnum][odornum] - airspikes
        if exc and inh: # mixed response
            mixed_list.append(response)
            mixedax.plot(response)
        elif exc: # purely excitatory in every bin
            excitatory_list.append(response)
            excax.plot(response)
        elif inh: # purely inhibitory in every bin
            inhibitory_list.append(response)
            inhax.plot(response)
        else:
            noresp_list.append(response)

excax.set_title(str(len(excitatory_list))+\
    ' Pure-Exc differentials (each bin\'s odor spikes >= resp spikes+-'+str(cutoff)+')')
inhax.set_title(str(len(inhibitory_list))+\
    ' Pure-Inh differentials (each bin\'s odor spikes <= resp spikes+-'+str(cutoff)+')')
mixedax.set_title(str(len(mixed_list))+\
    ' Mixed differentials (each bin\'s odor spikes >or< resp spikes+-'+str(cutoff)+')')
print "No response in any bin =",len(noresp_list)

show()
