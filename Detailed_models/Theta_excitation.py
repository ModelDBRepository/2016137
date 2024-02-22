import matplotlib.pyplot as plt
import numpy
from importlib import reload
import sys
import os
import importlib
import types

original_packs = sys.modules.copy()

import math

from neuron import h
import neuron

import json
import collections


#Modellke

def modellke(model_name, num_sim, synapse_number, noise, parameters, distance, width):

        neuron.load_mechanisms('./modfiles/')
        h.load_file('./HOC/' + model_name +'_skinner.hoc')

        v = parameters

        for sec in h.somatic:
            sec.gbar_Ikdrf = v[0]
            sec.gbar_Ikdrs = v[1]
            sec.gbar_Ika = v[2]
            sec.gkhbar_Ih = v[3]
            sec.g_passsd = v[6]
            sec.gkl_Kleaksd = v[8]
            sec.gna_Nasoma = v[10]   
            sec.gbar_IM = v[13]   
        for sec in h.basal:
            sec.gkhbar_Ih = v[3]
            sec.gbar_Ikdrf = v[0]
            sec.gbar_Ikdrs = v[1]
            sec.gbar_Ika = v[2]
            sec.gbar_cat = v[4]
            sec.gcalbar_cal = 10* v[4]
            sec.gkbar_kca = v[5]
            sec.g_passsd = v[6]
            sec.gkl_Kleaksd = v[8]
            sec.gna_Nadend = v[11]
            sec.gbar_IM = v[13]  
        for sec in h.axonal:  
            sec.gbar_Ikdrfaxon = v[0]
            sec.gbar_Ikdrsaxon = v[1]
            sec.g_passaxon = v[7]
            sec.gkl_Kleakaxon = v[9]
            sec.gna_Naaxon = v[12]

        h.finitialize(-80)
        h.fcurrent()


        """ Stting up synapses """

        stim_file = "stim_PSP_attenuation_test_hs.json" # shape and amplitude of the EPSC originated from the HippoUnit's PSP attenuation test. If you want to use this script with an OLM cell make sure to change this file accordingly
        with open(stim_file, 'r') as f:
            config = json.load(f, object_pairs_hook=collections.OrderedDict)

        distances = distance
        tolerance = width
        dist_range = [min(distances) - tolerance, max(distances)

        tau1 = config['tau_rise']
        tau2 = config['tau_decay']
        EPSC_amp = config['EPSC_amplitude']
        basal= []
        for sec in h.basal:
            basal.append(sec)
        basal = list(basal)

        num = synapse_number
        num_gaba = synapse_number_gaba

        AMPA_weight = 0.0001  #an approximate value that originates from HippoUnit's PSP attenuation test

        we_want_all = False #if we want to use all possible dendritic locations to recieve a synapse - typically not used


        "Getting the locations of the synapses"

        def random_locations_excit(we_want_all, distances, tolerance, basal, num):
            all_dend_loc_0_100 = []
            all_dend_loc_100_200 = []
            all_dend_loc_distal = []

            all_location_distances = {}

            locations = []
            location_distances = {}

            kumm_length_list = []
            kumm_length = 0
            num_of_secs = 0
            seed = 1

            for sec in h.basal:
                num_of_secs += sec.nseg
                kumm_length += sec.L
                kumm_length_list.append(kumm_length)

            if we_want_all == True:
                for sec in h.basal:
                    for seg in sec:
                        if h.distance(seg.x, sec=sec) <= 100:
                                all_dend_loc_0_100.append([sec.name, seg.x])
                                all_locations_distances[sec.name(), seg.x] = h.distance(seg.x, sec=sec)
                        if h.distance(seg.x, sec=sec) > 100 and h.distance(seg.x, sec=sec) <= 200:
                                all_dend_loc_100_200.append([sec.name, seg.x])
                                all_locations_distances[sec.name(), seg.x] = h.distance(seg.x, sec=sec)
                        if h.distance(seg.x, sec=sec) > 200:
                                all_dend_loc_dist.append([sec.name, seg.x])
                                all_locations_distances[sec.name(), seg.x] = h.distance(seg.x, sec=sec)

            else:   # the chance of a dendritic section to get selected is proportionate to the length of the dendritic section

                norm_kumm_length_list = [i/kumm_length_list[-1] for i in kumm_length_list]

                import random

                _num_ = num  # _num_ will be changed
                num_iterations = 0

                while len(locations) < num: # and num_iterations < 50 
                    random.seed(seed)
                    rand_list = [random.random() for j in range(_num_)]

                    for rand in rand_list:

                        for i in range(len(norm_kumm_length_list)):
                            if rand <= norm_kumm_length_list[i] and (rand > norm_kumm_length_list[i-1] or i==0):

                                seg_loc = (rand - norm_kumm_length_list[i-1]) / (norm_kumm_length_list[i] - norm_kumm_length_list[i-1])

                                segs = [seg.x for seg in basal[i]]
                                d_seg = [abs(seg.x - seg_loc) for seg in basal[i]]
                                min_d_seg = numpy.argmin(d_seg)
                                segment = segs[min_d_seg]
                                h.distance(sec=h.soma)
                                h('access ' + basal[i].name())
                                if h.distance(segment) >= dist_range[0] and h.distance(segment) < dist_range[1]:
                                    locations.append([basal[i].name(), segment])
                                    location_distances[basal[i].name(), segment] = h.distance(segment)
                    _num_ = num - len(locations)
                    seed += 10
                    num_iterations += 1

            return locations, location_distances


        ''' Excitation '''

        def frange(start, stop, step):
            """
            Generates range of real values.

            :param start: beginning of range
            :param stop: end of range
            :param step: step size between values

            """
            r = start
            while r < stop:
                yield r
                r += step

        def create_amp_vec(amplitude):
        
            time_np =numpy.array([])
            t_gen = frange(0, 40000 + 0.025, 0.025) #should be more than the length of the recording
            for n in t_gen:
                time_np = numpy.append(time_np, n)

            amp_np = numpy.array([])
            
            for t in time_np:
                    amp_np = numpy.append(amp_np, amplitude)

            return h.Vector(amp_np), h.Vector(time_np)


        def create_dc_vec(rise_start, decay_start, slope, baseline_dc, back_to_base):
            time_np =numpy.array([])
            t_gen = frange(0, 40000 + 0.025, 0.025) #should be more than the length of the recording
            for n in t_gen:
                time_np = numpy.append(time_np, n)
            dc_vec = numpy.array([])

            for t in time_np:
                if t <= rise_start:
                    dc_vec = numpy.append(dc_vec, baseline_dc)
                if t > rise_start and t<= decay_start:
                    dc_vec = numpy.append(dc_vec, (slope*(t-rise_start)+baseline_dc))
                if t > decay_start and t <= back_to_base:
                    dc_vec = numpy.append(dc_vec, (-slope*(t-decay_start)+((slope*(decay_start-rise_start))+baseline_dc)))
                if t > back_to_base:
                    dc_vec = numpy.append(dc_vec, baseline_dc)

            plt.figure()
            plt.plot(time_np,rate_function_2)
            dc_vec=h.Vector(dc_vec)
            time_np = h.Vector(time_np)
            return dc_vec, time_np


        synapse_lists = {}

        def set_ampa_nmda_multiple_loc_theta(dend_loc, AMPA_weight, dc_vec, amp_vec, time_dc, time_amp):

            for i in range(len(dend_loc)):

                ndend, xloc = dend_loc[i]

                exec("dend=h." + ndend)

                synapse_lists['ampa_list'][i] = h.Exp2Syn(xloc, sec=dend)
                synapse_lists['ampa_list'][i].tau1 = tau1
                synapse_lists['ampa_list'][i].tau2 = tau2


                synapse_lists['ns_list'][i] = h.NetStimVariable()
                synapse_lists['ns_list'][i].dur = 30000
                synapse_lists['ns_list'][i].delay = 500
                synapse_lists['ns_list'][i].noise = 1
                dc_vec.play(synapse_lists['ns_list'][i]._ref_dc_comp, time_dc, True)
                amp_vec.play(synapse_lists['ns_list'][i]._ref_amp, time_amp, True)
                synapse_lists['ns_list'][i].freq = 0.007
                synapse_lists['ns_list'][i].phase = 0
                synapse_lists['ns_list'][i].minFrequency = 1e-06
                synapse_lists['ampa_nc_list'][i] = h.NetCon(synapse_lists['ns_list'][i], synapse_lists['ampa_list'][i])
                synapse_lists['ampa_nc_list'][i].weight[0] = AMPA_weight

            return synapse_lists

        def activate_theta_stimuli(dend_loc, AMPA_weight, dc_vec, amp_vec, time_dc, time_amp):
            """From PathwayInteractionTest"""


            synapse_lists.update({'ampa_list' : [None] * len(dend_loc),
                                'ampa_nc_list' : [None]*len(dend_loc),
                                'ns_list' : [None]*len(dend_loc) 
                                })

            set_ampa_nmda_multiple_loc_theta(dend_loc, AMPA_weight, dc_vec, amp_vec, time_dc, time_amp)

        def run_simulation(dend_loc, recording_loc):

            (rec_ndend, xloc), distance = recording_loc

            exec("dendrite=h." + rec_ndend)

            exec("sect_loc=h.soma" + "("+str(0.5)+")")


            # initiate recording
            rec_t = h.Vector()
            rec_t.record(h._ref_t)

            rec_v = h.Vector()
            rec_v.record(sect_loc._ref_v)

            rec_v_dend = h.Vector()
            rec_v_dend.record(dendrite(xloc)._ref_v)
            
            
            v_stim = []
            dend_loc_rec =[]
            
            for i in range(len(dend_loc)):
                exec("dend_loc_rec.append(h." + str(dend_loc[i][0])+"("+str(dend_loc[i][1])+"))")
                v_stim.append(h.Vector())

            for i in range(len(dend_loc_rec)):
                v_stim[i].record(dend_loc_rec[i]._ref_v)


            h.stdinit()
            dt = 0.025
            h.dt = dt
            h.steps_per_ms = 1 / dt
            h.v_init = -60
            h.celsius = 24
            h.init()
            h.tstop = 30500
            h.run()

            t = numpy.array(rec_t)
            v = numpy.array(rec_v)
            v_dend = numpy.array(rec_v_dend)

            
            v_stim_locs = collections.OrderedDict()
            for i in range(len(dend_loc)):
                loc_key = (dend_loc[i][0],dend_loc[i][1])
                v_stim_locs[loc_key] = numpy.array(v_stim[i])

            return t, v, v_dend, v_stim_locs

        dend_loc, dend_loc_distance = random_locations_excit(we_want_all, distances, tolerance, basal, num)

        dc_vec, time_dc = create_dc_vec(0, 1000, 0, 0.011, 40000)  #rise_start, decay_start, slope, baseline_dc, back_to_base - slope is 0 so there will be no rise and decay only baseline, but decay start have to be bigger than rise start
        amp_vec, time_amp = create_amp_vec(0.01)

        recording_loc = min(dend_loc_distance.items(), key=lambda kv : abs(kv[1] - distances[0]))
        activate_theta_stimuli(dend_loc, AMPA_weight, dc_vec, amp_vec, time_dc, time_amp)
        t, v, v_dend, v_stim_locs = run_simulation(dend_loc, recording_loc)

        """ Saving the data for further evalutations """

        somatic_trace = numpy.column_stack([t, v])
        if not os.path.exists('./theta_things/' + model_name + '/'):
            os.makedirs('./theta_things/' + model_name + '/')
        numpy.savetxt('./theta_things/' + model_name + '/'+ model_name + '_' + str(num_sim) + '_' + str(distance) + '_noise_' + str(noise) +'_soma_rec_theta_excit_5_percent_long.dat', somatic_trace)

        """ Create plots to see theta responses"""

        plt.figure()
        plt.plot(t,v)


"PARAMETERS FOR THE MODEL AND THE SIMULATION"


model_name = 'HS_0731'

"the parameter set is coming from the csv file corresponding the given model name. The csv files can be found in the Parameters_from_neuroptimus folder"

parameters = [0.103144490247811, 0.000106775331027, 0.006819469725367, 2.12635184953167E-05, 0.008998494911798, 0.000212037977726, 6.0890455155313E-06, 8.917984284895E-05, 1.10037600392373E-05, 0.000126283068342, 0.020327589012284, 0.024489392592666, 0.086833274421172, 1.63896480048089E-05]
noise = [1] #all synapses are activated at the same time: 0, synapses are activated in a poisson fashion: 1

start_num_syn = int(10621 * 0.05) #we want to activate 5% of the total number of synapses - the total number is different for each cells, see paper.

end_num_syn = int(start_num_syn +10)
increase_syn_by = 10 #these two lines exist for technical reason, don't change it

stim_distances = [[200]] #here you can specify the location distances. If you want to stimulate the cell in the first 100 um AND the from 300 to 400 um you have to type [[50], [350]] with width = 50. Currently [[200]] with width = 200 setup the script will select dendritic location in the first 400 um of the dendritic arbor
width = 200

for z in noise:
    for j in stim_distances:
        for i in range(start_num_syn, end_num_syn, increase_syn_by):

            modellke(model_name, i, i, z, parameters, j, width)

            plt.savefig('./theta_things/' + model_name + '/'+ model_name + '_' + str(i) + '_' + str(j) + '_noise_' + str(z) +'_soma_rec_theta_excit_5_percent_long.svg')

