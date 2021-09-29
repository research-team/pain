from neuron import h
from neuron.units import ms, mV
h.load_file("stdgui.hoc")
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import math
#neuron.load_mechanisms("./mod")
import logging
logging.basicConfig(filename='logs.log',
                        filemode='a',
                        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.DEBUG)
logging.info("let's get it started")
from adelta import adelta
from cfiber import cfiber
import random

from adelta2 import adelta2
from interneuron import interneuron
from interneuron import interneuron
from WDR_neuron import WDR_model
from HT_neuron import HT_model
import numpy as np

recs = []
diffs = []
gaba_ncs = []
stims = []

def set_recording_vectors(compartment):
    ''' recording voltage
    Parameters
    ----------
    compartment: NEURON section
        compartment for recording
    Returns
    -------
    v_vec: h.Vector()
        recorded voltage
    t_vec: h.Vector()
        recorded time
    '''
    v_vec = h.Vector()   # Membrane potential vector at compartment
    t_vec = h.Vector()        # Time stamp vector

    v_vec.record(compartment(0.5)._ref_v)
    t_vec.record(h._ref_t)
    return v_vec, t_vec

def balance(cell, vinit=-55):
    ''' voltage balance
    Parameters
    ----------
    cell: NEURON cell
        cell for balance
    vinit: int (mV)
        initialized voltage
    '''
    for sec in cell.all_secs:
        if ((-(sec.ina_nattxs + sec.ina_nav1p8 + sec.ina_Nav1_3  + sec.ina_nakpump + sec.ina_nav1p9) / (vinit - sec.ena)) < 0):
            sec.pumpina_extrapump = -(sec.ina_nattxs + sec.ina_nav1p8  + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_nav1p9)
        else:
            sec.gnaleak_leak = -(sec.ina_nattxs + sec.ina_nav1p8  + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_nav1p9) / (vinit - sec.ena)

        if ((-(sec.ik_kdr + sec.ik_nakpump + sec.ik_kap + sec.ik_kad + sec.ik_iKCa) / (vinit - sec.ek)) < 0):
            sec.pumpik_extrapump = -(sec.ik_kdr + sec.ik_nakpump + sec.ik_kad + sec.ik_kap+ sec.ik_iKCa)
        else:
            sec.gkleak_leak = -(sec.ik_kdr + sec.ik_nakpump + sec.ik_kap + sec.ik_kad  + sec.ik_iKCa) / (vinit - sec.ek)

def simulate(cells, tstop=1000, vinit=-55):
    ''' simulation control
    Parameters
    ----------
    cell: NEURON cell
        cell for simulation
    tstop: int (ms)
        simulation time
    vinit: int (mV)
        initialized voltage
    '''
    h.finitialize(vinit * mV)
    for cell in cells:
        balance(cell)
    if h.cvode.active():
        h.cvode.active()
    else:
        h.fcurrent()
    h.frecord_init()
    h.tstop = tstop
    h.v_init = vinit
    h.continuerun(tstop * ms)

def spike_time_rec(section, th):
    vec = h.Vector()
    netcon = h.NetCon(section(0.5)._ref_v, None, sec=section)
    netcon.threshold = th
    netcon.record(vec)
    return vec

def add_gaba_receptors(synapses, pre, g):
    for syn in synapses:
        nc = h.NetCon(pre(1)._ref_v, syn, sec = pre)
        nc.threshold = 0
        nc.delay = random.gauss(1, 0.2)
        nc.weight[0] = random.gauss(g, g / 10)
        gaba_ncs.append(nc)

def gen_connect(synapses, freq):
    ees = h.NetStim()
    ees.number = 1000000
    ees.start = 10
    ees.interval = int(1000 / freq)
    stims.append(ees)

    for syn in synapses:
        nc = h.NetCon(ees, syn)
        nc.delay = random.gauss(1, 0.01)
        nc.weight[0] = random.gauss(10.5, 1.2)
        recs.append(nc)

def add_5HT_receptors(compartment, dist, g):
    for sec in compartment:
        diff = h.slow_5HT(sec(0.5))
        diff.h = dist#self.distances.get(compartment)
        diff.tx1 = 0 + (diff.h/50)*10#00
        diff.c0cleft = 2
        diff.a = 10
        rec = h.r5ht3a(sec(0.5))
        rec.gmax = random.gauss(g, g / 10)
        h.setpointer(diff._ref_serotonin, 'serotonin', rec)
        diffs.append(diff)
        recs.append(rec)

def pool_recording(pool, th):
    v_vecs = []
    spikes_vec = []
    for cell in pool:
        v_vec, t_vec = set_recording_vectors(list(cell.all_secs)[0])
        # print(list(cell.all_secs)[0])
        spike_vec = spike_time_rec(list(cell.all_secs)[0], th)
        spikes_vec.append(spike_vec)
        v_vecs.append(v_vec)
    return v_vecs, spikes_vec, t_vec

def set_syn_vectors(syn):
    v_vec = h.Vector()
    v_vec.record(syn._ref_dw)
    return v_vec

def syn_recording(pool):
    syn_vecs = []
    for cell in pool:
        # for i in range(len(cell.dend)-1):
        v_vec = set_syn_vectors(cell.synliststpd[0])
        syn_vecs.append(v_vec)
    return syn_vecs

def spike_time_plot(pool_spikes, ymax, color, file_name):
    with open(file_name,'w') as spk_file: # Append
        for spike in pool_spikes:
        # print(list(spike))
            for x in list(spike):
                pyplot.vlines(x, ymax-1, ymax, colors=color)
                spk_file.write(str(x)+"\n")
            pyplot.xlim(0, h.tstop)


def add_glu_receptors(cell, pre, g1, g2, withSTDP=False):
    diff = h.GrC_Gludif3(pre(0.5))
    h.setpointer(pre(0.5)._ref_v, 'affv', diff)
    diffs.append(diff)
    for i in range(len(cell.dend)-1):
        # print(section)
        rec = h.ampa_rec(cell.dend[i](0.5))
        rec.gmax=random.gauss(g1, g1 / 10)
        h.setpointer(diff._ref_glu, 'glu_m', rec)        # vc = h.IClamp(0.5, sec=cell.axon1.node[5])
        recs.append(rec)
        rec2 = h.nmda_rec(cell.dend[i](0.5))
        rec2.gmax=random.gauss(g2, g2 / 10)
        # print(g2)
        h.setpointer(diff._ref_glu, 'glu_m', rec2)        # vc = h.IClamp(0.5, sec=cell.axon1.node[5])
        recs.append(rec2)

        if withSTDP:
            syn = cell.synliststpd[i]
            nc = h.NetCon(pre(0.5)._ref_v, syn, sec=pre)
            presyn = h.NetCon(pre(0.5)._ref_v, syn, sec=pre)
            presyn.weight[0] = 1
            presyn.threshold = -20

            postsyn = h.NetCon(cell.dend[i](0.5)._ref_v, syn, sec=cell.dend[i])
            postsyn.weight[0] = -1
            postsyn.threshold = -10

            h.setpointer(rec2._ref_gmax, 'wsyn', syn)

            nc.weight[0] = syn.wsyn
            nc.delay = 1
            recs.append(nc)
            recs.append(presyn)
            recs.append(postsyn)

        # return stdp_recs

def show_output(v_vec, t_vec, label_v):
    ''' show graphs
    Parameters
    ----------
    v_vec: h.Vector()
        recorded voltage
    t_vec: h.Vector()
        recorded time
    '''
    outavg = []
    for i in v_vec:
        outavg.append(list(i))
        # pyplot.plot(t_vec, i)
    outavg = np.mean(np.array(outavg), axis = 0, dtype=np.float32)

    pyplot.plot(t_vec, outavg, label = label_v)

    pyplot.legend()
    pyplot.xlabel('time (ms)')
    pyplot.ylabel('mV')
    # pyplot.savefig(f"./results/ad_q_cur2_{dt}.pdf", format="pdf")

    with open(f'./results/{label_v}_res.txt', 'w') as f:
        for v in list(outavg):
            f.write(str(v)+"\n")

if __name__ == '__main__':
    # numofmodel = int(sys.argv[3])
    # if numofmodel < 1 or numofmodel > 14:
    #     print("ERROR! Please input model number in range 1...14")
    # else:
    # for i in range(5):
    # diffs = []
    IN2s = []
    IN1s = []
    LCNs = []
    WDRs = []
    HTs = []
    a_deltas = []
    c_fibers = []


    for i in range(10):
        WDRs.append(WDR_model())
        c_fibers.append(cfiber(250, 0.25, 0, 15000, True, 14))
        a_deltas.append(adelta2(10, True))

        HTs.append(HT_model())
        IN2s.append(interneuron())
        IN1s.append(interneuron())
        LCNs.append(interneuron())

    #
    logging.info("create circuit")

    # # add_glu_receptors(IN2.dend, a_delta.axon1.node[0], 10, 1)

    for i in range(20):
        add_glu_receptors(HTs[random.randint(0, len(HTs)-1)], a_deltas[random.randint(0, len(a_deltas)-1)].axon2.node[5], 150, 250, False)
        add_glu_receptors(IN2s[random.randint(0, len(IN2s)-1)], a_deltas[random.randint(0, len(a_deltas)-1)].axon3.node[5], 20, 10)
        add_glu_receptors(IN1s[random.randint(0, len(IN1s)-1)], c_fibers[random.randint(0, len(c_fibers)-1)].branch, 15, 10)
        add_glu_receptors(LCNs[random.randint(0, len(WDRs)-1)], c_fibers[random.randint(0, len(c_fibers)-1)].branch, 150, 300)
        add_glu_receptors(WDRs[random.randint(0, len(WDRs)-1)], LCNs[random.randint(0, len(LCNs)-1)].axon, 150, 400, False)
        add_gaba_receptors(HTs[random.randint(0, len(HTs)-1)].synlistinh, IN2s[random.randint(0, len(IN2s)-1)].axon, 1.5)
        add_gaba_receptors(WDRs[random.randint(0, len(WDRs)-1)].synlistinh, IN2s[random.randint(0, len(IN2s)-1)].axon, 1.5)
        # add_glu_receptors(IN1.dend, c_fiber.branch, 10, 1)
        add_gaba_receptors(LCNs[random.randint(0, len(LCNs)-1)].synlistinh, IN1s[random.randint(0, len(IN1s)-1)].axon, 1.5)
        add_gaba_receptors(c_fibers[random.randint(0, len(c_fibers)-1)].synlistinh, IN1s[random.randint(0, len(IN1s)-1)].axon, 0.5)

        # gen_connect(a_deltas[random.randint(0, len(a_deltas)-1)].synlistex, 50)
        # gen_connect(c_fibers[random.randint(0, len(c_fibers)-1)].synlistex, 50)

    logging.info("add connections")
    #
    # add_5HT_receptors(IN2s[random.randint(0, len(IN2s)-1)].dend, 100, 30)
    # add_5HT_receptors(IN1s[random.randint(0, len(IN1s)-1)].dend, 100, 30)
    #
    # add_5HT_receptors(a_deltas[random.randint(0, len(a_deltas)-1)].axon1.node, 100, 50)
    # add_5HT_receptors(c_fibers[random.randint(0, len(c_fibers)-1)].stimsec, 100, 25)

    # for sec in h.allsec():
    #     h.psection(sec=sec) #show parameters of each section
    wdr_syn_vecs = syn_recording(WDRs)
    # ht_syn_vecs = syn_recording(HTs)

    v_vecs, spikes_vec, t_vec = pool_recording(a_deltas, -20)
    v_vecs_cf, spikes_vecs_cf, t_vec = pool_recording(c_fibers, 0)
    v_vecs_in2, spikes_vecs_in2, t_vec = pool_recording(IN2s, 0)
    v_vecs_in1, spikes_vecs_in1, t_vec = pool_recording(IN1s, 0)
    v_vecs_lcn, spikes_vecs_lcn, t_vec = pool_recording(LCNs, 0)
    v_vecs_wdr, spikes_vecs_wdr, t_vec = pool_recording(WDRs, -10)
    v_vecs_ht, spikes_vecs_ht, t_vec = pool_recording(HTs, -5)

    logging.info("add recorders")

    simulate(c_fibers)

    logging.info("simulated")

    # show_output(v_vecs, t_vec, "a_delta")
    # show_output(v_vecs_cf, t_vec, "c_fiber")
    # show_output(v_vecs_in2, t_vec, "IN2")
    #
    # show_output(v_vecs_lcn, t_vec, "LCN")
    #
    # show_output(ht_syn_vecs, t_vec, "HT")
    # show_output(v_vecs_in1, t_vec, "IN1")
    # pyplot.show()
    #
    # show_output(v_vecs_wdr, t_vec, "WDR")
    # pyplot.show()
    spike_time_plot(spikes_vec, 1, 'k', "./results/a_delta_spk.txt")
    spike_time_plot(spikes_vecs_cf, 2, 'c', "./results/c_fiber_spk.txt")
    spike_time_plot(spikes_vecs_wdr, 3, 'm', "./results/WDR_spk.txt")
    spike_time_plot(spikes_vecs_lcn, 4, 'c', "./results/lcn_spk.txt")
    spike_time_plot(spikes_vecs_ht, 5, 'r', "./results/HT_spk.txt")

    logging.info("recorded")

    # spike_time_plot(spikes_vecs_in1, 5, 'k')
    # spike_time_plot(spikes_vecs_in2, 6, 'm')

    pyplot.show()
