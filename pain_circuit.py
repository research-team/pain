from neuron import h, gui
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import math
#neuron.load_mechanisms("./mod")
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
    v_vec11 = h.Vector()
    v_vec13 = h.Vector()
    v_vec16 = h.Vector()
    v_vec17 = h.Vector()
    v_vec18 = h.Vector()
    v_vecka = h.Vector()
    v_veckd = h.Vector()
    v_veckca = h.Vector()
    t_vec = h.Vector()        # Time stamp vector
    #
    # v_vec11.record(compartment(0.5)._ref_v)
    # v_vec13.record(compartment(0.5)._ref_hinf_nav11_wt_R1648H)
    # v_vec16.record(compartment(0.5)._ref_sinf_nav11_wt_R1648H)
    # v_vec17.record(compartment(0.5)._ref_rinf_nav11_wt_R1648H)
    # v_vec18.record(compartment(0.5)._ref_ik_kv1)
    # v_vecka.record(compartment(0.5)._ref_ik_kv3)
    # v_veckd.record(compartment(0.5)._ref_ik_kv4)
    v_vec.record(compartment(0.5)._ref_v)
    # v_veckca.record(compartment(0.5)._ref_ik_iKCa)

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

def simulate(cells, tstop=500, vinit=-55):
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
    h.finitialize(vinit)
    for cell in cells:
        balance(cell)
    if h.cvode.active():
        h.cvode.active()
    else:
        h.fcurrent()
    h.frecord_init()
    h.tstop = tstop
    h.v_init = vinit
    h.run()

def spike_time_rec(section):
    vec = h.Vector()
    netcon = h.NetCon(section(0.5)._ref_v, None, sec=section)
    netcon.record(vec)
    return vec


def add_gaba_receptors(synapses, pre, g):

    for syn in synapses:
        nc = h.NetCon(pre(1)._ref_v, syn, sec = pre)
        nc.threshold = 0
        nc.delay = random.gauss(1, 0.2)
        nc.weight[0] = random.gauss(g, g / 10)
        gaba_ncs.append(nc)

def add_5HT_receptors(compartment, dist, g):
    for sec in compartment:
        diff = h.slow_5HT(sec(0.5))
        diff.h = dist#self.distances.get(compartment)
        diff.tx1 = 0 #+ (diff.h/50)*10#00
        diff.c0cleft = 2
        diff.a = 30
        rec = h.r5ht3a(sec(0.5))
        rec.gmax = random.gauss(g, g / 10)
        h.setpointer(diff._ref_serotonin, 'serotonin', rec)
        diffs.append(diff)
        recs.append(rec)

def add_glu_receptors(compartment, pre, g1, g2):
    '''
    Adds 5HT receptors
    Parameters
    ----------
    compartment: section of NEURON cell
        part of neuron
    x: int
        x - coordinate of serotonin application
    time: int (ms)
        time of serotonin application
    g: float
        receptor conductance
    '''
    diff = h.GrC_Gludif3(pre(0.5))
    h.setpointer(pre(0.5)._ref_v, 'affv', diff)
    diffs.append(diff)

    for sec in compartment:
        # print(compartment)
        rec = h.ampa_rec(sec(0.5))
        rec.gmax=random.gauss(g1, g1 / 10)
        h.setpointer(diff._ref_glu, 'glu_m', rec)        # vc = h.IClamp(0.5, sec=cell.axon1.node[5])
        recs.append(rec)
        rec2 = h.nmda_rec(sec(0.5))
        rec2.gmax=random.gauss(g2, g2 / 10)
        print(g2)
        h.setpointer(diff._ref_glu, 'glu_m', rec2)        # vc = h.IClamp(0.5, sec=cell.axon1.node[5])
        recs.append(rec2)

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
        # pyplot.plot(t_vec, i, label = label_v)
    outavg = np.mean(np.array(outavg), axis = 0, dtype=np.float32)

    pyplot.plot(t_vec, outavg, label = label_v)
    # pyplot.plot(v_vec, v_vec13, label = 'h')
    # pyplot.plot(v_vec, v_vec16, label = 's')
    # pyplot.plot(v_vec, v_vec17, label = 'r')
    # pyplot.plot(t_vec, v_vec18, label = 'Kv1')
    # pyplot.plot(t_vec, v_vecka, label = 'Kv3')
    # pyplot.plot(t_vec, v_veckd, label = 'Kv4')
    # pyplot.clf()
    # pyplot.plot(t_vec,  v_vec)
    # pyplot.plot(t_vec, v_veckca, label = 'K_Ca')

    pyplot.legend()
    pyplot.xlabel('time (ms)')
    pyplot.ylabel('mV')
    # pyplot.savefig(f"./results/ad_q_cur2_{dt}.pdf", format="pdf")


    f = open('./res.txt', 'w')
    for v in list(v_vec):
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
        IN2s.append(interneuron())
        IN1s.append(interneuron())
        LCNs.append(interneuron())
        WDRs.append(WDR_model())
        HTs.append(HT_model())
        a_deltas.append(adelta2(10, True))
        c_fibers.append(cfiber(250, 0.25, 0, 15000, True, 13))
    #
    #
    # # add_glu_receptors(IN2.dend, a_delta.axon1.node[0], 10, 1)
    for i in range(15):
        add_glu_receptors(HTs[random.randint(0, len(HTs)-1)].dend, a_deltas[random.randint(0, len(a_deltas)-1)].axon2.node[5], 150, 300)
        add_glu_receptors(IN2s[random.randint(0, len(IN2s)-1)].dend, a_deltas[random.randint(0, len(a_deltas)-1)].axon3.node[5], 15, 10)
        add_glu_receptors(IN1s[random.randint(0, len(IN2s)-1)].dend, c_fibers[random.randint(0, len(a_deltas)-1)].branch, 15, 10)
        add_glu_receptors(LCNs[random.randint(0, len(LCNs)-1)].dend, c_fibers[random.randint(0, len(c_fibers)-1)].branch, 150, 30)
        add_glu_receptors(WDRs[random.randint(0, len(WDRs)-1)].dend, LCNs[random.randint(0, len(LCNs)-1)].axon, 150, 100)
        add_gaba_receptors(HTs[random.randint(0, len(HTs)-1)].synlistinh, IN2s[random.randint(0, len(IN2s)-1)].axon, 0.005)
        add_gaba_receptors(WDRs[random.randint(0, len(WDRs)-1)].synlistinh, IN2s[random.randint(0, len(IN2s)-1)].axon, 0.005)
        # add_glu_receptors(IN1.dend, c_fiber.branch, 10, 1)
        add_gaba_receptors(LCNs[random.randint(0, len(LCNs)-1)].synlistinh, IN1s[random.randint(0, len(IN1s)-1)].axon, 0.0025)
        add_gaba_receptors(c_fibers[random.randint(0, len(c_fibers)-1)].synlistinh, IN1s[random.randint(0, len(IN1s)-1)].axon, 0.0035)


    # add_5HT_receptors(IN2s[random.randint(0, len(IN2s)-1)].dend, 100, 50)
    # add_5HT_receptors(IN1s[random.randint(0, len(IN1s)-1)].dend, 100, 50)
    #
    # add_5HT_receptors(a_deltas[random.randint(0, len(a_deltas)-1)].axon1.node, 100, 10)
    # add_5HT_receptors(c_fibers[random.randint(0, len(c_fibers)-1)].stimsec, 100, 20)

    # for sec in h.allsec():
    #     h.psection(sec=sec) #show parameters of each section
    v_vecs = []
    v_vecs_cf = []
    v_vecs_in2 = []
    v_vecs_in1 = []
    v_vecs_lcn = []
    v_vecs_wdr = []
    v_vecs_ht = []
    wdr_spikes = []

    for i in range(0, len(c_fibers)):
        v_vec, t_vec = set_recording_vectors(a_deltas[i].axon1.node[5])
        v_vec_cf, t_vec = set_recording_vectors(c_fibers[i].branch)
        v_vecs.append(v_vec)
        v_vecs_cf.append(v_vec_cf)
        v_vec_in2, t_vec = set_recording_vectors(IN2s[i].soma)
        v_vec_in1, t_vec = set_recording_vectors(IN1s[i].soma)

        v_vec_lcn, t_vec = set_recording_vectors(LCNs[i].soma)
        v_vec_wdr, t_vec = set_recording_vectors(WDRs[i].soma)
        wdr_spike_vec = spike_time_rec(WDRs[i].soma)
        wdr_spikes.append(wdr_spike_vec)

        v_vec_ht, t_vec = set_recording_vectors(HTs[i].soma)

        v_vecs_in2.append(v_vec_in2)
        v_vecs_in1.append(v_vec_in1)
        v_vecs_lcn.append(v_vec_lcn)
        v_vecs_wdr.append(v_vec_wdr)
        v_vecs_ht.append(v_vec_ht)


    simulate(c_fibers)
    show_output(v_vecs, t_vec, "a_delta")
    show_output(v_vecs_cf, t_vec, "c_fiber")
    show_output(v_vecs_in2, t_vec, "IN2")

    show_output(v_vecs_lcn, t_vec, "LCN")
    show_output(v_vecs_wdr, t_vec, "WDR")
    show_output(v_vecs_ht, t_vec, "HT")
    show_output(v_vecs_in1, t_vec, "IN1")

    for spike in wdr_spikes:
        print(list(spike))

    pyplot.show()

        # pyplot.show()
