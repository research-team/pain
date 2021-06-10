from neuron import h, gui
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import math
#neuron.load_mechanisms("./mod")
from adelta import adelta
from cfiber import cfiber

from adelta2 import adelta2
from interneuron import interneuron
from interneuron import interneuron
from WDR_neuron import WDR_model
from HT_neuron import HT_model


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
    for sec in cell.all:
        if ((-(sec.ina_nattxs + sec.ina_nav1p8 + sec.ina_Nav1_3  + sec.ina_nakpump + sec.ina_nav1p9) / (vinit - sec.ena)) < 0):
            sec.pumpina_extrapump = -(sec.ina_nattxs + sec.ina_nav1p8  + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_nav1p9)
        else:
            sec.gnaleak_leak = -(sec.ina_nattxs + sec.ina_nav1p8  + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_nav1p9) / (vinit - sec.ena)

        if ((-(sec.ik_kdr + sec.ik_nakpump + sec.ik_kv1 + sec.ik_kv4 + sec.ik_kv2 + sec.ik_iKCa) / (vinit - sec.ek)) < 0):
            sec.pumpik_extrapump = -(sec.ik_kdr + sec.ik_nakpump + sec.ik_kv1 + sec.ik_kv4 + sec.ik_kv2 + sec.ik_iKCa)
        else:
            sec.gkleak_leak = -(sec.ik_kdr + sec.ik_nakpump + sec.ik_kv1 + sec.ik_kv4 + sec.ik_kv2 + sec.ik_iKCa) / (vinit - sec.ek)

def simulate(cell, tstop=500, vinit=-60):
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
    balance(cell)
    if h.cvode.active():
        h.cvode.active()
    else:
        h.fcurrent()
    h.frecord_init()
    h.tstop = tstop
    h.v_init = vinit
    h.run()

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
    recs =[]
    diffs = []
    diff = h.GrC_Gludif3(pre(0.5))
    h.setpointer(pre(0.5)._ref_v, 'affv', diff)
    diffs.append(diff)

    for sec in compartment:
        rec = h.ampa_rec(sec(0.5))
        rec.gmax=g1
        h.setpointer(diff._ref_glu, 'glu_m', rec)        # vc = h.IClamp(0.5, sec=cell.axon1.node[5])
        recs.append(rec)
        rec2 = h.nmda_rec(sec(0.5))
        rec2.gmax=g2
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
    pyplot.plot(t_vec, v_vec, label = label_v)
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
    a_delta = adelta2(10, True)
    c_fiber = cfiber(250, 0.25, 0, 15020, True, 8)

    # diff = h.GrC_Gludif3(a_delta.axon1.node[0](0.5))
    # h.setpointer(a_delta.axon1.node[5](0.5)._ref_v, 'affv', diff)
    # diffs.append(diff)
    #
    # diff2 = h.GrC_Gludif3(c_fiber.branch(0.5))
    # h.setpointer(c_fiber.branch(0.5)._ref_v, 'affv', diff2)
    # diffs.append(diff2)

    IN1 = interneuron()
    IN2 = interneuron()
    LCN = interneuron()
    WDR =  WDR_model()
    HT = HT_model()

    add_glu_receptors(IN2.dend, a_delta.axon1.node[0], 150, 100)
    add_glu_receptors(HT.dend, a_delta.axon1.node[0], 150, 100)
    add_glu_receptors(LCN.dend, c_fiber.branch, 150, 100)
    add_glu_receptors(WDR.dend, LCN.axon, 150, 100)

    # recs =[]
    # for sec in cell2.dend:
    #     rec = h.ampa_rec(sec(0.5))
    #     rec.gmax=15000
    #     h.setpointer(diff._ref_glu, 'glu_m', rec)        # vc = h.IClamp(0.5, sec=cell.axon1.node[5])
    #     recs.append(rec)
    #     rec2 = h.nmda_rec(sec(0.5))
    #     rec2.gmax=100
    #     h.setpointer(diff._ref_glu, 'glu_m', rec2)        # vc = h.IClamp(0.5, sec=cell.axon1.node[5])
    #     recs.append(rec2)
        # vc.delay = 10.0
        # vc.dur = 80.0
        # vc.dur[2] = 0.0
        # vc.amp[0] = -70
        # vc.amp = 0.471 + i*0.005
        # vc.amp[2] = -45
        # for sec in cell.all:
        #
        #     vc = h.VClamp(sec(0.5))
        #     vc.dur[0] = 100.0
        #     vc.dur[1] = 1000.0
        #     vc.dur[2] = 1.0
        #     vc.amp[0] = -150
        #     vc.amp[1] = -10
        #     vc.amp[2] = 0

    for sec in h.allsec():
        h.psection(sec=sec) #show parameters of each section
    v_vec, t_vec = set_recording_vectors(a_delta.axon1.node[5])
    v_vec_cf, t_vec = set_recording_vectors(c_fiber.branch)
    v_vec_in2, t_vec = set_recording_vectors(IN2.soma)
    v_vec_lcn, t_vec = set_recording_vectors(LCN.soma)
    v_vec_wdr, t_vec = set_recording_vectors(WDR.soma)
    v_vec_ht, t_vec = set_recording_vectors(HT.soma)

    # print("Number of model - ",cell.numofmodel)
    # v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec2, v_veckca, t_vec = set_recording_vectors(cell2.soma)
    # v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec3, v_veckca, t_vec = set_recording_vectors(cell.axon2.node[0])

    simulate(c_fiber)
    show_output(v_vec, t_vec, "a_delta")
    show_output(v_vec_cf, t_vec, "c_fiber")
    show_output(v_vec_in2, t_vec, "IN2")
    show_output(v_vec_lcn, t_vec, "LCN")
    show_output(v_vec_wdr, t_vec, "WDR")
    show_output(v_vec_ht, t_vec, "HT")

    # show_output(v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec2, v_veckca, t_vec)
    # show_output(v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec3, v_veckca, t_vec)
    pyplot.show()

        # pyplot.show()
