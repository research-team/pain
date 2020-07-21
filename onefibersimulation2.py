from neuron import h, gui
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import math
#neuron.load_mechanisms("./mod")
from adelta import adelta
from adelta2 import adelta2

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

    # v_vec11.record(compartment(0.5)._ref_ina_nav1p1)
    # v_vec13.record(compartment(0.5)._ref_ina_Nav1_3)
    # v_vec16.record(compartment(0.5)._ref_ina_nav1p6)
    # v_vec17.record(compartment(0.5)._ref_ina_nattxs)
    # v_vec18.record(compartment(0.5)._ref_ina_navv1p8)
    # v_vecka.record(compartment(0.5)._ref_ik_kv3)
    # v_veckd.record(compartment(0.5)._ref_ik_kv4)
    v_vec.record(compartment(0.5)._ref_v)
    # v_veckca.record(compartment(0.5)._ref_ik_iKCa)

    t_vec.record(h._ref_t)
    return v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec, v_veckca, t_vec

def balance(cell, vinit=-70):
    ''' voltage balance
    Parameters
    ----------
    cell: NEURON cell
        cell for balance
    vinit: int (mV)
        initialized voltage
    '''
    for sec in cell.node:
        if ((-(sec.ina_nattxs + sec.ina_navv1p8 + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_na11a + sec.ina_na16a) / (vinit - sec.ena)) < 0):
            sec.pumpina_extrapump = -(sec.ina_nattxs + sec.ina_navv1p8 + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_na11a + sec.ina_na16a)
        else:
            sec.gnaleak_leak = -(sec.ina_nattxs + sec.ina_navv1p8 + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_na11a + sec.ina_na16a) / (vinit - sec.ena)

        if ((-(sec.ik_kdr  + sec.ik_nakpump + sec.ik_kap + sec.ik_kad) / (vinit - sec.ek)) < 0):
            sec.pumpik_extrapump = -(sec.ik_kdr + sec.ik_nakpump + sec.ik_kap + sec.ik_kad)
        else:
            sec.gkleak_leak = -(sec.ik_kdr + sec.ik_nakpump  + sec.ik_kap + sec.ik_kad) / (vinit - sec.ek)

def simulate(cell, tstop=120, vinit=-70):
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
    # h.finitialize(vinit)
    # balance(cell)
    # if h.cvode.active():
    #     h.cvode.active()
    # else:
    #     h.fcurrent()
    # h.frecord_init()
    h.tstop = tstop
    h.v_init = vinit
    h.run()

def show_output(v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec, v_veckca, t_vec, dt):
    ''' show graphs
    Parameters
    ----------
    v_vec: h.Vector()
        recorded voltage
    t_vec: h.Vector()
        recorded time
    '''
    # pyplot.plot(t_vec, v_vec11, label = 'Nav1.1')
    # pyplot.plot(t_vec, v_vec13, label = 'Nav1.3')
    # pyplot.plot(t_vec, v_vec16, label = 'Nav1.6')
    # pyplot.plot(t_vec, v_vec17, label = 'Nav1.7')
    # pyplot.plot(t_vec, v_vec18, label = 'Nav1.8')
    # pyplot.plot(t_vec, v_vecka, label = 'Kv3')
    # pyplot.plot(t_vec, v_veckd, label = 'Kv4')
    pyplot.clf()
    pyplot.plot(t_vec, v_vec)
    # pyplot.plot(t_vec, v_veckca, label = 'K_Ca')



    # f = open('./res.txt', 'w')
    # for v in list(v_vec):
    #     f.write(str(v)+"\n")
    pyplot.legend()
    pyplot.xlabel('time (ms)')
    pyplot.ylabel('mV')
    # pyplot.savefig(f"./results/ad_dt_{dt}.pdf", format="pdf")

if __name__ == '__main__':
    # numofmodel = int(sys.argv[3])
    # if numofmodel < 1 or numofmodel > 14:
    #     print("ERROR! Please input model number in range 1...14")
    # else:
    for i in range(1):
        cell = adelta2(i*15)
        # stim = h.IClamp(cell.branch(1))
        # stim.delay = 5
        # stim.dur = 1
        # stim.amp = 0.1
        # for sec in h.allsec():
        #     h.psection(sec=sec) #show parameters of each section
        v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec, v_veckca, t_vec = set_recording_vectors(cell.axon1.node[2])
        # print("Number of model - ",cell.numofmodel)
        simulate(cell)
        show_output(v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec, v_veckca, t_vec, i*5)
        pyplot.show()
