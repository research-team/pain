from neuron import h, gui
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import math
#neuron.load_mechanisms("./mod")
from cfiber2 import cfiber2
from adelta import adelta

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
    v_vec.record(compartment(0.5)._ref_vext[0])
    t_vec.record(h._ref_t)
    return v_vec, t_vec

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
        if ((-(sec.ina_nattxs + sec.ina_navv1p8 + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_na11a) / (vinit - sec.ena)) < 0):
            sec.pumpina_extrapump = -(sec.ina_nattxs + sec.ina_navv1p8 + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_na11a)
        else:
            sec.gnaleak_leak = -(sec.ina_nattxs + sec.ina_navv1p8 + sec.ina_Nav1_3 + sec.ina_nakpump + sec.ina_na11a) / (vinit - sec.ena)

        if ((-(sec.ik_kdr  + sec.ik_nakpump + sec.ik_kap + sec.ik_kad) / (vinit - sec.ek)) < 0):
            sec.pumpik_extrapump = -(sec.ik_kdr + sec.ik_nakpump + sec.ik_kap + sec.ik_kad)
        else:
            sec.gkleak_leak = -(sec.ik_kdr + sec.ik_nakpump  + sec.ik_kap + sec.ik_kad) / (vinit - sec.ek)

def simulate(cell, tstop=500, vinit=-70):
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

def show_output(v_vec, t_vec):
    ''' show graphs
    Parameters
    ----------
    v_vec: h.Vector()
        recorded voltage
    t_vec: h.Vector()
        recorded time
    '''
    dend_plot = pyplot.plot(t_vec, v_vec)
    f = open('./res.txt', 'w')
    for v in list(v_vec):
        f.write(str(v)+"\n")
    pyplot.xlabel('time (ms)')
    pyplot.ylabel('mV')

if __name__ == '__main__':
    # numofmodel = int(sys.argv[3])
    # if numofmodel < 1 or numofmodel > 14:
    #     print("ERROR! Please input model number in range 1...14")
    # else:
    cell = adelta(15)
    # stim = h.IClamp(cell.branch(1))
    # stim.delay = 5
    # stim.dur = 1
    # stim.amp = 0.1
    for sec in h.allsec():
        h.psection(sec=sec) #show parameters of each section
    branch_vec, t_vec = set_recording_vectors(cell.node[4])
    # print("Number of model - ",cell.numofmodel)
    simulate(cell)
    show_output(branch_vec, t_vec)
    pyplot.show()
