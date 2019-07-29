from neuron import h, gui
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import math
#neuron.load_mechanisms("./mod")
from cfiber import cfiber

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
    for sec in cell.all:
        if ((-(sec.ina_nattxs + sec.ina_navv1p8 + sec.ina_Nav1_3 + sec.ina_nakpump) / (vinit - sec.ena)) < 0):
            sec.pumpina_extrapump = -(sec.ina_nattxs + sec.ina_navv1p8 + sec.ina_Nav1_3 + sec.ina_nakpump)
        else:
            sec.gnaleak_leak = -(sec.ina_nattxs + sec.ina_navv1p8 + sec.ina_Nav1_3 + sec.ina_nakpump) / (vinit - sec.ena)

        if ((-(sec.ik_kdr + sec.ik_nakpump + sec.ik_kap + sec.ik_kad) / (vinit - sec.ek)) < 0):
            sec.pumpik_extrapump = -(sec.ik_kdr + sec.ik_nakpump + sec.ik_kap + sec.ik_kad)
        else:
            sec.gkleak_leak = -(sec.ik_kdr + sec.ik_nakpump + sec.ik_kap + sec.ik_kad) / (vinit - sec.ek)

def simulate(cell, tstop=300, vinit=-55):
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
    if cell.numofmodel == 9 or cell.numofmodel == 10 or cell.numofmodel == 12:
        running_ = 1
        if cell.numofmodel == 9:
            dl = 0
            dt = 120000
        elif cell.numofmodel == 10:
            dl = 1000
            dt = 40
        else:
            dl = 800
            dt = 50
        h.stdinit()
        for n in range(2):
            cell.x_application = cell.x_application + dl
            cell.distance()
            for item in cell.diffs:
                item.tx1 = h.t + 1 
                item.initial = item.atp
                item.c0cleft = item.c0cleft
                item.h = cell.distances.get(cell.diffusions.get(item))
            h.continuerun(h.t+dt)
        h.continuerun(h.t+500)

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
    pyplot.xlabel('time (ms)')
    pyplot.ylabel('mV')

if __name__ == '__main__':
    numofmodel = int(sys.argv[3])
    if numofmodel < 1 or numofmodel > 14:
        print("ERROR! Please input model number in range 1...14")
    else:
        cell = cfiber(250, 0.25, 0, 15000, True, numofmodel)
        for sec in h.allsec():
            h.psection(sec=sec) #show parameters of each section
        branch_vec, t_vec = set_recording_vectors(cell.branch)
        print("Number of model - ",cell.numofmodel)
        simulate(cell)
        show_output(branch_vec, t_vec)
        pyplot.show()