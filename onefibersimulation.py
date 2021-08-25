from neuron import h, gui
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import math
import sys
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
    v_vec11 = h.Vector()
    v_vec13 = h.Vector()
    v_vec16 = h.Vector()
    v_vec17 = h.Vector()
    v_vec18 = h.Vector()
    v_vecka = h.Vector()
    v_veckd = h.Vector()
    v_veckca = h.Vector()
    t_vec = h.Vector()        # Time stamp vector

    v_vec11.record(compartment(0.5)._ref_ina_nav1p9)
    v_vec13.record(compartment(0.5)._ref_ica_iCaAN)
    v_vec16.record(compartment(0.5)._ref_ik_kdr)
    v_vec17.record(compartment(0.5)._ref_ina_nattxs)
    v_vec18.record(compartment(0.5)._ref_ina_nav1p8)
    v_vecka.record(compartment(0.5)._ref_ik_kap)
    v_veckd.record(compartment(0.5)._ref_ik_kad)
    v_vec.record(compartment(0.5)._ref_vext[0])
    v_veckca.record(compartment(0.5)._ref_ik_iKCa)

    t_vec.record(h._ref_t)
    return v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec, v_veckca, t_vec

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

def simulate(cell, tstop=350, vinit=-55):
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
            d_t = 60000
        elif cell.numofmodel == 10:
            dl = 1000
            d_t = 40
        else:
            dl = 800
            d_t = 10
        h.stdinit()
        for n in range(3):
            cell.x_application = cell.x_application + dl
            if n == 0:
                cell.dl = 5000
                cell.x_application = 0
                print(cell.x_application)
            else:
                cell.dl = 30050
                cell.x_application = -400
            cell.distance()
            for item in cell.diffs:
                item.tx1 = h.t + 5
                # item.initial = 0#item.atp
                # if n == 0:
                #     item.c0cleft = 0.01#item.c0cleft
                # else:
                item.c0cleft = item.c0cleft
                item.h = cell.distances.get(cell.diffusions.get(item))
            h.continuerun(h.t+d_t)
        h.continuerun(h.t+500)

def show_output(v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec, v_veckca, t_vec, dt):
    ''' show graphs
    Parameters
    ----------
    v_vec: h.Vector()
        recorded voltage
    t_vec: h.Vector()
        recorded time
    '''
    pyplot.plot(t_vec, v_vec11, label = 'Nav1.9')
    pyplot.plot(t_vec, v_vec13, label = 'iCan')
    pyplot.plot(t_vec, v_vec16, label = 'KDr')
    pyplot.plot(t_vec, v_vec17, label = 'Nav1.7')
    pyplot.plot(t_vec, v_vec18, label = 'Nav1.8')
    pyplot.plot(t_vec, v_vecka, label = 'Kv2')
    pyplot.plot(t_vec, v_veckd, label = 'Kv4')
    pyplot.clf()
    pyplot.plot(t_vec, v_vec, label = 'V')
    print(f'max - {max(v_vec)}')
    print(f'min - {min(v_vec)}')

    # pyplot.plot(t_vec, v_veckca, label = 'K_Ca')



    f = open('./res.txt', 'w')
    for v in list(v_vec):
        f.write(str(v)+"\n")
    pyplot.legend()
    pyplot.xlabel('time (ms)')
    pyplot.ylabel('mV')

if __name__ == '__main__':
    numofmodel = int(sys.argv[3])
    if numofmodel < 1 or numofmodel > 14:
        print("ERROR! Please input model number in range 1...14")
    else:
        cell = cfiber(250, 0.25, 0, 15000, True, numofmodel)
        # for sec in h.allsec():
        #     h.psection(sec=sec) #show parameters of each section
        # branch_vec, t_vec = set_recording_vectors(cell.stimsec[9])
        v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec, v_veckca, t_vec = set_recording_vectors(cell.branch)

        # vc = h.VClamp(0.5, sec=cell.stimsec[9])
        # vc.dur[0] = 0.0
        # vc.dur[1] = 90.0
        # vc.dur[2] = 5.0
        # vc.amp[0] = -55
        # vc.amp[1] = -50
        # vc.amp[2] = -45
        # branch_vec1, t_vec1 = set_recording_vectors(cell.stimsec[1])
        # branch_vec2, t_vec2 = set_recording_vectors(cell.stimsec[4])
        print("Number of model - ",cell.numofmodel)
        simulate(cell)
        show_output(v_vec11, v_vec13, v_vec16, v_vec17, v_vec18, v_vecka, v_veckd, v_vec, v_veckca, t_vec, 5)
        # show_output(branch_vec1, t_vec1)
        # show_output(branch_vec2, t_vec2)

        pyplot.show()
