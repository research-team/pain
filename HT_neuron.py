from neuron import h, gui
import random

import random

class HT_model(object):
  '''
  Interneuron class with parameters:
    delay: bool
      Does it have 5ht receptors?
      -Yes: True
      -No: False
    soma: NEURON Section (creates by topol())
    dend: NEURON Section (creates by topol())
    axon: NEURON Section (creates by topol())
    synlistinh: list (creates by synapses())
      list of inhibitory synapses
    synlistex: list (creates by synapses())
      list of excitatory synapses
    synlistees: list (creates by synapses())
      list of excitatory synapses for connection with generators
    x, y, z: int
      3D coordinates (isn't used)
    diffs: list
      list of diffusion mechanisms (NEURON staff)
    recs: list
      list of receptors mechanisms (NEURON staff)
  '''
  def __init__(self):
    self.diffs = []
    self.recs = []
    self.topol()
    # self.subsets()
    self.geom()
    # self.geom_nseg()
    self.biophys()
    self.synlistinh = []
    self.synliststpd = []
    self.synapses()
    self.x = self.y = self.z = 0.

    def __del__(self):
    #print 'delete ', self
      pass

  def topol(self):
    '''
    Creates sections soma, dend, axon and connects them
    if it's delay creates section dend[]: array
    '''
    self.soma = h.Section(name='soma', cell=self)
    self.axon = h.Section(name='axon', cell=self)
    self.hillock = h.Section(name='hillock', cell= self)
    self.dend = [h.Section(name='dend[%d]' % i) for i in range(random.randint(5,10))]
    for sec in self.dend:
      sec.connect(self.soma(0))
    self.hillock.connect(self.soma(1))
    self.axon.connect(self.hillock(1))
    self.all_secs = h.SectionList()
    # for sec in self.branch:
    self.all_secs.append(sec=self.soma)
    self.all_secs.append(sec=self.axon)
    self.all_secs.append(sec=self.hillock)
    for sec in self.dend:
        self.all_secs.append(sec=sec)

  def subsets(self):
    '''
    NEURON staff
    adds sections in NEURON SectionList
    '''
    self.all = h.SectionList()
    for sec in h.allsec():
      self.all.append(sec=sec)

  def geom(self):
    '''
    Adds length and diameter to sections
    '''
    self.soma.L = self.soma.diam = 15#random.randint(5, 15) # microns
    self.soma.nseg = 3

    self.hillock.L = 9 # microns
    self.hillock.diam = 2 # microns
    self.hillock.nseg = 3

    self.axon.L = 1000 # microns
    self.axon.diam = 1 # microns
    self.axon.nseg = 5

    for sec in self.dend:
      sec.L = 250 # microns
      sec.diam = 2.5
      sec.nseg = 5#random.gauss(1, 0.1) # microns

  # def geom_nseg(self):
  #   '''
  #   Calculates numder of segments in section
  #   '''
  #   for sec in self.all:
  #     sec.nseg = int((sec.L/(0.1*h.lambda_f(100)) + .9)/2.)*2 + 1

  def biophys(self):
    '''
    Adds channels and their parameters
    if delay is true, adds 5ht receptors
    '''
    # for sec in self.all:
    #   sec.cm = 0.8#random.gauss(1, 0.01) # cm uf/cm2 - membrane capacitance
    #   # sec.ena = 55
    #   # sec.ek = -70
    #   sec.Ra = 220

    self.soma.insert('fastchannels')
    self.soma.gnabar_fastchannels = 0.08
    self.soma.gkbar_fastchannels = 0.02
    self.soma.gl_fastchannels = 0.0004
    self.soma.el_fastchannels = -60
    self.soma.insert('iKCa')
    self.soma.insert('iCaL')
    self.soma.insert('Kv7M')
    self.soma.gbar_Kv7M = 50
    self.soma.cm = 0.8
    self.soma.Ra = 220

    # self.soma.insert('iNaP')
    # self.soma.insert('iCaAN')
    self.soma.insert('CaIntraCellDyn')
    self.soma.gbar_iKCa = 0.0001
    self.soma.depth_CaIntraCellDyn = 0.1
    self.soma.cai_tau_CaIntraCellDyn = 1.0
    self.soma.cai_inf_CaIntraCellDyn = 50.0e-6
    self.soma.pcabar_iCaL = 0.00001
    # self.soma.gnabar_iNaP = 0.0001

    self.soma.insert('extracellular') #adds extracellular mechanism for recording extracellular potential

    for sec in self.dend:
      sec.insert('pas')
      sec.g_pas = 0.0005
      sec.e_pas = -65
      sec.insert('iKCa')
      sec.insert('iCaL')
      sec.insert('CaIntraCellDyn')
      # sec.insert('iCaAN')
      sec.gbar_iKCa = 0.002
      sec.depth_CaIntraCellDyn = 0.1
      sec.cai_tau_CaIntraCellDyn = 2.0
      sec.cai_inf_CaIntraCellDyn = 50.0e-6
      sec.pcabar_iCaL = 3e-5
      sec.cm = 0.8#random.gauss(1, 0.01) # cm uf/cm2 - membrane capacitance
      sec.Ra = 220
      # sec.gbar_iCaAN =  0.00007

    self.hillock.insert('fastchannels')
    self.hillock.gnabar_fastchannels = 2.45
    self.hillock.gkbar_fastchannels = 0.076
    self.hillock.gl_fastchannels = 0.000042
    self.hillock.el_fastchannels = -65

    self.axon.insert('fastchannels')
    self.axon.gnabar_fastchannels = 0.01
    self.axon.gkbar_fastchannels = 0.04
    self.axon.gl_fastchannels = 0.00001
    self.axon.el_fastchannels = -65

    # self.axon.Ra = 50
    # self.axon.insert('hh')

  def add_5HTreceptors(self, compartment, time, g):
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
    diff = h.slow_5HT(compartment(0.5))
    diff.h = random.uniform(10, 2500)
    diff.tx1 = time + 0 + (diff.h/50)*10#00
    diff.c0cleft = 3
    diff.a = 0.1
    rec = h.r5ht3a(compartment(0.5))
    rec.gmax = g
    h.setpointer(diff._ref_serotonin, 'serotonin', rec)
    self.diffs.append(diff)
    self.recs.append(rec)

  def position(self, x, y, z):
    '''
    NEURON staff
    Adds 3D position
    '''
    soma.push()
    for i in range(h.n3d()):
      h.pt3dchange(i, x-self.x+h.x3d(i), y-self.y+h.y3d(i), z-self.z+h.z3d(i), h.diam3d(i))
    self.x = x; self.y = y; self.z = z
    h.pop_section()

  def connect2target(self, target):
    '''
    NEURON staff
    Adds presynapses
    Parameters
    ----------
    target: NEURON cell
        target neuron
    Returns
    -------
    nc: NEURON NetCon
        connection between neurons
    '''
    nc = h.NetCon(self.soma(1)._ref_v, target, sec = self.soma)
    nc.threshold = 0
    return nc

  def synapses(self):
    '''
    Adds synapses
    '''
    for i in range(20):
        s = h.GABAa_DynSyn(self.soma(0.5)) # Inhibitory
        # s.tau1 = 0.5
        # s.tau2 = 3.5
        # s.e = -80
        self.synlistinh.append(s)

    for sec in self.dend:
        for i in range(2):
            syn = h.StdwaSA(sec(0.5))
            self.synliststpd.append(syn)


  def is_art(self):
    return 0
