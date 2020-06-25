from neuron import h
h.load_file('stdlib.hoc') #for h.lambda_f

import random

class axon(object):
  '''
  bio-axon class with parameters:
    axon parameters from: https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=3810&file=/MRGaxon/MRGaxon.hoc#tabs-2
  number:
    number of nodes of Ranvier
  '''
  def __init__(self, number):
    PI = 3.14
    #topological parameters
    self.number = number
    self.axonnodes = number + 1
    self.paranodes1 = 2*number
    self.paranodes2 = 2*number
    self.axoninter = 6*number
    #morphological parameters
    self.fiberD = 5.0
    self.paralength1 = 5
    self.nodelength = 3
    space_p1 = 0.002
    space_p2 = 0.004
    space_i = 0.004
    self.g = 0.690
    self.axonD = 3
    self.nodeD = 2
    self.paraD1 = 2
    self.paraD2 = 3
    deltax = 1150
    self.paralength2 = 50
    self.nl = 120
    self.interlength=(deltax-self.nodelength-(2*self.paralength1)-(2*self.paralength2))/2
    #electrical parameters
    self.rhoa=6e5 #Ohm-um
    self.mycm=0.1 #uF/cm2/lamella membrane
    self.mygm=0.001 #S/cm2/lamella membrane
    self.Rpn0=(self.rhoa*0.01)/(PI*((((self.nodeD/2)+space_p1)**2)-((self.nodeD/2)**2)))
    self.Rpn1=(self.rhoa*0.01)/(PI*((((self.paraD1/2)+space_p1)**2)-((self.paraD1/2)**2)))
    self.Rpn2=(self.rhoa*0.01)/(PI*((((self.paraD2/2)+space_p2)**2)-((self.paraD2/2)**2)))
    self.Rpx=(self.rhoa*0.01)/(PI*((((self.axonD/2)+space_i)**2)-((self.axonD/2)**2)))
    self.topol()
    self.subsets()
    self.geom()
    self.biophys()

  def __del__(self):
    #print 'delete ', self
    pass

  def topol(self):
    '''
    Creates sections soma, dend, axon and connects them
    '''
    self.MYSA = [h.Section(name='MYSA[%d]' % i, cell=self) for i in range(self.paranodes1)]
    self.FLUT = [h.Section(name='FLUT[%d]' % i, cell=self) for i in range(self.paranodes2)]
    self.STIN = [h.Section(name='STIN[%d]' % i, cell=self) for i in range(self.axoninter)]
    self.node = [h.Section(name='node[%d]' % i, cell=self) for i in range(self.axonnodes)]

    for i in range(self.number):
      self.MYSA[2*i].connect(self.node[i](1))
      self.FLUT[2*i].connect(self.MYSA[2*i](1))
      self.STIN[6*i].connect(self.FLUT[2*i](1))
      self.STIN[6*i+1].connect(self.STIN[6*i](1))
      self.STIN[6*i+2].connect(self.STIN[6*i+1](1))
      self.STIN[6*i+3].connect(self.STIN[6*i+2](1))
      self.STIN[6*i+4].connect(self.STIN[6*i+3](1))
      self.STIN[6*i+5].connect(self.STIN[6*i+4](1))
      self.FLUT[2*i+1].connect(self.STIN[6*i+5](1))
      self.MYSA[2*i+1].connect(self.FLUT[2*i+1](1))
      self.node[i+1].connect(self.MYSA[2*i+1](1))
    #self.basic_shape()

  def subsets(self):
    '''
    NEURON staff
    adds sections in NEURON SectionList
    '''
    self.all = h.SectionList()
    for sec in h.allsec():
      self.all.append(sec=sec)

  def geom(self):
    for sec in self.node:
      sec.L = self.nodelength   # microns
      sec.diam = self.nodeD  # microns
      sec.nseg = 1

    for sec in self.MYSA:
      sec.L = self.paralength1   # microns
      sec.diam = self.fiberD  # microns
      sec.nseg = 1

    for sec in self.FLUT:
      sec.L = self.paralength2   # microns
      sec.diam = self.fiberD  # microns
      sec.nseg = 1

    for sec in self.STIN:
      sec.L = self.interlength   # microns
      sec.diam = self.fiberD  # microns
      sec.nseg = 1

    h.define_shape()


  def biophys(self):
    '''
    Adds channels and their parameters
    '''
    for sec in self.node:
      sec.Ra = self.rhoa/20000
      sec.cm = 2
      sec.insert('extracellular')
      sec.xraxial[1] = self.Rpn0
      sec.xg[1] = 1e10
      sec.xc[1] = 0
      sec.insert('navv1p8')
      # sec.insert('extrapump')
      # sec.insert('koi')
      # sec.insert('naoi')
      # sec.insert('nakpump')
      sec.insert('nattxs')
      sec.insert('kdr')
      sec.insert('kv1')
      sec.insert('kv3')
      sec.insert('kv4')
      # sec.insert('leak')
      sec.insert('Nav1_3')
      sec.insert('nav1p1')
      sec.insert('nav1p6')
      sec.insert('pas')
      sec.insert('iKCa')
      sec.insert('iCaL')
      sec.insert('CaIntraCellDyn')

      sec.gbar_navv1p8 = 0.1
      sec.gnabar_nav1p1 = 0.4
      sec.gnabar_nav1p6 = 0.4
      sec.gbar_kdr = 0.0#1
      sec.gkbar_kv1 = 0.1
      sec.gkbar_kv3 = 0.2
      sec.gkbar_kv4 = 0.15
      sec.gbar_nattxs = 0.2
      sec.gbar_Nav1_3 = 0.2
      sec.g_pas = 0.002
      sec.e_pas = -70
      sec.gbar_iKCa = 0.015
      sec.depth_CaIntraCellDyn = 0.1
      sec.cai_tau_CaIntraCellDyn = 2.0
      sec.cai_inf_CaIntraCellDyn = 50.0e-6
      sec.pcabar_iCaL = 0.0001

      # sec.smalla_nakpump = -0.0047891
      # sec.theta_naoi = 0.029
      # sec.theta_koi = 0.029
      sec.celsiusT_nattxs = 37
      sec.celsiusT_navv1p8 = 37
      # sec.celsiusT_nakpump = 37

      # self.add_5HTreceptors(sec, 10, 1)

    for sec in self.MYSA:
      sec.Ra = self.rhoa*(1/((self.paraD1/self.fiberD)**2))/20000
      sec.cm = 2*self.paraD1/self.fiberD
      sec.insert('pas')
      sec.g_pas = 0.001*self.paraD1/self.fiberD
      sec.e_pas = -70
      sec.insert('extracellular')
      sec.xraxial[1] = self.Rpn1
      sec.xg[1] = self.mygm/(self.nl*2)
      sec.xc[1] = self.mycm/(self.nl*2)

    for sec in self.FLUT:
      sec.Ra = self.rhoa*(1/((self.paraD2/self.fiberD)**2))/20000
      sec.cm = 2*self.paraD2/self.fiberD
      sec.insert('pas')
      sec.g_pas = 0.0001*self.paraD2/self.fiberD
      sec.e_pas = -70
      sec.insert('extracellular')
      sec.xraxial[1] = self.Rpn2
      sec.xg[1] = self.mygm/(self.nl*2)
      sec.xc[1] = self.mycm/(self.nl*2)

    for sec in self.STIN:
      sec.Ra = self.rhoa*(1/((self.axonD/self.fiberD)**2))/20000
      sec.cm = 2*self.axonD/self.fiberD
      sec.insert('pas')
      sec.g_pas = 0.0001*self.axonD/self.fiberD
      sec.e_pas = -70
      sec.insert('extracellular')
      sec.xraxial[1] = self.Rpx
      sec.xg[1] = self.mygm/(self.nl*2)
      sec.xc[1] = self.mycm/(self.nl*2)
