from neuron import h
h.load_file('stdlib.hoc') #for h.lambda_f

import random
import math

from axon import axon

class adelta2(object):
  '''
  bio-axon class with parameters:
    axon parameters from: https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=3810&file=/MRGaxon/MRGaxon.hoc#tabs-2
  number:
    number of nodes of Ranvier
  '''
  def __init__(self, dt, diff_type):
    self.dt = dt
    self.coordinates = dict()
    self.distances = dict()
    self.diffusions = dict()
    self.diffs = []
    self.recs = []
    self.axons = []
    self.synlistinh = []
    self.synlistex = []
    self.axon1 = axon(17)
    self.axon2 = axon(10)
    self.axon3 = axon(10)
    self.synapses()
    # self.axon4 = axon(10)
    self.axons.append(self.axon2)
    self.axons.append(self.axon3)
    self.axons.append(self.axon1)
    self.x_application = 5600
    self.fast_diff = diff_type
    self.build_subsets()
    self.add_receptors()

  def build_subsets(self):
      '''
      adds sections in NEURON SectionList
      '''
      self.soma = h.Section(name='soma', cell=self)
      self.all_secs = h.SectionList()
      # for sec in self.branch:
      for axon in self.axons:
          for sec in axon.all_secs:
              self.all_secs.append(sec=sec)
      self.all_secs.append(sec=self.soma)

      # self.all = h.SectionList()
      # for sec in h.allsec():
      #   self.all.append(sec=sec)

  def position(self, axon, last_x, plus_y):
      '''
      Adds 3D position
      '''
      i = 0
      z = 0#random.randint(0, 2500)
      for sec in axon.node:
        # h.pt3dclear()
        # h.pt3dadd((last_x + axon.interlength*i),  i*plus_y, z, axon.nodeD)
        # h.pt3dadd((last_x + axon.interlength*i + axon.nodelength), i*2*plus_y, z, axon.nodeD)
        xyz = dict(x=(last_x + axon.interlength*i + axon.nodelength), y=i*2*plus_y, z=z)
        self.coordinates.update({sec: xyz})
        i+=1

  def distance(self, compartment, x_app, y_app, z_app):
      '''
      Adds distances from application for every compartment
      '''
      #self.distances.clear()
      distance = math.sqrt((x_app-self.coordinates.get(compartment).get('x'))**2 + (y_app-self.coordinates.get(compartment).get('y'))**2 + (z_app-self.coordinates.get(compartment).get('z'))**2)
      return distance

  def __del__(self):
    #print 'delete ', self
    pass

  def add_receptors(self):
      # self.axon1.node[0].connect(self.soma(0))
      self.axon2.node[0].connect(self.axon1.MYSA[len(self.axon1.MYSA)-1](0))
      self.axon3.node[0].connect(self.axon1.MYSA[len(self.axon1.MYSA)-1](0))
      # self.axon4.node[0].connect(self.axon1.MYSA[len(self.axon1.MYSA)-1](1))
      self.position(self.axon1, 0, 0)
      # self.distance(self.axon1)
      # print(self.distances.get(self.axon1.node[len(self.axon1.node)-1]))
      self.position(self.axon2, self.coordinates.get(self.axon1.node[len(self.axon1.node)-1]).get('x'), 100)
      self.position(self.axon3, self.coordinates.get(self.axon1.node[len(self.axon1.node)-1]).get('x'), -100)
      # self.position(self.axon4, self.coordinates.get(self.axon1.node[len(self.axon1.node)-1]).get('x'), 10)
      # self.distance(self.axon2)
      # self.distance(self.axon3)

      self.soma.Ra = 100
      self.soma.nseg = 10
      self.soma.L = self.soma.diam = 20
      # self.soma.insert('nav1p8')
      # self.soma.insert('nattxs')
      # self.soma.insert('nav11_L263V')
      # self.soma.insert('nav1p6')
      # self.soma.insert('pas')
      # self.soma.insert('kv1')
      # self.soma.insert('kv3')
      # self.soma.insert('kv4')
      #
      # self.soma.gbar_nav1p8 = 0.0005
      # self.soma.gnabar_nav1p6 = 0.1#random.uniform(0.35, 0.5)
      # self.soma.gnabar_nav11_L263V = 0.1#random.uniform(0.35, 0.5)
      #
      # self.soma.gkbar_kv1 = 0.002#random.uniform(0.02, 0.06)
      # self.soma.gkbar_kv3 = 0.002#random.uniform(0.02, 0.06)
      # self.soma.gkbar_kv4 = 0.0001
      # self.soma.gbar_nattxs = 0.05# random.uniform(0.3, 0.5)
      #
      # self.soma.g_pas = 0.002
      # self.soma.e_pas = -60
      # self.soma.ena = 55
      # self.soma.ek = -90
      #
      # self.soma.celsiusT_nattxs = 37
      # self.soma.celsiusT_nav1p8 = 37

      # print(self.coordinates)

      # self.add_5HTreceptors(sec, 10, 1)
      # for sec in self.axon1.node:
      #     self.add_5HTreceptors(sec, 10, 15)
          # self.add_P2Xreceptors(sec, 10, 15)
      for sec in self.axon2.node:
          self.add_5HTreceptors(sec, random.randint(10, 20), 29)
          # self.add_P2Xreceptors(sec, 10, 20)
      for sec in self.axon3.node:
          self.add_5HTreceptors(sec, random.randint(10, 20), 29)
          # self.add_P2Xreceptors(sec, 10, 20)


  def add_P2Xreceptors(self, compartment, time, g):
      '''
      Adds P2X3 receptors
      Parameters
      ----------
      compartment: section of NEURON cell
          part of neuron
      x: int
          x - coordinate of ATP application
      time: int (ms)
          time of ATP application
      g: float
          receptor conductance
      '''
      x = [13485, 13485]
      y = [1800, -1800]
      z = [0, 0, 0]

      if self.fast_diff:
          for i in range(len(x)):
              diff = h.AtP_42(compartment(0.5))
              # if i > 0:
              #     y = y*(-1)
              diff.h = self.distance(compartment, x[i], y[i], z[i])
              # print(compartment)
              # print(diff.h)
              diff.tx1 = time + i*self.dt
              diff.Deff = 0.8
              diff.c0cleft = 1
              diff.k = 0.001
              rec = h.p2x3(compartment(0.5))
              rec.gmax = random.gauss(g, g / 10)
              rec.Ev = 5
              rec2 = h.p2x2(compartment(0.5))
              rec2.gmax = 0
              rec2.Ev = -7
              h.setpointer(diff._ref_atp, 'patp', rec)
              self.recs.append(rec)
              h.setpointer(diff._ref_atp, 'patp', rec2)
              self.recs.append(rec2)
              self.diffs.append(diff)
      else:
          diff = h.AtP_slow(compartment(0.5))
          diff.h = self.distance(compartment, x[0], y[0], z[0])
          diff.tx1 = time + 0 #+ (diff.h/1250)*1000
          diff.c0cleft = 0.200
      # self.diffusions.update({diff: compartment})
          rec = h.p2x3(compartment(0.5))
          rec.gmax = random.gauss(g, g / 10)
          rec.Ev = 5
          h.setpointer(diff._ref_c0cleft, 'patp', rec)
          self.recs.append(rec)
          self.diffs.append(diff)

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
      x = [13485, 13485]
      y = [1800, -1800]
      z = [0, 0, 0]
      if self.fast_diff:
          for i in range(len(x)):
            diff = h.diff_5HT(compartment(0.5))
            diff.h = self.distance(compartment, x[i], y[i], z[i])
            # print(compartment)
            # print(diff.h)
            diff.tx1 = time + i*self.dt
            diff.a = 1
            diff.Deff = 0.4
            diff.c0cleft = 2
            rec = h.r5ht3a(compartment(0.5))
            rec.gmax = random.gauss(g, g / 10)
            h.setpointer(diff._ref_serotonin, 'serotonin', rec)
            self.recs.append(rec)
            self.diffs.append(diff)
      else:
          diff = h.slow_5HT(compartment(0.5))
          diff.h = self.distance(compartment, x[0], y[0], z[0])
          diff.tx1 = time + 0 + (diff.h/50)*10#00
          diff.c0cleft = 3
          diff.a = 10
          rec = h.r5ht3a(compartment(0.5))
          rec.gmax = random.gauss(g, g / 10)
          h.setpointer(diff._ref_serotonin, 'serotonin', rec)
          self.diffs.append(diff)
          self.recs.append(rec)

  def synapses(self):
    '''
    Adds synapses
    '''
    for i in range(10):
        s = h.GABAa_DynSyn(self.axon2.node[0](0.5)) # Inhibitory
        self.synlistinh.append(s)
        s = h.ExpSyn(self.axon2.node[0](0.5)) # Excitatory
        s.tau = 0.15
        s.e = 50
        self.synlistex.append(s)

  def connect2target(self, target):
      '''
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
      nc = h.NetCon(self.axon1.node[len(self.axon1.node)-1](1)._ref_v, target, sec = self.axon1.node[len(self.axon1.node)-1])
      nc.threshold = 10
      return nc
