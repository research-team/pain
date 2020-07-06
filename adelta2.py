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
  def __init__(self):
    self.coordinates = dict()
    self.distances = dict()
    self.diffusions = dict()
    self.diffs = []
    self.recs = []
    self.axons = []
    self.axon1 = axon(17)
    self.axon2 = axon(10)
    self.axon3 = axon(10)
    self.axons.append(self.axon1)
    self.axons.append(self.axon2)
    self.axons.append(self.axon3)
    self.x_application = 5600
    self.fast_diff = True
    self.build_subsets()
    self.add_receptors()

  def build_subsets(self):
      '''
      adds sections in NEURON SectionList
      '''
      self.all = h.SectionList()
      for sec in h.allsec():
        self.all.append(sec=sec)

  def position(self, axon, last_x, plus_y):
      '''
      Adds 3D position
      '''
      i = 0
      for sec in axon.node:
        h.pt3dclear()
        h.pt3dadd((last_x + axon.interlength*i),  i*plus_y, 0, axon.nodeD)
        h.pt3dadd((last_x + axon.interlength*i + axon.nodelength), i*2*plus_y, 0, axon.nodeD)
        xyz = dict(x=(last_x + axon.interlength*i + axon.nodelength), y=i*2*plus_y, z=0)
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
      self.axon2.node[0].connect(self.axon1.MYSA[len(self.axon1.MYSA)-1](1))
      self.axon3.node[0].connect(self.axon1.MYSA[len(self.axon1.MYSA)-1](1))
      self.position(self.axon1, 0, 0)
      # self.distance(self.axon1)
      # print(self.distances.get(self.axon1.node[len(self.axon1.node)-1]))
      self.position(self.axon2, self.coordinates.get(self.axon1.node[len(self.axon1.node)-1]).get('x'), 100)
      self.position(self.axon3, self.coordinates.get(self.axon1.node[len(self.axon1.node)-1]).get('x'), -100)

      # self.distance(self.axon2)
      # self.distance(self.axon3)

      print(self.coordinates)

      # self.add_5HTreceptors(sec, 10, 1)
      for sec in self.axon1.node:
          self.add_P2Xreceptors(sec, 20, 10)
      for sec in self.axon2.node:
          self.add_P2Xreceptors(sec, 20, 10)
      for sec in self.axon3.node:
          self.add_P2Xreceptors(sec, 20, 10)

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
      x = [13500, 13500]
      y = [1800, -1800]
      z = [0, 0]

      if self.fast_diff:
          for i in range(2):
              diff = h.AtP_42(compartment(0.5))
              # if i > 0:
              #     y = y*(-1)
              diff.h = self.distance(compartment, x[i], y[i], z[i])
              print(diff.h)
              diff.tx1 = time + i*35
              diff.Deff = 0.8
              diff.c0cleft = 25
              diff.k = 0.0#1
              rec = h.p2x3(compartment(0.5))
              rec.gmax = g
              rec.Ev = 2
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
          diff.h = self.distances.get(compartment)
          diff.tx1 = time + 0 + (diff.h/1250)*1000
          diff.c0cleft = 100
      # self.diffusions.update({diff: compartment})
          rec = h.p2x3(compartment(0.5))
          rec.gmax = g
          rec.Ev = 5
          h.setpointer(diff._ref_atp, 'patp', rec)
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
      if self.fast_diff:
          diff = h.diff_5HT(compartment(0.5))
          diff.h = self.distances.get(compartment)
          diff.tx1 = time
          diff.a = 0.25
          diff.Deff = 0.2
          diff.c0cleft = 3
      else:
          diff = h.slow_5HT(compartment(0.5))
          diff.h = self.distances.get(compartment)
          diff.tx1 = time + 0 + (diff.h/50)*10#00
          diff.c0cleft = 3
      rec = h.r5ht3a(compartment(0.5))
      rec.gmax = g
      h.setpointer(diff._ref_serotonin, 'serotonin', rec)
      self.diffs.append(diff)
      self.recs.append(rec)

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
      nc = h.NetCon(self.branch(1)._ref_v, target, sec = self.branch)
      nc.threshold = 10
      return nc
