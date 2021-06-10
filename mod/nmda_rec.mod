NEURON {
	POINT_PROCESS nmda_rec
		POINTER glu_m
	RANGE g, gmax, Ev
	NONSPECIFIC_CURRENT i}

UNITS{
	(nA) = (nanoamp)
	(molar) = (1/liter)
	(uM) = (micromolar)
	(mV) = (millivolt)
	(pS) = (picosiemens)
}

PARAMETER {

	kon1 = 12 (/uM /s)
	k2 = 15  (/s)
	kon3 = 6 (/uM /s)
	k4 = 30 (/s)
	k5 = 205 (/s)
	k6 =	130 (/s)
	k7 = 485 (/s)
	k8 = 2140 (/s)
	k9 = 1270 (/s)
	k10 = 275 (/s)
	k11 = 3.5 (/s)
	k12 = 0.7 (/s)
	k13 = 7.9 (/uM /s)
	k14 = 23 (/s)

	gmax = 32.4 (pS)	: conductance
	Ev = 0 (mV)

}

ASSIGNED {
	v (mV)	: voltage
	i (nA)	: current
	g  (pS)	: conductance
    glu_m (uM) : concentration
    k1 (/s)   : binding
    k3 (/s)   : binding

  B :  Mg blocking function
}

STATE {
	Ro
	R1
	C3
	C2
	C1
	O
	C5
  C4
}

INITIAL {
	Ro=1
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
  B = 1/(1+exp(-0.093*v)*(1/3.57)) :  Mg blocking function
	g = gmax*O
	i = g * (v - Ev)*B*0.001
}

KINETIC kstates{

	k1 = kon1*glu_m
  k3 = kon3*glu_m

	~ Ro <-> R1 (k1, k2)
	~ R1 <-> C3 (k3, k4)
	~ C3 <-> C2 (k5, k6)
	~ C2 <-> C1 (k7, k8)
	~ C1 <-> O (k9, k10)
	~ C3 <-> C5 (k11, k12)
	~ C2 <-> C4 (k13, k14)


	CONSERVE Ro+R1+C3+C2+C1+O+C4+C5=1
}
