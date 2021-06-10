NEURON {
	POINT_PROCESS ampa_rec
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

	kon1 = 1.8 (/uM /s)
	k2 = 2400  (/s)
	kon3 = 10 (/uM /s)
	k4 = 10000 (/s)
	k5 = 16000 (/s)
	k6 =	5000 (/s)
	k7 = 700 (/s)
	k8 = 150 (/s)
	k9 = 100 (/s)
	k10 = 2.1 (/s)
	k11 = 300 (/s)
	k12 = 15 (/s)
	kon13 = 10 (/uM /s)
	k14 = 1000 (/s)
	k15 = 16000 (/s)
	k16 = 12000 (/s)
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
    k13 (/s)   : binding
}

STATE {
	C
	CA
	CA2
	OA2
	D1
	D2
	D3
}

INITIAL {
	C=1
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
	g = gmax*OA2
	i = g * (v - Ev)*0.001

}

KINETIC kstates{

	k1 = kon1*glu_m
  k3 = kon3*glu_m
  k13 = kon13*glu_m

	~ C <-> CA (k1, k2)
	~ CA <-> CA2 (k3, k4)
	~ CA2 <-> OA2 (k5, k6)
	~ CA <-> D1 (k7, k8)
	~ CA2 <-> D2 (k9, k10)
	~ D1 <-> D2 (k13, k14)
	~ D2 <-> D3 (k15, k16)
	~ OA2 <-> D3 (k11, k12)


	CONSERVE C+CA+CA2+OA2+D1+D2+D3=1
}
