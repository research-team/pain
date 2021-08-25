NEURON {
	POINT_PROCESS p2x3
		POINTER patp
	RANGE K1, L1, K2, L2, K3, L3, K4, L4, R4, D4, R3, D3, R5, D5, R2, D2, M4, M4, M3, N3, M2, N2, M1, N1
	RANGE Re, AR, A2R, A3R, Ro, AD, A2D, A3D, A3Df, D
	RANGE g, gmax, Ev, i, celsiusT
	NONSPECIFIC_CURRENT i}

UNITS{
	(nA) = (nanoamp)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(pS) = (picosiemens)
}

PARAMETER {

	K1 = 120000 (/mM /s)
	L1 = 20 (/s)
	K2 = 80000 (/mM /s)
	L2 = 40 (/s)
	K3 = 40000 (/mM /s)
	L3 = 60 (/s)
	K4 = 70 (/s)
	L4 = 1 (/s)
	R4 = 0.00001 (/s)
	D4 = 0.00001 (/s)
	R3 = 0.00001 (/s)
	D3 = 0.00001 (/s)
	R2 = 0.00001 (/s)
	D1 = 0.00001 (/s)
	R1 = 0.25 (/s)
	D2 = 0.2 (/s)
	R5 = 0.001 (/s)
	D5 = 23 (/s)
	N4 = 1.0 (/s)
	M4 = 0.0001 (/mM /s)
	N3 = 0.0255 (/s)
	M3 = 8000 (/mM /s)
	N2 = 0.017 (/s)
	M2 = 16000 (/mM /s)
	N1 = 0.085 (/s)
	M1 = 24000 (/mM /s)

	celsiusT = 37
	gmax = 32.4 (pS)	: conductance
	Ev = 0 (mV)

}

ASSIGNED {
	v (mV)	: voltage
	i (nA)	: current
	g  (pS)	: conductance
  patp (mM) : concentration
  k1 (/s)   : binding
  k2 (/s)   : binding
  k3 (/s)   : binding
  m1 (/s)
  m2 (/s)
  m3 (/s)
  m4 (/s)
	r1 (/s)
	r2 (/s)
	r3 (/s)
	r4 (/s)
	r5 (/s)


	kvot_qt
}

STATE {
	Re
	AR
	A2R
	A3R
	Ro
	AD
	A2D
	A3D
	A3Df
	D
}

INITIAL {
	Re=1
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
	g = gmax*Ro
	i = g * (v - Ev)
}

KINETIC kstates{

	k1 = K1*patp
  k2 = K2*patp
  k3 = K3*patp

  m1 = M1*patp
  m2 = M2*patp
  m3 = M3*patp
  m4 = M4*patp

  kvot_qt = 1/((9^((celsiusT-20)/10)))
	r1 = R1 * kvot_qt
	r2 = R2 * kvot_qt
	r3 = R3 * kvot_qt
	r4 = R4 * kvot_qt
	r5 = R5 * kvot_qt

	~ Re <-> AR (k1, L1)
	~ AR <-> A2R (k2, L2)
	~ AR <-> AD (D2, r2)
	~ A2R <-> A3R (k3, L3)
	~ A2R <-> A2D (D3, r3)
	~ Re <-> D (D1, r1)
	~ A3R <-> Ro (K4, L4)
	~ A3R <-> A3D (D4, r4)
	~ Ro <-> A3Df (D5, r5)
	~ A3Df <-> A3D (N4, M4)
	~ A3D <-> A2D (N3, m3)
	~ A2D <-> AD (N2, m2)
	~ AD <-> D (N1, m1)


	CONSERVE Re+AR+A2R+A3R+Ro+AD+A2D+A3D+A3Df+D=1
}
