NEURON	{
	SUFFIX nav11_L1670W
	USEION na READ ena WRITE ina
	RANGE gmax, gbar, ina, hTau
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gmax = 0.00001 (S/cm2)
	vHminf = -8.4
	vHhinf = -28
	kminf = 7.7
	khinf = 5.1
}

ASSIGNED	{
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gbar	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
	hInf
	hTau
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gbar = gmax*m*m*m*h
	ina = gbar*(v-ena)
}

DERIVATIVE states	{
	rates(v)
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates(v)
	m = mInf
	h = hInf
}

PROCEDURE rates(v (mV)){
	UNITSOFF
	TABLE mInf, mTau, hInf, hTau
	FROM -100 TO 100 WITH 200

		mInf = 1/(1+exp(-(v-vHminf)/kminf))
		mTau = 0.3 + (0.27 * exp(-v/11))
		hInf = 1/(1+exp((v-vHhinf)/khinf))
		hTau = 0.40 + (0.48 * exp(-v/11))
	UNITSON
}
