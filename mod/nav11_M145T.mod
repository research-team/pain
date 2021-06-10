NEURON	{
	SUFFIX nav11_M145T
	USEION na READ ena WRITE ina
	RANGE gmax, gbar, ina, hTau, m, h, hInf, mInf
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gmax = 0.00001 (S/cm2)
	vHminf = -11.7
	vHhinf = -65
	kminf = 7.1
	khinf = 9.1
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
		mTau = 0.8 + (0.27 * exp(-v/11))
		hInf = 1/(1+exp((v-vHhinf)/khinf))
		if (v > -50){
		  hTau = 10.6 :* (1/exp(1.5*((v+65)/18.8)^2))
		}else{
		  hTau = 0.6 :* (1/exp(0.85*((v+65)/19)^2))
		}
		: hTau = 0.1 + (0.48 * exp(-v/11))
	UNITSON
}
