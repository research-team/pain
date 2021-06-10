
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (uF) = (microfarad)
    (molar) = (1/liter)
    (nA) = (nanoamp)
    (mM) = (millimolar)
    (um) = (micron)
    (S) = (siemens)
    FARADAY = 96520 (coul)
    R = 8.3134  (joule/degC)

}


NEURON {
    SUFFIX nav11_WT_Q1478K
    USEION na READ ena WRITE ina
    RANGE gnabar, ina
    RANGE minf, mtau, hinf, htau, sinf, stau, rinf, rtau, m, h, s, r
}


INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}


PARAMETER {

    celsiusT = 37 (degC)
    dt (ms)
    enat  (mV)
    gnabar (mho/cm2)
    gl (mho/cm2)
    el (mV)
}


ASSIGNED {

    v (mV)
    gnat (mho/cm2)
    ena	(mV)
    ina	(mA/cm2)
    minf hinf sinf rinf
    mtau (ms) htau (ms) stau (ms) rtau (ms)
    mexp hexp sexp
    s2
}


STATE {
    m h s r
}


BREAKPOINT {
    SOLVE states METHOD cnexp
    gnat = gnabar*m*m*m*h*s*r
    ina = gnat*(v - ena)
}


UNITSOFF

INITIAL{
    rates(v)
    m = minf
    h = hinf
    s = sinf
    r = rinf
}

DERIVATIVE states{
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau
    s' = (sinf-s)/stau
    r' = (rinf-r)/rtau
    : s2 = 0.3*(1-exp(-t/2732))
}

PROCEDURE rates(v (mV)) {   :Computes rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, sum, q10
    q10 = 3^((celsiusT - 23)/10)

  TABLE minf, mtau, hinf, htau, sinf, stau, rinf, rtau
  FROM -100 TO 100 WITH 200

    :"m" sodium activation system
    minf = 1/(1+exp(-(v+24.5)/7.9))
    mtau = 0.01 + 0.035*(1/(1+exp((v+12)/10.7)))

    :"h" sodium fast inactivation system
    hinf = 1/(1+exp((v+65.1)/6.0))
    if (v > -70){
      htau = 40.6 * (1/exp(1.5*((v+65)/18.8)^2))
    }else{
      htau = 40.6 * (1/exp(0.85*((v+65)/19)^2))
    }
    :"s" sodium slow inactivation system
    sinf = 0.06 + 0.94/(1+exp((v+75.8)/6.2))
    stau = 1258*exp(110*((v+11)/1649)^2)

    :"r" recovery sodium slow inactivation system
    rinf = 1/(1+exp((v+60.4)/8.1))
    rtau = 3425*exp(3.75*((v+12)/548.65)^2)

    mtau=mtau*(1/q10)
		htau=htau*(1/q10)
    stau=stau*(1/q10)
    rtau=rtau*(1/q10)
}

UNITSON
