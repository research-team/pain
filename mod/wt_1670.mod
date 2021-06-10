
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
    SUFFIX nav11_wt_L1670W
    USEION na READ ena WRITE ina
    RANGE gnabar, ina
    RANGE minf, mtau, hinf, htau, sinf, stau, m, h, s, rtau, rinf, r
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

}

PROCEDURE rates(v (mV)) {   :Computes rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, sum, q10
    q10 = 3^((celsiusT - 30)/10)

  TABLE minf, mtau, hinf, htau, sinf, stau, rinf, rtau
  FROM -100 TO 100 WITH 200

    :"m" sodium activation system
    minf = 1/(1+exp(-(v+21.6)/6.6))
    mtau = 23.5*exp(0.008*((v-35.6)/15.5)^2)

    :"h" sodium fast inactivation system
    hinf = 0.04+0.96/(1+exp((v+55.6)/4.9))
    htau = 10.75*exp(0.05*((v-70)/15.232)^2)

    :"s" sodium slow inactivation system
    sinf = 1/(1+exp((v+53.5)/7.6)) + 0.14
    stau = 1310 :+ (7.5*exp(((v+68.7)/19.6)^2))

    :"r" recovery sodium slow inactivation system
    rinf = 1/(1+exp((v+66.9)/7.4)) + 0.15
    rtau = 3252 :+ (7.5*exp(((v+68.7)/19.6)^2))

    mtau=mtau*(1/q10)
    htau=htau*(1/q10)
    stau=stau*(1/q10)
    rtau=rtau*(1/q10)

}

UNITSON
