
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
    SUFFIX nav11_L263V
    USEION na READ ena WRITE ina
    RANGE gnabar, ina
    RANGE minf, mtau, hinf, htau, sinf, rinf, rtau, stau, m, h, s, r
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
    alphar betar
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

    : printf("v: %g,  h: %g,  r: %g, rinf: %g, rtau: %g\n", v, h, r, rinf, rtau)

}

PROCEDURE rates(v (mV)) {   :Computes rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, sum, q10
    q10 = 3^((celsiusT - 30)/10)

  TABLE minf, mtau, hinf, htau, sinf, stau, rinf, rtau
  FROM -100 TO 100 WITH 200

    :"m" sodium activation system
    minf = 1/(1+exp(-(v+26.5)/7.2))
    mtau = 0.01 :+ 0.05*(1/(1+exp((v+15)/9.7)))

    :"h" sodium fast inactivation system
    hinf = 1/(1+exp((v+55.3)/6.7))
    htau = 5.75*exp(0.12*((v-38)/15.232)^2)

    :"s" sodium slow inactivation system
    sinf = 1/(1+exp((v+55)/5.1))
    stau = 1525*exp(106*((v+146)/319)^2)

    alphar = 0.0028*exp((v-16)/144)
    betar = 0.01*exp((v+10)/10)

    rinf = alphar/(alphar + betar)
    rtau = 1/(alphar + betar)

    mtau=mtau*(1/q10)
    htau=htau*(1/q10)
    stau=stau*(1/q10)
    rtau=rtau*(1/q10)

}

UNITSON
