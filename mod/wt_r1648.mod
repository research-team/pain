
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
    SUFFIX nav11_wt_R1648H
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

}

PROCEDURE rates(v (mV)) {   :Computes rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, sum, q10
    q10 = 3^((celsiusT - 30)/10)

  TABLE minf, mtau, hinf, htau, sinf, stau, rinf, rtau
  FROM -100 TO 100 WITH 200

    :"m" sodium activation system
    minf = 1/(1+exp(-(v+19.4)/7.9))
    mtau = 0.2 :23.5*exp(0.008*((v-35.6)/15.5)^2)

    :"h" sodium fast inactivation system
    hinf = 0.04+0.96/(1+exp((v+62.7)/6.9))
    htau = 15.049*(0.0145*exp(0.0924*((v-35.5)/13.07)^2))

    :"s" sodium slow inactivation system
    sinf = 1/(1+exp((v+66.7)/7.8)) :+ 0.14
    stau = 51769*(0.0174*exp(10.69*((v+5.79)/174.4)^2))

    :"r" recovery sodium slow inactivation system
    alphar = 0.00143*exp((v-16)/144)
    betar = 0.01*exp((v+8)/10)

    rinf = alphar/(alphar + betar)
    rtau = 1/(alphar + betar)


    mtau=mtau*(1/q10)
    htau=htau*(1/q10)
    stau=stau*(1/q10)
    rtau=rtau*(1/q10)

}

UNITSON
