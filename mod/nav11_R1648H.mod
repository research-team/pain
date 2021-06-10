
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
    SUFFIX nav11_R1648H
    USEION na READ ena WRITE ina
    RANGE gnabar, ina
    RANGE minf, mtau, hinf, htau, sinf, stau, m, h, s
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
    minf hinf sinf
    mtau (ms) htau (ms) stau (ms)
    mexp hexp sexp
    q10
}


STATE {
    m h s
}


BREAKPOINT {
    SOLVE states METHOD cnexp
    gnat = gnabar*m*m*m*h*s
    ina = gnat*(v - ena)
}


UNITSOFF

INITIAL{
    rates(v)
    m = minf
    h = hinf
    s = sinf
}

DERIVATIVE states{
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau
    s' = (sinf-s)/stau
}

PROCEDURE rates(v (mV)) {   :Computes rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, sum
    q10 = 3^((celsiusT - 6.3)/10)


    :"m" sodium activation system
    minf = 1/(1+exp(-(v+21.2)/4.9))
    mtau = 0.15 : 0.05*(1/(1+exp((v+5)/9.7)))

    :"h" sodium fast inactivation system
    hinf = 1/(1+exp((v+39.7)/7.7))
    htau = 11.8*exp(-0.5*((v+57.4)/28.8)^2)

    :"s" sodium slow inactivation system
    sinf = 1/(1+exp((v+46.1)/5.4))
    stau = 1000 + (106.7*exp(-0.5*((v+52.7)/18.3)^2))

    mtau=mtau*(1/q10)
    htau=htau*(1/q10)
    stau=stau*(1/q10)
}

UNITSON
