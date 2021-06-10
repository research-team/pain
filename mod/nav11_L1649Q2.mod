
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
    SUFFIX nav11_L1649Q2
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
    : printf("v: %g, s: %g, sinf: %g, stau: %g \n", v, s, sinf, stau)

}

PROCEDURE rates(v (mV)) {   :Computes rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, sum, q10
    q10 = 3^((celsiusT - 30)/10)

  TABLE minf, mtau, hinf, htau, sinf, stau
  FROM -100 TO 100 WITH 200

    :"m" sodium activation system
    minf = 1/(1+exp(-(v+23.3)/7))
    mtau = 0.23 :01 + 0.05*(1/(1+exp((v+5)/9.7)))


    :"h" sodium fast inactivation system
    hinf = 0.08 + (0.92/(1+exp((v+35.3)/7.5)))
    htau = 1 + (1/((5.16*(exp(v/9.43))) + (1/((exp((v+104.42)/14.12))-1.31))))


    :"s" sodium slow inactivation system
    sinf = 1/(1+exp((v+55)/6.2))
    stau = 992.7 + (180.2*(31.53*exp(-0.53*((v+73.23)/29.6)^2)))

    mtau=mtau*(1/q10)
    htau=htau*(1/q10)
    stau=stau*(1/q10)
}

UNITSON
