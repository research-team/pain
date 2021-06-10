
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
    SUFFIX nav11_WT_L1649Q
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
    q10 = 3^((celsiusT - 30)/10)

  TABLE minf, mtau, hinf, htau, sinf, stau, rinf, rtau
  FROM -100 TO 100 WITH 200

    :"m" sodium activation system
    minf = 1/(1+exp(-(v+21.0)/6.6))
    mtau = 0.1193*exp(0.125*((v-16.338)/26.66)^2)

    :"h" sodium fast inactivation system
    hinf = 1/(1+exp((v+54.2)/5.5))
    if (v < -70){
      htau = 155.63*(0.438*exp(0.0148*((v+112)/20.63)^2))
    }else{
      htau = 42.66*(0.405*exp(-1.104*((v+57.17)/16.23)^2))
    }
    :"s" sodium slow inactivation system
    sinf = 0.01 + 0.99/(1+exp((v+59.7)/9.7))
    stau = 1327*exp(110*((v+8)/1649)^2)

    :"r" recovery sodium slow inactivation system
    rinf = 1/(1+exp((v+77.2)/8.3))
    rtau = 436*exp(44.8*((v-34)/449.28)^2)

    mtau=mtau*(1/q10)
		htau=htau*(1/q10)
    stau=stau*(1/q10)
    rtau=rtau*(1/q10)
}

UNITSON
