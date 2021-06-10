TITLE KM Kv7 KCNQ channel
: Jose Ambros-Ingerson jose@kiubo.net 2006
: Modified by Bill Holmes with new tau calculation and for M
NEURON {
  SUFFIX 	Kv7M
  USEION 	k READ ek WRITE ik
  RANGE		gbar, g, i

  GLOBAL 	inf_n, tau_n
}

PARAMETER {
  ek 			(mV)
  gbar		= 30 	(pS/um2)

  : inf_n parameters
  gates_n	= 1
  vhalf_n	= -41.0	(mV)
  slope_n	= -7	(mV)
  needAdj = 1
  vhalfA_n 	= 0 (mV)		: adjusted for ngate power; Set in ngate_adjust()
  slopeA_n	= 0 (mV)
  v5_adj	= 0 (mV)		: for return values in ngate_adjust
  slp_adj	= 0 (mV)

  : tau_n parameters
  tauA_n	= 54	(ms)
  tauDv_n	= 0	(mV)	: Delta to vhalf_n
  tauG_n	= 0.5		: Left-right bias. range is (0,1)
  tauF_n	= 0		: Up-Down bias. range is ~ -3.5(cup-shape), -3(flat), 0(from k), 1(sharper)
  tau0_n	= 1.5 	(ms)	: minimum tau

}

STATE {
  n
}

ASSIGNED {
  v		(mV)
  celsius	(degC)
  ik 		(mA/cm2)
  g		(pS/um2)
  i		(mA/cm2)
  inf_n
  tau_n		(ms)
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  g = 0
  if( n >=0 ){	: make sure no domain error for pow. Cvode may test it
    g 	= gbar * n^gates_n
  }
  i	= g * ( v - ek ) * (1e-4)
  ik	= i
}

INITIAL {
  needAdj = 1
  rates( v )
  n = inf_n
}

UNITS {
  (mA)	= (milliamp)
  (mV)	= (millivolt)
  (pS)	= (picosiemens)
  (um)	= (micrometer)
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
  rates( v )
  n' = ( inf_n - n )/ tau_n
}

PROCEDURE rates( v (mV)){
  if( needAdj > 0 ){  : gives given values when gates_n=1
    needAdj = 0
    ngate_adjust( gates_n, vhalf_n, slope_n )
    vhalfA_n = v5_adj
    slopeA_n = slp_adj
  }
  inf_n = Boltzmann( v, vhalfA_n, slopeA_n )
  tau_n = BorgMod_tau( v, vhalfA_n, slopeA_n, tau0_n, tauA_n, tauG_n, tauF_n, tauDv_n )

}

FUNCTION Boltzmann( v (mV), v5 (mV), s (mV) ){
  Boltzmann = 1 / (1 + exp( (v - v5) / s ))
}

FUNCTION BorgMod_tau( v (mV), v5 (mV), s (mV), tau0 (ms), tauA (ms), tauG, tauF, tauDv (mV) ) (ms) {
  LOCAL kc, kr, Dv, wr, kf

:  kr = 1000
:  wr = 1000
  Dv = (v - ( v5 + tauDv ))
:  kc = kr * 10^tauF / s * 1(mV)
  kf =  10^tauF

  BorgMod_tau = tau0 + tauA * 4 * sqrt( tauG * (1-tauG))
 	        / ( exp( - Dv *(1-tauG)*kf/s ) + exp( Dv *tauG*kf/s ))
}

: Boltzmann's inverse
FUNCTION Boltz_m1( x, v5 (mV), s (mV) ) (mV) {
  Boltz_m1 = s * log( 1/x - 1 ) + v5
}

: Find parameters for a Boltzmann eq that when taken to the ngate power matches one with a single power
: return result in v5_adj and slp_adj
: We solve for exact match on two points
PROCEDURE ngate_adjust( ng, vh (mV), slp (mV) ) {
  LOCAL x1, x2, v1, v2
  x1 = 0.3
  x2 = 0.7
  v1 = Boltz_m1( x1, vh, slp )
  v2 = Boltz_m1( x2, vh, slp )
  slp_adj = (v2 - v1)/( log( (1/x2)^(1/ng) - 1 ) - log( (1/x1)^(1/ng) - 1 ) )
  v5_adj = v1 - slp_adj * log( 1 / x1^(1/ng) - 1 )
}
