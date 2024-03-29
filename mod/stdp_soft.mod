COMMENT
Spike Timing Dependent Weight Adjuster
based on Song and Abbott, 2001.
Andrew Davison, UNIC, CNRS, 2003-2004
ENDCOMMENT

NEURON {
	POINT_PROCESS StdwaSA
	RANGE interval, tlast_pre, tlast_post, M, P
	RANGE deltaw, wmax, aLTP, aLTD, dw
	GLOBAL tauLTP, tauLTD, on
	POINTER wsyn
}

ASSIGNED {
	interval	(ms)	: since last spike of the other kind
	tlast_pre	(ms)	: time of last presynaptic spike
	tlast_post	(ms)	: time of last postsynaptic spike
	M			: LTD function
	P			: LTP function
	deltaw			: change in weight
	wsyn			: weight of the synapse
  dw
}

INITIAL {
	interval = 0
	tlast_pre = 0
	tlast_post = 0
	M = 0
	P = 0
	deltaw = 0
}

PARAMETER {
	tauLTP  = 20	(ms)    : decay time for LTP part ( values from           )
	tauLTD  = 20	(ms)    : decay time for LTD part ( Song and Abbott, 2001 )
	wmax    = 1000		: min and max values of synaptic weight
	aLTP    = 0.00005		: amplitude of LTP steps
	aLTD    = 0.000049	: amplitude of LTD steps
	on	= 1		: allows learning to be turned on and off globally
}

NET_RECEIVE (w) {
	if (w >= 0) {				: this is a pre-synaptic spike
		P = P*exp((tlast_pre-t)/tauLTP) + aLTP
		interval = tlast_post - t	: interval is negative
		tlast_pre = t
		deltaw = wmax * M * exp(interval/tauLTD)
    : printf("entry w=%g tlast_pre=%g tlast_post=%g deltaw=%g interval=%g \n", w, tlast_pre, tlast_post, deltaw, interval)

	} else {				: this is a post-synaptic spike
		M = M*exp((tlast_post-t)/tauLTD) - aLTD
		interval = t - tlast_pre	: interval is positive
		tlast_post = t
		deltaw = wmax * P * exp(-interval/tauLTP)
		: printf("entry w=%g tlast_pre=%g tlast_post=%g deltaw=%g interval=%g \n", w, tlast_pre, tlast_post, deltaw, interval)
	}
	if (on) {
		: printf("before wsyn=%g \n", wsyn)
    dw = wsyn

		wsyn = wsyn + deltaw
    dw = wsyn
		: printf("after wsyn=%g deltaw=%g \n", wsyn, deltaw)

		if (wsyn > wmax) {
			wsyn = wmax
      dw = wsyn
		}
		if (wsyn < 0) {
			wsyn = 0
      dw = wsyn
		}
	}
}
