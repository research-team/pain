TITLE 3D diffusion without boundary
COMMENT
Author: Elena Saftenku, 2003
ENDCOMMENT
NEURON{
  POINT_PROCESS GrC_Gludif3
    POINTER affv
  RANGE    glu,rPSD,h,nu,gludir,gluspill
  RANGE Deff,meandist,rabs,alpha,h,Rmf
  RANGE inclugludir,inclugluspill, Popeak,alpha,Podir,Pospill
  RANGE ts1,td1,tm1
}
UNITS{
  (molar)=(1/liter)
  (mM)=(millimolar)
  (um)=(micron)
  (nA)=(nanoamp)
  PI=(pi) (1)
}
PARAMETER {
  nu=2(/um2): density of release sites
  rabs=0(um):radius of absorbing boundary
  h=0.02 (um):cleft width
  alpha=5 : 1/extracellular volume fraction
  Deff=0.043 (um2/ms):effective diffusion coefficient
  c0cleft = 8.769 (mM):initial [glu] in the cleft after release
  rPSD=0.11 (um) :  radius of post-synaptic density
  meandist=0.29 (um) : lowest limit of integration
  Rmf=2.9(um) :radius of mossy fiber terminal
  Popeak=0.678: adjusted peak open probability of AMPA receptors
  inclugludir=1 : inclusion of direct component
  inclugluspill=1 : inclusion of spillover component
  tm1=0 (ms) : 0.09 (ms) , shift of experimental mEPSC
  td1=0 (ms) : 0.16(ms) , shift of direct EPSC
  ts1=0 (ms) : 0.16 (ms) , shift of spillover EPSC
  vth = -40
 }
ASSIGNED{
  tx1(ms)
  gludir (mM)
  gluspill(mM)
  glu (mM)
  Podir
  Pospill
  Spike_On
  R_On
  affv (mv)
}
INITIAL {
  tx1=10000000
  glu=0
  gludir=0
  gluspill=0
}
BREAKPOINT
{
  SPK_DETECT (affv, t)
  if (t <= tx1){
    glu = 0
    gludir = 0
    gluspill = 0
    Podir = 0
    Pospill = 0
  }
  if(t > tx1) {
  UNITSOFF
    gludir = (2*c0cleft*h*alpha/sqrt(4*PI*Deff*(t-tx1)))*(1-exp(rPSD*rPSD/(4*Deff*(tx1-t))))
    if(gludir>c0cleft){gludir=c0cleft}
    gluspill = 2*nu*c0cleft*h*rPSD*rPSD*PI*alpha*(1/sqrt(4*PI*Deff*(t-tx1)))*(exp(meandist*meandist/(4*Deff*(tx1-t)))-exp(Rmf*Rmf/(4*Deff*(tx1-t))))

    glu= inclugluspill*gluspill + inclugludir*gludir
    :glu=(2*c0cleft*PI*h*alpha*rPSD*rPSD*
    :exp(0.46(um)*0.46(um)/(4*Deff*(tx1-t))))/
    :sqrt(4*4*4*PI*PI*PI*Deff*Deff*Deff*(t-tx1)*(t-tx1)*(t-tx1))

    : Experimental waveforms
    Podir=(0.94*exp((tx1-t)/0.37(ms))+0.06*exp((tx1-t)/2.2(ms))-exp((tx1-t)/0.199(ms)))/0.249*(0.43/0.484)*Popeak
    Pospill=(0.39*exp((tx1-t)/2.0(ms))+0.61*exp((tx1-t)/9.1(ms))-exp((tx1-t)/0.44(ms)))/0.682*(0.125/0.484)*Popeak
  }
}

PROCEDURE SPK_DETECT (affv (mv), t (ms)) {
	if (Spike_On == 0 && affv > vth) {
	Spike_On = 1
	tx1 = t + 1
	R_On = 1
	} else if (affv < vth) {
	Spike_On = 0
	}
}
NET_RECEIVE (weight)
{
  tx1=t
}
