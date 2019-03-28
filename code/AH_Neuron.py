from brian2 import *

# Setup Parameters
duration = 520*ms
neuron_pop = 2
defaultclock.dt = 50*usecond #100*usecond

# Neuron Model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Neuron Parameters
area = 20000 * umetre**2
C_N = 1 * ufarad * cm**-2
Erest_N = -17.0 * mV
EK_N = -72 * mV
ENa = 55.0 * mV
EA = -75 * mV
grest_n = 0.3*msiemens/cm**2
gNa = 120.0*msiemens/cm**2
gK_N = 20*msiemens/cm**2
gA = 47.7*msiemens/cm**2
E_EPSP = 40 * mV
E_IPSP = -80 * mV
tauEPSP = 0.5 * ms
tauIPSP = 0.5 * ms
gEPSPmax_N = 2 * usiemens
refract = 1 * ms

T_on = 0.3 * second
T_off = 3 * second
E_AH = -89 * mV


# General Neuron Equations
eqs_N = '''
Irest = grest_n * (Erest_N - v) : amp / meter**2
IK = gK_N * n_n**4 * (EK_N - v) : amp / meter**2
INa= gNa * m_n**3 * h_n * (ENa - v) :  amp / meter**2
iA = gA * a**3 * b * (EA-v) :  amp/meter**2


dgAH/dt = (3 * nsiemens - gAH) / T_off + zAH : siemens
dzAH/dt = -zAH/T_on : siemens * Hz

on_time = exp(-(t-tspike)/T_on) : 1
off_time = exp(-(t-tspike)/T_off) : 1

Iepsp = gEPSP / area * (E_EPSP - v) : amp / meter**2
Iipsp = gIPSP / area * (E_IPSP - v) : amp / meter**2
Istim : amp

dn_n/dt = alpha_n * (1 - n_n) - beta_n * n_n : 1
alpha_n = (0.02 * (v / mV + 45.7)) / (1 - exp(-0.1 * (v / mV + 45.7))) / ms : Hz
beta_n = 0.25 * exp(-0.0125 * (v / mV + 55.7)) / ms : Hz

dm_n/dt = alpham_n * (1 - m_n) - betam_n * m_n : 1
dh_n/dt = alphah_n * (1 - h_n) - betah_n * h_n : 1
alpham_n = 0.38 * (v / mV + 29.7)/(1 - exp(-0.1 * (v / mV + 29.7))) / ms : Hz
betam_n = 15.2 * exp(-0.0556 * (v / mV + 54.7)) / ms : Hz
alphah_n = 0.266 * exp(-0.05*(v / mV + 48)) / ms : Hz
betah_n = 3.8 / (1 + exp(-0.1 * (v / mV + 18))) / ms : Hz

da/dt = (ainf-a)/taua : 1
ainf = ((0.0761*exp(0.0314/mV*(v+94.22*mV)))/(1+exp(0.0346/mV*(v+1.17*mV))))**(1/3) : 1
taua = 0.3622*ms+1.158*ms/(1+exp(0.0497/mV*(v+55.96*mV))) : second
db/dt = (binf-b)/taub : 1
binf = (1/(1+exp(0.0688/mV*(v+53.3*mV))))**4 : 1
taub = 1.24*ms+2.678*ms/(1+exp(0.0624/mV*(v+50*mV))) : second

dgEPSP/dt = -gEPSP/(tauEPSP)+zEPSP / ms : siemens
dzEPSP/dt = -zEPSP/(tauEPSP) : siemens
EPalpha = (1 / (1 + exp(-(v / mV + 30))) + 4) / 5 : 1

Istretch = gEPSPmax_N / area * SensoryStrength * DSTND * (E_EPSP - v) : amp / meter**2
DSTND : 1
SensoryStrength = (1 / (exp(-(Pdiam/mm-5))+49)/0.02) + 0.6 : 1

dgIPSP/dt = -gIPSP/(tauIPSP)+zIPSP / ms : siemens
dzIPSP/dt = -zIPSP/(tauIPSP) : siemens
IPalpha = (1 / (1 + exp(-(v / mV + 30))) + 4) / 5 : 1

Isum = (Irest + IK + INa + iA + Iepsp + Iipsp + Istretch) * area + Istim : amp

Pdiam : meter

Cm = Cm_area * area : farad
Cm_area : farad / meter**2
tspike : second
'''

eqs_withoutAH = '''
I_AH = 0 * nA : amp
dv/dt = (Isum) / Cm : volt
change = (Isum) / Cm : volt * Hz
'''

eqs_withAH = '''
I_AH = gAH * (E_AH - v) : amp
dv/dt = Isum / Cm + I_AH / Cm : volt
change = (Isum) / Cm + I_AH / Cm : volt * Hz
'''

eqs_N1 = eqs_N + eqs_withoutAH
eqs_N2 = eqs_N + eqs_withAH

# Neuron Implementation
N1 = NeuronGroup(1,eqs_N1,
                  threshold='v > 30*mV and change < 0 * mV * Hz',
                  reset='zAH=50*nsiemens/second',
                  refractory=refract,
                  method='exponential_euler')
N1.v = -65 * mV
mN1 = StateMonitor(N1,('v','Irest','IK','INa','iA','I_AH','tspike','Isum','Cm','gAH','zAH'),record=True)
N1.Cm_area = C_N
N1.Istim = 0 * nA
N1.tspike=-100*ms
N1.gAH = 2.4 * nsiemens
N1.zAH = 0 * nsiemens/ms
N1.Pdiam = 7 * mm
N1.DSTND = 0

N2 = NeuronGroup(1,eqs_N2,
                  threshold='v > 30*mV and change < 0 * mV * Hz',
                  reset='zAH=50*nsiemens/second',
                  refractory=refract,
                  method='exponential_euler')
N2.v = -65 * mV
mN2 = StateMonitor(N2,('v','Irest','IK','INa','iA','I_AH','tspike','Isum','Cm','gAH','zAH'),record=True)
N2.Cm_area = C_N
N2.Istim = 0 * nA
N2.tspike=-100*ms
N2.gAH = 2.4 * nsiemens
N2.zAH = 0 * nsiemens/ms
N2.Pdiam = 7 * mm
N2.DSTND = 0


# Run and Plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pre=0.5*second
stimmag =3 * nA
run(pre, report = 'text')

N1.DSTND = 0
N2.DSTND = 0
run(10*ms, report = 'text')


N1.DSTND = 0.008
N2.DSTND = 0.008

run(500*ms, report = 'text')

N1.DSTND = 0
N2.DSTND = 0
run(10*ms, report = 'text')

figure(1)
title('Cell without AH current')
subplot(2,1,1)
plot((mN1.t - pre)/ms,mN1.v[0]/mV,'-b')
ylim([-80, 60])
xlim([0, duration/ms])
ylabel('Transmembrane \n Potential (mV)')
xlabel('Time (ms)')
subplot(2,1,2)
plot((mN1.t - pre)/ms,mN1.I_AH[0]/nA,'-k',label='AH')
xlim([0, duration/ms])
ylim([-3, 1])
ylabel('AH Current (nA)')
xlabel('Time (ms)')

figure(2)
title('Cell with AH current')
subplot(2,1,1)
plot((mN2.t - pre)/ms,mN2.v[0]/mV,'-b')
ylim([-80, 60])
xlim([0, duration/ms])
ylabel('Transmembrane \n Potential (mV)')
xlabel('Time (ms)')
subplot(2,1,2)
plot((mN2.t - pre)/ms,mN2.I_AH[0]/nA,'-k',label='AH')
xlim([0, duration/ms])
ylim([-3, 1])
ylabel('AH Current (nA)')
xlabel('Time (ms)')


show()