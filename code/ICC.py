from brian2 import *
import matplotlib


def MPofICC(conductance = 81):
    # Setup Parameters
    duration = 30*second
    neuron_pop = 12               # Neuron cells / group
    muscle_pop = 40               # Muscle cells / group
    MP_layer = 400 * um
    CM_layer = 600 * um
    LM_layer = 0 * mm
    MEA_Tissue_dist = -0.75 * mm
    length = 10 * cm              # Model section length
    Plength = 1 * cm              # Pellet length
    Pdiam = 7 * mm       # Pellet diameter
    offset = 1 * cm               # X-direction offset
    defaultclock.dt = 50*usecond #100*usecond

    #moderate intensity
    SNrefract = 1 * ms
    SNthreshold = 'v > 30 * mV and norm < 0 * volt/second'
    eqs_AH = '''
    I_AH = gAH * (E_AH - v) : amp
    SensoryStrength = (1 / (exp(-(Pdiam/mm-5))+49)/0.02) + 0.6 : 1
    '''
    gEJPmax = 12 * usiemens
    gIJPmax = 18 * usiemens
    refract = 7 * ms
    scaleFactor = 0.1

    # ICC Parameters
    area = 20000 * umetre**2
    C_ICC = 4 * nfarad
    Erest_ICC = -65 * mV
    Grest_ICC = conductance * nsiemens
    Eprim = -20 * mV
    Gprim = 35 * usiemens
    tauAicc = 0.15 * second

    # ICC Equations
    eqs_ICC = '''
        dv/dt = (Isum) / Cm + Stim : volt
        Irest = Grest_ICC * (Erest_ICC - v) : amp
        Iprim = gprim * (Eprim - v) : amp
        Igap_icm : amp
        Igap_ilm : amp
        I_ejp : amp
        I_ijp : amp
        Stim : volt / second
        gprim = Gprim * Aicc * Nicc : siemens
        dAicc/dt = (Aicc_inf - Aicc) / tauAicc : 1
        Aicc_inf = 1 / (1 + exp(-55.8 - v / mV)) : 1
        dNicc/dt = (Nicc_inf - Nicc) / tauNicc : 1
        Nicc_inf = 1/(1 + exp(v / mV + 61.3)) : 1
        tauNicc = 0.25 * second + 25.75 * second / (1 + exp(2 * (v / mV + 56.3))) : second
        x : meter
        y : meter
        Isum = Irest + Iprim + Igap_icm + Igap_ilm + I_ejp + I_ijp : amp
        Cm : farad
        Cm_area = Cm / area : farad / meter**2
        '''

    # ICC Implementation
    ICC = NeuronGroup(muscle_pop, eqs_ICC)
    mICC = StateMonitor(ICC,('v','Irest','Iprim','Igap_ilm','Isum','Cm','Aicc','Nicc'),record=True)
    ICC.v = -65 * mV
    ICC.Aicc = 0.1
    ICC.Nicc = 0.01
    for i in range(muscle_pop):
        ICC.x[i] = i * length / muscle_pop
    ICC.x += offset
    ICC.y = MP_layer
    ICC.Cm = C_ICC

    pre=10*second
    stimmag =3 * nA
    run(pre, report = 'text')
    run(10*second, report = 'text')
    run(10*second, report = 'text')
    run(10*second, report = 'text')

    A = figure()
    title('ICC without Stimulation')
    plot((mICC.t - pre)/second,mICC.v[0]/mV,'-b')
    ylim([-80, 60])
    xlim([0, duration/second])
    ylabel('Transmembrane \n Potential (mV)')
    xlabel('Time (s)')

    return A

conductance = [50, 81, 100, 200]
for i in range(len(conductance)):
    MPofICC(conductance[i]).show()
