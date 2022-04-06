from brian2 import *
import matplotlib as plt


def ExtracellStim(pattern, freq=0.5, amp=-1, pw=200, phase=0, dur=50):
    dtime = defaultclock.dt
    frequency = float(freq) * Hz
    amplitude = float(amp) * mA
    pulseWidth = float(pw) * us
    duration = float(dur) * second
    shift = (1 / frequency) * float(phase) / (2 * pi) / dtime

    if pattern == "Sine":
        G = arange(0, duration / dtime) * dtime
        H = sin(frequency * 2 * pi * G + shift)
        stimulus = TimedArray(H, dt=dtime)
        string = pattern + '_' + str(frequency / Hz) + 'Hz'

    elif pattern == "Pulse":
        base = zeros(int(1 / frequency / dtime))
        stim = zeros(int(duration / dtime))
        for i in range(int(0 + shift), int(pulseWidth / dtime + shift) + int(1)):
            base[i] = 1
        for k in range(int(duration * frequency)):
            for i in range(len(base)):
                stim[i + k * len(base)] = base[i]
        H = stim
        stimulus = TimedArray(H, dt=dtime)
        string = pattern + '_' + str(frequency / Hz) + 'Hz'

    else:
        H = zeros(int(1 / frequency / dtime))
        stimulus = TimedArray(H, dt=dtime)
        string = pattern

    return {'Amp': amplitude, 'Stim': stimulus, 'Duration': duration, 'Name': string}
    # ===============================================================================

def MotilityModel_CES(stimPackage,conductance = 81, intensity="moderate", mode="deterministic", pelletDiam=7):
    # Motility model with colonic electrical stimulation
    # _______________________________________________________________________________

    # Intensity parameter select
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if intensity != "low":
        SNrefract = 1 * ms
        SNthreshold = 'v > 30 * mV and norm < 0 * volt/second'
        eqs_AH = '''
        I_AH = gAH * (E_AH - v) : amp
        SensoryStrength = (1 / (exp(-(Pdiam/mm-5))+49)/0.02) + 0.6 : 1
        '''
        if intensity == "high":
            gEJPmax = 20 * usiemens
            gIJPmax = 30 * usiemens
            refract = 1 * ms
            scaleFactor = 1
        else:  # intensity == "moderate"
            gEJPmax = 12 * usiemens
            gIJPmax = 18 * usiemens
            refract = 7 * ms
            scaleFactor = 0.1
    else:  # intensity == "low"
        mode = "deterministic"
        gEJPmax = 10 * usiemens
        gIJPmax = 15 * usiemens
        refract = 10 * ms
        SNrefract = refract
        scaleFactor = 1
        SNthreshold = 'v > 30 * mV'
        eqs_AH = '''
        I_AH = 0 * amp : amp
        SensoryStrength = 1 : 1
        '''

    # Setup Parameters
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    duration = stimPackage['Duration']
    neuron_pop = 12  # Neuron cells / group
    muscle_pop = 40  # Muscle cells / group
    MP_layer = 400 * um
    CM_layer = 600 * um
    LM_layer = 0 * mm
    MEA_Tissue_dist = -0.75 * mm
    length = 10 * cm  # Model section length
    Plength = 1 * cm  # Pellet length
    Pdiam = float(pelletDiam) * mm  # Pellet diameter
    offset = 1 * cm  # X-direction offset

    # ICC Model
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # ICC Parameters
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
    ICC.v = -65 * mV
    ICC.Aicc = 0.1
    ICC.Nicc = 0.01
    for i in range(muscle_pop):
        ICC.x[i] = i * length / muscle_pop
    ICC.x += offset
    ICC.y = MP_layer
    ICC.Cm = C_ICC

    # Muscle Model
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Muscle Parameters
    Cthresh = 0 * mV  # contract threshold
    C_M = 80 * nfarad
    E_EJP = 40 * mV
    tauEJP = 2 * ms
    E_IJP = -65 * mV
    tauIJP = 5 * ms
    tau_f1 = 12 * ms
    gCaL = 0.8 * nsiemens / pfarad * C_M
    ECaL = 45 * mV
    EK = -65 * mV
    gK = 0.8 * nsiemens / pfarad * C_M
    C_CM = 80 * nfarad
    Erest_CM = -65 * mV
    Grest_CM = 300 * nsiemens

    # Smooth Muscle (circular muscle) Equations
    eqs_CM = '''
    dv/dt = (Isum) / Cm + Stim : volt
    Irest = Grest_CM * (Erest_CM - v) : amp
    Igap_icm : amp
    Isum = Irest + Igap_icm + ICaL + IK + I_ejp + I_ijp : amp
    ICaL = gCaL * d**2 * (0.8 * f1 + 0.2 * f2) * (ECaL - v) : amp
    IK = gK * q**2 * (0.38 * r1 + 0.62 * r2) * (EK - v) : amp
    I_ejp = gEJP * (E_EJP - v) : amp
    I_ijp = gIJP * (E_IJP - v) : amp
    dinf = 1 / (1 + exp(-(v / mV + 22) / 7)) : 1
    dd/dt = (dinf - d) / tau_d : 1
    tau_d = 1 * (2.29 * ms + 5.7 / (1 + ((v / mV + 29.97) / 9)**2) * ms) : second
    df1/dt = (finf - f1) / tau_f1 : 1
    finf = 1 / (1 + exp((v / mV + 38) / 7)) : 1
    df2/dt = (finf - f2) / tau_f2 : 1
    tau_f2 = 1.5 * (90.97 * (1 - (1 / ((1 + exp((v / mV + 13.96) / 45.38)) * (1 + exp(-(v / mV + 9.5) / 3.39)))))) * ms : second
    qinf = 1 / (1 + exp(-(v / mV + 18.67) / 26.66)) : 1
    rinf = 1 / (1 + exp((v / mV + 21) / 6.3)) : 1 # 21 should be 63
    tau_q = 500 / (1 + ((v / mV + 60.71) / 15.79)**2) * ms : second
    tau_r1 = 5 / (1 + ((v / mV + 62.71) / 35.86)**2) * ms : second # 5e4
    tau_r2 = 3 * ms + 2.2e5 / (1 + exp((v / mV + 22) / 4)) * ms : second #3e4
    dq/dt = (qinf - q) / tau_q : 1
    dr1/dt = (rinf - r1) / tau_r1 : 1
    dr2/dt = (rinf - r2) / tau_r2 : 1
    dgEJP/dt = -gEJP/(tauEJP) + zEJP / ms : siemens
    dzEJP/dt = -zEJP/(tauEJP) : siemens
    EJPalpha = (99 / (1 + exp(-(v / mV + 54))) + 1) / 100 : 1
    dgIJP/dt = -gIJP/(tauIJP)+zIJP / ms : siemens
    dzIJP/dt = -zIJP/(tauIJP) : siemens
    IJPalpha = (99 / (1 + exp(-(v / mV + 54))) + 1) / 100 : 1
    x : meter
    y : meter
    Stim : volt / second
    Cm : farad
    Cm_area = Cm / area : farad / meter**2
    T = 1 / (1 + exp(-(t + 6 * ms) / ms * ((v / mV - Cthresh / mV)/abs(v / mV - Cthresh / mV)))) : 1
    '''

    # CM Implementation
    CM = NeuronGroup(muscle_pop, eqs_CM, threshold='v>=Cthresh')
    CM.v = -65 * mV
    CM.gEJP = 0 * usiemens
    CM.zEJP = 0 * usiemens
    for i in range(muscle_pop):
        CM.x[i] = i * length / muscle_pop
    sCM = SpikeMonitor(CM)
    mT = StateMonitor(CM, 'T', record=True)
    CM.x += offset
    CM.y = CM_layer
    CM.Cm = C_CM

    # Neuron Model
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Neuron Parameters
    area = 20000 * umetre ** 2
    C_N = 1 * ufarad * cm ** -2
    Erest_N = -17.0 * mV
    EK_N = -72 * mV
    ENa = 55.0 * mV
    EA = -75 * mV
    grest_n = 0.3 * msiemens / cm ** 2
    gNa = 120.0 * msiemens / cm ** 2
    gK_N = 20 * msiemens / cm ** 2
    gA = 47.7 * msiemens / cm ** 2
    E_EPSP = 40 * mV
    E_IPSP = -80 * mV
    tauEPSP = 0.5 * ms
    tauIPSP = 0.5 * ms
    gEPSPmax_N = 2 * usiemens
    E_AH = -89 * mV
    T_on = scaleFactor * 0.3 * second
    T_off = scaleFactor * 3 * second

    # General Neuron Equations
    eqs_GN = '''
    Irest = grest_n * (Erest_N - v) : amp / meter**2
    IK = gK_N * n_n**4 * (EK_N - v) : amp / meter**2
    INa= gNa * m_n**3 * h_n * (ENa - v) :  amp / meter**2
    iA = gA * a**3 * b * (EA-v) :  amp/meter**2
    Iepsp = gEPSP / area * (E_EPSP - v) : amp / meter**2
    Iipsp = gIPSP / area * (E_IPSP - v) : amp / meter**2
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
    dgIPSP/dt = -gIPSP/(tauIPSP)+zIPSP / ms : siemens
    dzIPSP/dt = -zIPSP/(tauIPSP) : siemens
    IPalpha = (1 / (1 + exp(-(v / mV + 30))) + 4) / 5 : 1
    x : meter
    y : meter
    '''

    eqs_N = '''
    dv/dt = (Isum) / Cm + Stim : volt
    Isum = (Irest + IK + INa + iA + Iepsp + Iipsp) * area : amp
    Stim : volt / second
    Cm = Cm_area * area : farad
    Cm_area : farad / meter**2
    '''
    eqs_N += eqs_GN

    # Sensory-specific Equations
    eqs_SN = '''
    dv/dt = (Isum) / Cm + Stim : volt
    Istretch = gEPSPmax_N / area * SensoryStrength * DSTND * (E_EPSP - v) : amp / meter**2
    DSTND : 1
    Isum = (Irest + IK + INa + iA + Iepsp + Iipsp + Istretch) * area + I_AH: amp
    dgAH/dt = (3 * nsiemens - gAH) / T_off + zAH : siemens
    dzAH/dt = -zAH/T_on : siemens * Hz
    Stim : volt / second
    Cm = Cm_area * area : farad
    Cm_area : farad / meter**2
    norm = (Isum) / Cm : volt / second
    '''
    eqs_SN += eqs_GN
    eqs_SN += eqs_AH

    # ECMN Implementation
    ECMN = NeuronGroup(neuron_pop, eqs_N,
                       threshold='v > 30*mV',
                       refractory=refract,
                       method='exponential_euler')
    ECMN.v = -65 * mV
    for i in range(neuron_pop):
        ECMN.x[i] = i * length / neuron_pop
    ECMN.x += offset
    ECMN.y = MP_layer
    ECMN.Cm_area = C_N

    # ICMN Implementation
    ICMN = NeuronGroup(neuron_pop, eqs_N,
                       threshold='v > 30*mV',
                       refractory=refract,
                       method='exponential_euler')
    ICMN.v = -65 * mV
    for i in range(neuron_pop):
        ICMN.x[i] = i * length / neuron_pop
    ICMN.x += offset
    ICMN.y = MP_layer
    ICMN.Cm_area = C_N

    # AI Implementation
    AI = NeuronGroup(neuron_pop, eqs_N,
                     threshold='v > 30 * mV',
                     refractory=refract,
                     method='exponential_euler')
    AI.v = -65 * mV
    for i in range(neuron_pop):
        AI.x[i] = i * length / neuron_pop
    AI.x += offset
    AI.y = MP_layer
    AI.Cm_area = C_N

    # DI Implementation
    DI = NeuronGroup(neuron_pop, eqs_N,
                     threshold='v > 30 * mV',
                     refractory=refract,
                     method='exponential_euler')
    DI.v = -65 * mV
    for i in range(neuron_pop):
        DI.x[i] = i * length / neuron_pop
    DI.x += offset
    DI.y = MP_layer
    DI.Cm_area = C_N

    # SN Implementation
    SN = NeuronGroup(neuron_pop, eqs_SN,
                     threshold=SNthreshold,
                     reset='zAH+=50*nsiemens/second',
                     refractory=SNrefract,
                     method='exponential_euler')
    SN.v = -65 * mV
    for i in range(neuron_pop):
        SN.x[i] = i * length / neuron_pop
    sSN = SpikeMonitor(SN)
    SN.x += offset
    SN.y = MP_layer
    SN.Cm_area = C_N
    SN.gAH = 2.4 * nsiemens
    SN.zAH = 0 * nsiemens / ms

    # Pellet Model
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Pellet Parameters
    fSnot = 0.0005 * cm * second ** -1  # Static friction maximum
    fDnot = 0.0000025 * cm * second ** -1  # Dynamic friction maximum

    # Pellet Equations
    eqs_P = '''
    dvel/dt = forceNet / ms : meter / second
    forceNet = forceA - forceS  - forceD : meter / second
    forceA = (cos(b2) - cos(b1)) * cm * second**-1 : meter / second
    forceS = (1 - abs(mov)) * (forceA / (1 + exp(200 * ((forceA - fSnot) / (um * ms**-1)))) + fSnot / (1 + exp(-200 * ((forceA - fSnot) / (um * ms**-1) + 0.00001)))) : meter / second
    forceD = mov * fDnot : meter / second
    mov =  2 / (1 + exp(-9e9 * (vel * second / cm))) - 1 : 1
    dx/dt = vel : meter
    b1 : 1
    b2 : 1
    '''

    # Pellet Implementation
    P = NeuronGroup(1, eqs_P, threshold='x>=length')
    P.x = Plength / 2
    P.vel = 0 * um / ms
    P.b1 = 0
    P.b2 = 0
    mP = StateMonitor(P, ('x'), record=True)
    P.x += offset
    sP = SpikeMonitor(P)

    # Serosal Electrode Model
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Electrode Parameters
    num_electrodes = 1
    sigma = 0.3 * siemens / meter  # Resistivity of extracellular field
    ri = 100 * ohm
    stimCurrent = stimPackage['Stim']
    mag = -1 * mA

    # Electrode Equations
    eqs_MEA = '''
    iStim = mag * stimCurrent(t) : amp
    x : meter
    y : meter
    '''

    # Electrode Implementation
    MEA = NeuronGroup(num_electrodes, eqs_MEA)
    mMEA = StateMonitor(MEA, 'iStim', record=True)
    MEA.x = length / 2
    MEA.y = MEA_Tissue_dist
    MEA.x += offset

    # Synapses
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Gap Junctions
    GapICM = Synapses(ICC, CM, '''gGapICC = 200 * nsiemens : siemens
                                #Igap_icm_pre = gGapICC * (v_post - v_pre) : amp (summed)
                                Igap_icm_post = gGapICC * (v_pre - v_post) : amp (summed) ''')
    GapICM.connect(condition='x_pre == x_post')

    # Neuromuscular Junctions
    EJP_CM = Synapses(ECMN, CM, on_pre='zEJP_post += EJPalpha * gEJPmax')
    EJP_CM.connect(condition='abs(x_pre - x_post) <= length / (2 * muscle_pop)')

    IJP_CM = Synapses(ICMN, CM, on_pre='zIJP_post += IJPalpha * gIJPmax')
    IJP_CM.connect(condition='abs(x_pre - x_post) <= length / (2 * muscle_pop)')

    # Ascending Interneuron - Circular Motor Neuron
    AI_ECMN = Synapses(AI, ECMN, on_pre='zEPSP_post += EPalpha * gEPSPmax_N')
    AI_ECMN.connect(condition='abs((x_pre - length / neuron_pop) - x_post) <= length / (2 * neuron_pop)')

    # Descending Interneuron - Circular Motor Neuron
    DI_ICMN = Synapses(DI, ICMN, on_pre='zEPSP_post += EPalpha * gEPSPmax_N')
    DI_ICMN.connect(condition='abs((x_pre + length / neuron_pop) - x_post) <=length / (2 * neuron_pop)')

    # Sensory Neuron - Ascending Interneuron
    SN_AI = Synapses(SN, AI, on_pre='zEPSP_post += EPalpha * gEPSPmax_N')
    SN_AI.connect(condition='x_pre == x_post')

    # Sensory Neuron - Ascending Interneuron
    SN_DI = Synapses(SN, DI, on_pre='zEPSP_post += EPalpha * gEPSPmax_N')
    SN_DI.connect(condition='x_pre == x_post')

    # Circular Mucscle - Pellet Propulsion
    CM_PP = Synapses(CM, P, '''
                    #T_post = 1 / (1 + exp(-(t + 6 * ms) / ms * ((v_pre / mV - Cthresh / mV)/abs(v_pre / mV - Cthresh / mV)))) : 1 (summed)
                    b1_post = pi / 2 / (1 + exp(-(t + 6 * ms) / ms * ((v_pre / mV - Cthresh / mV)/abs(v_pre / mV - Cthresh / mV)))) / exp(15 * abs(x_pre - (x_post - Plength / 2)) / cm) / muscle_pop : 1 (summed)
                    b2_post = pi / 2 / (1 + exp(-(t + 6 * ms) / ms * ((v_pre / mV - Cthresh / mV)/abs(v_pre / mV - Cthresh / mV)))) / exp(15 * abs(x_pre - (x_post + Plength / 2)) / cm)  / muscle_pop : 1 (summed)
           	''')
    CM_PP.connect(True)

    # Pellet - Sensory Neuron
    PP_SN = Synapses(P, SN, '''
                    DSTND_post = ((Plength / 2 - abs(x_post - x_pre))/abs(Plength / 2 - abs(x_post - x_pre)) + 1) / (100 * (Plength / length) * neuron_pop) : 1 (summed)
                    ''')
    PP_SN.connect(True)

    # Serosal MEA Stimulation
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    eqs_Extracellular = '''
    deltaX = x_pre - x_post : meter
    deltaY = y_pre - y_post : meter
    Cm_syn = Cm_area_post : farad * meter**-2
    w = 1 / (4 * pi * sigma) * (2 * deltaX**2 - deltaY**2) / (deltaX**2 + deltaY**2)**(5/2) : 1 / siemens
    Stim_post = (1 / (ri * Cm_syn)) * w * iStim_pre : volt / second (summed)   
    '''

    # MEA - ICC
    MEA_ICC = Synapses(MEA, ICC, eqs_Extracellular)
    MEA_ICC.connect()

    # MEA - CM
    MEA_CM = Synapses(MEA, CM, eqs_Extracellular)
    MEA_CM.connect()

    # MEA - ECMN
    MEA_ECMN = Synapses(MEA, ECMN, eqs_Extracellular)
    MEA_ECMN.connect()

    # MEA - ICMN
    MEA_ICMN = Synapses(MEA, ICMN, eqs_Extracellular)
    MEA_ICMN.connect()

    # MEA - AI
    MEA_AI = Synapses(MEA, AI, eqs_Extracellular)
    MEA_AI.connect()

    # MEA - DI
    MEA_DI = Synapses(MEA, DI, eqs_Extracellular)
    MEA_DI.connect()

    # MEA - SN
    MEA_SN = Synapses(MEA, SN, eqs_Extracellular)
    MEA_SN.connect()

    # Stochastic and more complicated "intense" mechanisms
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if mode == "stochastic":
        sEMN = PoissonGroup(muscle_pop, scaleFactor * 1 * Hz)
        sIMN = PoissonGroup(muscle_pop, scaleFactor * 3 * Hz)
        sEJP = Synapses(sEMN, CM, on_pre='zEJP_post += EJPalpha * scaleFactor * gEJPmax')
        sIJP = Synapses(sIMN, CM, on_pre='zIJP_post += IJPalpha * scaleFactor * gIJPmax')
        sEJP.connect(j='i')
        sIJP.connect(j='i')

    if intensity != "low":
        # Extrinsic feedback implementation
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Afferent
        Afferent = NeuronGroup(neuron_pop, eqs_SN,
                               threshold='v > 30 * mV',
                               refractory=SNrefract,
                               method='exponential_euler')
        Afferent.v = -65 * mV
        for i in range(neuron_pop):
            Afferent.x[i] = i * length / neuron_pop
        sAfferent = SpikeMonitor(Afferent)
        Afferent.x += offset
        Afferent.y = MP_layer
        Afferent.Cm_area = C_N

        # Efferent
        Efferent = NeuronGroup(neuron_pop, eqs_N,
                               threshold='v > 30 * mV',
                               refractory=refract,
                               method='exponential_euler')
        Efferent.v = -65 * mV
        for i in range(neuron_pop):
            Efferent.x[i] = i * length / neuron_pop
        Efferent.x += offset
        Efferent.y = MP_layer
        Efferent.Cm_area = C_N

        # Pellet - Afferent
        PP_Afferent = Synapses(P, Afferent, '''
                        DSTND_post = ((Plength / 2 - abs(x_post - x_pre))/abs(Plength / 2 - abs(x_post - x_pre)) + 1) / (100 * (Plength / length) * neuron_pop) : 1 (summed)
                        ''')
        PP_SN.connect(True)

        # Afferent - Efferent
        Aff_Eff = Synapses(Afferent, Efferent, on_pre='zEPSP_post += EPalpha * scaleFactor * gEPSPmax_N')
        Aff_Eff.connect(condition='x_pre == x_post')
        Aff_Eff.delay = '500 * ms'

        # Efferent - ICMN - ascending inhibiton
        Eff_ICMN = Synapses(Efferent, ICMN, on_pre='zEPSP_post += EPalpha * gEPSPmax_N')
        Eff_ICMN.connect(condition='abs((x_pre - length / neuron_pop) - x_post) <= length / (2 * neuron_pop)')
        Eff_ICMN.delay = '500 * ms'

        # Descending Excitation Implementation
        DI_ECMN = Synapses(DI, ECMN, on_pre='zEPSP_post += EPalpha * gEPSPmax_N')
        DI_ECMN.connect(condition='abs((x_pre + length / neuron_pop) - x_post) <=length / (2 * neuron_pop)')
        DI_ECMN.delay = '500 * ms'

        # MEA - Afferent
        MEA_Afferent = Synapses(MEA, Afferent, eqs_Extracellular)
        MEA_Afferent.connect()

        # MEA - Efferent
        MEA_Efferent = Synapses(MEA, Efferent, eqs_Extracellular)
        MEA_Efferent.connect()

        if intensity == "high":
            EJP_CM.delay = '5*ms'
            IJP_CM.delay = '5*ms'
            AI_ECMN.delay = '100*ms'
            DI_ICMN.delay = '100*ms'
            SN_AI.delay = '5*ms'
            SN_DI.delay = '5*ms'
        else:  # intensity == "moderate"
            EJP_CM.delay = '1*ms'
            IJP_CM.delay = '1*ms'
            AI_ECMN.delay = '20*ms'
            DI_ICMN.delay = '20*ms'
            SN_AI.delay = '1*ms'
            SN_DI.delay = '1*ms'

    # Run and Plot
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    run(duration, report='text')

    if max(mP.x[0]) / cm >= length / cm:
        transittime = sP.t[0] / second
    if max(mP.x[0]) / cm < length / cm:
        transittime = "Incomplete"

    A = figure()
    plot(sSN.t / second, sSN.i * length / cm / neuron_pop, '*k')
    plot(mP.t / second, (mP.x[0] - offset) / cm, '-b')
    plot(sCM.t / second, sCM.i * length / cm / muscle_pop, '.r')
    xlim([0, duration / second])
    ylim([(length) / cm, 0])
    title(stimPackage['Name'])
    xlabel('Time (s)')
    ylabel('Position (cm)')
    # plt.close()
    return {'Fig': A, 'Transit_time': transittime, 'Name': stimPackage['Name']}

if __name__ == "__main__":
    controlStim = ExtracellStim("Control")
    conductance = [50, 81, 100, 200]
    for i in range(len(conductance)):
        control = MotilityModel_CES(controlStim,conductance[i])
        control['Fig'].show()
