def CombustionSolver(products,a_products,b,Ho,Po):
    'Solve for the Conditions at Combustion (aka Injector)'
    import EquilibriumFunctions
    P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound = EquilibriumFunctions.PHsolver(products,a_products,b,Ho,Po)
    '''print("\n\nTHERMODYNAMIC PROPERTIES - COMBUSTION")
    print("------------------------")
    print(f"Pressure:    {P} Bar")
    print(f"Temperature: {T} KG/M^3")
    print(f"Density:     {rho} K")
    print(f"H:           {H_0} KJ/KG")
    print(f"G:           {G_0} KJ/KG")
    print(f"S:           {S_0} K\n")
    print(f"M, (1/Ng):   {M} KG/KMOL")
    print(f"Cp:          {C_p_comb} KJ/(KG K)")
    print(f"GAMMAs:      {gamma}")
    print(f"SON VEL:     {a_sound} M/S")'''
    return P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound

def ThroatSolver(P_c,T_c,gammas,Ng,Ho,S0,products,a_products,b,wc=0):
    'Solve for the Conditions at the Throat'
    import EquilibriumFunctions
    import math
    Po = P_c/math.pow((gammas + 1)/2,(gammas/(gammas-1)))
    abs_error = 1
    T = T_c
    iter = 0
    while abs_error > 1e-12:
        P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound = EquilibriumFunctions.PSsolver(Po,T,Ng,S0,products,a_products,b)
        usqrd = 2*(Ho - H_0)*1000 + wc**2
        Machsqrd = usqrd/(a_sound*a_sound)
        Po = Po*(1+gamma*abs(Machsqrd))/(1+gamma)
        iter = iter + 1
        Ng = 1/M
        abs_error = abs((usqrd - a_sound*a_sound)/usqrd)
    Mach = math.sqrt(Machsqrd)
    '''print("\n\nTHERMODYNAMIC PROPERTIES - THROAT")
    print("---------------------------------")
    print(f"Pressure:    {P} Bar")
    print(f"Temperature: {T} KG/M^3")
    print(f"Density:     {rho} K")
    print(f"H:           {H_0} KJ/KG")
    print(f"G:           {G_0} KJ/KG")
    print(f"S:           {S_0} K\n")
    print(f"M, (1/Ng):   {M} KG/KMOL")
    print(f"Cp:          {C_p_comb} KJ/(KG K)")
    print(f"GAMMAs:      {gamma}")
    print(f"SON VEL:     {a_sound} M/S")
    print(f"MACH:        {round(Mach,3)}")'''
    Tpackage = [P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound,Mach,1]
    return Tpackage

def assigned_subsonic(AeAt,Pinf,Ng,T,Ho,S_0,products,a_products,b,THROAT_VALUES):
    'Give the Conditions for a Subsonic Area Ratio'
    import math
    import EquilibriumFunctions
    import CompFlow
    AeAtDes = AeAt
    gamEst = 1.225
    Mest = CompFlow.get_mach_num(AeAt,gamEst,'sub')
    PR = CompFlow.PresRatio(gamEst,Mest)
    Po = Pinf/PR
    abs_error = 1
    while abs_error > 5e-5:
        P, T, rho, H_0, G_0, S0, M, C_p_comb, gamma, a_sound = EquilibriumFunctions.PSsolver(Po,T,Ng,S_0,products,a_products,b)
        
        Ue = math.sqrt(2*(Ho - H_0)*1000)
        diff = abs(S0 - S_0)
        Mach = Ue/a_sound
        AeAt = THROAT_VALUES[2]*THROAT_VALUES[9]/(rho*Ue)
        update = (AeAt - AeAtDes)/AeAtDes*0.2
        nextP = Po - update*1e5
        Po = nextP
        abs_error = abs((AeAt - AeAtDes)/AeAtDes)
    Entry = [P, T, rho, H_0, G_0, S0, M, C_p_comb, gamma, a_sound, Mach,AeAt]
    
    ############################################################################################################################################## EntropyEquilibrium Function END
    '''print("\n\nTHERMODYNAMIC PROPERTIES - EXIT")
    print("-------------------------------")
    print(f"Pressure:    {P} Bar")
    print(f"Temperature: {T} KG/M^3")
    print(f"Density:     {rho} K")
    print(f"H:           {H_0} KJ/KG")
    print(f"G:           {G_0} KJ/KG")
    print(f"S:           {S0} K\n")
    print(f"M, (1/Ng):   {M} KG/KMOL")
    print(f"Cp:          {C_p_comb} KJ/(KG K)")
    print(f"GAMMAs:      {gamma}")
    print(f"SON VEL:     {a_sound} M/S")
    print(f"MACH:        {Mach}")
    print(f"Ae/At:       {round(AeAt,3)}")'''
    return Entry

def assigned_supersonic(AeAt,Pinf,Ng,T,Ho,S_0,products,a_products,b,THROAT_VALUES,wc=0):
    'Give the Exit conditions for an Area Ratio'
    import math
    import EquilibriumFunctions
    import CompFlow
    AeAtDes = AeAt
    gamEst = 1.225
    Mest = CompFlow.get_mach_num(AeAt,gamEst,'sup')
    PR = CompFlow.PresRatio(gamEst,Mest)
    Po = Pinf/PR
    abs_error = 1
    while abs_error > 5e-5:
        P, T, rho, H_0, G_0, S0, M, C_p_comb, gamma, a_sound = EquilibriumFunctions.PSsolver(Po,T,Ng,S_0,products,a_products,b)
        
        Ue = math.sqrt(2*(Ho - H_0)*1000 + wc**2)
        diff = abs(S0 - S_0)
        Mach = Ue/a_sound
        AeAt = THROAT_VALUES[2]*THROAT_VALUES[9]/(rho*Ue)
        update = (AeAt - AeAtDes)/AeAtDes
        nextP = Po*(1+update)
        Po = nextP
        abs_error = abs((AeAt - AeAtDes)/AeAtDes)
    Exit = [P, T, rho, H_0, G_0, S0, M, C_p_comb, gamma, a_sound,Mach,AeAt]
    
    ############################################################################################################################################## EntropyEquilibrium Function END
    '''print("\n\nTHERMODYNAMIC PROPERTIES - EXIT")
    print("-------------------------------")
    print(f"Pressure:    {P} Bar")
    print(f"Temperature: {T} KG/M^3")
    print(f"Density:     {rho} K")
    print(f"H:           {H_0} KJ/KG")
    print(f"G:           {G_0} KJ/KG")
    print(f"S:           {S0} K\n")
    print(f"M, (1/Ng):   {M} KG/KMOL")
    print(f"Cp:          {C_p_comb} KJ/(KG K)")
    print(f"GAMMAs:      {gamma}")
    print(f"SON VEL:     {a_sound} M/S")
    print(f"MACH:        {Mach}")
    print(f"Ae/At:       {round(AeAt,3)}")'''
    return Exit
def assigned_pressure(Pe,Ng,T,Ho,S_0,products,a_products,b,THROAT_VALUES):
    'Give Exit Conditions for Assigned Pressure'
    import math
    import EquilibriumFunctions
    Po = Pe
    P, T, rho, H_0, G_0, S0, M, C_p_comb, gamma, a_sound = EquilibriumFunctions.PSsolver(Po,T,Ng,S_0,products,a_products,b)   
    Ue = math.sqrt(2*(Ho - H_0)*1000)
    diff = abs(S0 - S_0)
    Mach = Ue/a_sound
    AeAt = THROAT_VALUES[2]*THROAT_VALUES[9]/(rho*Ue)
    Exit = [P, T, rho, H_0, G_0, S0, M, C_p_comb, gamma, a_sound,Mach,AeAt]
    ############################################################################################################################################## EntropyEquilibrium Function END
    '''print("\n\nTHERMODYNAMIC PROPERTIES - EXIT")
    print("-------------------------------")
    print(f"Pressure:    {P} Bar")
    print(f"Temperature: {T} KG/M^3")
    print(f"Density:     {rho} K")
    print(f"H:           {H_0} KJ/KG")
    print(f"G:           {G_0} KJ/KG")
    print(f"S:           {S0} K\n")
    print(f"M, (1/Ng):   {M} KG/KMOL")
    print(f"Cp:          {C_p_comb} KJ/(KG K)")
    print(f"GAMMAs:      {gamma}")
    print(f"SON VEL:     {a_sound} M/S")
    print(f"MACH:        {Mach}")
    print(f"Ae/At:       {round(AeAt,3)}")'''
    return Exit

def FACsolver(Po,products,a_products,b,Ho,AcAt=6.0,ExitConditions=[]):
    'Solve for conditions using Finite Area Combustor model'
    import EquilibriumFunctions
    import CompFlow
    from math import sqrt, pi
    P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound =EquilibriumFunctions.PHsolver(products,a_products,b,Ho,Po)
    Comb_Conditions = [P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound,]
    EstMACH = CompFlow.get_mach_num(AcAt,gamma,'sub')
    Pc = Po/CompFlow.PresRatio(gamma,EstMACH)
    # Begin iteration to find conditions at nozzle inlet
    P_inj = Comb_Conditions[0]*1e5
    P_inj_bar = P_inj
    abs_error = 1
    iteration = 0
    while abs_error > 1e-8:
        Pc = Pc*P_inj/P_inj_bar
        wc = sqrt((Po - Pc)/rho)
        Hc = Ho - (wc**2)/2/1000
        P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound = EquilibriumFunctions.PHsolver(products,a_products,b,Hc,Pc)
        Nozzle_Start_CONDITIONS = [P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound]
        THROAT_CONDITIONS = ThroatSolver(Pc,T,gamma,1/M,Ho,S_0,products,a_products,b)
        wc_bar = (THROAT_CONDITIONS[2]*THROAT_CONDITIONS[9]/rho)*1/AcAt
        P_inj_bar = Pc + rho*wc_bar*wc_bar
        abs_error = abs((P_inj - P_inj_bar)/P_inj)
        Mach = wc/a_sound
        Nozzle_Start_CONDITIONS.append(Mach)
        Nozzle_Start_CONDITIONS.append(AcAt)
        iteration += 1
        if iteration > 20:
            print("Iteration past 20, breaking")
            break
    Pinf = P*(1 + 0.5*(gamma-1)*Mach*Mach)**(gamma/(gamma-1))
    Ng = 1/THROAT_CONDITIONS[6]
    Comb_Conditions.append(0)
    Comb_Conditions.append("COMB")
    THROAT_CONDITIONS.append("THROAT")
    Nozzle_Start_CONDITIONS.append("COMB END")
    ENGINE = [Comb_Conditions, Nozzle_Start_CONDITIONS, THROAT_CONDITIONS]
    for N in range(0,len(ExitConditions),2):
        if ExitConditions[N] == 'Pressure':
            Exit = assigned_pressure(ExitConditions[N+1],Ng,T,Ho,S_0,products,a_products,b,THROAT_CONDITIONS)
            Exit.append("EXIT")
        elif ExitConditions[N] == 'Area':
            Exit = assigned_supersonic(ExitConditions[N+1],Po,Ng,T,Ho,S_0,products,a_products,b,THROAT_CONDITIONS)
            Exit.append("EXIT")
        elif ExitConditions[N] == 'Separation':
            Palt = ExitConditions[N+1]
            Exit =  assigned_pressure(0.4e5,Ng,T,Ho,S_0,products,a_products,b,THROAT_CONDITIONS)
            Pexit = Exit[0]*1e5
            Mach = Exit[10]
            Psep = pi/(3*Mach) * Palt
            while Pexit > Psep:
                Exit =  assigned_pressure(Psep,Ng,T,Ho,S_0,products,a_products,b,THROAT_CONDITIONS)
                Pexit = Exit[0]*1e5
                Mach = Exit[10]
                Psep = pi/(3*Mach) * Palt
            Exit =  assigned_pressure(Psep,Ng,T,Ho,S_0,products,a_products,b,THROAT_CONDITIONS)
            Pexit = Exit[0]*1e5
            Mach = Exit[10]
            Psep = pi/(3*Mach) * Palt
            Exit.append("EXIT")
        ENGINE.append(Exit)
    ENGINE.append(Pinf)
    return ENGINE

def IACsolver(Po,products,a_products,b,Ho,ExitConditions=[],EntryConditions=[]):
    'Solve for conditions assuming Infinite Area Model'
    from math import pi
    P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound = CombustionSolver(products,a_products,b,Ho,Po)
    Comb_Conditions = [P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound]
    Ng = 1/M
    THROAT_VALUES = ThroatSolver(Po,T,gamma,1/M,Ho,S_0,products,a_products,b)
    THROAT_VALUES.append("THROAT")
    Comb_Conditions.append(0)
    Comb_Conditions.append("COMB")
    ENGINE = [Comb_Conditions]
    if len(EntryConditions) > 0:
        for N in range(0,len(EntryConditions)):
            Entry = assigned_subsonic(EntryConditions[N],Po,Ng,T,Ho,S_0,products,a_products,b,THROAT_VALUES)
            Entry.append("ENTRY")
        ENGINE.append(Entry)
    ENGINE.append(THROAT_VALUES)
    for N in range(0,len(ExitConditions),2):
        if ExitConditions[N] == 'Pressure':
           Exit =  assigned_pressure(ExitConditions[N+1],Ng,T,Ho,S_0,products,a_products,b,THROAT_VALUES)
           Exit.append("EXIT")
        elif ExitConditions[N] == 'Area':
            Exit = assigned_supersonic(ExitConditions[N+1],Po,Ng,T,Ho,S_0,products,a_products,b,THROAT_VALUES)
            Exit.append("EXIT")
        elif ExitConditions[N] == 'Separation':
            Palt = ExitConditions[N+1]
            Exit =  assigned_pressure(0.4e5,Ng,T,Ho,S_0,products,a_products,b,THROAT_VALUES)
            Pexit = Exit[0]*1e5
            Mach = Exit[10]
            Psep = pi/(3*Mach) * Palt
            while Pexit > Psep:
                Exit =  assigned_pressure(Psep,Ng,T,Ho,S_0,products,a_products,b,THROAT_VALUES)
                Pexit = Exit[0]*1e5
                Mach = Exit[10]
                Psep = pi/(3*Mach) * Palt
            Exit =  assigned_pressure(Psep,Ng,T,Ho,S_0,products,a_products,b,THROAT_VALUES)
            Pexit = Exit[0]*1e5
            Mach = Exit[10]
            Psep = pi/(3*Mach) * Palt
            Exit.append("EXIT")
        ENGINE.append(Exit)
    ENGINE.append(P)
    return ENGINE


def EnginePerformance(ENGINE,OF,viewResults=True):
    'Return the values for the performance parameters of the Engine'
    Pinf = ENGINE[len(ENGINE)-1]*1e5 #Pinf (Pa)
    Pinj = ENGINE[0][0]
    P = []
    T = []
    rho = []
    H_0 = []
    G_0 = []
    S_0 =[]
    M =[]
    C_p =[]
    gamma = []
    a_sound =[]
    MACH = []
    aeat = []
    NAME = []
    # Pereformance Parameters
    A_mdot = []
    Cstar = []
    Cf_SL = []
    Cf_opt = []
    Cf_vac = []
    Isp_SL = []
    Isp_opt = []
    Isp_vac = []

    for N in range(0,len(ENGINE)-1):
        if ENGINE[N][len(ENGINE[N])-1] == 'THROAT':
            A_mdot_t = 1/(ENGINE[N][2]*ENGINE[N][10]*ENGINE[N][9])
            CstarVal = Pinf*A_mdot_t
            break



    
    for NUM in range(0,len(ENGINE)-1):
        if NUM == 0: #Chamber conditions
            P.append(ENGINE[NUM][0])
            T.append(ENGINE[NUM][1])
            rho.append(ENGINE[NUM][2])
            H_0.append(ENGINE[NUM][3])
            G_0.append(ENGINE[NUM][4])
            S_0.append(ENGINE[NUM][5])
            M.append(ENGINE[NUM][6])
            C_p.append(ENGINE[NUM][7])
            gamma.append(ENGINE[NUM][8])
            a_sound.append(ENGINE[NUM][9])
            MACH.append(ENGINE[NUM][10])
            aeat.append('')
            NAME.append(ENGINE[NUM][len(ENGINE[NUM])-1])
        else:
            P.append(ENGINE[NUM][0])
            T.append(ENGINE[NUM][1])
            rho.append(ENGINE[NUM][2])
            H_0.append(ENGINE[NUM][3])
            G_0.append(ENGINE[NUM][4])
            S_0.append(ENGINE[NUM][5])
            M.append(ENGINE[NUM][6])
            C_p.append(ENGINE[NUM][7])
            gamma.append(ENGINE[NUM][8])
            a_sound.append(ENGINE[NUM][9])
            MACH.append(ENGINE[NUM][10])
            aeat.append(ENGINE[NUM][11])
            NAME.append(ENGINE[NUM][len(ENGINE[NUM])-1])
            A_mdot.append(1/(ENGINE[NUM][2]*ENGINE[NUM][10]*ENGINE[NUM][9]))
            Cstar.append(CstarVal)
            Cf_opt.append(ENGINE[NUM][10]*ENGINE[NUM][9]/Cstar[NUM-1])
            Cf_SL.append(Cf_opt[NUM-1] + (ENGINE[NUM][0] - 1)*1e5/Pinf*aeat[NUM])
            Cf_vac.append(Cf_opt[NUM-1] + (ENGINE[NUM][0])*1e5/Pinf*aeat[NUM])
            Isp_opt.append(Cf_opt[NUM-1]*Cstar[NUM-1]/9.81)
            Isp_SL.append(Cf_SL[NUM-1]*Cstar[NUM-1]/9.81)
            Isp_vac.append(Cf_vac[NUM-1]*Cstar[NUM-1]/9.81)

    if viewResults == True:
        str0 = ''
        strline = ''
        strpinf = ''
        str1 = ''
        str2 = ''
        str3 = ''
        str4 = ''
        str5 = ''
        str6 = ''
        str7 = ''
        str8 = ''
        str9 = ''
        str10 = ''
        str11 = ''
        str12 = ''
        str13 = ''
        #length = 22*len(NAME)
        
        for NUM in range(len(NAME)):
            str0 = str0 + '{:^18}'.format(NAME[NUM])
            strline = strline + '----------------------'
            strpinf = strpinf + '{:<18.4f}'.format(Pinj/P[NUM])
            str1 = str1 + '{:<18.4f}'.format(P[NUM])
            str2 = str2 + '{:<18.4f}'.format(T[NUM])
            str3 = str3 + '{:<18.4f}'.format(rho[NUM])
            str4 = str4 + '{:<18.4f}'.format(H_0[NUM])
            str5 = str5 + '{:<18.4f}'.format(G_0[NUM])
            str6 = str6 + '{:<18.4f}'.format(S_0[NUM])
            str7 = str7 + '{:<18.4f}'.format(M[NUM])
            str8 = str8 + '{:<18.4f}'.format(C_p[NUM])
            str9 = str9 + '{:<18.4f}'.format(gamma[NUM])
            str10 = str10 + '{:<18.4f}'.format(a_sound[NUM])
            str11 = str11 + '{:<18.4f}'.format(MACH[NUM])
            if NUM == 0:
                str12 = str12 + '{:<13}'.format(aeat[NUM])    
            else:
                str12 = str12 + '{:<13.3f}'.format(aeat[NUM])


        strings = [str0,strline,strpinf,str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str13]
        
        strnam0 = '{:15}'.format(' ')
        strnameline=''
        strnpinf = '{:18}'.format('Pinj/P:')
        strnam1 = '{:18}'.format('P, BAR:')
        strnam2 = '{:18}'.format('T, K:')
        strnam3 = '{:18}'.format('rho, kg/m^3:')
        strnam4 = '{:18}'.format('H_0, kJ/kg:')
        strnam5 = '{:18}'.format('G_0, kJ/kg:')
        strnam6 = '{:18}'.format('S_0, kJ/(kg)(K):')
        strnam7 = '{:18}'.format('M, kg/kmol:')
        strnam8 = '{:18}'.format('Cp, kJ/(kg)(K):')
        strnam9 = '{:18}'.format('GAMMAs:')
        strnam10 = '{:18}'.format('SON VEL, (m/s):')
        strnam11 = '{:18}'.format('MACH:')
        #strnam12 = '{:18}'.format('Ae/At:')
        strnam13 = '{:18}'.format('\nPERFORMANCE CHARACTERISTICS')
        strname = [strnam0,strnameline,strnpinf,strnam1,strnam2,strnam3,strnam4,strnam5,strnam6,strnam7,strnam8,strnam9,strnam10,strnam11,strnam13]
        

        
        strnam1 = '{:18}'.format('Ae/At:')
        strnam2 = '{:18}'.format('C* (m/s):')
        strnam3 = '{:18}'.format('CF_SL:')
        strnam4 = '{:18}'.format('ISP_SL (s):')
        strnam5 = '{:18}'.format('CF_opt:')
        strnam6 = '{:18}'.format('ISP_opt (s):')
        strnam7 = '{:18}'.format('CL_vac:')
        strnam8 = '{:18}'.format('ISP_vac (s):')


        str1 = ''
        str2 = ''
        str3 = ''
        str4 = ''
        str5 = ''
        str6 = ''
        str7 = ''
        str8 = ''
        for NUM in range(len(NAME)):
            if NUM == 0:
                str1 = str1 + '{:<18}'.format('') 
                str2 = str2 + '{:<18}'.format('') 
                str3 = str3 + '{:<18}'.format('') 
                str4 = str4 + '{:<18}'.format('') 
                str5 = str5 + '{:<18}'.format('') 
                str6 = str6 + '{:<18}'.format('')
                str7 = str7 + '{:<18}'.format('')
                str8 = str8 + '{:<18}'.format('')
            else:
                str1 = str1 + '{:<18.3f}'.format(aeat[NUM])
                str2 = str2 + '{:<18.3f}'.format(Cstar[NUM-1])
                str3 = str3 + '{:<18.3f}'.format(Cf_SL[NUM-1])
                str4 = str4 + '{:<18.3f}'.format(Isp_SL[NUM-1])
                str5 = str5 + '{:<18.3f}'.format(Cf_opt[NUM-1])
                str6 = str6 + '{:<18.3f}'.format(Isp_opt[NUM-1])
                str7 = str7 + '{:<18.3f}'.format(Cf_vac[NUM-1])
                str8 = str8 + '{:<18.3f}'.format(Isp_vac[NUM-1])
        strval1 = [str1,str2,str3,str5,str7,str4,str6,str8]
        strname1 = [strnam1,strnam2,strnam3,strnam5,strnam7,strnam4,strnam6,strnam8]
        
        for NUM in range(len(strval1)):
            strname.append(strname1[NUM])
            strings.append(strval1[NUM])
            

    else:
        strname = 'name'
        strings = 'strings'    
    performance = [strname,strings]
    return CstarVal,Cf_SL,Cf_opt,Cf_vac,Pinf,performance


def OptimizeOF(Po,InitialReactants,OF1,OF2):
    import EquilibriumFunctions
    'Used to optimize the OF ratio for maximum efficiency'
    # Uses IAC model to optimize efficiency since it's fast
    # No exit conditions
    # Compares values of Cstar
    # varies OF on a range of 0 to 8 by 0.1 intervals
    OFs = []
    cstars = []
    length = int((OF2-OF1)/0.1)
    for NUM in range(0,length):
        OF = OF1 + round(NUM*0.1,2)
        products, a_products, b, Ho = EquilibriumFunctions.INITIAL_CONDITIONS(Po,InitialReactants,OF)
        ENGINE = IACsolver(Po,products,a_products,b,Ho)
        CstarVal,Cf_SL,Cf_opt,Cf_vac,Pinf,performance = EnginePerformance(ENGINE,OF,viewResults=False)
        OFs.append(OF)
        cstars.append(CstarVal)
        #print(f"OF: {OF}    Cstar: {CstarVal}")
        if NUM > 0:
            if CstarVal < cstars[NUM-1]:
                break
    maxCstar = max(cstars)
    for num in range(0,len(cstars)):
        if maxCstar == cstars[num]:
            maxOF = OFs[num]
            break
    #print(f"Maximum O/F Ratio at: {maxOF}")
    return maxOF


def FACconditionAcMdot(Po,products,a_products,b,Ho,AcMdot,ExitConditions=[]):
    'Solve FAC conditions using Ac/mdot values'
    desAcMdot = AcMdot
    ENGINE = FACsolver(Po,products,a_products,b,Ho)
    rho = ENGINE[1][2]
    vel = ENGINE[1][9] * ENGINE[1][10]
    AcMdot = 1/rho/vel
    rhot = ENGINE[2][2]
    velt = ENGINE[2][9] * ENGINE[2][10]
    AcAt = rhot*velt/rho/vel
    while (abs(AcMdot - desAcMdot)/desAcMdot > 1e-5):
        ENGINE = FACsolver(Po,products,a_products,b,Ho,AcAt=AcAt)
        rho = ENGINE[1][2]
        vel = ENGINE[1][9] * ENGINE[1][10]
        AcMdot = 1/rho/vel
        rhot = ENGINE[2][2]
        velt = ENGINE[2][9] * ENGINE[2][10]
        AcAt = rhot*velt*desAcMdot
    ENGINE = FACsolver(Po,products,a_products,b,Ho,AcAt=AcAt,ExitConditions=ExitConditions)
    return ENGINE

def EngineSize(Rt,AeAt,ThetaConv,ThetaN,ThetaE,AcAt=6.0,Lpercent=80,Lstar=1):
    'Rt=Throat Radius, AeAt=Exit_Area, ThetaConv - Convergent Angle, Theta_N and Theta_E are set for Rao Nozzles, Lpercent is the percentage of a conical length, Lstar is the L* value for the Chamber'
    from math import pi, sqrt, tan, cos, sin
    
    'Returns the X,Y values of the Engine (chamber included) - Units must be in meters'
    # THROAT AT VALUE (0,0)
    X = []
    Y = []
    At = Rt*Rt*pi
    Ae = At*AeAt
    Ac = At*AcAt
    Rc = sqrt(Ac/pi)
    m1 = tan(ThetaN*pi/180)
    m2 = tan(ThetaE*pi/180)
    Nx = 0.382*Rt*sin(ThetaN*pi/180)
    Ny = 1.382*Rt - 0.382*Rt*cos(ThetaN*pi/180)
    Ex = Lpercent/100*(sqrt(AeAt)-1)*Rt/tan(15*pi/180)
    Ey = Rt*sqrt(AeAt)
    C1 = Ny - m1*Nx
    C2 = Ey - m2*Ex
    Qx = (C2-C1)/(m1-m2)
    Qy = (m1*C2 - m2*C1)/(m1 - m2)
    Nlist = [Nx,Ny]
    Q = [Qx,Qy]
    E = [Ex,Ey]
    
    #Chamber Size Calulations
    Vchamber = At*Lstar
    Lcone = (Rt*(sqrt(AcAt)-1) + 1.5*Rt*(1/cos(ThetaConv*pi/180) -1))/tan(ThetaConv*pi/180)
    Vcone = pi/3*Lcone*((Rt**2) + Rc**2 + Rt*Rc)
    Vcomb = Vchamber - Vcone
    Lcombustion = Vcomb/Ac

    #Lcombustion            #Lcone(convergent only)
    #|----------------------|  |     _   
    #                       \       /
    #                        \     /
    #                         \   /
    #                          \_/
    
    # Calculate Combustion Chamber Data Points
    # Straight Section - 10 points
    
    Xcombstraight = -Lcombustion-Lcone
    #Ycombstraight = Rc
    XcombstraightEnd = -Lcone
    YcombstraightEnd = Rc
    
    for N in range(0,11):#  Points [A,B]
        #N = 10-N
        Xdiff = (XcombstraightEnd-Xcombstraight)/10*N + Xcombstraight
        Ydiff = Rc
        X.append(Xdiff)
        Y.append(Ydiff)
    # Calculate the Curves for the entry and exit points
    # Entry Curve - Circle with deg ThetaComb
    # R=1.5Rt
    XcurveEnd = -1.5*Rt*sin(ThetaConv*pi/180)
    YcurveEnd = 1.5*Rt - 1.5*Rt*cos(ThetaConv*pi/180) + Rt

    #Convergent Cone Portion
    # Connect points between XcurveEnd,YcurveEnd and XcombstraightEnd,YcombstraightEnd - 10 Points
    m = (YcombstraightEnd-YcurveEnd)/(XcombstraightEnd-XcurveEnd)
    Xselect = (XcurveEnd - XcombstraightEnd)/10
    Yselect = YcombstraightEnd + m*Xselect
    Xend = Xselect*10
    for N in range(1,10):# Points(B,C)
        Xint = Xselect*N
        Yint = YcombstraightEnd + m*Xint
        Xdiff = Xint + XcombstraightEnd
        Y.append(Yint)
        X.append(Xdiff)

    # Points [C,T]
    #Entrance Curve Start - Circle with deg ThetaConv
    for N in range(0,11):
        Theta = ThetaConv*pi/180 - ThetaConv*pi/180/10*N
        X.append(-1.5*Rt*sin(Theta))
        Y.append(1.5*Rt - 1.5*Rt*cos(Theta) + Rt)
    
    #Exit Curve Start - Circle with deg ThetaN
    #XEXITcurveEnd = .382*Rt*sin(ThetaN*pi/180)
    #YEXITcurveEnd = .382*Rt - .382*Rt*cos(ThetaN*pi/180) + Rt

    # Points (T,N)
    for N in range(1,10):
        Theta = ThetaN*pi/180/10*N
        X.append(.382*Rt*sin(Theta))
        Y.append(.382*Rt - .382*Rt*cos(Theta) + Rt)
    #Lastly, the Parabola
    
    for N in range(0,21):
        t = 0.05*N
        Xparabola = (1-t)*Nx + 2*(1-t)*t*Qx + Ex*t**2
        Yparabola = (1-t)*Ny + 2*(1-t)*t*Qy + Ey*t**2
        X.append(Xparabola)
        Y.append(Yparabola)
    '''
    print("\nDimensions:\n")
    print(f"L_chamber:    {Lcombustion} m      {Lcombustion*39.36} in")
    print(f"InjToThroat:  {Lcombustion+Lcone} m      {(Lcombustion+Lcone)*39.36} in")
    print(f"L_Engine:     {Xparabola} m     {Xparabola*39.36} in")
    print(f"Total Length: {Xparabola - Xcombstraight} m     {(Xparabola - Xcombstraight)*39.36} in")
    print(f"\nN: {Nlist}")
    print(f"Q: {Q}")
    print(f"E: {E}")'''

    return X,Y,Lcombustion,Lcone,Xparabola,Xcombstraight,Nlist,Q,E