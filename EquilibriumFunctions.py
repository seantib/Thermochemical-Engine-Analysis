def PHsolver(products,a_products,b,Ho,Po):
    import Thermo
    import math
    import numpy
    P = Po/100000 
    num_iterations = 0
    abs_max_diff = 1
    T = 3800
    Ng = 0.1
    while abs_max_diff > 9e-2:
        if (num_iterations == 0 and len(products[0])==4):
            # upon the first iteration, create the 'B' matrix
            products,a_products,pi_vals,del_ln_Ng,del_ln_T,del_ln_n_j,lamda = Thermo.H_calculate_initial_variables(products,a_products,b,Ho,P,Ng,T,initial=True)
        else:
            products,a_products,pi_vals,del_ln_Ng,del_ln_T,del_ln_n_j,lamda = Thermo.H_calculate_initial_variables(products,a_products,b,Ho,P,Ng,T)
        
        for NUM in range(len(products)):
            ln_num_mole = math.log(products[NUM][9]) + lamda*del_ln_n_j[NUM]
            products[NUM][9] = math.exp(ln_num_mole)
        
        ln_Ng = math.log(Ng) + lamda*del_ln_Ng
        Ng = math.exp(ln_Ng)

        ln_T = math.log(T) + lamda*del_ln_T
        abs_max_diff = abs(T - math.exp(ln_T))
        T = math.exp(ln_T)

        num_iterations = num_iterations + 1

    products,a_products,pi_vals,del_ln_Ng,del_ln_T,del_ln_n_j,lamda = Thermo.H_calculate_initial_variables(products,a_products,b,Ho,P,Ng,T)
    # After convergence, Calculate Thermodynamic DATA
    H_0,S_0,G_0,C_p_comb,C_v_comb,gamma,a_sound = Thermo.Evaluate_Derivatives(products,a_products,Ng,T,P)
    rho = Po/Ng/8314/T #kg/m^3
    M = 1/Ng
    '''print("\nTHERMODYNAMIC PROPERTIES\n")
    print("\n------------------------")
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

    # Output all calculated values
    return P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound

def PSsolver(Po,T,Ng,S_0,products,a_products,b):
    import Thermo
    import math
    import numpy
    P = Po/100000 
    num_iterations = 0
    abs_max_diff = 1
    T = T
    Ng = Ng
    while abs_max_diff > 9e-2:
        if num_iterations == 0:
            # upon the first iteration, create the 'B' matrix
            products,a_products,pi_vals,del_ln_Ng,del_ln_T,del_ln_n_j,lamda = Thermo.S_calculate_initial_variables(products,a_products,b,S_0,P,Ng,T,initial=True)
        else:
            products,a_products,pi_vals,del_ln_Ng,del_ln_T,del_ln_n_j,lamda = Thermo.S_calculate_initial_variables(products,a_products,b,S_0,P,Ng,T)
        
        for NUM in range(len(products)):
            ln_num_mole = math.log(products[NUM][9]) + lamda*del_ln_n_j[NUM]
            products[NUM][9] = math.exp(ln_num_mole)
        
        ln_Ng = math.log(Ng) + lamda*del_ln_Ng
        Ng = math.exp(ln_Ng)

        ln_T = math.log(T) + lamda*del_ln_T
        abs_max_diff = abs(T - math.exp(ln_T))
        T = math.exp(ln_T)

        num_iterations = num_iterations + 1

    products,a_products,pi_vals,del_ln_Ng,del_ln_T,del_ln_n_j,lamda = Thermo.S_calculate_initial_variables(products,a_products,b,S_0,P,Ng,T)
    # After convergence, Calculate Thermodynamic DATA
    H_0,S_0,G_0,C_p_comb,C_v_comb,gamma,a_sound = Thermo.Evaluate_Derivatives(products,a_products,Ng,T,P)
    rho = Po/Ng/8314/T #kg/m^3
    M = 1/Ng

    # Output all calculated values
    return P, T, rho, H_0, G_0, S_0, M, C_p_comb, gamma, a_sound,



def INITIAL_CONDITIONS(Po,InitReactants,OF):
    'Generate the products, a_products, b, and Ho required'
    import Thermo
    #from ThermoFileSavedData import AVAILABLE_PRODUCTS
    'Po - chamber pressure (Bar), fuel(name), oxidizer(name), OF - Oxidizer to Fuel Ratio, AcAt-contract ratio, exitType either "Pressure" or "Area", Exit = # matches exitType'
    # if Exit == NONE for either type, only calculate to Throat Conditions!!!
    numConditions = len(InitReactants)

    ############### Calculate the Reactant Conditions (Enthalpy, Mass, Entropy) ############################################
    #P = Po/100000
    # First Things First, Lets pull the data on the fuel and oxidizer
    

    # Calculate Properties of Reactants, such as ni, bi, Ho
    Master_element_list = []
    a_react = []
    n = [[],[]]
    b_react = []
    b0 = []
    h_react = []
    h_j = [[],[]]
    #h0 = []    

    # Count number of fuel, number of oxidizers
    fuels = 0
    oxids = 0
    for NUM in range(0,numConditions):
        if InitReactants[NUM][3] == 'fuel':
            fuels += 1
        else:
            oxids += 1

    # Generate a Master Element List
    for NUM in range(0,numConditions):
        if InitReactants[NUM][3] == 'fuel':
            aval = 0
            name,elements,elenum,MW,hform,T,isGas = Thermo.getThermVar(InitReactants[NUM][0])
            
            for I in elements:
                if I not in Master_element_list:
                    Master_element_list.append(I)
            nj = InitReactants[NUM][2]/MW # kilogram moles per kilogram of total fuel
            nval = n[aval]
            nval.append(nj)
            C,h,S = Thermo.CalcThermoData(InitReactants[NUM][0],InitReactants[NUM][1])
            h_j[aval].append(h)
        elif InitReactants[NUM][3] == 'oxid':
            aval = 1
            name,elements,elenum,MW,hform,T,isGas = Thermo.getThermVar(InitReactants[NUM][0])
            for I in elements:
                if I not in Master_element_list:
                    Master_element_list.append(I)
            nj = InitReactants[NUM][2]/MW # kilogram moles per kilogram of total fuel
            nval = n[aval]
            nval.append(nj)
            C,h,S = Thermo.CalcThermoData(InitReactants[NUM][0],InitReactants[NUM][1])
            h_j[aval].append(h)
    





    # generate "a" matrix for the reactants - should be a list organized as a_react[ELEMENT][SPECIES]
    for NUM in range(len(Master_element_list)):
        a = []
        b = []
        a_react.append(a)
        b_react.append(b)
    # Fuels, then Oxidizers
    for NUM in range(0,numConditions):
        name,elements,elenum,MW,h,T,isGas = Thermo.getThermVar(InitReactants[NUM][0]) #Pull up data on the reactant
        for I in range(len(Master_element_list)):   #Select each element list individually
            J = Master_element_list[I]
            a_element = a_react[I]      #Selecting the correct element
            if J in elements:           #If the master element is present in the reactant
                for i in range(0,len(elements)):
                    if J == elements[i]:    #Find that location of the master element and add the number of said element
                        a_element.append(elenum[i])
            else:
                a_element.append(0)         #Otherwise, just append a zero
    
    # Generate initial b matrix
    # fuels, then oxidizers
    for NUM in range(0,2):
        for I in range(0,len(b_react)):
            b = 0
            if NUM == 0: # Calculate b for fuels
                for J in range(0,fuels):
                    b = b + (a_react[I][J] * n[NUM][J])
                b_react[I].append(b)
            else:        # Calculate b for oxid
                for J in range(fuels,fuels+oxids):
                    b = b + (a_react[I][J] * n[NUM][J-fuels])
                b_react[I].append(b)
    

    # Calculate b0
    for NUM in range(0,len(b_react)):
        bnot = (b_react[NUM][0] + OF*b_react[NUM][1])/(1+OF)
        b0.append(bnot)

    # Calculate Ho
    # fuels, then oxidizers
    for NUM in range(0,2):
        ho = 0
        if NUM == 0: # Calculate b for fuels
            for J in range(0,fuels):
                ho = ho + (h_j[NUM][J] * n[NUM][J])
            h_react.append(ho)
        else:        # Calculate b for oxid
            for J in range(0,oxids):
                ho = ho + (h_j[NUM][J] * n[NUM][J])
            h_react.append(ho)

    Ho = (h_react[0] + OF*h_react[1])/(1+OF)
    #this list will be replaced with the findProducts function
    products = Thermo.findProducts(Master_element_list)
    

    # Now that all products (gaseous) are in a list, we may start setting up the iteration equations

    # We have some things to create, first, a matrix of the number of elements per mole of each species for each element, organized by a_products = [element][species]
    a_products = []
    for i in Master_element_list:
        a = []
        for num in products:
            if i in num[1]:
                for L in range(0,len(num[1])):
                    if num[1][L] == i:
                        a.append(num[2][L])
            else:
                a.append(0)
        a_products.append(a)
    

    return products, a_products, b0, Ho