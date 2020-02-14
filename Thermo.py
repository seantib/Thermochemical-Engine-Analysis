# used to hold the structure that the rest of the programs use
# Needs to import the thermo.txt file in NASA's CEA program, create elements and structures based off that
# Each element needs it's own structure, to define and work with the programs easier.
# Would be a lot easier to read from the thermo.txt file used by CEA
# Perform a search in the file, parse it, and output it as a set of variables
# REQUIREMENTS: List of Elements, Number of Elements (total), Enthalpies, T-constants

# FORMAT
#############################################################################################
## LINE                     Constants                               Columns          Parsed?
##  1                 Species Name or formula                         1-24              %

##  2                 Number of T-intervals                           2                 %
##                    Optional Identification Code                    4-9               %
##                    Chemicals Formulas, Symbols, Numbers            11-50             %
##                    Zero for gas, Non-Zero for condensed phases     52                %
##                    Molecular Weight                                53-65             %
##                    Heat of Formation at 298.15K, J/mol             66-80             %

##  3                 Temperature Range                               2-21              %
##                    Number of coefficients for Cp^o R               23                %
##                    T exponents for empirical equation for Cp^o R   24-63             %
##                    H^o (298.15) - H^o (0), J/mol                   66-80

##  4                 First Five coefficients for Cp^o R              1-80              %

##  5                 Last Three Coefficients for Cp^o R              1-48              %
##                    Integration constants b1 and b2                 49-80             %

## Repeat 3, 4, and 5 for each interval
#################################################################################################
import os
import math
import numpy
#from ThermoFileSavedData import AVAILABLE_PRODUCTS
from ALLELEMENTS import ALL_SPECIES_PROPERTIES
import json # Used to import the ALL_ELEMENT_SPECIES file - Will need to adjust getThermVar and findProducts


def getThermVar(name):
    'output name, elements, element numbers, MW, heat of formation, Temperature values (including exponents and coefficients) and isGas'
    
    for NUM in range(len(ALL_SPECIES_PROPERTIES[0])):
        if name == ALL_SPECIES_PROPERTIES[0][NUM]:
            speciesname = ALL_SPECIES_PROPERTIES[0][NUM]
            elements = ALL_SPECIES_PROPERTIES[1][NUM]
            elementnum = ALL_SPECIES_PROPERTIES[2][NUM]
            MW = ALL_SPECIES_PROPERTIES[3][NUM]
            Hform = ALL_SPECIES_PROPERTIES[4][NUM]
            T = ALL_SPECIES_PROPERTIES[5][NUM]
            zeroISgas = ALL_SPECIES_PROPERTIES[6][NUM]
            break           

    return speciesname, elements, elementnum, MW, Hform, T, zeroISgas

def findProducts(Master_Element_List,AVAILABLE_PRODUCTS=ALL_SPECIES_PROPERTIES):
    'find all possible products for the elements provided, list names'
    Names = AVAILABLE_PRODUCTS[0]
    Ele = AVAILABLE_PRODUCTS[1]
    Elenum = AVAILABLE_PRODUCTS[2]
    isgas = AVAILABLE_PRODUCTS[6]
    ProductNames = []
    ProductISgas = []
    elementlist = Master_Element_List
    # search the file between line , and point which says "END PRODUCTS", including fuel and oxidizer, just in case not all is combusted
    for n in range(0,len(Names)):
        if (len(Ele[n])<=len(elementlist) and Names[n] != 'error, not in list'): #if (Names[n] != fspeciesname and Names[n] != ospeciesname and len(Ele[n])<=len(elementlist) and Names[n] != 'error, not in list'):
            check = []
            for i in range(0,len(Ele[n])):
                if Ele[n][i] in elementlist:
                    check.append(1)
            if len(check) == len(Ele[n]):
                ProductNames.append(Names[n])
                ProductISgas.append(isgas[n])
    products = []
    for n in range(0,len(ProductNames)):
        i = ProductNames[n]
        if ProductISgas[n] == 0:
            name,elements,elenum,MW,h,T,isGas = getThermVar(i)
            add = [name,elements,elenum,MW]
            products.append(add)
    

    return products

def CalcThermoData(name,Temp):
    'Calculate Cp, H, and S, IF liquid, return only H (all units in J/mol)'
    speciesname, elements, elementnum, MW, Hform, T, zeroISgas = getThermVar(name)
    # check to see if the Temp is in the range of the variables
    R = 8.314
    flag = 0
    for n in range(0,len(T)):
        if (Temp >= T[n][0] and Temp <= T[n][1]):
            rng = n
            flag = 1
            break
        elif (len(T) == 1 and Temp == T[0][0]):
            rng = 0
            flag = 1
            break
    if (flag == 0):
        print(f"Temp Range does not exist for {name}, try something in at {T[0][0]} K")
        print("Exiting Program")
        quit()
    # Output the Constant Values
    if len(T) > 1:
        a1 = T[rng][11]
        a2 = T[rng][12]
        a3 = T[rng][13]
        a4 = T[rng][14]
        a5 = T[rng][15]
        a6 = T[rng][16]
        a7 = T[rng][17]
        #a8 = T[rng][18]
        b1 = T[rng][19]
        b2 = T[rng][20]

        # Output the Exponents
        e1 = T[rng][3]
        e2 = T[rng][4]
        e3 = T[rng][5]
        e4 = T[rng][6]
        e5 = T[rng][7]
        e6 = T[rng][8]
        e7 = T[rng][9]
        #e8 = T[rng][10]

        Cp = (a1*math.pow(Temp,e1) + a2*math.pow(Temp,e2) + a3*math.pow(Temp,e3) + a4*math.pow(Temp,e4) + a5*math.pow(Temp,e5) + a6*math.pow(Temp,e6) + a7*math.pow(Temp,e7))*R
        H = (-a1*math.pow(Temp,e1) + a2*math.pow(Temp,e2)*math.log(Temp) + a3*math.pow(Temp,e3) + a4*math.pow(Temp,e4)/2 + a5*math.pow(Temp,e5)/3 + a6*math.pow(Temp,e6)/4 + a7*math.pow(Temp,e7)/5 + b1/Temp)*Temp*R
        S = (-a1/2*math.pow(Temp,e1) - a2*math.pow(Temp,e2) + a3*math.log(Temp) + a4*math.pow(Temp,e4) + a5/2*math.pow(Temp,e5) + a6*math.pow(Temp,e6)/3 + a7*math.pow(Temp,e7)/4 + b2)*R
        return Cp, H, S
    else:
        # for liquids!
        H = Hform
        Cp = 'Error'
        S = 'Error'
        return Cp, H, S


def H_calculate_initial_variables(products,a_products,b,Ho,P,Ng,T,initial=False):
    'Calculates pi_1 ... pi_n, delta_LN_Ng, delta_T, and delta_LN_n_j'
    #Calculate the values of the species, including Gibbs Free Energy
    R = 8.314 #kJ/mol K
    if (initial == True):
        num_mole = Ng/len(products)
        for n in range(0,len(products)):
            i = products[n]
            #Cpf,Hf,Sf = CalcThermoData(i[0],298.15)
            Cp_j,H_j,S0_j = CalcThermoData(i[0],T)
            S_j = S0_j - R*math.log(num_mole/Ng) - R*math.log(P)
            Gibbs0_j = H_j- T*S0_j
            Gibbs_j = (Gibbs0_j/R/T + math.log(num_mole) - math.log(Ng) + math.log(P))*R*T
            i.append(Cp_j)
            i.append(H_j)
            i.append(Gibbs0_j)
            i.append(Gibbs_j)
            i.append(S_j)
            i.append(num_mole)
            i.append(S0_j)
    
    else:
        for n in range(0,len(products)):
            i = products[n]
            #Cpf,Hf,Sf = CalcThermoData(i[0],298.15)
            Cp_j,H_j,S0_j = CalcThermoData(i[0],T)
            S_j = S0_j - R*math.log(i[9]/Ng) - R*math.log(P)
            Gibbs0_j = H_j- T*S0_j
            Gibbs_j = (Gibbs0_j/R/T + math.log(i[9]) - math.log(Ng) + math.log(P))*R*T
            i[4] = Cp_j
            i[5] = H_j
            i[6] = Gibbs0_j
            i[7] = Gibbs_j
            i[8] = S_j
            i[10] = S0_j
    #step 1: set up matrix to solve the four equations w/4 unknowns
    len_square_matrix = len(b) + 2
    A = []
    for I in range(len_square_matrix):
        anum = []
        for J in range(len_square_matrix):
            if I < len(b):
                if J < len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        val = a_products[I][NUM]*a_products[J][NUM]*products[NUM][9]
                        sumval = sumval + val
                    anum.append(sumval)
                elif (J==len(b)):
                    sumval = 0
                    for NUM in range(len(products)):
                        val = a_products[I][NUM]*products[NUM][9]
                        sumval = sumval + val
                    anum.append(sumval)
                else:
                    sumval = 0
                    for NUM in range(len(products)):
                        val = a_products[I][NUM]*products[NUM][9]*products[NUM][5]/R/T
                        sumval = sumval + val
                    anum.append(sumval)
            elif (I >= len(b) and I < len_square_matrix-1):
                if J < len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        val = a_products[J][NUM]*products[NUM][9]
                        sumval = sumval + val
                    anum.append(sumval)
                elif J == len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        val = products[NUM][9]
                        sumval = sumval + val
                    value = sumval - Ng
                    anum.append(value)
                else:
                    sumval = 0
                    for NUM in range(len(products)):
                        val = products[NUM][9]*products[NUM][5]/R/T
                        sumval = sumval + val
                    anum.append(sumval)
            else:
                if J < len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        sumval = sumval + a_products[J][NUM]*products[NUM][9]*products[NUM][5]/R/T
                    anum.append(sumval)
                elif J == len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        sumval = sumval + products[NUM][9]*products[NUM][5]/R/T
                    anum.append(sumval)
                else:
                    sumval1 = 0
                    sumval2 = 0
                    for NUM in range(len(products)):
                        sumval1 = sumval1 + products[NUM][9]*products[NUM][4]/R
                        sumval2 = sumval2 + products[NUM][9]*products[NUM][5]/R/T*products[NUM][5]/R/T
                    value = sumval1 + sumval2
                    anum.append(value)
        A.append(anum)
    
    #step 2: set up the solution to this matrix (i.e. Ax = 'B')
    B = []
    for I in range(len_square_matrix):
        if I < len(b):
            bval = 0
            sum1 = 0
            for NUM in range(len(products)):
                sum1 = sum1 + a_products[I][NUM]*products[NUM][9]*products[NUM][7]/R/T
                bval = bval + a_products[I][NUM]*products[NUM][9]
            Bvalue = b[I] - bval + sum1
            B.append(Bvalue)
        elif I == len(b):
            njval = 0
            sum1 = 0
            for NUM in range(len(products)):
                sum1 = sum1 + products[NUM][9]*products[NUM][7]/R/T
                njval = njval + products[NUM][9]
            Bvalue = Ng - njval + sum1
            B.append(Bvalue)
        else:
            sum1 = 0
            sum2 = 0
            for NUM in range(len(products)):
                sum1 = sum1 + products[NUM][9]*products[NUM][5]
                sum2 = sum2 + products[NUM][9]*products[NUM][5]/R/T*products[NUM][7]/R/T
            val = Ho/R/T - sum1/R/T + sum2
            B.append(val)

    solution_vector = numpy.linalg.solve(A,B)    
    pi_vals = []
    for num in range(len(b)):
        pi_vals.append(solution_vector[num])
    del_ln_Ng = solution_vector[len(b)]
    del_ln_T = solution_vector[len(solution_vector)-1]
    #step 3: create new matrix using equations in CEA
    C = []
    for I in range(len(products)):
        cnum = []
        for J in range(len(products)):
            if I == J:
                cnum.append(1)
            else:
                cnum.append(0)
        C.append(cnum)
    
    #step 4: create new solution to this matrix
    D = []
    for I in range(len(products)):
        dval = 0
        for J in range(len(b)):
            dval = dval + a_products[J][I]*pi_vals[J]
        dval = -products[I][7]/R/T + dval + del_ln_Ng + products[I][5]/R/T*del_ln_T
        D.append(dval)
    
    del_ln_n_j = numpy.linalg.solve(C,D)

    C1 = []
    C2 = []
    check = 0
    for NUM in range(len(del_ln_n_j)):
        C1.append((del_ln_n_j[NUM]))
        c2val = abs((-math.log(products[NUM][9]/Ng) - 9.2103404)/(del_ln_n_j[NUM] - del_ln_Ng))
        C2.append(c2val)
        checkval = math.log(products[NUM][9]/Ng)
        if checkval <= -18.420681:
            check = 0
    maxC1 = max(C1)

    L1 = 2/max(5*abs(del_ln_T),5*abs(del_ln_Ng),maxC1)
    L2 = min(C2)
    if check == 1:
        lamda = min(1,L1,L2)
    else:
        lamda = min(1,L1)
    
    # An added step before the products are finished; This may have to go in Evaluate Derivatives
    # Should search through the products, remove any that have a n_j less than 10^-18
    newProducts = []
    newa_products = []
    for num in range(0,len(a_products)):
        a = []
        for n in range(0,len(products)):
            if products[n][9] > 1e-999:
                a.append(a_products[num][n])
        newa_products.append(a)
    for species in products:
        if species[9] > 1e-999:
            newProducts.append(species)
        else:
            print(species)
    products = newProducts
    a_products = newa_products
    return products,a_products,pi_vals,del_ln_Ng,del_ln_T,del_ln_n_j,lamda

def S_calculate_initial_variables(products,a_products,b,So,P,Ng,T,initial=False):
    'Calculates pi_1 ... pi_n, delta_LN_Ng, delta_T, and delta_LN_n_j'
    #Calculate the values of the species, including Gibbs Free Energy
    R = 8.314 #kJ/mol K
    if (initial == True):
        for n in range(0,len(products)):
            i = products[n]
            #Cpf,Hf,Sf = CalcThermoData(i[0],298.15)
            Cp_j,H_j,S0_j = CalcThermoData(i[0],T)
            S_j = S0_j - R*math.log(i[9]/Ng) - R*math.log(P)
            Gibbs0_j = H_j- T*S0_j
            Gibbs_j = (Gibbs0_j/R/T + math.log(i[9]) - math.log(Ng) + math.log(P))*R*T
            i[4] = Cp_j
            i[5] = H_j
            i[6] = Gibbs0_j
            i[7] = Gibbs_j
            i[8] = S_j
            i[10] = S0_j
    
    else:
        for n in range(0,len(products)):
            i = products[n]
            #Cpf,Hf,Sf = CalcThermoData(i[0],298.15)
            Cp_j,H_j,S0_j = CalcThermoData(i[0],T)
            Gibbs0_j = H_j- T*S0_j
            Gibbs_j = (Gibbs0_j/R/T + math.log(i[9]) - math.log(Ng) + math.log(P))*R*T
            S_j = S0_j - R*math.log(i[9]/Ng) - R*math.log(P)
            i[4] = Cp_j
            i[5] = H_j
            i[6] = Gibbs0_j
            i[7] = Gibbs_j
            i[8] = S_j
            i[10] = S0_j
    #step 1: set up matrix to solve the four equations w/4 unknowns
    len_square_matrix = len(b) + 2
    A = []
    for I in range(len_square_matrix):
        anum = []
        for J in range(len_square_matrix):
            if I < len(b):
                if J < len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        val = a_products[I][NUM]*a_products[J][NUM]*products[NUM][9]
                        sumval = sumval + val
                    anum.append(sumval)
                elif (J==len(b)):
                    sumval = 0
                    for NUM in range(len(products)):
                        val = a_products[I][NUM]*products[NUM][9]
                        sumval = sumval + val
                    anum.append(sumval)
                else:
                    sumval = 0
                    for NUM in range(len(products)):
                        val = a_products[I][NUM]*products[NUM][9]*products[NUM][5]/R/T
                        sumval = sumval + val
                    anum.append(sumval)
            elif (I >= len(b) and I < len_square_matrix-1):
                if J < len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        val = a_products[J][NUM]*products[NUM][9]
                        sumval = sumval + val
                    anum.append(sumval)
                elif J == len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        val = products[NUM][9]
                        sumval = sumval + val
                    value = sumval - Ng
                    anum.append(value)
                else:
                    sumval = 0
                    for NUM in range(len(products)):
                        val = products[NUM][9]*products[NUM][5]/R/T
                        sumval = sumval + val
                    anum.append(sumval)
            else:
                if J < len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        sumval = sumval + a_products[J][NUM]*products[NUM][9]*products[NUM][8]/R
                    anum.append(sumval)
                elif J == len(b):
                    sumval = 0
                    for NUM in range(len(products)):
                        sumval = sumval + products[NUM][9]*products[NUM][8]/R
                    anum.append(sumval)
                else:
                    sumval1 = 0
                    sumval2 = 0
                    for NUM in range(len(products)):
                        sumval1 = sumval1 + products[NUM][9]*products[NUM][4]/R
                        sumval2 = sumval2 + products[NUM][9]*products[NUM][5]/R/T*products[NUM][8]/R
                    value = sumval1 + sumval2
                    anum.append(value)
        A.append(anum)
    
    #step 2: set up the solution to this matrix (i.e. Ax = 'B')
    B = []
    for I in range(len_square_matrix):
        if I < len(b):
            bval = 0
            sum1 = 0
            for NUM in range(len(products)):
                sum1 = sum1 + a_products[I][NUM]*products[NUM][9]*products[NUM][7]/R/T
                bval = bval + a_products[I][NUM]*products[NUM][9]
            Bvalue = b[I] - bval + sum1
            B.append(Bvalue)
        elif I == len(b):
            njval = 0
            sum1 = 0
            for NUM in range(len(products)):
                sum1 = sum1 + products[NUM][9]*products[NUM][7]/R/T
                njval = njval + products[NUM][9]
            Bvalue = Ng - njval + sum1
            B.append(Bvalue)
        else:
            s = 0
            sum2 = 0
            sumn = 0
            for NUM in range(len(products)):
                s = s + products[NUM][9]*products[NUM][8]
                sum2 = sum2 + products[NUM][9]*products[NUM][8]/R*products[NUM][7]/R/T
                sumn = sumn + products[NUM][9]
            val = So/R - s/R + Ng - sumn + sum2
            B.append(val)

    solution_vector = numpy.linalg.solve(A,B)    
    pi_vals = []
    for num in range(len(b)):
        pi_vals.append(solution_vector[num])
    del_ln_Ng = solution_vector[len(b)]
    del_ln_T = solution_vector[len(solution_vector)-1]
    #step 3: create new matrix using equations in CEA
    C = []
    for I in range(len(products)):
        cnum = []
        for J in range(len(products)):
            if I == J:
                cnum.append(1)
            else:
                cnum.append(0)
        C.append(cnum)
    
    #step 4: create new solution to this matrix
    D = []
    for I in range(len(products)):
        dval = 0
        for J in range(len(b)):
            dval = dval + a_products[J][I]*pi_vals[J]
        dval = -products[I][7]/R/T + dval + del_ln_Ng + products[I][5]/R/T*del_ln_T
        D.append(dval)
    
    del_ln_n_j = numpy.linalg.solve(C,D)

    C1 = []
    C2 = []
    check = 0
    for NUM in range(len(del_ln_n_j)):
        C1.append((del_ln_n_j[NUM]))
        c2val = abs((-math.log(products[NUM][9]/Ng) - 9.2103404)/(del_ln_n_j[NUM] - del_ln_Ng))
        C2.append(c2val)
        checkval = math.log(products[NUM][9]/Ng)
        if checkval <= -18.420681:
            check = 0
    maxC1 = max(C1)

    L1 = 2/max(5*abs(del_ln_T),5*abs(del_ln_Ng),maxC1)
    L2 = min(C2)
    if check == 1:
        lamda = min(1,L1,L2)
    else:
        lamda = min(1,L1)
    
    # An added step before the products are finished; This may have to go in Evaluate Derivatives
    # Should search through the products, remove any that have a n_j less than 10^-18
    newProducts = []
    newa_products = []
    for num in range(0,len(a_products)):
        a = []
        for n in range(0,len(products)):
            if products[n][9] > 0:
                a.append(a_products[num][n])
        newa_products.append(a)
    for species in products:
        if species[9] > 0:
            newProducts.append(species)
    products = newProducts
    a_products = newa_products
    return products,a_products,pi_vals,del_ln_Ng,del_ln_T,del_ln_n_j,lamda

def Evaluate_Derivatives(products,a_products,Ng,T,P):
    'Evaluate the Properties of the Combustion'
    R = 8.314
    # Constant Pressure
    A = []
    B = []
    C = []
    len_matrix = len(a_products)+1
    for I in range(len_matrix):
        anum = []
        for J in range(len_matrix):
            if I < len(a_products):
                if J < len(a_products):
                    sumval = 0
                    for NUM in range(len(products)):
                        sumval = sumval + a_products[I][NUM]*a_products[J][NUM]*products[NUM][9]
                    anum.append(sumval)
                elif J == len(a_products):
                    sumval = 0
                    for NUM in range(len(products)):
                        sumval = sumval + a_products[I][NUM]*products[NUM][9]
                    anum.append(sumval)
            if I == len(a_products):
                if J < len(a_products):
                    sumval = 0
                    for NUM in range(len(products)):
                        sumval = sumval + a_products[J][NUM]*products[NUM][9]
                    anum.append(sumval)
                else:
                    anum.append(0)
        A.append(anum)
    
    # Matrix B
    for I in range(len_matrix):
        if I < len(a_products):
            sumval = 0
            for NUM in range(len(products)):
                sumval = sumval + a_products[I][NUM]*products[NUM][9]*products[NUM][5]/R/T
            B.append(-sumval)
        else:
            sumval = 0
            for NUM in range(len(products)):
                sumval = sumval + products[NUM][9]*products[NUM][5]/R/T
            B.append(-sumval)

    # Constant Temperature
    
    # Matrix C
    for I in range(len_matrix):
        if I < len(a_products):
            sumval = 0
            for NUM in range(len(products)):
                sumval = sumval + a_products[I][NUM]*products[NUM][9]
            C.append(sumval)
        else:
            sumval = 0
            for NUM in range(len(products)):
                sumval = sumval + products[NUM][9]
            C.append(sumval)

    # Calculating Constant Pressure Variables

    Deriv_ConstPress_Values = numpy.linalg.solve(A,B)
    Deriv_Pi_Temp_constP = []
    for NUM in range(len(a_products)):
        Deriv_Pi_Temp_constP.append(Deriv_ConstPress_Values[NUM])
    Deriv_Ng_Temp_constP = Deriv_ConstPress_Values[len(Deriv_ConstPress_Values)-1]

    # Calculating Constant Temperature Variables
    Deriv_ConstTemp_Values = numpy.linalg.solve(A,C)
    Deriv_Pi_Pres_constT = []
    for NUM in range(len(a_products)):
        Deriv_Pi_Pres_constT.append(Deriv_ConstTemp_Values[NUM])
    Deriv_Ng_Pres_constT = Deriv_ConstTemp_Values[len(Deriv_ConstTemp_Values)-1]

    # Calculating C_p/R

    len_CP = 4
    ans_part = 0
    for N in range(len_CP):
        if N == 0:
            sum2 = 0
            for I in range(len(a_products)):
                sum1 = 0
                for J in range(len(products)):
                    sum1 = sum1 + a_products[I][J]*products[J][9]*products[J][5]/R/T
                sum2 = sum2 + sum1*Deriv_Pi_Temp_constP[I]
            ans_part = ans_part + sum2
        if N == 1:
            sum1 = 0
            for I in range(len(products)):
                sum1 = sum1 + products[I][9]*products[I][5]/R/T
            ans_part = ans_part + sum1*Deriv_Ng_Temp_constP
        if N == 2:
            sum1 = 0
            for I in range(len(products)):
                sum1 = sum1 + products[I][9]*products[I][4]/R
            ans_part = ans_part + sum1
        if N == 3:
            sum1 = 0
            for I in range(len(products)):
                sum1 = sum1 + products[I][9]*products[I][5]*products[I][5]/R/R/T/T
            ans_part = ans_part + sum1
    
    C_p_comb = ans_part*R

    deriv_lnV_lnT_constP = 1 + Deriv_Ng_Temp_constP

    deriv_lnV_lnP_constT = -1 + Deriv_Ng_Pres_constT

    C_v_comb = C_p_comb + Ng*R*(deriv_lnV_lnT_constP*deriv_lnV_lnT_constP)/deriv_lnV_lnP_constT

    gamma = C_p_comb/C_v_comb

    gamma_s = -gamma/deriv_lnV_lnP_constT

    a_sound = math.sqrt(Ng*R*T*gamma_s*1000)

    deriv_lnP_lnT_s = C_p_comb/Ng/R/deriv_lnV_lnT_constP
    H_0 = 0
    S_0 = 0
    G_0 = 0
    for NUM in range(len(products)):
        S_0 = S_0 + products[NUM][9]*(products[NUM][10] - R*math.log(products[NUM][9]/Ng) - R*math.log(P))
        H_0 = H_0 + products[NUM][5]*products[NUM][9]
        G_0 = G_0 + products[NUM][7]*products[NUM][9]

    return H_0,S_0,G_0,C_p_comb,C_v_comb,gamma_s,a_sound