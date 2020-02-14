#AUTHOR: SEAN TIBBETTS
import Thermo                
import math
import EquilibriumFunctions
import RocketEngineCalculations
import matplotlib.pyplot as plt
import json_load_species as JS

############################################################## FUNCTIONS #########################################################################
def SaveandPlot(X,Y):
    'Save and Plot the Engine Contour'
    plt.figure(1,figsize=(10,2))
    plt.plot(X,Y)
    plt.title("Engine Contour")
    plt.xlabel("Length (m)")
    plt.ylabel("Radius (m)")
    plt.grid()
    plt.savefig('EngineContour.png')



########################################################## COMBUSTION PARAMETERS ###################################################################
#                                       Units                 #Inputs
Po = 28e5                                #Pa             # Chamber Pressure
fuel = ('CH4',298.15,1,'fuel')
oxidizer = ('O2(L)',90.17,1,'oxid')
Exit = ['Area',8]

errval = JS.parseThermoJSON()
######################################################################################################################### Create Initial Conditions
InitialReactants = [fuel,oxidizer]
OF = 2.7 #RocketEngineCalculations.OptimizeOF(Po,InitialReactants,2.5,3)    #Optional Solver for Maximizing O/F ratio
products, a_products, b, Ho = EquilibriumFunctions.INITIAL_CONDITIONS(Po,InitialReactants,OF)


############################################################## SOLVER #######################################################################################
print("Running Solver...")
ENGINE = RocketEngineCalculations.FACsolver(Po,products,a_products,b,Ho,ExitConditions=Exit,AcAt=6.606)
#ENGINE = RocketEngineCalculations.FACconditionAcMdot(Po,products,a_products,b,Ho,0.00444625,ExitConditions=['Separation',1e5,'Area',10.9])

CstarVal,Cf_SL,Cf_opt,Cf_vac,Pinf,performance = RocketEngineCalculations.EnginePerformance(ENGINE,OF)
inputstring = ''
for NUM in range(len(performance[0])):
    inputstring = inputstring + performance[0][NUM] + performance[1][NUM] + '\n'
print(inputstring)
Thrust = CstarVal*Cf_SL[len(Cf_SL)-1]*1.96
Topt = CstarVal*Cf_opt[len(Cf_opt)-1]*1.96
mdot = Thrust/CstarVal/Cf_SL[len(Cf_SL)-1]
Astar = Thrust/Pinf/Cf_SL[len(Cf_SL)-1]
Rstar = math.sqrt(Astar/math.pi)*100

print("\n\nDesired Properties:")
print(f"Thrust(SL): {Thrust/1000} kN")
print(f"Thrust(opt): {Topt/1000} kN")
print(f"Mass Flow Rate: {mdot} kg/s")
print(f"Throat Radius: {Rstar} cm")


X,Y,Lcombustion,Lcone,Xparabola,Xcombstraight,Nlist,Q,E = RocketEngineCalculations.EngineSize(Rstar/100,10.9,30,28,11)



SaveandPlot(X,Y)


'''
DataPoints = open('Data.txt','w')
for N in range(0,len(X)):
    string = f'{X[N]}\t{Y[N]}\n'
    DataPoints.write(string)
DataPoints.close()'''
#____________________________________________________________TO DO LIST __________________________________________________________________
# Need to overhaul my regenerative cooling code with updated values - specifically the pressure drops
# Add options for IAC, FAC (AcAt or AcMdot) to GUI code
# Output tables for GUI
# SQLite for the ALLELEMENTS data - Also, need to re-parse the data, some species (C10 species) are parsing incorrectly
