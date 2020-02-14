import PySimpleGUI as sg
from guiClass import Gui, SaveandPlot, WriteXY              
import math
import EquilibriumFunctions
import RocketEngineCalculations
import json_load_species as JS

JS.parseThermoJSON() # Load the Thermo_JSON.json file
g = Gui()   # Load the GUI Class outlined in guiClass.py
InitialReactants = []   # List of the initial reactants for the chemical equation
while True: #Continuosly run the script
    try:
        event,values = g.window.Read()  #read the GUI - the following lines depict the events
        if event is None or event == 'Quit':
            break
        # Inputting Values
        Po = float(values['Pres'])*1e5 #Pascals
        OF = float(values['OFRatio'])  #O/F Ratio
        
        #Set the solver type
        if values['IACsolver'] == True:     
            solverType = 'IAC'
        elif values['FACsolver'] == True:
            solverType = 'FAC'
        
        #Entrance Conditions
        if values['noneEntrance'] == True:      
            Entry = []
            if solverType == 'FAC':
                AcAt = 6.0
        elif values['AcAt'] == True:
            if solverType == 'IAC':
                Entry = ['Area',float(values['ENTRANCECONDITIONS'])]
            elif solverType == 'FAC':
                AcAt = float(values['ENTRANCECONDITIONS'])
        elif values['MdotAc'] == True:
            mdotAc = float(values['ENTRANCECONDITIONS'])
        
        #Adding propellants/Clearing propellants
        if event == 'Add':
            name = values['PropName']
            Temp = float(values['PropTemp'])
            wtpct = float(values['WeightPerc'])
            if values['Fuel'] == True:
                typeProp = 'fuel'
            elif values['Oxid'] == True:
                typeProp = 'oxid'
            val = (name,Temp,wtpct,typeProp)
            InitialReactants.append(val)
            g.window['AddReactant'].update(InitialReactants)
        elif event == 'Clear':
            Temp = []
            for L in range(0,len(InitialReactants)-1):
                Temp.append(InitialReactants[L])
            InitialReactants = Temp 
            g.window['AddReactant'].update(InitialReactants)
        elif event == 'ClearAll':
            InitialReactants = []
            g.window['AddReactant'].update(InitialReactants)
        
        #Exit Conditions
        if values['ExitNone'] == True:
            ExitConditions = []
        elif values['Pressure'] == True:
            ExitConditions = ['Pressure',float(values['ExitNum'])]
        elif values['AreaRatio'] == True:
            ExitConditions = ['Area',float(values['ExitNum'])]
        elif values['Separation'] == True:
            ExitConditions = ['Separation',float(values['ExitNum'])]
        
        # Run the Calculations
        if event == 'RunCalculation':
            products, a_products, b, Ho = EquilibriumFunctions.INITIAL_CONDITIONS(Po,InitialReactants,OF)
            print("Running Solver...")
            if solverType == 'FAC':
                if values['MdotAc'] == False:
                    ENGINE = RocketEngineCalculations.FACsolver(Po,products,a_products,b,Ho,ExitConditions=ExitConditions,AcAt=AcAt)
                else:
                    ENGINE = RocketEngineCalculations.FACconditionAcMdot(Po,products,a_products,b,Ho,mdotAc,ExitConditions=ExitConditions)
            elif solverType == 'IAC':
                ENGINE = RocketEngineCalculations.IACsolver(Po,products,a_products,b,Ho,ExitConditions=ExitConditions,EntryConditions=Entry)
            CstarVal,Cf_SL,Cf_opt,Cf_vac,Pinf,performance = RocketEngineCalculations.EnginePerformance(ENGINE,OF)
            
            # This will need to have a separate tab to create the Thrust Values, mdot values, etc.
            Thrust = CstarVal*Cf_SL[len(Cf_SL)-1]*1.96
            Topt = CstarVal*Cf_opt[len(Cf_opt)-1]*1.96
            mdot = Thrust/CstarVal/Cf_SL[len(Cf_SL)-1]
            Astar = Thrust/Pinf/Cf_SL[len(Cf_SL)-1]
            Rstar = math.sqrt(Astar/math.pi)*100
            
            # Output the Calculations to the 'RUN_calculation' Tab 
            inputstring = ''
            for NUM in range(len(performance[0])):
                inputstring = inputstring + performance[0][NUM] + performance[1][NUM] + '\n'
            g.window['RUN_CALCULATION'].update(inputstring)
            
        # Optimize OF Ratio - SOMETIMES CRASHES
        elif event == 'Optimize':
            OF = RocketEngineCalculations.OptimizeOF(Po,InitialReactants,OF,5)
            g.window['tabgroup'].Widget.select(0)
            printval = str(InitialReactants) + '\n' + 'Optimum OF Ratio is: {:.2f}'.format(OF)
            g.window['AddReactant'].update(printval)
            g.window['OFRatio'].update(str(OF))

        # Clear the screen on the Run Calculation tab
        elif event == 'ClrScrn':
            g.window['RUN_CALCULATION'].update('')

        # Engine Sizing Tab - generates a graph of the engine, and allows you to save it as a text
        elif event == 'RESIZE':
            Theta_C = float(values['ThetaC'])
            Theta_N = float(values['ThetaN'])
            Theta_E = float(values['ThetaE'])
            if values['MdotAc'] == False and values['noneEntrance'] == False:
                AcAt = float(values['ENTRANCECONDITIONS'])
            else: AcAt = 6.0
            X,Y,Lcombustion,Lcone,Xparabola,Xcombstraight,Nlist,Q,E = RocketEngineCalculations.EngineSize(Rstar/100,ENGINE[len(ENGINE)-2][11],Theta_C,Theta_N,Theta_E,AcAt=AcAt)
            SaveandPlot(X,Y)
            g.window['ENGINECONTOUR'].update('EngineContour.png')
            strval =str(InitialReactants) + '\n' f'L_Chamber: {Lcombustion}m\nL_Engine: {Xparabola-Xcombstraight} m'
            g.window['AddReactant'].update(strval)
        elif event == 'SAVE':
            WriteXY(X,Y)

    # Error Checking
    except Exception as e:
        e = str(e)
        InitialReactants = []
        if event == 'RunCalculation' or event == 'Optimize':
            g.window['tabgroup'].Widget.select(0)
        if e == "local variable 'speciesname' referenced before assignment":
           string = 'Choose Propellant from Thermo.txt'
        elif "could not convert string to float" in e:
            string = 'Check Pressure or O/F Ratio is not empty' 
        else:
            string = e
        g.window['AddReactant'].update(string)