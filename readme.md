Attached is a compressed file with a set of python files. To run this code from the set (i.e. without downloading any packages). If you run into an issue executing the code, you will most likely need to install the numpy package (uses it for linear algebra calculations).

Quick setup:
1. Open main_NoGUI.py (this holds the backbone of the gui version) - in it, right at the top, you can see all the imported files (very similar to matlab)
2. You can see the function SaveandPlot() - This requires matplotlib to run. If you don't want to use this, you may delete lines 6, (10-18), and 63. (you may also need to delete line 59
3. If you'd like to run the gui verison, run "python main.py" in the terminal - REQUIRES matplotlib and PySimpleGUI to run!


Importing Fuels and Oxidizers:
These are stored as tuples (special kind of list/array). The order goes (name,temperature,wt%,type), where the 'name' must be a string object (hence the quotations). For a list of available fuels/oxidizers, check out the Thermo_json.json file attached in the folder.

1. you'll need a fuel name, and the temperature at which it is being injected (or before combustion) - for some reactants (usually the liquids), this is only defined at one temperature, and will throw an error if it is not set at that temperature.
2. The weight percentage (wt%) is the percentage of the total fuel or oxidizer (separate of course) that this particular reactant takes up. This functionality allows for multiple fuels and oxidizers to be entered, as in Aerozene-50!
3. Lastly, the type will need to be specified. This is ALSO a string!

Directly underneath, there is an array containing [fuel,oxidizer]. If you wish to use more than one fuel, simply use [fuel,fuel2,fuel3,fuel...,oxidizer1,oxidizer2,....]

Setting the O/F Ratio:
There is a variable called OF, and a commented out function which allows the optimization between two O/F ratios! (This function does work!)

Exit Conditions:
THIS MUST BE IN THE BRACKETS
It allows for multiple types of exit conditions to be calculated simultaneously
Currently, you have three options 'Pressure', 'Area', or 'Separation'
Pressure - Enter in a numerical pressure value in PACALS
Area - Enter the desired area ratio
Separation - Enter the pressure at which you'd like separation conditions to be calculated

Currently it is calling the FAC (Finite Area Combustor Solver). It has a default value of Ac/At = 6.0. If you'd like to use a separate value, simply replace the line with:

ENGINE = RocketEngineCalculations.FACsolver(Po,products,a_products,b,Ho,ExitConditions=Exit,AcAt=##.##)

Where your ###.## are any value you'd like.

If you want to call the Infinite Area Model, simply use:
ENGINE = RocketEngineCalculations.IACsolver(Po,products,a_products,b,Ho,ExitConditions=Exit)


Underneath is the line:
CstarVal,Cf_SL,Cf_opt,Cf_vac,Pinf,performance = RocketEngineCalculations.EnginePerformance(ENGINE,OF)

This returns all these values from the engine and O/F ratio chosen.


This set of lines prints all those nice values you get with the CEA output to the terminal
for NUM in range(len(performance[0])):
    inputstring = inputstring + performance[0][NUM] + performance[1][NUM] + '\n'
print(inputstring)

Lastly, these lines calculate Sea-Level Thrust, Optimal Thrust, mdot, Astar, and Rstar.

To update these values:
1. Comment out the Right side of the Thrust Equation, replace with your desired Sea-Level Thrust.
2. move Topt below mdot, and replace 1.96 in the Topt equation with mdot.

Thrust = CstarVal*Cf_SL[len(Cf_SL)-1]*1.96
Topt = CstarVal*Cf_opt[len(Cf_opt)-1]*1.96
mdot = Thrust/CstarVal/Cf_SL[len(Cf_SL)-1]
Astar = Thrust/Pinf/Cf_SL[len(Cf_SL)-1]
Rstar = math.sqrt(Astar/math.pi)*100


The Last function is a Rao Nozzle Creation Function. This actually returns a text file called Data.txt which contains the X,Y points for the entire combustion chamber and nozzle curve.
