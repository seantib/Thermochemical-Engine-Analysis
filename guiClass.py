import PySimpleGUI as sg
sg.change_look_and_feel('Dark')
class Gui:
    def __init__(self):
        self.Inital_TAB = [[sg.Text('Chamber Pressure (BAR):'),sg.Input(focus=True,key='Pres',size=(10,4))],
                        [sg.Text('Propellant          Temperature         Wt% (0-1)')],
                        [sg.Input(key='PropName',size=(12,4)),
                        sg.Input('298.15',key='PropTemp',size=(12,4)),
                        sg.Input('1',key='WeightPerc',size=(12,4)),
                        sg.Radio('Fuel',group_id='type',key='Fuel',default=True),
                        sg.Radio('Oxidizer',group_id='type',key='Oxid')],
                        [sg.Text("O/F Ratio:"), sg.Input('2.0',key='OFRatio',size=(4,4))],
                        [sg.Text("Exit Condition:"),sg.Input(size=(6,4),key='ExitNum'),sg.Radio('None',key='ExitNone',default=True,group_id='ExitCond'),
                        sg.Radio('Pressure',key='Pressure',group_id='ExitCond'), sg.Radio('Area Ratio',key='AreaRatio',group_id='ExitCond'),
                        sg.Radio('Separation',key='Separation',group_id='ExitCond')],
                        [sg.Text('Solver Type'), sg.Radio('IAC',group_id='SolverType',key='IACsolver'), sg.Radio('FAC',group_id='SolverType',key='FACsolver',default=True)],
                        [sg.Text('Entrance Conditions'), sg.Input(size=(6,4), key='ENTRANCECONDITIONS'),sg.Radio('None',group_id='inputType',key='noneEntrance',default=True), sg.Radio('AcAt',group_id='inputType',key='AcAt'),sg.Radio('Ac/mdot',group_id='inputType',key='MdotAc')],
                        [sg.Button('Add',bind_return_key=True), sg.Button('Clear',key='Clear'), sg.Button('Clear All',key='ClearAll')],
                        [sg.Output(key='AddReactant',size=(100,4))]]
        self.PERFORMANCE_TAB = [[sg.Button('Run Calculation',key='RunCalculation'),sg.Button('Optimize O/F',key='Optimize'),sg.Button('Clear Screen',key='ClrScrn')],
                                [sg.Output(key='RUN_CALCULATION',size=(100,35))]]
                                
        self.SIZING_TAB = [[sg.Button('ReSize',key='RESIZE'), sg.Text('Theta_C:'),sg.Input('30',size=(6,4),key='ThetaC'),
                            sg.Text('Theta_N:'),sg.Input('25',size=(6,4),key='ThetaN'),
                            sg.Text('Theta_E:'),sg.Input('10',size=(6,4),key="ThetaE"),sg.Button('Save Contour',key='SAVE')],[sg.Image('EngineContour.png',key="ENGINECONTOUR")]]
        self.SIZING_OUTPUT = [[sg.Text('Engine Length:'),sg.Text(key='EngineLength',size=(10,4))],
                              [sg.Text('Engine Length:'),sg.Text(key='EngineSize',size=(10,4))]]
        self.layout = [[sg.TabGroup([[sg.Tab('Initial Conditions',self.Inital_TAB),
                        sg.Tab('Performance',self.PERFORMANCE_TAB),
                        sg.Tab('Graphing',self.SIZING_TAB),
                        sg.Tab('Data',self.SIZING_OUTPUT)]],key='tabgroup')]]
        self.window = sg.Window('Thermodynamic Engine Analysis',resizable=True,default_element_size=(150,40)).layout(self.layout)

def SaveandPlot(X,Y):
    import matplotlib.pyplot as plt
    'Save and Plot the Engine Contour'
    plt.cla()
    plt.figure()
    plt.ylim((0,0.2))
    plt.plot(X,Y)
    plt.title("Engine Contour")
    plt.xlabel("Length (m)")
    plt.ylabel("Radius (m)")
    plt.grid()
    plt.savefig('EngineContour.png')

def WriteXY(X,Y):
    DataPoints = open('Data.txt','w')
    for N in range(0,len(X)):
        string = f'{X[N]}\t{Y[N]}\n'
        DataPoints.write(string)
    DataPoints.close()
    return





