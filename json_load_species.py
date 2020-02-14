# Used to load the Thermo_JSON file to a 2D List (similar to ALLELEMENTS.py)
import json #load the json module
from ALLELEMENTS import ALL_SPECIES_PROPERTIES

def parseThermoJSON(ALL_SPECIES_PROPERTIES=ALL_SPECIES_PROPERTIES):
    'parse the Thermo_JSON file to a series of lists and update ALLELEMENTS - returns zero'
    SPECIES_NAMES = ALL_SPECIES_PROPERTIES[0]
    SPECIES_ELEMENTS = ALL_SPECIES_PROPERTIES[1]
    SPECIES_ELENUM = ALL_SPECIES_PROPERTIES[2]
    SPECIES_MW = ALL_SPECIES_PROPERTIES[3]
    SPECIES_Hform = ALL_SPECIES_PROPERTIES[4]
    SPECIES_T = ALL_SPECIES_PROPERTIES[5]
    SPECIES_zeroISgas = ALL_SPECIES_PROPERTIES[6]

    with open('Thermo_JSON.json') as json_file:
        data = json.load(json_file)
        for value in data['species']:
            SPECIES_NAMES.append(value['name'])
            SPECIES_ELEMENTS.append(value['elements'])
            SPECIES_ELENUM.append(value['element_num'])
            SPECIES_MW.append(value['mw'])
            SPECIES_Hform.append(value['hform'])
            SPECIES_T.append(value['T'])
            SPECIES_zeroISgas.append(value['zeroISgas'])
        
    return 0

def addThermoJSON(name,elements,element_num,mw,Hform,T=[[298.15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]], zeroISgas=1,ALL_SPECIES_PROPERTIES=ALL_SPECIES_PROPERTIES):
    'Add a new Species to the Thermo_JSON file & appends the values to the current ALL_SPECIES_Properties list such that a restart is not requried'
    SPECIES_NAMES = ALL_SPECIES_PROPERTIES[0]
    SPECIES_ELEMENTS = ALL_SPECIES_PROPERTIES[1]
    SPECIES_ELENUM = ALL_SPECIES_PROPERTIES[2]
    SPECIES_MW = ALL_SPECIES_PROPERTIES[3]
    SPECIES_Hform = ALL_SPECIES_PROPERTIES[4]
    SPECIES_T = ALL_SPECIES_PROPERTIES[5]
    SPECIES_zeroISgas = ALL_SPECIES_PROPERTIES[6]

    SPECIES_NAMES.append(name)
    SPECIES_ELEMENTS.append(elements)
    SPECIES_ELENUM.append(element_num)
    SPECIES_MW.append(mw)
    SPECIES_Hform.append(Hform)
    SPECIES_T.append(T)
    SPECIES_zeroISgas.append(zeroISgas)

    data = {}
    data['species'] = []
    for NUM in range(len(ALL_SPECIES_PROPERTIES[0])):
        name = ALL_SPECIES_PROPERTIES[0][NUM]
        elements = ALL_SPECIES_PROPERTIES[1][NUM]
        element_num = ALL_SPECIES_PROPERTIES[2][NUM]
        mw = ALL_SPECIES_PROPERTIES[3][NUM]
        Hform = ALL_SPECIES_PROPERTIES[4][NUM]
        T = ALL_SPECIES_PROPERTIES[5][NUM]
        zeroISgas = ALL_SPECIES_PROPERTIES[6][NUM]
        data['species'].append({
            'name':name,
            'elements': elements,
            'element_num': element_num,
            'mw': mw,
            'hform': Hform,
            'T': T,
            'zeroISgas': zeroISgas
        })
    with open('Thermo_JSON.json','w') as outfile:
        json.dump(data,outfile,indent=4)
        
    return 0

def removeThermoJSON(nameRemove,ALL_SPECIES_PROPERTIES=ALL_SPECIES_PROPERTIES):
    'Remove the specified species from both ALL_SPECIES_PROPERTIES and the Thermo_JSON.json file'
    if ALL_SPECIES_PROPERTIES[0] == []:
        return 1
    data = {}
    data['species'] = []
    for NUM in range(len(ALL_SPECIES_PROPERTIES[0])):
        if ALL_SPECIES_PROPERTIES[0][NUM] != nameRemove:
            name = ALL_SPECIES_PROPERTIES[0][NUM]
            elements = ALL_SPECIES_PROPERTIES[1][NUM]
            element_num = ALL_SPECIES_PROPERTIES[2][NUM]
            mw = ALL_SPECIES_PROPERTIES[3][NUM]
            Hform = ALL_SPECIES_PROPERTIES[4][NUM]
            T = ALL_SPECIES_PROPERTIES[5][NUM]
            zeroISgas = ALL_SPECIES_PROPERTIES[6][NUM]
            data['species'].append({
                'name':name,
                'elements': elements,
                'element_num': element_num,
                'mw': mw,
                'hform': Hform,
                'T': T,
                'zeroISgas': zeroISgas
            })
    with open('Thermo_JSON.json','w') as outfile:
        json.dump(data,outfile,indent=4)

    return 0