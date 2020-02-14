def AreaRatio(gam,M):
    from math import sqrt
    'Return the Area ratio for a mach number and gamma'
    gp1 = gam + 1
    gm1 = gam - 1
    epsilon = 1/M * sqrt(((1 + gm1/2*M*M)/(gp1/2))**(gp1/gm1))
    return epsilon

def FAC_AreaRatio(gam,M):
    gp1 = gam + 1
    gm1 = gam - 1
    eps = (gp1/2)**(1/gm1) * M * (1 + gm1/gp1*M*M)**(1/gm1)
    PIfmach = (1 + gm1/gp1*M*M)**(1/gm1)
    Fc = 1/eps
    return Fc,PIfmach
def TempRatio(gam,M):
    TR = 1 + (gam-1)/2*M*M
    return TR

def PresRatio(gam,M):
    'return pressure ratio'
    TR = TempRatio(gam,M)
    PR = (TR)**(gam/(gam-1))
    return PR

def get_mach_num(exp,gam,type):
    'Return the Mach number for an Area Ratio'
    diff = 1
    if type == 'sup':
        X1 = 5
        X0 = 1
    elif type == 'sub':
        X1 = 0.1
        X0 = 0.01
    while diff > 0.0001:
        Y = exp
        Y0 = AreaRatio(gam,X0)
        Y1 = AreaRatio(gam,X1)
        m = (Y1 - Y0)/(X1 - X0)
        X2 = (Y - Y0)/m + X0
        X0 = X1
        X1 = X2
        diff = abs(Y - Y0)
    MACH = X2
    return MACH

def get_mach_num_FAC(exp,gam,type):
    'Return the Mach number for an Area Ratio'
    diff = 1
    if type == 'sup':
        X1 = 5
        X0 = 1
    elif type == 'sub':
        X1 = 0.1
        X0 = 0.01
    while diff > 0.0001:
        Y = exp
        Y0,Pifmach = FAC_AreaRatio(gam,X0)
        Y1,Pifmach = FAC_AreaRatio(gam,X1)
        m = (Y1 - Y0)/(X1 - X0)
        X2 = (Y - Y0)/m + X0
        X0 = X1
        X1 = X2
        diff = abs(Y - Y0)
    MACH = X2
    return MACH