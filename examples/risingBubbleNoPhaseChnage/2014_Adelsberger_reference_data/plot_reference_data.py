import numpy as np

def plotAdelsberger(axs):
    folder = './2014_Adelsberger_reference_data/'
    ######################### load Drops data ##################################
    t_centerOfMass_DROPS = list()
    dat_centerOfMass_DROPS = list()
    t_riseVelocity_DROPS = list()
    dat_riseVelocity_DROPS = list()
    t_sphericity_DROPS = list()
    dat_sphericity_DROPS = list()
    #read center of mass:
    input = open(folder+'TC1_centerOfMass_DROPS.txt', 'r')
    lines = input.readlines()
    for line in lines[:]:
        numbers = line.split()
        t_centerOfMass_DROPS.append(float(numbers[0]))
        dat_centerOfMass_DROPS.append(float(numbers[1]))
    t_centerOfMass_DROPS   = np.array(t_centerOfMass_DROPS)/np.sqrt(0.5/abs(-9.8*0.1))
    dat_centerOfMass_DROPS = np.array(dat_centerOfMass_DROPS)
    #read rise velocity:
    input = open(folder+'TC1_riseVelocity_DROPS.txt', 'r')
    lines = input.readlines()
    for line in lines[:]:
        numbers = line.split()
        t_riseVelocity_DROPS.append(float(numbers[0]))
        dat_riseVelocity_DROPS.append(float(numbers[1]))
    t_riseVelocity_DROPS   = np.array(t_riseVelocity_DROPS)/np.sqrt(0.5/abs(-9.8*0.1))
    dat_riseVelocity_DROPS = np.array(dat_riseVelocity_DROPS)/np.sqrt(0.5*abs(-0.98))
    #read sphereicity:
    input = open(folder+'TC1_sphericity_DROPS.txt', 'r')
    lines = input.readlines()
    for line in lines[:]:
        numbers = line.split()
        t_sphericity_DROPS.append(float(numbers[0]))
        dat_sphericity_DROPS.append(float(numbers[1]))
    t_sphericity_DROPS   = np.array(t_sphericity_DROPS)/np.sqrt(0.5/abs(-9.8*0.1))
    dat_sphericity_DROPS = np.array(dat_sphericity_DROPS)
    ######################### load NaSt3D data ##################################
    t_centerOfMass_NaSt3D = list()
    dat_centerOfMass_NaSt3D = list()
    t_riseVelocity_NaSt3D = list()
    dat_riseVelocity_NaSt3D = list()
    t_sphericity_NaSt3D = list()
    dat_sphericity_NaSt3D = list()
    #read center of mass:
    input = open(folder+'TC1_centerOfMass_NaSt3D.txt', 'r')
    lines = input.readlines()
    for line in lines[:]:
        numbers = line.split()
        t_centerOfMass_NaSt3D.append(float(numbers[0]))
        dat_centerOfMass_NaSt3D.append(float(numbers[1]))
    t_centerOfMass_NaSt3D   = np.array(t_centerOfMass_NaSt3D)/np.sqrt(0.5/abs(-9.8*0.1))
    dat_centerOfMass_NaSt3D = np.array(dat_centerOfMass_NaSt3D)
    #read rise velocity:
    input = open(folder+'TC1_riseVelocity_NaSt3D.txt', 'r')
    lines = input.readlines()
    for line in lines[:]:
        numbers = line.split()
        t_riseVelocity_NaSt3D.append(float(numbers[0]))
        dat_riseVelocity_NaSt3D.append(float(numbers[1]))
    t_riseVelocity_NaSt3D   = np.array(t_riseVelocity_NaSt3D)/np.sqrt(0.5/abs(-9.8*0.1))
    dat_riseVelocity_NaSt3D = np.array(dat_riseVelocity_NaSt3D)/np.sqrt(0.5*abs(-0.98))
    #read sphereicity:
    input = open(folder+'TC1_sphericity_NaSt3D.txt', 'r')
    lines = input.readlines()
    for line in lines[:]:
        numbers = line.split()
        t_sphericity_NaSt3D.append(float(numbers[0]))
        dat_sphericity_NaSt3D.append(float(numbers[1]))
    t_sphericity_NaSt3D   = np.array(t_sphericity_NaSt3D)/np.sqrt(0.5/abs(-9.8*0.1))
    dat_sphericity_NaSt3D = np.array(dat_sphericity_NaSt3D)
    ######################### load OpenFOAM data ##################################
    t_centerOfMass_OpenFOAM = list()
    dat_centerOfMass_OpenFOAM = list()
    t_riseVelocity_OpenFOAM = list()
    dat_riseVelocity_OpenFOAM = list()
    t_sphericity_OpenFOAM = list()
    dat_sphericity_OpenFOAM = list()
    #read center of mass:
    input = open(folder+'TC1_centerOfMass_OpenFOAM.txt', 'r')
    lines = input.readlines()
    for line in lines[:]:
        numbers = line.split()
        t_centerOfMass_OpenFOAM.append(float(numbers[0]))
        dat_centerOfMass_OpenFOAM.append(float(numbers[1]))
    t_centerOfMass_OpenFOAM   = np.array(t_centerOfMass_OpenFOAM)/np.sqrt(0.5/abs(-9.8*0.1))
    dat_centerOfMass_OpenFOAM = np.array(dat_centerOfMass_OpenFOAM)
    #read rise velocity:
    input = open(folder+'TC1_riseVelocity_OpenFOAM.txt', 'r')
    lines = input.readlines()
    for line in lines[:]:
        numbers = line.split()
        t_riseVelocity_OpenFOAM.append(float(numbers[0]))
        dat_riseVelocity_OpenFOAM.append(float(numbers[1]))
    t_riseVelocity_OpenFOAM   = np.array(t_riseVelocity_OpenFOAM)/np.sqrt(0.5/abs(-9.8*0.1))
    dat_riseVelocity_OpenFOAM = np.array(dat_riseVelocity_OpenFOAM)/np.sqrt(0.5*abs(-0.98))
    #read sphereicity:
    input = open(folder+'TC1_sphericity_OpenFOAM.txt', 'r')
    lines = input.readlines()
    for line in lines[:]:
        numbers = line.split()
        t_sphericity_OpenFOAM.append(float(numbers[0]))
        dat_sphericity_OpenFOAM.append(float(numbers[1]))
    t_sphericity_OpenFOAM   = np.array(t_sphericity_OpenFOAM)/np.sqrt(0.5/abs(-9.8*0.1))
    dat_sphericity_OpenFOAM = np.array(dat_sphericity_OpenFOAM)
    
    axs[0].plot(t_centerOfMass_DROPS,dat_centerOfMass_DROPS,c='0',linestyle='solid')
    axs[0].plot(t_centerOfMass_NaSt3D,dat_centerOfMass_NaSt3D,c='0',linestyle='dashed')
    axs[0].plot(t_centerOfMass_OpenFOAM,dat_centerOfMass_OpenFOAM,c='0',linestyle='dotted')
    axs[2].plot(t_riseVelocity_DROPS,dat_riseVelocity_DROPS,c='0',linestyle='solid',label='DROPS')
    axs[2].plot(t_riseVelocity_NaSt3D,dat_riseVelocity_NaSt3D,c='0',linestyle='dashed',label='NaSt3D')
    axs[2].plot(t_riseVelocity_OpenFOAM,dat_riseVelocity_OpenFOAM,c='0',linestyle='dotted',label='OpenFOAM')
    axs[1].plot(t_sphericity_DROPS,dat_sphericity_DROPS,c='0',linestyle='solid')
    axs[1].plot(t_sphericity_NaSt3D,dat_sphericity_NaSt3D,c='0',linestyle='dashed')
    axs[1].plot(t_sphericity_OpenFOAM,dat_sphericity_OpenFOAM,c='0',linestyle='dotted')
    return axs
