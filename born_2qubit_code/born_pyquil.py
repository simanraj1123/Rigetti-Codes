import numpy as np, matplotlib.pyplot as plt, random, time, os, datetime
from functools import lru_cache
from pyquil import Program, get_qc
from pyquil.gates import *
import os
from pyquil.quilatom import quil_sin, quil_cos, Parameter
from pyquil.quilbase import DefGate
from pyquil.latex import display, to_latex
# import Peres_helpers as hf
import pickle
from collections import Counter
from datetime import date

N_SHOTS = 10_000
Z_SCORE = 3

e=0

# Generate random parameters for random states.
def params_complex():
	'''
	Generates parameters to prepare COMPLEX quantum states.
    
    Returns:
        A list of three tuples. The first element of each tuple is the value of theta 
        and the second element is the value of phi.
    
	'''
	theta = np.arccos(np.cos(e) - 2 * np.array([random.uniform(0,1) for _ in range(2)]))
	phi = np.array([2*np.pi*random.uniform(0,1) for _ in range(2)])
	params = zip(theta, phi)
	return list(params)

# The construction of the parametric circuit.
def circuit_123(a,b):
    circ = Program()
    
    theta = circ.declare('theta', 'REAL', 2)
    phi = circ.declare('phi', 'REAL', 2)
    t = circ.declare('t', 'REAL', 2)
    c = circ.declare('ro', 'BIT', 2)
    
    circ += RY(theta[0], a)
    circ += RZ(phi[0], a)
    
    circ += RY(theta[1], b)
    circ += RZ(phi[1]/2, b)
    
    circ += CNOT(a, b)
    
    circ += RZ(-phi[1]/2, b)
    circ += RY((t[1]-theta[1])/2, b)
    
    circ += CNOT(a,b)
    
    circ += RY((-theta[1] - t[1])/2, b)
    
    circ += RY(-t[0], a)
    
    circ += X(a)
    
    circ += MEASURE(a, c[0])
    circ += MEASURE(b, c[1])
    
    circ.wrap_in_numshots_loop(N_SHOTS)
    
    return circ


def f(u):
    '''
    Compile and run the circuit given the parameters. The list of outputs is returned.
    '''
    theta, phi = [], []
    outcomes = {}
    for i in range(len(u)):
        theta.append(u[i][0])
        phi.append(u[i][1])
    t = [2 * np.arccos(1/np.sqrt(3)), 2 * np.arccos(1/np.sqrt(2))]
    s123 = qc.run(exe, memory_map={'theta': theta, 'phi': phi, 't': t}) # Stores the output of the circuit run.
    counts_s123 = Counter([''.join(list(map(str, elem))) for elem in s123])
    
    t = [np.pi/2, 0]
    s12 = qc.run(exe, memory_map={'theta': theta, 'phi': phi, 't': t}) # Stores the output of the circuit run.
    counts_s12 = Counter([''.join(list(map(str, elem))) for elem in s12])
    
    t = [np.pi, np.pi/2]
    s23 = qc.run(exe, memory_map={'theta': theta, 'phi': phi, 't': t}) # Stores the output of the circuit run.
    counts_s23 = Counter([''.join(list(map(str, elem))) for elem in s23])
    
    t = [np.pi/2, np.pi]
    s31 = qc.run(exe, memory_map={'theta': theta, 'phi': phi, 't': t}) # Stores the output of the circuit run.
    counts_s31 = Counter([''.join(list(map(str, elem))) for elem in s31])
    
    t = [0, 0]
    ssingles = qc.run(exe, memory_map={'theta': theta, 'phi': phi, 't': t}) # Stores the output of the circuit run.
    counts_ssingles = Counter([''.join(list(map(str, elem))) for elem in ssingles])
    
#     outcomes['Kappa'] = (3*outcomes['P_123'] - 2*(outcomes['P_12']+outcomes['P_23']+outcomes['P_31']) + outcomes['Singles']['00'] + outcomes['Singles']['10'] + outcomes['Singles']['01']) / N_SHOTS
    
    return {'Clicks_123': s123, 'Counts_123': counts_s123, 'Clicks_12': s12, 'Counts_12': counts_s12, 'Clicks_23': s23, 'Counts_23': counts_s23, 'Clicks_31': s31, 'Counts_31': counts_s31, 'Clicks_singles': ssingles, 'Counts_singles': counts_ssingles}





def run_born(q1, q2, trial, engine, states):
    global qc
    global exe
    result_list = []
    print(f'Engine requested: {engine}')
    if engine == 'qvm':
        qc = get_qc('Aspen-9', as_qvm=True) # Initialise QPU.
    elif engine == 'Aspen':
        qc = get_qc('Aspen-9')
    else:
        qc = get_qc('2q-qvm')
    # qc = get_qc('2q-qvm')

    circ = circuit_123(q1,q2)
    exe = qc.compile(circ)
    
    print('Running Bell-state measurements')
    for i in range(len(states)):
        data = states[i]
#         if engine == 'Aspen-9':
#             data['Device_details'] = qc.device.get_isa()

        data[f'Counts_{engine}'] = f(data['State_params'])

#         states.append(data)

        print(f'Done with iteration {i}', end='\r')
    
    folder = f'born_{engine}_{datetime.date.today()}_{q1}_{q2}_bits_{N_SHOTS}_shots_trial_{trial}'
    print(f'Creating folder {folder}')
    os.system(f'mkdir -p {folder}')
    with open(f'{folder}/result_list_trial_{trial}', 'wb') as file:
        pickle.dump(states, file)
    print(f'Results saved in file {folder}/result_list_trial_{trial}')
    
    print()
    print('Completed.')
    return states


def get_kappa(outcomes):
    counts_123 = Counter([''.join(list(map(str, elem))) for elem in outcomes['Clicks_123']])['10']
    counts_12 = Counter([''.join(list(map(str, elem))) for elem in outcomes['Clicks_12']])['10']
    counts_23 = Counter([''.join(list(map(str, elem))) for elem in outcomes['Clicks_23']])['10']
    counts_31 = Counter([''.join(list(map(str, elem))) for elem in outcomes['Clicks_31']])['10']
    counts_singles = Counter([''.join(list(map(str, elem))) for elem in outcomes['Clicks_singles']])
    
    kappa = 3 * counts_123 - 2 * (counts_12 + counts_23 + counts_31) + counts_singles['00'] + counts_singles['01'] + counts_singles['10']
    return kappa/N_SHOTS


def compute_kappa(result_list, engine):
    for data in result_list:
        data[f'Kappa_{engine}'] = get_kappa(data[f'Counts_{engine}'])
        
    return result_list