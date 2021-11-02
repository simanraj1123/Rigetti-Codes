import numpy as np, matplotlib.pyplot as plt, random, time, datetime
from functools import lru_cache
from pyquil import Program, get_qc
from pyquil.gates import *
import os, sys
from pyquil.quilatom import quil_sin, quil_cos, Parameter
from pyquil.quilbase import DefGate
from pyquil.latex import display, to_latex
# import Peres_helpers as hf
import pickle
from collections import Counter
from scipy.optimize import curve_fit as cf
sys.path.append('binomial_cython')
from binomial import binomial_dist

Z_SCORE = 3
N_SHOTS = 10_000

# The circuit for constructing the product state and measuring in Bell basis.
def circuit_bell(qubit1, qubit2, config):
    circ = Program()
    
    c = circ.declare('ro', 'BIT', 2)
    theta = circ.declare('theta', 'REAL', 2)
    phi = circ.declare('phi', 'REAL', 2)
    
    # Preparation of states.
    circ += RY(theta[0], qubit1)
    circ += RZ(phi[0], qubit1)
    
    circ += RY(theta[1], qubit2)
    circ += RZ(phi[1], qubit2)
    
    # Measuring in psi+ basis
    # REMOVED FOR TEST OF SIMPLER CIRCUIT
#     circ += CNOT(qubit1, qubit2)
#     circ += H(qubit1)
    
    # Adding conditionals for different projections.
    if config == '1100':
        circ += H(qubit2)
    elif config == '1010':
        circ += H(qubit1)

    elif config == '0110':
#         circ += X(qubit2)
#         circ += H(qubit1)
        circ += CNOT(qubit1, qubit2)
        circ += H(qubit1)
        circ += X(qubit2)

    elif config == '1110':
#         circ += RY(2*np.arccos(1/np.sqrt(3)), qubit1)
#         circ += RY(np.arccos(1/np.sqrt(2)), qubit2)
#         circ += CNOT(qubit1, qubit2)
#         circ += RY(-np.arccos(1/np.sqrt(2)),qubit2)
        circ += CNOT(qubit1, qubit2)
        circ += RY(np.arccos(1/np.sqrt(2)),qubit2)
        circ += CNOT(qubit1, qubit2)
        circ += RY(-np.arccos(1/np.sqrt(2)), qubit2)
        circ += RY(-2*np.arccos(1/np.sqrt(3)), qubit1)
        
    circ += MEASURE(qubit1, c[0])
    circ += MEASURE(qubit2, c[1])
    
    circ.wrap_in_numshots_loop(N_SHOTS)
    
    return circ

# Declaring variables to store quantum circuit and the executable.
qc = None
exe = None
def run_peres(q1, q2, trial, engine, states):
    global qc
    global exe
    result_list = []
    print(f'Engine requested: {engine}')
    if engine == 'qvm':
        qc = get_qc('Aspen-9', as_qvm=True) # Initialise QPU.
        
    elif engine == 'Aspen':
        qc = get_qc('Aspen-9')
        
    elif engine == 'noisy-Aspen':
        qc = get_qc('Aspen-9', as_qvm=True, noisy=True)
        
    else:
        qc = get_qc('2q-qvm')
    
    configs = ['1100', '0110', '1010', '1110', '0000']
    
    for state in states:
        print(f'State: {state}\r')
        for config in configs:
            circ = circuit_bell(q1,q2, config)
            exe = qc.compile(circ)
#             print(f'Running projection {config}')
            data = state
#         if engine == 'Aspen-9':
#             data['Device_details'] = qc.device.get_isa()
            params = list(zip(*data['State_params'])) # Unpack parameters
            theta, phi = params[0], params[1] # Store thetas and phis in seperate tuples.
#             print(theta)
#             print(phi)

            s12 = qc.run(exe, memory_map={'theta': theta, 'phi': phi}) # Stores the output of the circuit run.
            counts_s12 = Counter([''.join(list(map(str, elem))) for elem in s12])

            res = {'Clicks': s12, 'Counts': counts_s12}

            data[f'Counts_{config}_{engine}'] = res

#         states.append(data)

#         print(f'Done with iteration {i}', end='\r')
    
    
    
#     for data in result_list:
#         data['Gamma'] = get_gammas(data['Counts_bell'], data['Counts_comp'])

#         data = gamma_theory(data)
    
    folder = f'product_peres_{engine}_{datetime.date.today()}_{q1}_{q2}_bits_{N_SHOTS}_shots_trial_{trial}'
    print(f'Creating folder {folder}')
    os.system(f'mkdir -p {folder}')
    with open(f'{folder}/result_list_trial_{trial}', 'wb') as file:
        pickle.dump(states, file)
    print(f'Results saved in file {folder}/result_list_trial_{trial}')
    
    print()
    print('Completed.')
    return states

def get_good_qbits(cz_fid, meas_fid):
    qc = get_qc('Aspen-9')
    details = qc.device.get_isa()

    good_edges = []
    for edge in details.edges:
        gates = edge.gates
        for gate in gates:
            if gate.operator == 'CZ':
                try:
                    if gate.fidelity > cz_fid:
                        good_edges.append(edge)
                except TypeError:
                    continue
    good_qubits = {}
    for edge in good_edges[:]:
        q1, q2 = edge.targets
        for q in details.qubits:
            if q.id in [q1,q2]:
                gates = q.gates
                for gate in gates:
                    if gate.operator == 'MEASURE':
                        if gate.fidelity < meas_fid:
                            try:
                                good_edges.remove(edge)
                            except ValueError:
#                                 print('Edge already removed')
                                continue
    good_qubits = []
    for edge in good_edges:
        q1, q2 = edge.targets
        good_qubits.append({'Edge': edge})
        for q in details.qubits:
            if q.id == q1:
                good_qubits[-1]['Qubit1'] = q
            elif q.id == q2:
                good_qubits[-1]['Qubit2'] = q
    return good_qubits