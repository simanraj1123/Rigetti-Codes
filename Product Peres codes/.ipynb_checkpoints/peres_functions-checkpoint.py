import numpy as np, matplotlib.pyplot as plt, random, time, datetime
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
from scipy.optimize import curve_fit as cf

Z_SCORE = 3
N_SHOTS = 10_000




########################
def damping_channel(damp_prob=.1):
    """
    Generate the Kraus operators corresponding to an amplitude damping
    noise channel.

    :params float damp_prob: The one-step damping probability.
    :return: A list [k1, k2] of the Kraus operators that parametrize the map.
    :rtype: list
    """
    damping_op = np.sqrt(damp_prob) * np.array([[0, 1],
                                                [0, 0]])

    residual_kraus = np.diag([1, np.sqrt(1-damp_prob)])
    return [residual_kraus, damping_op]

def append_kraus_to_gate(kraus_ops, g):
    """
    Follow a gate `g` by a Kraus map described by `kraus_ops`.

    :param list kraus_ops: The Kraus operators.
    :param numpy.ndarray g: The unitary gate.
    :return: A list of transformed Kraus operators.
    """
    return [kj.dot(g) for kj in kraus_ops]


def append_damping_to_gate(gate, damp_prob=.1):
    """
    Generate the Kraus operators corresponding to a given unitary
    single qubit gate followed by an amplitude damping noise channel.

    :params np.ndarray|list gate: The 2x2 unitary gate matrix.
    :params float damp_prob: The one-step damping probability.
    :return: A list [k1, k2] of the Kraus operators that parametrize the map.
    :rtype: list
    """
    return append_kraus_to_gate(damping_channel(damp_prob), gate)

def dephasing_kraus_map(p=.1):
    """
    Generate the Kraus operators corresponding to a dephasing channel.

    :params float p: The one-step dephasing probability.
    :return: A list [k1, k2] of the Kraus operators that parametrize the map.
    :rtype: list
    """
    return [np.sqrt(1-p)*np.eye(2), np.sqrt(p)*np.diag([1, -1])]

def tensor_kraus_maps(k1, k2):
    """
    Generate the Kraus map corresponding to the composition
    of two maps on different qubits.

    :param list k1: The Kraus operators for the first qubit.
    :param list k2: The Kraus operators for the second qubit.
    :return: A list of tensored Kraus operators.
    """
    return [np.kron(k1j, k2l) for k1j in k1 for k2l in k2]
########################

def params_real():
	'''
	Generates parameters to prepare random REAL quantum states.
	'''
	theta = np.arccos(1 - 2 * np.array([random.uniform(0,1) for _ in range(3)]))
	phi = np.array([(np.pi)*random.randint(0,1) for _ in range(3)])
	params = zip(theta, phi)
	return list(params)
def params_complex():
	'''
	Generates parameters to prepare COMPLEX quantum states.
	'''
	theta = np.arccos(1 - 2 * np.array([random.uniform(0,1) for _ in range(3)]))
	phi = np.array([2*np.pi*random.uniform(0,1) for _ in range(3)])
	params = zip(theta, phi)
	return list(params)

# Gammas for different pairs of states.
def g(u):
    '''
    Calls the sigma function with different values of parameters correponding to the configurations, |ψ12>, |ψ1> and |ψ2>. Returns a
    dictionary with configurations as keys and output as values (which are lists).
    '''
    params = list(zip(*u)) # Unpack parameters
    theta, phi = params[0], params[1] # Store thetas and phis in seperate tuples.
    
    s12 = qc.run(exe, memory_map={'theta': theta, 'phi': phi}) # Stores the output of the circuit run.
    counts_s12 = Counter([''.join(list(map(str, elem))) for elem in s12])
    
    return {'Clicks': s12, 'Counts': counts_s12}

# Computing all the three gammas.
# @dura # To calculate the time taken for all the circuits to run.
def f(u):
	'''
	Calls the g function to run the circuit for different configurations and returns a dictionary with 'a', 'b', 'c' as keys and the corresponding 
	outputs of the three configurations. This marks the end of what the Quantum computer must be used for. After this it is all about post-
	processing the data.
	'''
	alpha = g([u[0], u[1]]) # Running for alpha
	beta = g([u[1], u[2]]) # Running for beta
	gamma = g([u[2], u[0]]) # Running for gamma

	res = {'a': alpha, 'b': beta, 'c': gamma}

	return res

def get_gammas(counts_bell, counts_comp, N):
    res = {}
    for gamma in counts_bell.keys():
        counts_12 = Counter([''.join(list(map(str, elem))) for elem in counts_bell[gamma]['Clicks'][:N]])['01']
        counts_1 = Counter([''.join(list(map(str, elem))) for elem in counts_comp[gamma]['Clicks'][:N]])['01']
        counts_2 = Counter([''.join(list(map(str, elem))) for elem in counts_comp[gamma]['Clicks'][:N]])['10']
        #Counter([''.join(list(map(str, elem))) for elem in s2])
        g = (2*counts_12 - counts_1 - counts_2) / (2 * np.sqrt(counts_1*counts_2))
        res[gamma] = g

    res['F'] = res['a']**2 + res['b']**2 + res['c']**2 - 2 * res['a'] * res['b'] * res['c']
    return res
p=0.1
def circuit_bell(qubit1, qubit2):
    circ = Program()
    corrupted_CZ = append_kraus_to_gate(
    tensor_kraus_maps(
        dephasing_kraus_map(p),
        dephasing_kraus_map(p)
    ),
    np.diag([1, 1, 1, -1]))
    
    c = circ.declare('ro', 'BIT', 2)
    theta = circ.declare('theta', 'REAL', 2)
    phi = circ.declare('phi', 'REAL', 2)
    
    # Preparation of states.
    circ += RY(theta[0], qubit1)
    circ += RZ(phi[0], qubit1)
    
    circ += RY(theta[1], qubit2)
    circ += RZ(phi[1], qubit2)
    
#     circ += X(qubit1)
    
    # Measuring in psi+ basis
    circ += CNOT(qubit1, qubit2)
    circ += H(qubit1)
#     circ += H(qubit1)

    circ += I(qubit1)
    circ += I(qubit2)
    
    
    circ += MEASURE(qubit1, c[0])
    circ += MEASURE(qubit2, c[1])
    
    circ.define_noisy_gate("I", [qubit1], append_damping_to_gate(np.eye(2), 0.2))
    circ.define_noisy_gate("I", [qubit2], append_damping_to_gate(np.eye(2), 0.2))
    circ.define_noisy_gate("CZ", [qubit1, qubit2], corrupted_CZ)
    
    
    
    circ.wrap_in_numshots_loop(N_SHOTS)
    
    return circ

def circuit_comp(qubit1, qubit2):
    circ = Program()
    corrupted_CZ = append_kraus_to_gate(
    tensor_kraus_maps(
        dephasing_kraus_map(p),
        dephasing_kraus_map(p)
    ),
    np.diag([1, 1, 1, -1]))
    
    c = circ.declare('ro', 'BIT', 2)
    theta = circ.declare('theta', 'REAL', 2)
    phi = circ.declare('phi', 'REAL', 2)
    
    # Preparation of states.
    circ += RY(theta[0], qubit1)
    circ += RZ(phi[0], qubit1)
    
    circ += RY(theta[1], qubit2)
    circ += RZ(phi[1], qubit2)
    
    circ += I(qubit1)
    circ += I(qubit2)
    
    
    circ += MEASURE(qubit1, c[0])
    circ += MEASURE(qubit2, c[1])
    
    circ.define_noisy_gate("I", [qubit1], append_damping_to_gate(np.eye(2), 0.2))
    circ.define_noisy_gate("I", [qubit2], append_damping_to_gate(np.eye(2), 0.2))
    circ.define_noisy_gate("CZ", [qubit1, qubit2], corrupted_CZ)
    
    
    
    circ.wrap_in_numshots_loop(N_SHOTS)
    
    return circ

def gamma_theory(data):
    u = data['State_params']
    g12 = np.cos(u[1][1] - u[0][1])
    g23 = np.cos(u[2][1] - u[1][1])
    g31 = np.cos(u[0][1] - u[2][1])
    
    data['Gammas_theory'] = {'a': g12, 'b': g23, 'c': g31}
    
    f = g12**2 + g23**2 + g31**2 - 2 * g12 * g23 * g31
    
    data['Gammas_theory']['F'] = f
    
    return data

qc = None
exe = None
def run_peres(q1, q2, trial, engine, iters, specified_params='None'):
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

    circ = circuit_bell(q1,q2)
    exe = qc.compile(circ)
    
    print('Running Bell-state measurements')
    for i in range(iters):
        data = {}
        data['State_params'] = params_complex() if specified_params == 'None' else specified_params[i]

        data['Counts_bell'] = f(data['State_params'])

        result_list.append(data)

        print(f'Done with iteration {i}', end='\r')
    
    print('\n')
    print('Running Computational measurements')
    circ = circuit_comp(q1,q2)
#     global exe
    exe = qc.compile(circ)
    i=0
    for data in result_list:
        data['Counts_comp'] = f(data['State_params'])
        print(f'Done with iteration {i}', end='\r')
        i += 1
    
    folder = f'product_peres_{engine}_{datetime.date.today()}_{q1}_{q2}_bits_{N_SHOTS}_shots_trial_{trial}'
    print(f'Creating folder {folder}')
    os.system(f'mkdir -p {folder}')
    with open(f'{folder}/result_list_trial_{trial}', 'wb') as file:
        pickle.dump(result_list, file)
    print(f'Results saved in file {folder}/result_list_trial_{trial}')
    
    print()
    print('Completed.')
    return result_list

def compute_gammas(result_list, sub_N, cfs, cfs_as):
    fig = plt.figure(figsize=(17, 10))
    for sub_N in [sub_N]:
        all_a = []
        all_b = []
        all_c = []
        all_F = []

        all_a_theory = []
        all_b_theory = []
        all_c_theory = []
        all_F_theory = []

        all_a_errors = []
        all_b_errors = []
        all_c_errors = []
        all_F_errors = []

        for data in result_list:
            data['Gamma'] = get_gammas(data['Counts_bell'], data['Counts_comp'], sub_N)
            all_a.append(data['Gamma']['a'])
            all_b.append(data['Gamma']['b'])
            all_c.append(data['Gamma']['c'])
            all_F.append(data['Gamma']['F'])
            
            data = gamma_theory(data)
            all_a_theory.append(data['Gammas_theory']['a'])
            all_b_theory.append(data['Gammas_theory']['b'])
            all_c_theory.append(data['Gammas_theory']['c'])
            all_F_theory.append(data['Gammas_theory']['F'])

        x = np.array(list(range(1,len(result_list)+1)))
    #     fig = plt.figure(figsize=(17,10))
    
        #### Extra code ####
#         with open('all_cfs.pkl', 'rb') as file:
#             cfs = np.array(pickle.load(file))
        cf12 = cfs[:,0]
        cf23 = cfs[:,1]
        cf31 = cfs[:,2]
        cfF = cfs[:,3]

        m12 = list(map(np.mean, cf12))
        m23 = list(map(np.mean, cf23))
        m31 = list(map(np.mean, cf31))
        mF = list(map(np.mean, cfF))
        
        cf12_as = cfs_as[:,0]
        cf23_as = cfs_as[:,1]
        cf31_as = cfs_as[:,2]
        cfF_as = cfs_as[:,3]

        m12_as = list(map(np.mean, cf12_as))
        m23_as = list(map(np.mean, cf23_as))
        m31_as = list(map(np.mean, cf31_as))
        mF_as = list(map(np.mean, cfF_as))
        #######################
        plt.subplot(2,2,1)
        plt.plot(x, all_a, 'o', color='blue', markersize=4)
#         plt.plot(x, all_a_theory, '*', color='red', markersize=4)
        plt.errorbar(x, m12, yerr=[m12 - cf12[:,0], cf12[:,1] - m12], fmt='.', marker='', capsize=4, color='red')
        plt.errorbar(x, m12_as, yerr=[m12_as - cf12_as[:,0], cf12_as[:,1] - m12_as], fmt='.', marker='', capsize=4, color='blue')
        plt.axhline(y=1, ls='dashed', alpha=0.5)
        plt.axhline(y=-1, ls='dashed', alpha=0.5)
        plt.xticks(x)
        plt.xlabel('$\\gamma_{12}$', size=14)

        plt.subplot(2,2,2)
        plt.plot(x, all_b, 'o', color='blue', markersize=4)
#         plt.plot(x, all_b_theory, '*', color='blue', markersize=4)
        plt.errorbar(x, m23, yerr=[m23 - cf23[:,0], cf23[:,1] - m23], fmt='.', marker='', capsize=4, color='red')
        plt.errorbar(x, m23_as, yerr=[m23_as - cf23_as[:,0], cf23_as[:,1] - m23_as], fmt='.', marker='', capsize=4, color='blue')
        plt.axhline(y=1, ls='dashed', alpha=0.5)
        plt.axhline(y=-1, ls='dashed', alpha=0.5)
        plt.xticks(x)
        plt.xlabel('$\\gamma_{23}$', size=14)

        plt.subplot(2,2,3)
        plt.plot(x, all_c, 'o', color='blue', markersize=4)
#         plt.plot(x, all_c_theory, '*', color='darkgreen', markersize=4)
        plt.errorbar(x, m31, yerr=[m31 - cf31[:,0], cf31[:,1] - m31], fmt='.', marker='', capsize=4, color='red')
        plt.errorbar(x, m31_as, yerr=[m31_as - cf31_as[:,0], cf31_as[:,1] - m31_as], fmt='.', marker='', capsize=4, color='blue')
        plt.axhline(y=1, ls='dashed', alpha=0.5)
        plt.axhline(y=-1, ls='dashed', alpha=0.5)
        plt.xticks(x)
        plt.xlabel('$\\gamma_{31}$', size=14)

        plt.subplot(2,2,4)
        plt.plot(x, all_F, 'o', color='blue', markersize=4)
#         plt.plot(x, all_F_theory, '*', color='maroon', markersize=4)
        plt.errorbar(x, mF, yerr=[mF - cfF[:,0], cfF[:,1] - mF], fmt='.', marker='', capsize=4, color='red')
        plt.errorbar(x, mF_as, yerr=[mF_as - cfF_as[:,0], cfF_as[:,1] - mF_as], fmt='.', marker='', capsize=4, color='blue')
        plt.axhline(y=1, ls='dashed', alpha=0.5)
        plt.axhline(y=-1, ls='dashed', alpha=0.5)
        plt.xticks(x)
        plt.xlabel('$F$', size=14)
        fit = lambda x, a, b: a*x + b
        xdata = np.arange(1,len(result_list)+1)
        ydata = [result_list[i]['Gamma']['F'] for i in range(len(result_list))]
        popt, pcov = cf(fit, xdata, ydata)
        plt.plot(xdata, fit(xdata, *popt), ls='dashed', alpha=0.5)
        print(popt)

    return fig




