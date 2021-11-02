# import numpy as np, matplotlib.pyplot as plt, random, time, datetime
# from functools import lru_cache
# from pyquil import Program, get_qc
# from pyquil.gates import *
# import os, sys
# from pyquil.quilatom import quil_sin, quil_cos, Parameter
# from pyquil.quilbase import DefGate
# from pyquil.latex import display, to_latex
# # import Peres_helpers as hf
# import pickle
# from collections import Counter
# from scipy.optimize import curve_fit as cf
# sys.path.append('binomial_cython')
# from binomial import binomial_dist


# The z-score and number of shots for circuits
Z_SCORE = 3
N_SHOTS = 10_000

# # Generarating parameters theta and phi for randomly generating states on the Bloch-sphere.
# def params_real():
# 	'''
# 	Generates parameters to prepare random REAL quantum states.
# 	'''
# 	theta = np.arccos(1 - 2 * np.array([random.uniform(0,1) for _ in range(3)]))
# 	phi = np.array([(np.pi)*random.randint(0,1) for _ in range(3)])
# 	params = zip(theta, phi)
# 	return list(params)
# def params_complex():
# 	'''
# 	Generates parameters to prepare COMPLEX quantum states.
# 	'''
# 	theta = np.arccos(1 - 2 * np.array([random.uniform(0,1) for _ in range(3)]))
# 	phi = np.array([2*np.pi*random.uniform(0,1) for _ in range(3)])
# 	params = zip(theta, phi)
# 	return list(params)

# # Gammas for different pairs of states. This function returns the click-counts for the different configs.
# def g(u):
#     '''
#     Calls the sigma function with different values of parameters correponding to the configurations, |ψ12>, |ψ1> and |ψ2>. Returns a
#     dictionary with configurations as keys and output as values (which are lists).
#     '''
#     params = list(zip(*u)) # Unpack parameters
#     theta, phi = params[0], params[1] # Store thetas and phis in seperate tuples.
    
#     s12 = qc.run(exe, memory_map={'theta': theta, 'phi': phi}) # Stores the output of the circuit run.
#     counts_s12 = Counter([''.join(list(map(str, elem))) for elem in s12])
    
#     return {'Clicks': s12, 'Counts': counts_s12}

# # Computing all the three gammas. This function returns the click-counts for all the configs.
# def f(u):
# 	'''
# 	Calls the g function to run the circuit for different configurations and returns a dictionary with 'a', 'b', 'c' as keys and the corresponding 
# 	outputs of the three configurations. This marks the end of what the Quantum computer must be used for. After this it is all about post-
# 	processing the data.
# 	'''
# 	alpha = g([u[0], u[1]]) # Running for alpha
# 	beta = g([u[1], u[2]]) # Running for beta
# 	gamma = g([u[2], u[0]]) # Running for gamma

# 	res = {'a': alpha, 'b': beta, 'c': gamma}

# 	return res

# # Calculate the gamma values and the F value.
# def get_gammas(counts_bell, counts_comp):
#     res = {}
#     for gamma in counts_bell.keys():
#         counts_12 = Counter([''.join(list(map(str, elem))) for elem in counts_bell[gamma]['Clicks']])['01']
#         counts_1 = Counter([''.join(list(map(str, elem))) for elem in counts_comp[gamma]['Clicks']])['01']
#         counts_2 = Counter([''.join(list(map(str, elem))) for elem in counts_comp[gamma]['Clicks']])['10']
#         #Counter([''.join(list(map(str, elem))) for elem in s2])
#         g = (2*counts_12 - counts_1 - counts_2) / (2 * np.sqrt(counts_1*counts_2))
#         res[gamma] = g

#     res['F'] = res['a']**2 + res['b']**2 + res['c']**2 - 2 * res['a'] * res['b'] * res['c']
#     return res


p=0.1
# The circuit for constructing the product state and measuring in Bell basis.
# def circuit_bell(qubit1, qubit2):
#     circ = Program()
    
#     c = circ.declare('ro', 'BIT', 2)
#     theta = circ.declare('theta', 'REAL', 2)
#     phi = circ.declare('phi', 'REAL', 2)
    
#     # Preparation of states.
#     circ += RY(theta[0], qubit1)
#     circ += RZ(phi[0], qubit1)
    
#     circ += RY(theta[1], qubit2)
#     circ += RZ(phi[1], qubit2)
    
#     # Measuring in psi+ basis
#     circ += CNOT(qubit1, qubit2)
#     circ += H(qubit1)

#     circ += MEASURE(qubit1, c[0])
#     circ += MEASURE(qubit2, c[1])
    
#     circ.wrap_in_numshots_loop(N_SHOTS)
    
#     return circ

# # The circuit for constructing the product state and measuring in Computational basis.
# def circuit_comp(qubit1, qubit2):
#     circ = Program()
    
#     c = circ.declare('ro', 'BIT', 2)
#     theta = circ.declare('theta', 'REAL', 2)
#     phi = circ.declare('phi', 'REAL', 2)
    
#     # Preparation of states.
#     circ += RY(theta[0], qubit1)
#     circ += RZ(phi[0], qubit1)
    
#     circ += RY(theta[1], qubit2)
#     circ += RZ(phi[1], qubit2)
    
#     circ += MEASURE(qubit1, c[0])
#     circ += MEASURE(qubit2, c[1])
    
#     circ.wrap_in_numshots_loop(N_SHOTS)
    
#     return circ

# # Theoretical value of gammas and F.
# def gamma_theory(data):
#     u = data['State_params']
#     g12 = np.cos(u[1][1] - u[0][1])
#     g23 = np.cos(u[2][1] - u[1][1])
#     g31 = np.cos(u[0][1] - u[2][1])
    
#     data['Gammas_theory'] = {'a': g12, 'b': g23, 'c': g31}
    
#     f = g12**2 + g23**2 + g31**2 - 2 * g12 * g23 * g31
    
#     data['Gammas_theory']['F'] = f
    
#     return data

# # Declaring variables to store quantum circuit and the executable.
# qc = None
# exe = None
# def run_peres(q1, q2, trial, engine, iters, specified_params='None'):
#     global qc
#     global exe
#     result_list = []
#     print(f'Engine requested: {engine}')
#     if engine == 'qvm':
#         qc = get_qc('Aspen-9', as_qvm=True) # Initialise QPU.
#     elif engine == 'Aspen':
#         qc = get_qc('Aspen-9')
#     else:
#         qc = get_qc('2q-qvm')
#     # qc = get_qc('2q-qvm')

#     circ = circuit_bell(q1,q2)
#     exe = qc.compile(circ)
    
#     print('Running Bell-state measurements')
#     for i in range(iters):
#         data = {}
#         data['State_params'] = params_complex() if specified_params == 'None' else specified_params[i]

#         data['Counts_bell'] = f(data['State_params'])

#         result_list.append(data)

#         print(f'Done with iteration {i}', end='\r')
    
#     print('\n')
#     print('Running Computational measurements')
#     circ = circuit_comp(q1,q2)
# #     global exe
#     exe = qc.compile(circ)
#     i=0
#     for data in result_list:
#         data['Counts_comp'] = f(data['State_params'])
#         print(f'Done with iteration {i}', end='\r')
#         i += 1
    
# #     for data in result_list:
# #         data['Gamma'] = get_gammas(data['Counts_bell'], data['Counts_comp'])

# #         data = gamma_theory(data)
    
#     folder = f'product_peres_{engine}_{datetime.date.today()}_{q1}_{q2}_bits_{N_SHOTS}_shots_trial_{trial}'
#     print(f'Creating folder {folder}')
#     os.system(f'mkdir -p {folder}')
#     with open(f'{folder}/result_list_trial_{trial}', 'wb') as file:
#         pickle.dump(result_list, file)
#     print(f'Results saved in file {folder}/result_list_trial_{trial}')
    
#     print()
#     print('Completed.')
#     return result_list


# def compute_gammas(result_list):
#     for data in result_list:
#         data['Gamma'] = get_gammas(data['Counts_bell'], data['Counts_comp'])

#         data = gamma_theory(data)
        
#     return result_list
# 
# def plot_gammas(result_list):
#     fig = plt.figure(figsize=(17, 10))
#     all_a = []
#     all_b = []
#     all_c = []
#     all_F = []

#     all_a_theory = []
#     all_b_theory = []
#     all_c_theory = []
#     all_F_theory = []

#     all_a_errors = []
#     all_b_errors = []
#     all_c_errors = []
#     all_F_errors = []

#     for data in result_list:
#         all_a.append(data['Gamma']['a'])
#         all_b.append(data['Gamma']['b'])
#         all_c.append(data['Gamma']['c'])
#         all_F.append(data['Gamma']['F'])

#         all_a_theory.append(data['Gammas_theory']['a'])
#         all_b_theory.append(data['Gammas_theory']['b'])
#         all_c_theory.append(data['Gammas_theory']['c'])
#         all_F_theory.append(data['Gammas_theory']['F'])

#     x = np.array(list(range(1,len(result_list)+1)))
# #     fig = plt.figure(figsize=(17,10))

#     plt.subplot(2,2,1)
#     plt.plot(x, all_a, 'o', color='blue', markersize=4)
# #         plt.plot(x, all_a_theory, '*', color='red', markersize=4)
# #         plt.errorbar(x, m12, yerr=[m12 - cf12[:,0], cf12[:,1] - m12], fmt='.', marker='', capsize=4, color='red')
# #         plt.errorbar(x, m12_as, yerr=[m12_as - cf12_as[:,0], cf12_as[:,1] - m12_as], fmt='.', marker='', capsize=4, color='blue')
#     plt.axhline(y=1, ls='dashed', alpha=0.5)
#     plt.axhline(y=-1, ls='dashed', alpha=0.5)
#     plt.xticks(x)
#     plt.xlabel('$\\gamma_{12}$', size=14)

#     plt.subplot(2,2,2)
#     plt.plot(x, all_b, 'o', color='blue', markersize=4)
# #         plt.plot(x, all_b_theory, '*', color='blue', markersize=4)
# #         plt.errorbar(x, m23, yerr=[m23 - cf23[:,0], cf23[:,1] - m23], fmt='.', marker='', capsize=4, color='red')
# #         plt.errorbar(x, m23_as, yerr=[m23_as - cf23_as[:,0], cf23_as[:,1] - m23_as], fmt='.', marker='', capsize=4, color='blue')
#     plt.axhline(y=1, ls='dashed', alpha=0.5)
#     plt.axhline(y=-1, ls='dashed', alpha=0.5)
#     plt.xticks(x)
#     plt.xlabel('$\\gamma_{23}$', size=14)

#     plt.subplot(2,2,3)
#     plt.plot(x, all_c, 'o', color='blue', markersize=4)
# #         plt.plot(x, all_c_theory, '*', color='darkgreen', markersize=4)
# #         plt.errorbar(x, m31, yerr=[m31 - cf31[:,0], cf31[:,1] - m31], fmt='.', marker='', capsize=4, color='red')
# #         plt.errorbar(x, m31_as, yerr=[m31_as - cf31_as[:,0], cf31_as[:,1] - m31_as], fmt='.', marker='', capsize=4, color='blue')
#     plt.axhline(y=1, ls='dashed', alpha=0.5)
#     plt.axhline(y=-1, ls='dashed', alpha=0.5)
#     plt.xticks(x)
#     plt.xlabel('$\\gamma_{31}$', size=14)

#     plt.subplot(2,2,4)
#     plt.plot(x, all_F, 'o', color='blue', markersize=4)
# #         plt.plot(x, all_F_theory, '*', color='maroon', markersize=4)
# #         plt.errorbar(x, mF, yerr=[mF - cfF[:,0], cfF[:,1] - mF], fmt='.', marker='', capsize=4, color='red')
# #         plt.errorbar(x, mF_as, yerr=[mF_as - cfF_as[:,0], cfF_as[:,1] - mF_as], fmt='.', marker='', capsize=4, color='blue')
#     plt.axhline(y=1, ls='dashed', alpha=0.5)
#     plt.axhline(y=-1, ls='dashed', alpha=0.5)
#     plt.xticks(x)
#     plt.xlabel('$F$', size=14)
# #         fit = lambda x, a, b: a*x + b
# #         xdata = np.arange(1,len(result_list)+1)
# #         ydata = [result_list[i]['Gamma']['F'] for i in range(len(result_list))]
# #         popt, pcov = cf(fit, xdata, ydata)
# #         plt.plot(xdata, fit(xdata, *popt), ls='dashed', alpha=0.5)
# #         print(popt)

#     return fig


# def get_theory_cfs(result_list, q, proc_no): 
#     print(f'Started process no.: {proc_no}.')
#     st = time.time()
#     for res in result_list:
#         counts = binomial_dist(np.array(res['State_params'], dtype=float), 10_000, 10_000)
#         counts = np.array(counts)

#         all_cfs = []

#         for i in range(4):
#             data = counts[:,i]
#             data.sort()
#             ca, cb = data[int(0.005 * len(data))], data[int(0.995 * len(data))]
#             all_cfs.append([ca, cb])

#         res['cfs_theory'] = {}
#         res['cfs_theory']['a'] = all_cfs[0]
#         res['cfs_theory']['b'] = all_cfs[1]
#         res['cfs_theory']['c'] = all_cfs[2]
#         res['cfs_theory']['F'] = all_cfs[3]
        
#     q.put((proc_no, result_list))
#     print(f'Ended process no.: {proc_no}. Time taken: {np.round(time.time() - st,3)}')
#     return None