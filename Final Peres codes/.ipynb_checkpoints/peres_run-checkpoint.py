import peres_functions as pf
import pickle
import numpy as np, matplotlib.pyplot as plt
from collections import Counter
import multiprocessing as mp

q1, q2 = 12, 25
trial = 1
engine = '2q-qvm'
iters = 10

specific_states = [pf.params_complex() for _ in range(10)]

try:
    result_list = pf.run_peres(q1, q2, trial, engine, iters, specific_states)
except Exception as e:
    print(e)
    
result_list = pf.compute_gammas(result_list)



qout = mp.Queue()

if __name__ == '__main__':
    incre = 2
    procs = [mp.Process(target=pf.get_theory_cfs, args=(result_list[i*incre:(i+1)*incre], qout, i)) for i in range(5)]
    for proc in procs:
        proc.start()
    for proc in procs:
        proc.join()
        
    all_results = []
    for proc in procs:
        all_results.append(qout.get())
    
    print(all_results)