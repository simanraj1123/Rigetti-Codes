import peres_functions as pf
import os
import matplotlib.pyplot as plt
# import pyquil.api._errors.UserMessageError as ue
os.system('clear')

q1, q2 = 12, 25
trial = 1
engine = 'qvm'
iters = 20

try:
    result_list = pf.run_peres(q1, q2, trial, engine, iters)
    plot = pf.compute_gammas(result_list, 10000)
    plot.savefig('test_fig_2.pdf')
    print('Plot saved as "test_fig_2.pdf".')

except Exception as e:
    print(e)
    
