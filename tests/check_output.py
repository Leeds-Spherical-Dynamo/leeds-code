"""
Script to test output of Dynamo code run on 
Christensen 2001 benchmark.

Currently averages last 500 lines of output
on vel_energy.dat, mag_energy.dat and
rot_torque.dat. This maybe changed to ensure
that a set amount of duration has passed.
"""

import csv
import numpy as np

# Load dat file and calculate average of column data
def calc_average(case, datfile):
    columns=[]
    with open("./test-{}/{}".format(case, datfile)) as dat:
        for column in zip(*[line for line in csv.reader(dat, delimiter=" ", skipinitialspace=True)]):
            columns.append(np.array(column[1::]).astype(float))
    average = np.average(columns[1][-500::])
    return average

# Create Dictionary and add benchmark values
cases = {}
cases['case0'] = {'ke': {'bench':851.6181, 'tol':0.72994 }}
cases['case1'] = {'ke': {'bench': 11231.222, 'tol': 7.2994}}
cases['case1']['me'] = {'bench': 228620.5, 'tol': 145.98}
cases['case2'] = {'ke': {'bench': 15470.3, 'tol': 18.2485}}
cases['case2']['me'] = {'bench': 308619, 'tol': 145.98}
cases['case2']['rot'] = {'bench': -13.2975, 'tol': 0.3}

cases['case0_symmetry'] = {'ke': {'bench':851.6181, 'tol':0.72994 }}
cases['case1_symmetry'] = {'ke': {'bench': 11231.222, 'tol': 7.2994}}
cases['case1_symmetry']['me'] = {'bench': 228620.5, 'tol': 145.98}
cases['case2_symmetry'] = {'ke': {'bench': 15470.3, 'tol': 18.2485}}
cases['case2_symmetry']['me'] = {'bench': 308619, 'tol': 145.98}
cases['case2_symmetry']['rot'] = {'bench': -13.2975, 'tol': 0.3}

# Add simulation results to dict
cases['case0']['ke']['test'] = calc_average('case0', "vel_energy.dat")
cases['case1']['ke']['test'] = calc_average('case1', "vel_energy.dat")
cases['case1']['me']['test'] = calc_average('case1', "mag_energy.dat")
cases['case2']['ke']['test'] = calc_average('case2', "vel_energy.dat")
cases['case2']['me']['test'] = calc_average('case2', "mag_energy.dat")
cases['case2']['rot']['test'] = calc_average('case2', "rot_torque.dat")

cases['case0_symmetry']['ke']['test'] = calc_average('case0_symmetry', "vel_energy.dat")
cases['case1_symmetry']['ke']['test'] = calc_average('case1_symmetry', "vel_energy.dat")
cases['case1_symmetry']['me']['test'] = calc_average('case1_symmetry', "mag_energy.dat")
cases['case2_symmetry']['ke']['test'] = calc_average('case2_symmetry', "vel_energy.dat")
cases['case2_symmetry']['me']['test'] = calc_average('case2_symmetry', "mag_energy.dat")
cases['case2_symmetry']['rot']['test'] = calc_average('case2_symmetry', "rot_torque.dat")

# Check that results lie within benchmark tolerances
for case in cases:
    for test in cases[case]:
        print(case, test)
        print(cases[case][test])
        bench=cases[case][test]['bench']
        tol = cases[case][test]['tol']
        val = cases[case][test]['test']
        max = bench+tol
        min = bench - tol
        error = abs(val-bench)
        if val > min and val < max:
            print("Pass, error:{}\n".format(error))
        else:
            print("fail, error:{}\n".format(error))
