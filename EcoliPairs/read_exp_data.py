from __future__ import division
import sys
sys.path.append('../')
sys.path.append('results/')
import itertools
from tools.ancillary.importData import importData
import re

def read_exp_data():
    """
    Reads the experimental data and returns the results as python dicitonaries
    """
    # Initial cell concentration
    init_cell_conc = 7.5e6

    # Load the experimental growth data
    day1Rep1Inst = importData(inputFile = 'expData/day1Rep1.txt',delType = 'tab',dataType = 'float')
    day1Rep1 = day1Rep1Inst.run()

    day1Rep2Inst = importData(inputFile = 'expData/day1Rep2.txt',delType = 'tab',dataType = 'float')
    day1Rep2 = day1Rep2Inst.run()

    day2Rep1Inst = importData(inputFile = 'expData/day2Rep1.txt',delType = 'tab',dataType = 'float')
    day2Rep1 = day2Rep1Inst.run()

    day2Rep2Inst = importData(inputFile = 'expData/day2Rep2.txt',delType = 'tab',dataType = 'float')
    day2Rep2 = day2Rep2Inst.run()

    day3Rep1Inst = importData(inputFile = 'expData/day3Rep1.txt',delType = 'tab',dataType = 'float')
    day3Rep1 = day3Rep1Inst.run()

    day3Rep2Inst = importData(inputFile = 'expData/day3Rep2.txt',delType = 'tab',dataType = 'float')
    day3Rep2 = day3Rep2Inst.run()

    day4Rep1Inst = importData(inputFile = 'expData/day4Rep1.txt',delType = 'tab',dataType = 'float')
    day4Rep1 = day4Rep1Inst.run()

    day4Rep2Inst = importData(inputFile = 'expData/day4Rep2.txt',delType = 'tab',dataType = 'float')
    day4Rep2 = day4Rep2Inst.run()

    # List of all mutants in these pairs
    mutants_list = list(set([m for m1m2 in day1Rep1.keys() for m in m1m2]))

    # Names of the mutant pairs
    mutant_pairs = list(itertools.combinations([m for m in mutants_list],r=2))

    # Also include in the mutants list all mutant pairs where both mutants are the same
    for m in mutants_list:
        mutant_pairs.append((m,m))

    # Compute the measured fold change in growth of the pairs 
    # First find the average over replicates
    # Note that the data matrixes should actually be symmetric as (mutant1,mutant2) is the
    # same as (mutant2,mutant1), however, this is not always the case due to the experiiemntal
    # errors. As such we shouls also take average on the values of (mutant1,mutant2) and
    # (mutant2,mutant1). 
    day1Ave = dict([((m1,m2),(day1Rep1[(m1,m2)] + day1Rep1[(m2,m1)] + day1Rep2[(m1,m2)] + day1Rep2[(m2,m1)])/4) for (m1,m2) in mutant_pairs])

    day2Ave = dict([((m1,m2),(day2Rep1[(m1,m2)] + day2Rep1[(m2,m1)] + day2Rep2[(m1,m2)] + day2Rep2[(m2,m1)])/4) for (m1,m2) in mutant_pairs])

    day3Ave = dict([((m1,m2),(day3Rep1[(m1,m2)] + day3Rep1[(m2,m1)] + day3Rep2[(m1,m2)] + day3Rep2[(m2,m1)])/4) for (m1,m2) in mutant_pairs])

    day4Ave = dict([((m1,m2),(day4Rep1[(m1,m2)] + day4Rep1[(m2,m1)] + day4Rep2[(m1,m2)] + day4Rep2[(m2,m1)])/4) for (m1,m2) in mutant_pairs])

    return [day1Ave,day2Ave,day3Ave,day4Ave]

#-------------------------------
if __name__ == '__main__':
    [day1Ave,day2Ave,day3Ave,day4Ave] = read_exp_data()
