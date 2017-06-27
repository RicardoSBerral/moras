#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
from math import sqrt
from decimal import Decimal
import re
import csv

class DataInfo:
    def __init__(self, data, n_variables, number_objectives):
        self.data = data
        self.n_variables = n_variables
        self.number_objectives = number_objectives
        self.first_variable = n_variables - number_objectives

prog = "./nucleo/cons"
n_partitions = 10
min_size_forest = 10
max_size_forest = 100
step_size_forest = 10

data_wtr = DataInfo("wtr", 22, 5)
data_cosmo = DataInfo("music2_z0.00_500_1e12", 30, 7)

def call_program (*args):
    return subprocess.Popen(prog + " " + " ".join(str(x) for x in args), shell=True, stdout=subprocess.PIPE).stdout.read()

def get_num_nodes (str):
    m = re.search("Num nodos: (\d+) ", str)
    if m:
        return int(m.group(1))

def average_and_sd (list):
    average = sum(list)/float(len(list))
    average_of_squares = sum(x**2 for x in list)/float(len(list))
    sd = sqrt(average_of_squares - average**2)
    return average, sd

def average_with_sd (list):
    average, sd = average_and_sd (list)
    return "{:.2E} ± {:.3E}".format(Decimal(average), Decimal(sd))

def get_results_multi (dataInfo, ret, matrix):
    lines = ret.splitlines()[-dataInfo.number_objectives-1:-1]
    for i in xrange(0, dataInfo.number_objectives):
        m = re.match("\d+: ([\.\d]+)", lines[i])
        if m:
            matrix[i].append(float(m.group(1)))
    return matrix

def get_results_multi_variable_trees_number (dataInfo, ret, matrix, size):
    lines = ret.splitlines()[-dataInfo.number_objectives-1:-1]
    for i in xrange(0, dataInfo.number_objectives):
        m = re.match("\d+: ([\.\d]+)", lines[i])
        if m:
            matrix[i][size].append(float(m.group(1)))
    return matrix

def single_objective_tree(dataInfo):
    for i in xrange(dataInfo.first_variable, dataInfo.n_variables):
        partial_results = []
        partial_nodes = []
        for j in xrange(0, n_partitions):
            train_partition = dataInfo.data + "_train_cv" + str(j)
            test_partition = dataInfo.data + "_test_cv" + str(j)
            ret = call_program("testSingleObjectiveTree.mKo", train_partition, test_partition, i)
            partial_results.append(float(ret.splitlines()[-1]))
            partial_nodes.append(get_num_nodes(ret))
        print "{:d}º variable: {}, nodes: {}".format(i, average_with_sd(partial_results), average_with_sd(partial_nodes)) 

def multi_objective_tree(dataInfo):
    matrix = [[] for x in xrange(0, dataInfo.number_objectives)]
    nodes = []
    for j in xrange(0, n_partitions):
        train_partition = dataInfo.data + "_train_cv" + str(j)
        test_partition = dataInfo.data + "_test_cv" + str(j)
        ret = call_program("testMultiObjectiveTree.mKo", train_partition, test_partition, dataInfo.first_variable, dataInfo.number_objectives)
        matrix = get_results_multi(dataInfo, ret, matrix)
        nodes.append(get_num_nodes(ret))
    for i in xrange(0, dataInfo.number_objectives):
        print "{:d}º variable: {}".format(i + dataInfo.first_variable, average_with_sd(matrix[i]))
    print "Nodes: {}".format(average_with_sd(nodes))

def single_objective_forest_fixed_trees_number(dataInfo):
    for i in xrange(dataInfo.first_variable, dataInfo.n_variables):
        partial_results = []
        for j in xrange(0, n_partitions):
            train_partition = dataInfo.data + "_train_cv" + str(j)
            test_partition = dataInfo.data + "_test_cv" + str(j)
            ret = call_program("testSingleObjectiveForest.mKo", train_partition, test_partition, i, 10)
            partial_results.append(float(ret.splitlines()[-1]))
        print "{:d}º variable: {}".format(i, average_with_sd(partial_results))

def multi_objective_forest_fixed_trees_number(dataInfo):
    matrix = [[] for x in xrange(0, dataInfo.number_objectives)]
    for j in xrange(0, n_partitions):
        train_partition = dataInfo.data + "_train_cv" + str(j)
        test_partition = dataInfo.data + "_test_cv" + str(j)
        ret = call_program("testMultiObjectiveForest.mKo", train_partition, test_partition, dataInfo.first_variable, dataInfo.number_objectives, 10)
        matrix = get_results_multi(dataInfo, ret, matrix)
    for i in xrange(0, dataInfo.number_objectives):
        print "{:d}º variable: {}".format(i + dataInfo.first_variable, average_with_sd(matrix[i]))

def single_objective_forest_variable_trees_number(dataInfo):
    for i in xrange(dataInfo.first_variable, dataInfo.n_variables):
        with open("{}_{:d}_variable.csv".format(dataInfo.data, i), 'wb') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for s in xrange(min_size_forest, max_size_forest, step_size_forest):
                partial_results = []
                for j in xrange(0, n_partitions):
                    train_partition = dataInfo.data + "_train_cv" + str(j)
                    test_partition = dataInfo.data + "_test_cv" + str(j)
                    ret = call_program("testSingleObjectiveForest.mKo", train_partition, test_partition, i, s)
                    partial_results.append(float(ret.splitlines()[-1]))
                average, sd = average_and_sd(partial_results)
                csvwriter.writerow([s, average, sd])

def multi_objective_forest_variable_trees_number(dataInfo):
    matrix = [[[] for y in xrange(min_size_forest, max_size_forest, step_size_forest)] for x in xrange(0, dataInfo.number_objectives)]
    for s in xrange(min_size_forest, max_size_forest, step_size_forest):
        for j in xrange(0, n_partitions):
            train_partition = dataInfo.data + "_train_cv" + str(j)
            test_partition = dataInfo.data + "_test_cv" + str(j)
            ret = call_program("testMultiObjectiveForest.mKo", train_partition, test_partition, dataInfo.first_variable, dataInfo.number_objectives, s)
            matrix = get_results_multi_variable_trees_number(dataInfo, ret, matrix, s)
    for i in xrange(0, dataInfo.number_objectives):
        with open("{}_{:d}_variable.csv".format(dataInfo.data, i), 'wb') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for s in xrange(min_size_forest, max_size_forest, step_size_forest):
                average, sd = average_and_sd(matrix[i][s])
                csvwriter.writerow([s, average, sd])

def main():

    print "\n***** GENERATE PARTITIONS *****"

    call_program("generatePartitions.mKo", data_wtr.data, n_partitions)
    call_program("generatePartitions.mKo", data_cosmo.data, n_partitions)

    print "\n***** WTR *****"

    print "\nSINGLE-OBJECTIVE TREE"
    # single_objective_tree(data_wtr)

    print "\nMULTI-OBJECTIVE TREE"
    # multi_objective_tree(data_wtr)

    print "\nSINGLE-OBJECTIVE FOREST WITH 10 TREES"
    # single_objective_forest_fixed_trees_number(data_wtr)

    print "\nMULTI-OBJECTIVE FOREST WITH 10 TREES"
    # multi_objective_forest_fixed_trees_number(data_wtr)

    print "\nSINGLE-OBJECTIVE FOREST WITH VARIABLE TREES"
    single_objective_forest_variable_trees_number(data_wtr)

    print "\nMULTI-OBJECTIVE FOREST WITH 10 TREES"
    multi_objective_forest_variable_trees_number(data_wtr)

    print "\n***** COSMO *****"

    print "\nSINGLE-OBJECTIVE TREE"
    # single_objective_tree(data_cosmo)

    print "\nMULTI-OBJECTIVE TREE"
    # multi_objective_tree(data_cosmo)

    print "\nSINGLE-OBJECTIVE FOREST WITH 10 TREES"
    # single_objective_forest_fixed_trees_number(data_cosmo)

    print "\nMULTI-OBJECTIVE FOREST WITH 10 TREES"
    # multi_objective_forest_fixed_trees_number(data_cosmo)

    print "\nSINGLE-OBJECTIVE FOREST WITH VARIABLE TREES"
    single_objective_forest_variable_trees_number(data_cosmo)

    print "\nMULTI-OBJECTIVE FOREST WITH 10 TREES"
    multi_objective_forest_variable_trees_number(data_cosmo)

    # Para cada ejemplo, Classify
    # ErrSecClass

main()