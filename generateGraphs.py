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

def call_gnuplot (args_dictionary):
    args_string = ""
    for var,value in args_dictionary.iteritems():
        args_string += "{}='{}';".format(var,value)
    command = "gnuplot -e \"{}\" plotForestVariableSize.gpi".format(args_string)
    return subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()

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

def get_results_multi (data_info, ret, matrix):
    lines = ret.splitlines()[-data_info.number_objectives-1:-1]
    for i in xrange(0, data_info.number_objectives):
        m = re.match("\d+: ([\.\d]+)", lines[i])
        if m:
            matrix[i].append(float(m.group(1)))
    return matrix

def get_results_multi_variable_trees_number (data_info, ret, matrix, size_index):
    lines = ret.splitlines()[-data_info.number_objectives-1:-1]
    for i in xrange(0, data_info.number_objectives):
        m = re.match("\d+: ([\.\d]+)", lines[i])
        if m:
            matrix[i][size_index].append(float(m.group(1)))
    return matrix

def single_objective_tree(data_info):
    for i in xrange(data_info.first_variable, data_info.n_variables):
        partial_results = []
        partial_nodes = []
        for j in xrange(0, n_partitions):
            train_partition = data_info.data + "_train_cv" + str(j)
            test_partition = data_info.data + "_test_cv" + str(j)
            ret = call_program("testSingleObjectiveTree.mKo", train_partition, test_partition, i)
            partial_results.append(float(ret.splitlines()[-1]))
            partial_nodes.append(get_num_nodes(ret))
        print "{:d}º variable: {}, nodes: {}".format(i, average_with_sd(partial_results), average_with_sd(partial_nodes)) 

def multi_objective_tree(data_info):
    matrix = [[] for x in xrange(0, data_info.number_objectives)]
    nodes = []
    for j in xrange(0, n_partitions):
        train_partition = data_info.data + "_train_cv" + str(j)
        test_partition = data_info.data + "_test_cv" + str(j)
        ret = call_program("testMultiObjectiveTree.mKo", train_partition, test_partition, data_info.first_variable, data_info.number_objectives)
        matrix = get_results_multi(data_info, ret, matrix)
        nodes.append(get_num_nodes(ret))
    for i in xrange(0, data_info.number_objectives):
        print "{:d}º variable: {}".format(i + data_info.first_variable, average_with_sd(matrix[i]))
    print "Nodes: {}".format(average_with_sd(nodes))

def single_objective_forest_fixed_trees_number(data_info, size_forest = 10):
    for i in xrange(data_info.first_variable, data_info.n_variables):
        partial_results = []
        for j in xrange(0, n_partitions):
            train_partition = data_info.data + "_train_cv" + str(j)
            test_partition = data_info.data + "_test_cv" + str(j)
            ret = call_program("testSingleObjectiveForest.mKo", train_partition, test_partition, i, size_forest)
            partial_results.append(float(ret.splitlines()[-1]))
        print "{:d}º variable: {}".format(i, average_with_sd(partial_results))

def multi_objective_forest_fixed_trees_number(data_info, size_forest = 10):
    matrix = [[] for x in xrange(0, data_info.number_objectives)]
    for j in xrange(0, n_partitions):
        train_partition = data_info.data + "_train_cv" + str(j)
        test_partition = data_info.data + "_test_cv" + str(j)
        ret = call_program("testMultiObjectiveForest.mKo", train_partition, test_partition, data_info.first_variable, data_info.number_objectives, size_forest)
        matrix = get_results_multi(data_info, ret, matrix)
    for i in xrange(0, data_info.number_objectives):
        print "{:d}º variable: {}".format(i + data_info.first_variable, average_with_sd(matrix[i]))

def single_objective_forest_variable_trees_number(data_info):
    for i in xrange(data_info.first_variable, data_info.n_variables):
        with open("{}_{:d}_variable_single.csv".format(data_info.data, i), 'wb') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for s in xrange(min_size_forest, max_size_forest + 1, step_size_forest):
                partial_results = []
                for j in xrange(0, n_partitions):
                    train_partition = data_info.data + "_train_cv" + str(j)
                    test_partition = data_info.data + "_test_cv" + str(j)
                    ret = call_program("testSingleObjectiveForest.mKo", train_partition, test_partition, i, s)
                    partial_results.append(float(ret.splitlines()[-1]))
                average, sd = average_and_sd(partial_results)
                csvwriter.writerow([s, average, sd])

def multi_objective_forest_variable_trees_number(data_info):
    tree_size_index = list(enumerate(xrange(min_size_forest, max_size_forest + 1, step_size_forest)))
    matrix = [[[] for s in tree_size_index] for x in xrange(0, data_info.number_objectives)]
    for index_size, size in tree_size_index:
        for j in xrange(0, n_partitions):
            train_partition = data_info.data + "_train_cv" + str(j)
            test_partition = data_info.data + "_test_cv" + str(j)
            ret = call_program("testMultiObjectiveForest.mKo", train_partition, test_partition, data_info.first_variable, data_info.number_objectives, size)
            matrix = get_results_multi_variable_trees_number(data_info, ret, matrix, index_size)
    for i in xrange(0, data_info.number_objectives):
        with open("{}_{:d}_variable_multi.csv".format(data_info.data, i + data_info.first_variable), 'wb') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for index_size, size in tree_size_index:
                average, sd = average_and_sd(matrix[i][index_size])
                csvwriter.writerow([size, average, sd])

def generate_graphs(data_info):
    for i in xrange(data_info.first_variable, data_info.n_variables):
        call_gnuplot({
            "filename": "{}_{:d}_variable.png".format(data_info.data, i),
            "data": data_info.data,
            "variable": str(i),
        })

def main():

    print "\n***** GENERATE PARTITIONS *****"

    call_program("generatePartitions.mKo", data_wtr.data, n_partitions)
    call_program("generatePartitions.mKo", data_cosmo.data, n_partitions)

    print "\n***** WTR *****"

    print "\nSINGLE-OBJECTIVE TREE"
    single_objective_tree(data_wtr)

    print "\nMULTI-OBJECTIVE TREE"
    multi_objective_tree(data_wtr)

    print "\nSINGLE-OBJECTIVE FOREST WITH 10 TREES"
    single_objective_forest_fixed_trees_number(data_wtr)

    print "\nMULTI-OBJECTIVE FOREST WITH 10 TREES"
    multi_objective_forest_fixed_trees_number(data_wtr)

    print "\nSINGLE-OBJECTIVE FOREST WITH 100 TREES"
    single_objective_forest_fixed_trees_number(data_wtr, 100)

    print "\nMULTI-OBJECTIVE FOREST WITH 100 TREES"
    multi_objective_forest_fixed_trees_number(data_wtr, 100)

    print "\nSINGLE-OBJECTIVE FOREST WITH VARIABLE TREES"
    single_objective_forest_variable_trees_number(data_wtr)

    print "\nMULTI-OBJECTIVE FOREST WITH VARIABLE TREES"
    multi_objective_forest_variable_trees_number(data_wtr)

    print "\nGENERATE GRAPHS"
    generate_graphs(data_wtr)

    print "\n***** COSMO *****"

    print "\nSINGLE-OBJECTIVE TREE"
    single_objective_tree(data_cosmo)

    print "\nMULTI-OBJECTIVE TREE"
    multi_objective_tree(data_cosmo)

    print "\nSINGLE-OBJECTIVE FOREST WITH 10 TREES"
    single_objective_forest_fixed_trees_number(data_cosmo)

    print "\nMULTI-OBJECTIVE FOREST WITH 10 TREES"
    multi_objective_forest_fixed_trees_number(data_cosmo)

    print "\nSINGLE-OBJECTIVE FOREST WITH 100 TREES"
    single_objective_forest_fixed_trees_number(data_cosmo, 100)

    print "\nMULTI-OBJECTIVE FOREST WITH 100 TREES"
    multi_objective_forest_fixed_trees_number(data_cosmo, 100)

    print "\nSINGLE-OBJECTIVE FOREST WITH VARIABLE TREES"
    single_objective_forest_variable_trees_number(data_cosmo)

    print "\nMULTI-OBJECTIVE FOREST WITH VARIABLE TREES"
    multi_objective_forest_variable_trees_number(data_cosmo)

    print "\nGENERATE GRAPHS"
    generate_graphs(data_cosmo)

    # Para cada ejemplo, Classify
    # ErrSecClass

main()