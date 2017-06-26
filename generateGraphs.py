#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
from math import sqrt
import re

prog = "./nucleo/cons"
data_wtr = "wtr"
data_cosmo = "music2_z0.00_500_1e12"
n_partitions = 10
n_variables = 22
number_objectives = 5
first_variable= n_variables - number_objectives
min_size_forest = 50
max_size_forest = 100

def call_program (*args):
    return subprocess.Popen(prog + " " + " ".join(str(x) for x in args), shell=True, stdout=subprocess.PIPE).stdout.read()

def get_num_nodes (str):
    m = re.search("Num nodos: (\d+) ", str)
    if m:
        return int(m.group(1))

def average_with_sd (list):
    average = sum(list)/float(len(list))
    average_of_squares = sum(x**2 for x in list)/float(len(list))
    sd = sqrt(average_of_squares - average**2)
    return "{:.2f} ± {:.3f}".format(average, sd)

def get_results_multi (ret, matrix):
    lines = ret.splitlines()[-number_objectives-1:-1]
    for i in xrange(0, number_objectives):
        m = re.match("\d+: ([\.\d]+)", lines[i])
        if m:
            matrix[i].append(float(m.group(1)))
    return matrix

def single_objective_tree(data):
    for i in xrange(first_variable, n_variables):
        partial_results = []
        partial_nodes = []
        for j in xrange(0, n_partitions):
            train_partition = data + "_train_cv" + str(j)
            test_partition = data + "_test_cv" + str(j)
            ret = call_program("testSingleObjectiveTree.mKo", train_partition, test_partition, i)
            partial_results.append(float(ret.splitlines()[-1]))
            partial_nodes.append(get_num_nodes(ret))
        print "{:d}º variable: {}, nodes: {}".format(i, average_with_sd(partial_results), average_with_sd(partial_nodes)) 

def multi_objective_tree(data):
    matrix = [[] for x in xrange(0, number_objectives)]
    nodes = []
    for j in xrange(0, n_partitions):
        train_partition = data + "_train_cv" + str(j)
        test_partition = data + "_test_cv" + str(j)
        ret = call_program("testMultiObjectiveTree.mKo", train_partition, test_partition, first_variable, number_objectives)
        matrix = get_results_multi(ret, matrix)
        nodes.append(get_num_nodes(ret))
    for i in xrange(0, number_objectives):
        print "{:d}º variable: {}".format(i + first_variable, average_with_sd(matrix[i]))
    print "Nodes: {}".format(average_with_sd(nodes))

def single_objective_forest_fixed_trees_number(data):
    for i in xrange(first_variable, n_variables):
        partial_results = []
        partial_nodes = []
        for j in xrange(0, n_partitions):
            train_partition = data + "_train_cv" + str(j)
            test_partition = data + "_test_cv" + str(j)
            ret = call_program("testSingleObjectiveForest.mKo", train_partition, test_partition, i, 10)
            partial_results.append(float(ret.splitlines()[-1]))
        print "{:d}º variable: {}".format(i, average_with_sd(partial_results))

def multi_objective_forest_fixed_trees_number(data):
    matrix = [[] for x in xrange(0, number_objectives)]
    nodes = []
    for j in xrange(0, n_partitions):
        train_partition = data + "_train_cv" + str(j)
        test_partition = data + "_test_cv" + str(j)
        ret = call_program("testMultiObjectiveForest.mKo", train_partition, test_partition, first_variable, number_objectives, 10)
        print ret
        matrix = get_results_multi(ret, matrix)
    for i in xrange(0, number_objectives):
        print "{:d}º variable: {}".format(i + first_variable, average_with_sd(matrix[i]))

def main():
    call_program("generatePartitions.mKo", data_wtr, n_partitions)
    # call_program("generatePartitions.mKo", data_cosmo, n_partitions)

    print "\nSINGLE-OBJECTIVE TREE"
    # single_objective_tree(data_wtr)

    print "\nMULTI-OBJECTIVE TREE"
    # multi_objective_tree(data_wtr)

    print "\nSINGLE-OBJECTIVE FOREST WITH 10 TREES"
    # single_objective_forest_fixed_trees_number(data_wtr)

    print "\nMULTI-OBJECTIVE FOREST WITH 10 TREES"
    multi_objective_forest_fixed_trees_number(data_wtr)

main()

# # Para cada ejemplo, Classify
# # ErrSecClass
# echo "Single-objective forest"
# for size_forest in $(seq $min_size_forest 10 $max_size_forest)
# do
#     for j in $(seq 0 1 $n_partitions_zero_based)
#     do
#         train_partition=${data}_train_cv${j}
#         test_partition=${data}_test_cv${j}
#         partial_results[$j]=$($prog testSingleObjectiveForest.mKo $train_partition $test_partition $size_forest | tail -1)
#     done
#     results[$i]=$( average partial_results )
#     echo "With $size_forest trees: ${results[${i}]}"
# done

# echo "Multi-objective forest"
# for size_forest in $(seq $min_size_forest 10 $max_size_forest)
# do
#     for j in $(seq 0 1 $n_partitions_zero_based)
#     do
#         train_partition=${data}_train_cv${j}
#         test_partition=${data}_test_cv${j}
#         partial_results[$j]=$($prog testSingleObjectiveForest.mKo $train_partition $test_partition $size_forest | tail -1)
#     done
#     results[$i]=$( average partial_results )
#     echo "With $size_forest trees: ${results[${i}]}"
# done