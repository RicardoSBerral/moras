#!/bin/bash

prog=./nucleo/cons
data=wtr
n_partitions=10
n_partitions_zero_based=$(( $n_partitions - 1 ))
last_variable=21
number_objectives=5
first_variable=$(( $last_variable - $number_objectives + 1 ))
min_size_forest=50
max_size_forest=100

average () {
    local -n array=$1
    local comma_array=$(IFS=,; echo "${array[*]}")
    echo $(echo $comma_array | awk '{
        s = 0;
        split($0,array,",");
        size = length(array)
        for (j=1;j<=size;j++) {
            s += array[j];
        }
    } END{ print s/size; }')
}

$prog generatePartitions.mKo $data $n_partitions

echo "Single-objective tree"
for i in $(seq $first_variable 1 $last_variable)
do
    for j in $(seq 0 1 $n_partitions_zero_based)
    do
        train_partition=${data}_train_cv${j}
        test_partition=${data}_test_cv${j}
        partial_results[$j]=$($prog testSingleObjectiveTree.mKo $train_partition $test_partition $i | tail -1)
    done
    results[$i]=$( average partial_results )
    echo "${i}º variable: ${results[${i}]}"
done

echo "Multi-objective tree"
for j in $(seq 0 1 $n_partitions_zero_based)
do
    train_partition=${data}_train_cv${j}
    test_partition=${data}_test_cv${j}
    partial_results[$j]=$($prog testMultiObjectiveTree.mKo $train_partition $test_partition $first_variable $number_objectives | tail -1)
done
results[$i]=$( average partial_results )
echo "${i}º variable: ${results[${i}]}"

for i in $(seq $first_variable 1 $last_variable)
do
done

# Para cada ejemplo, Classify
# ErrSecClass
echo "Single-objective forest"
for size_forest in $(seq $min_size_forest 10 $max_size_forest)
do
    for j in $(seq 0 1 $n_partitions_zero_based)
    do
        train_partition=${data}_train_cv${j}
        test_partition=${data}_test_cv${j}
        partial_results[$j]=$($prog testSingleObjectiveForest.mKo $train_partition $test_partition $size_forest | tail -1)
    done
    results[$i]=$( average partial_results )
    echo "With $size_forest trees: ${results[${i}]}"
done

echo "Multi-objective forest"
for size_forest in $(seq $min_size_forest 10 $max_size_forest)
do
    for j in $(seq 0 1 $n_partitions_zero_based)
    do
        train_partition=${data}_train_cv${j}
        test_partition=${data}_test_cv${j}
        partial_results[$j]=$($prog testSingleObjectiveForest.mKo $train_partition $test_partition $size_forest | tail -1)
    done
    results[$i]=$( average partial_results )
    echo "With $size_forest trees: ${results[${i}]}"
done