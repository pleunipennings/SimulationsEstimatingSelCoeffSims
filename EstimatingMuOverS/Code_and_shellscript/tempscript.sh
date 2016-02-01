#!/bin/bash
mu=0.001
output_every_Xgen=5
numgen_inN=0.101
start_output=0.001
cost=1
for seed in 1
do
echo "
$seed
$mu
$cost
$output_every_Xgen
$numgen_inN
$start_output
" | ./Code_and_shellscript/HIVevolution_HIV1site >../Data/Data_s_300_T_10_cost_1.txt
done
