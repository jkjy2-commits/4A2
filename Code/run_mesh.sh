#!/bin/bash
../Code/solver.x < input_bump_ni13.txt > screen_bump_ni13.txt &
../Code/solver.x < input_bump_ni25.txt > screen_bump_ni25.txt &
../Code/solver.x < input_bump_ni49.txt > screen_bump_ni49.txt &
../Code/solver.x < input_bump_ni97.txt > screen_bump_ni97.txt &
../Code/solver.x < input_bump_ni193.txt > screen_bump_ni193.txt &
../Code/solver.x < input_bump_ni385.txt > screen_bump_ni385.txt &
../Code/solver.x < input_bump_ni769.txt > screen_bump_ni769.txt &
wait
