#!/bin/bash

sed -iE 's/cConvergenceTest = false/cConvergenceTest = true/' modelconst.cpp
sed -iE 's/cDt = 1.0/cDt = 0.5/' modelconst.cpp
make efld_main
wait
./efld_main
wait
sed -iE 's/cConvergenceTest = true/cConvergenceTest = false/' modelconst.cpp
sed -iE 's/cDt = 0.5/cDt = 1.0/' modelconst.cpp
make efld_main
wait
./efld_main
wait
./plot.py
