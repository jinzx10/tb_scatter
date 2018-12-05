#!/bin/bash

sed -iE 's/cConvergenceTest = true/cConvergenceTest = false/' efld_main.cpp
make efld_main
wait
./efld_main

