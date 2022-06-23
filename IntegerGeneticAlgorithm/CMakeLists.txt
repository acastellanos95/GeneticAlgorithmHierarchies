cmake_minimum_required(VERSION 3.22)
project(IntegerGeneticAlgorithm)

set(CMAKE_CXX_STANDARD 14)

add_executable(IntegerGeneticAlgorithm main.cpp Utils.h Individual.h lib/RandomGA.h)
