cmake_minimum_required(VERSION 3.22)
project(BinaryGeneticAlgorithm)

set(CMAKE_CXX_STANDARD 14)

add_executable(BinaryGeneticAlgorithm main.cpp Utils.h Individual.h lib/RandomGA.h)
