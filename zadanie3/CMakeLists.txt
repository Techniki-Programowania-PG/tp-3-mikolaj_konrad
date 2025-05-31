cmake_minimum_required(VERSION 3.14)
project(Projekt)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(external/pybind11)

pybind11_add_module(sygnaly zadanie3.cpp)
