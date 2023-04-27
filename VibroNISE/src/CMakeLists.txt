cmake_minimum_required(VERSION 3.3)
project(VibroNISE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    BasisMod.c
    BasisMod.h
    DirectRun.c
    GlobalsMod.h
    OverlapsMod.c
    OverlapsMod.h
    OverlapsMod_.c
    ParmsHandleMod.c
    ParmsHandleMod.h
    ParmsMod.h
    Pulse.c
    RandomMod.c
    RandomMod.h
    Spectra.c
    SubsMod.c
    SubsMod.h
    ToolsMod.c
    ToolsMod.h
    VibroNISE.c
    VibroNISE.h)

add_executable(VibroNISE ${SOURCE_FILES})