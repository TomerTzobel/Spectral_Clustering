cmake_minimum_required(VERSION 3.17)
project(TomerPart C)

set(CMAKE_C_STANDARD 99)

add_executable(TomerPart main.c utils.c ddg.c eigengap.c wam.c utils.h wam.h ddg.h)