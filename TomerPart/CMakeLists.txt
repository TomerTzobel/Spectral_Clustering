cmake_minimum_required(VERSION 3.17)
project(TomerPart C)

set(CMAKE_C_STANDARD 99)

#link_directories(D:/Programs/mingw/mingw64/lib/)
include_directories(C:/Users/adamk/AppData/Local/Programs/Python/Python38/libs/include)
#link_libraries(C:/Users/adamk/AppData/Local/Programs/Python/Python38/libs/python38.lib)

add_executable(TomerPart include spkmeans.c spkmeans.c spkmeans.h spkmeansmodule.c)