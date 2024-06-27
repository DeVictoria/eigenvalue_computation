cmake_minimum_required(VERSION 3.24)
project(backlog_itmo_c_eigenvalues_DeVictoria C)

set(CMAKE_C_STANDARD 99)

include_directories(.)

add_executable(backlog_itmo_c_eigenvalues_DeVictoria
        main.c
        return_codes.h)
