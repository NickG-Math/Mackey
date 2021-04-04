#pragma once

///	@file
///	@brief Contains open MP macros that enable and control multithreading


#if defined MACKEY_USE_OPEN_MP & defined _OPENMP
#include "omp.h"
///Macro that does openMP parallel for if MACKEY_USE_OPEN_MP and _OPENMP are both defined, and nothing otherwise
#define MACKEY_RUN_LOOP_PARALLEL _Pragma("omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic)")
#define MACKEY_RUN_BLOCK_SERIAL _Pragma("omp critical")
#pragma message ("mackey: openMP enabled!")
#else
///Macro that does ```openMP parallel for``` if \c MACKEY_USE_OPEN_MP and _OPENMP are both defined, and nothing otherwise
#pragma message ("mackey: openMP not enabled! To enable, define macro MACKEY_USE_OPEN_MP before including any header of this library and use the compile option for openMP")
#define MACKEY_RUN_LOOP_PARALLEL
#define MACKEY_RUN_BLOCK_SERIAL
#endif
