/*******************************************************************************
 * Copyright (c) Huawei Technologies Co., Ltd. 2020-2021. All rights reserved.
 * Description: Part of KML library
 * Author: KML
 * Create: 2020
 ******************************************************************************/

#ifndef KML_SOLVER_DEFS_H_INCLUDED
#define KML_SOLVER_DEFS_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

/** Error code returned by KML solver functions.
 *
 * @sa{KmlIssGetErrorText, KmlIssClearError}
 */
enum KML_SOLVER_ERROR {
    NO_ERROR = 0,
    NONZERO_INDEXING = 1,
    NO_DIAGONAL_ELEMENT = 2,
    ZERO_DIAGONAL_ELEMENT = 3,
    NO_MEMORY = 4,
    BAD_HANDLE = 5,
    DATA_SIZE_ERROR = 6,
    OUT_OF_EXECUTION_ORDER = 7,
    INCONSISTENT_PERMUTATION = 8,
    KML_SOLVER_INVALID_ARGUMENT = 9,
    NOT_IMPLEMENTED = -12,
};

enum KML_SOLVER_PARAM {
    FILL_IN = 0,
    PERM = 1,
    /* For ISS */
    THRESHOLD = 3,
    MAX_ITERATIONS_COUNT = 4,
    RESTART_PARAM = 5,
    ITERATION_COUNT = 6,
    TOLERANCE = 7
};

#if defined(__cplusplus)
}
#endif

#endif
