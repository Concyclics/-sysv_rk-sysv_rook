/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of ilaenv.
 * Author: CHEN Han
 * Create: 2022-06-26
 *******************************************************************************/

#include "ilaenv.h"

int ilaenv(const int* ispec,

           const char* name,

           const char* opts,

           const int* n1,

           const int* n2,

           const int* n3,

           const int* n4) {
    const int n = *n1;

    if (KmlLsame(*name, 'S')) {
        if (n <= 4000) {
            return 8;

        } else if (n <= 7000) {
            return 32;

        } else {
            return 64;

        } 
    } else if (KmlLsame(*name, 'D')) {
        if (n <= 4000) {
            return 8;

        } else if (n <= 10000) {
            return 16;

        } else {
            return 64;
        }
    } else if (KmlLsame(*name, 'C')) {
        if (n <= 4000) {
            return 8;

        } else if (n <= 10000) {
            return 32;

        } else {
            return 64;
        }
    } else if (KmlLsame(*name, 'Z')) {
        if (n <= 1000) {
            return 8;

        } else if (n <= 15000) {
            return 32;

        } else {
            return 64;
        }

        return 0;
    }
}
