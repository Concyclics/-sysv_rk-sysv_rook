/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of ilaenv.
 * Author: CHEN Han
 * Create: 2022-06-26
 *******************************************************************************/

#pragma once

int ilaenv(const int* ispec,
           const char* name,
           const char* opts,
           const int* n1,
           const int* n2,
           const int* n3,
           const int* n4);