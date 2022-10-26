/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of ilaenv.
 * Author: CHEN Han
 * Create: 2022-06-26
 *******************************************************************************/

#pragma once

#include <kml_service.h>

/**
 * ////////////////////////////////_FuncTest///////////////////////////////////////
 * @Brief Returns problem-dependent parameters for the local environment
 * @param[in]		ispec		Specifies the parameter to be returned
 * as the value of ILAENV.
 * @param[in]		name		The name of the calling subroutine.
 * @param[in]		opts		The character options to the subroutine
 * NAME, concatenated into a single character string.
 * @param[in]		n1
 * @param[in]		n2
 * @param[in]		n3
 * @param[in]		n4
 */
int ilaenv(const int* ispec,
           const char* name,
           const char* opts,
           const int* n1,
           const int* n2,
           const int* n3,
           const int* n4);