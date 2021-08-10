/*
 *  fem.h
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2015 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

static const char EulerName[] = "Euler";
static const char RKName[] = "Runge-Kutta";
static const char NoName[] = "Undefined scheme";
const char* femTemporalSchemeName(femTemporalSchemeType theScheme)
{
    if (theScheme == FEM_EULER) return EulerName;
    if (theScheme == FEM_RK)    return RKName;
    return NoName;
}


double femMin(double *x, int n) 
{
    double myMin = x[0];
    for (int i=1 ;i < n; i++) 
        if (x[i] < myMin) myMin = x[i];
    return myMin;
}

double femMax(double *x, int n) 
{
    double myMax = x[0];
    for (int i=1 ;i < n; i++) 
        if (x[i] > myMax) myMax = x[i];
    return myMax;
}

void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}
