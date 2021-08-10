/*
 *  fem.h
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2015 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

# ifndef _FEM_H_
# define _FEM_H_
 
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# define Error(a)   femError(a,  __LINE__, __FILE__)
# define Warning(a) femWarning(a,  __LINE__, __FILE__)

typedef enum {FEM_EULER,FEM_RK} femTemporalSchemeType;


typedef struct {
    double epsilon;
    double beta;
    double zeta;
    int size;
    int sizePlot;
    double *X;
    double *U;
    double *x;
    double *u;
    double *uregime;
} femProblem;


femProblem *advdiffNew(double epsilon, double beta, double zeta, int size, int sizePlot);
void		    advdiffFree(femProblem *myProblem);
void 	      advdiffReset(femProblem *myProblem);
void        advdiffSteadySolution(femProblem *myProblem);
void        advdiffSolution(femProblem *myProblem, double time);
double 		  advdiffComputeTimeStep(femProblem *myProblem);
void		    advdiffComputeRightHandSide(femProblem *myProblem, double *U, double *F);
void	      advdiffUpdateEuler(femProblem *myProblem, double dt);
void	      advdiffUpdateRungeKutta(femProblem *myProblem, double dt);

const char* femTemporalSchemeName(femTemporalSchemeType theScheme);
double      femMin(double *x, int n);
double      femMax(double *x, int n);
void        femError(char *text, int line, char *file);
void        femWarning(char *text, int line, char *file);

# endif