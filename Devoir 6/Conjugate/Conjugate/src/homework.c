#include"fem.h"
// Théo Denis 27411800


#ifndef NOASSEMBLE


void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    int i,j;
    // assemblage du résidu pour la 1ère itération
    if(mySolver->iter==0){
        for (i = 0; i < nLoc; i++){
            mySolver->R[map[i]] += Bloc[i];
            for(j=0; j<nLoc; j++){
                mySolver->R[map[i]] -= Aloc[i*nLoc+j]*Uloc[j];
            }
        }
    }

    // calcul de S pour les itérations
    for (i = 0; i < nLoc; i++){
        for(j = 0; j< nLoc; j++){
            mySolver->S[map[i]] += Aloc[i*nLoc+j] * mySolver->D[map[j]];
        }
    }
}
#endif
#ifndef NOITERATE




double *femIterativeSolverEliminate(femIterativeSolver *mySolver){
    mySolver->iter++;
    double error = 0.0; 
    int i;

    for (i=0; i<mySolver->size; i++){
        error += mySolver->R[i] * mySolver->R[i];
    }

    // cas particulier de la première itération x0 = 0, d0=r0
    if (mySolver->iter==1){
        for (i=0; i<mySolver->size;i++){
            mySolver->X[i] = 0;
            mySolver->D[i] = mySolver->R[i];
        }
    } 
    else{
        
        double alphaDen = 0.0;
        for (i=0; i < mySolver->size; i++){
            alphaDen += (mySolver->S[i])*(mySolver->R[i]);
        }
        // alpha_k = r.r/Ad.r
        double alpha = error/alphaDen;

        double betaNum = 0.0; 
        for (i=0; i < mySolver->size; i++){
            mySolver->R[i] = mySolver->R[i] - (alpha * mySolver->S[i]);
            betaNum += mySolver->R[i]*mySolver->R[i];
        }

        // beta_k = r_k+1 . r_k+1/r_k.r_k
        double beta = betaNum/error;

        for (i=0; i< mySolver->size; i++){
            // xk+1 - xk = alpha dk
            mySolver->X[i] = alpha*mySolver->D[i]; 
            // d_k+1 = r_k+1 + beta_k * d_k
            mySolver->D[i] = mySolver->R[i] + beta*mySolver->D[i];
            mySolver->S[i] = 0.0;
        }
    }

    mySolver->error = sqrt(error);
    return(mySolver->X);
}
#endif


