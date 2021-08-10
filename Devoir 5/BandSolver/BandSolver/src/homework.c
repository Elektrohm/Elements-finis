
#include"fem.h"


#ifndef NORENUMBER 

// besoin d'une global pour récupérer x et y
femMesh * MeshGlobal;

// fonction de comparaison sur x pour qsort
int sortX(const void* node0, const void* node1){
    int * n0 = (int *) node0;
    int * n1 = (int *) node1;
    
    // needs to be *n to get the first value of the pseudo pointer
    double diag = MeshGlobal->X[*n0] - MeshGlobal->X[*n1]; 
    return (diag > 0) - (diag < 0);
}

// fonction de comparaison sur x pour qsort
int sortY(const void* node0, const void* node1){
    int * n0 = (int *) node0;
    int * n1 = (int *) node1;
    
    // needs to be *n to get the first value of the pseudo pointer
    double diag = MeshGlobal->Y[*n0] - MeshGlobal->Y[*n1]; 
    return (diag > 0) - (diag < 0);
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;
    // il me cassait les pieds et mettre ça là à corriger :| 
    int *runner;
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nNode; i++) 
                theMesh->number[i] = i;
            break;

        case FEM_XNUM : 
            // allocation structure temporaire 
            runner = (int *) malloc(sizeof(int)*theMesh->nNode);
            for(i = 0 ; i < theMesh->nNode; i++){
                runner[i] = i;
            }

            // besoin du mesh global pour comparer
            MeshGlobal = theMesh; 
            qsort(runner, theMesh->nNode, sizeof(int), sortX);

            // retranscription des résultats
            for(i = 0; i < theMesh->nNode; i++){
                theMesh->number[runner[i]] = i;
            }
            free(runner);
            break;


        case FEM_YNUM : 
            // allocation structure temporaire 
            runner = malloc(sizeof(int)*theMesh->nNode);
            for(i =0 ; i < theMesh->nNode; i++){
                runner[i] = i;
            }

            // besoin du mesh global pour comparer
            MeshGlobal = theMesh; 
            qsort(runner, theMesh->nNode, sizeof(int), sortY);

            // retranscription des résultats
            for(i =0; i < theMesh->nNode; i++){
                theMesh->number[runner[i]] = i;
            }
            free(runner);
            break;

        default : Error("Unexpected renumbering option"); }
}

#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh)
{
    int myBand = 0;
    int nLoc = theMesh->nLocalNode;
    int i, j, k, min, max, corners[nLoc];

    // boucle sur tous les éléments et chacun de leur noeud
    for(i=0; i < theMesh->nElem; i++){
        for(j=0; j < nLoc; j++){
            corners[j] = theMesh->number[theMesh->elem[i*nLoc + j]];
        }
        // corners donne les 3 points du triangles ou 4 coin du carré
        min = corners[0];
        max = corners[0];

        // boucle sur les coins au sein d'un élément
        for(k=1; k<nLoc; k++){
            min = fmin(corners[k], min);
            max = fmax(corners[k], max);
        }
        if(myBand < (max-min)){
            // piège à con, faire le +1 (exemple : myband = 5+1, myband > 6 => nope or intéressant)
            myBand = max-min;
        }
    }
    return (++myBand);
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        for(j = 0; j < nLoc; j++) {
            if (map[j] >= map[i]) myBandSystem->A[map[i]][map[j]] += Aloc[i*nLoc+j]; }
    myBandSystem->B[map[i]] += Bloc[i]; }
}


#endif
#ifndef NOBANDELIMINATE


double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-8 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        jend = fmin(k+band, size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) {
                A[i][j] = A[i][j] - A[k][j] * factor;
            }
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i+band, size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }  


    return(myBand->B);
}


#endif