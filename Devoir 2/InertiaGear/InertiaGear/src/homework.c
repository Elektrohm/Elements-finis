#include "fem.h"

// Théo Denis 27411800

#ifndef NORHOSTEEL
double inertiaGearSteelRho()
{

//
// Modifier la valeur pour avoir la masse volumique de l'acier [kg/m3]
// Une tolerance de 10% sur la valeur est admise
//

    double rho = 7800;
    return rho;
}
#endif

#ifndef NOINERTIA
double inertiaGearInertia(femMesh *theMesh, femIntegration *theRule, double rho)
{

    double I = 0;
    double xLoc[theRule->n];
    double yLoc[theRule->n];

    
    for(int k=0; k<theMesh->nElem; k++){
        
        // indices où aller chercher les points dans X et Y (sommets des triangles elem)
        int idx0 = theMesh->elem[k*3];
        int idx1 = theMesh->elem[k*3+1];
        int idx2 = theMesh->elem[k*3+2];

        // calcul du jacobien
        double J = fabs((theMesh->X[idx1]-theMesh->X[idx0])*(theMesh->Y[idx2]-theMesh->Y[idx0]) - (theMesh->Y[idx1]-theMesh->Y[idx0])*(theMesh->X[idx2]-theMesh->X[idx0]));

        // calcul des nouveaux sommets
        for(int i=0; i<theRule->n; i++){
            xLoc[i] = theMesh->X[idx0]*(1-theRule->xsi[i]-theRule->eta[i]) + theMesh->X[idx1]*theRule->xsi[i] + theMesh->X[idx2]*theRule->eta[i];
            yLoc[i] = theMesh->Y[idx0]*(1-theRule->xsi[i]-theRule->eta[i]) + theMesh->Y[idx1]*theRule->xsi[i] + theMesh->Y[idx2]*theRule->eta[i];
        }
        
        // calcul de l'intégral
        for(int j=0; j<theRule->n; j++){
            I += J*theRule->weight[j]*(xLoc[j]*xLoc[j] + yLoc[j]*yLoc[j]);
        }

    }
    
    return I*rho*0.1*1e-12;
}
#endif

#ifndef NOREAD
femMesh *inertiaGearMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash;
    
    // ouverture du fichier en mode lecture
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    // lecture de la première ligne, récupère le int à la position %d
    // initialise nNode de la structure à cette valeur
    ErrorScan(fscanf(file, "Number of nodes %d \n", &theMesh->nNode));

    // allocation des vecteurs X,Y
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        // initialisation des valeurs de X et Y
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i])); }

    // lecture de la 1ère ligne du tableau d'appartenance
    ErrorScan(fscanf(file, "Number of triangles %d \n", &theMesh->nElem));


    // hardcoding du nLocalNode triangle = 3
    theMesh->nLocalNode = 3;

    // allocation du vecteur *elem, on a nElem*nLocalNode qui sont chacun des int
    theMesh->elem = malloc(sizeof(int)*theMesh->nElem * theMesh->nLocalNode);

    for (i=0; i < theMesh->nElem; i++){
        // Ajout des éléments dans elem
        ErrorScan(fscanf(file,"%d : %d %d %d \n", &trash, &theMesh->elem[i*3], &theMesh->elem[i*3+1], &theMesh->elem[i*3+2]));
    }

    fclose(file);
    return theMesh;
}
#endif

#ifndef NOFREE
void inertiaGearMeshFree(femMesh *theMesh)
{
    // free les tableaux dans la structure en sens inverse de leur allocation
    free(theMesh->elem);
    free(theMesh->Y);
    free(theMesh->X);

    // free la structure
    free(theMesh);

}
#endif