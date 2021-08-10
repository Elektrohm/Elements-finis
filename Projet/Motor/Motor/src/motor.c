#include "motor.h"

//
// ========= Projet � r�aliser ===================
//

// macro pour accélérer certaines opérations
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

// variable globale pour optimiser le programme
int firstIteration = 1; // 1 = True, 0 = False
double d = 0; // épaisseur dans le calcul de l'intégrale
femIntegration* theRule;
femDiscrete* theSpace;
femEdges* theEdges;



void motorAdaptMesh(motor *theMotor, double delta)
{
    motorMesh *theMesh = theMotor->mesh;
    int startInd = 0;

    double x, y;
    for (int i = 0; i < theMesh->nNode; ++i)
    {
        if (theMotor->movingNodes[i] == 1)
        {
            x = theMesh->X[i] * cos(delta) - theMesh->Y[i] * sin(delta);
            y = theMesh->X[i] * sin(delta) + theMesh->Y[i] * cos(delta);
            theMesh->X[i] = x;
            theMesh->Y[i] = y;
        }
    }
    theMotor->theta += delta;

    // calcule �  partir d'où commence les éléments de air_gap
    for(int i=0; i<11; i++){
        startInd += theMesh->nElemDomain[i];
    }
    // pointeur vers les éléments de air_gap
    int *elemband = &(theMesh->elem[3*startInd]);

    int n1, n2, cas, val;
    // boucle sur tous les éléments de air_gap
    for(int i=0; i<theMesh->nElemDomain[11]; i++){
        int mv1 = theMotor->movingNodes[elemband[3*i]];
        int mv2 = theMotor->movingNodes[elemband[3*i+1]];
        int mv3 = theMotor->movingNodes[elemband[3*i+2]];

        // récupère les 2noeuds avec la même valeur de moving node
        if(mv1==mv2){
            n1 = elemband[3*i];
            n2 = elemband[3*i+1];
            cas = 0;
            val = mv1;
        }

        if(mv1==mv3){
            n1 = elemband[3*i];
            n2 = elemband[3*i+2];
            cas = 1;
            val = mv1;
        }

        else{
            n1 = elemband[3*i+1];
            n2 = elemband[3*i+2];
            cas = 2;
            val = mv2;
        } 

        int candidat = mv1;
        double best = 0.0;
    
        // boucle sur tous les autres éléments de air_gap
        for(int j=0; j<theMesh->nElemDomain[11];j++){
            
            // boucle sur les 3 noeuds de l'éléments 
            for(int k=0; k<3; k++){
                int n3 = elemband[3*j+k];

                // si noeud a une valeur moving node différente des deux noeuds 
                if((theMotor->movingNodes[n3]!=val) & (n3!=n1) & (n3!=n2)){

                    double X1 = theMesh->X[n1];
                    double Y1 = theMesh->Y[n1];
                    double X2 = theMesh->X[n2];
                    double Y2 = theMesh->Y[n2];
                    double X3 = theMesh->X[n3];
                    double Y3 = theMesh->Y[n3];

                    // calcul du jacobien 
                    double jac = fabs((X2-X1)*(Y3-Y1) - (Y2-Y1)*(X3-X1));
                    
                    // calcul des côtés du triangles
                    double c1 = sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1));
                    double c2 = sqrt((X3-X1)*(X3-X1)+(Y3-Y1)*(Y3-Y1));
                    double c3 = sqrt((X3-X2)*(X3-X2)+(Y3-Y2)*(Y3-Y2));

                    // calcul du périmètre
                    double perim = c1 + c2 + c3;

                    // calcul de la métrique de comparaison
                    double ratio = jac/(perim*perim*perim);

                    // si meilleur ratio
                    if(ratio > best){
                        // conserve le noeud comme candidat
                        candidat = n3;
                        best = ratio;
                    }

                }
            }
                
        } 

        // remplace la valeur du mesh par le candidat
        switch (cas)
        {
        case 0:
            elemband[3*i+2] = candidat;
            break;
        
        case 1:
            elemband[3*i+1] = candidat;
            break;

        case 2:
            elemband[3*i] = candidat;
            break;

        default:
            printf("there was an error");
            break;
        }
        
    }
}

double motorComputeCouple(motor *theMotor)
{
    motorMesh *theMesh = theMotor->mesh;

    int startInd = 0;

    // d = rmax - rmin
    double rmin=__DBL_MAX__;
    double rmax=0.0;
    double mu0 = theMotor->mu[10];

    for(int i=0; i<10; i++){
        startInd += theMesh->nElemDomain[i];
    }

    // pointeur vers les éléments de rotor_gap
    int *elemrotgap = &(theMesh->elem[3*startInd]);

    double x[3], y[3], phi[3], dphidx[3], dphidy[3], dphidxsi[3], dphideta[3];
    double dphidr, dphidtheta;
    double Couple = 0.0;

    // boucle sur les éléments
    for(int i =0; i<theMesh->nElemDomain[10]; i++){
        
        // récupère valeur de x et y
        for(int p=0; p<3; p++){
            x[p] = theMesh->X[elemrotgap[3*i+p]];
            y[p] = theMesh->Y[elemrotgap[3*i+p]];
        }

        // calcul uniquement �  la 1ère itération d (constant)
        if(firstIteration==1){
            // récupère rayon max et min pour le calcul de d
            for(int p=0; p<3; p++){
                double r = sqrt(x[p]*x[p]+y[p]*y[p]);
                rmin = MIN(rmin,r);
                rmax = MAX(rmax,r);
                d = rmax-rmin;
            }
            firstIteration = 0;
        }
        
        // changement de variable
        for(int j=0; j<theRule->n; j++){
            // récupère paramètre intégration élément parents
            double xsi = theRule->xsi[j];
            double eta = theRule->eta[j];
            double weight = theRule->weight[j];
 
            // récupère les fonctions de formes
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            // initialise variable du changement de variable
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            double xloc = 0.0;
            double yloc = 0.0;

            // jacobien du changement de variable
            double jac = 0.0;

            for(int k=0; k<theSpace->n; k++){
                xloc += x[k]*phi[k];
                yloc += y[k]*phi[k];

                dxdxsi += x[k]*dphidxsi[k];
                dxdeta += x[k]*dphideta[k];

                dydxsi += y[k]*dphidxsi[k];
                dydeta += y[k]*dphideta[k];
            }

            jac = fabs(dxdxsi*dydeta - dydxsi*dxdeta);

            // applique le changement de variable
            for(int l=0; l<theSpace->n; l++){
                dphidx[l] = (1.0/jac)*(dydeta*dphidxsi[l]-dydxsi*dphideta[l]);
                dphidy[l] = (1.0/jac)*(dxdxsi*dphideta[l]-dxdeta*dphidxsi[l]);
            }

            // calcul r, cos et sin pour passage en polaire
            double r = sqrt(xloc*xloc+yloc*yloc);
            double cos = xloc/r;
            double sin = yloc/r;
            
            double dadr = 0.0;
            double dadtheta = 0.0;


            // chain rule 
            for(int m=0; m<theSpace->n; m++){
                dphidtheta = dphidy[m]*r*cos - dphidx[m]*r*sin;
                dphidr = dphidx[m]*cos + dphidy[m]*sin;
                dadr += theMotor->a[elemrotgap[3*i+m]]*dphidr;
                dadtheta += theMotor->a[elemrotgap[3*i+m]]*dphidtheta;
            }

            Couple += jac*weight*dadr*dadtheta;
        }
    } 
    Couple *= (-theMotor->L/(mu0*d));
    return Couple;
}


void motorComputeCurrent(motor *theMotor)
{   
    // déphasage �  appliquer pour avoir le meilleur switch
    double offset = 12*M_PI/180;
    double theta = fabs(theMotor->theta+offset);
    double js = 0.0;

    // mise �  0 des bobines + récupération de js
    for(int i=0; i<12; i++){
        js = MAX(theMotor->js[i],js);
        theMotor->js[i] = 0.0;
    }

    // récupère le reste de la position du moteur modulo pi/2
    double position = fmod(theta,M_PI_2);

    // de 0 �  30° allumage des bobines 3 et 4
    if(position >=0 && position<M_PI/6){
        theMotor->js[3] = js;
        theMotor->js[4] = -js;
    }
    // de 30 �  60° allumage des bobines 1 et 2
    if(position >=M_PI/6 && position<M_PI/3){
        theMotor->js[1] = js;
        theMotor->js[2] = -js;
    }
    // de 60 �  90° allumage des bobines 5 et 6
    if(position >=M_PI/3 && position<M_PI_2){
        theMotor->js[5] = js;
        theMotor->js[6] = -js;
    }

    return;
}



void myFemDiffusionMeshLocal(motor *theMotor, const int iElem, int *map,  int *ctr, double *x, double *y, double *u, int *number,int* dirichlet)
{
    // permet de récupérer les valeurs X, Y local de l'élément 
    int j,nLocal = theMotor->mesh->nLocalNode;
    
    for (j=0; j < nLocal; ++j) {
        map[j] = theMotor->mesh->elem[iElem*nLocal+j];
        x[j]   = theMotor->mesh->X[map[j]];
        y[j]   = theMotor->mesh->Y[map[j]]; 
        u[j]   = theMotor->a[map[j]];
        ctr[j] = dirichlet[map[j]];
        map[j] = number[map[j]];
        }   
}

femMesh *myMeshConverter(motorMesh *motorMesh, int *number)
{
    // permet de recréer un femMesh �  partir de motorMesh
    femMesh *theFemMesh = malloc(sizeof(femMesh));
    theFemMesh->elem = motorMesh->elem;
    theFemMesh->X = motorMesh->X;
    theFemMesh->Y = motorMesh->Y;
    theFemMesh->nElem = motorMesh->nElem;
    theFemMesh->nNode = motorMesh->nNode;
    theFemMesh->nLocalNode = motorMesh->nLocalNode;
    theFemMesh->number = number;
    return theFemMesh;
}


void motorComputeMagneticPotential(motor *theMotor)
{   
    motorMesh* theMotorMesh = theMotor->mesh;      
    
    int *number = malloc(theMotorMesh->nNode * sizeof(int));

    for (int i = 0; i < theMotorMesh->nNode; i++) 
        number[i] = i; 
    femMesh *theFemMesh = myMeshConverter(theMotorMesh, number);
    femMeshRenumber(theFemMesh,FEM_YNUM);

    // malloc les structures utiles
    if(firstIteration==1){
        theEdges = femEdgesCreate(theFemMesh);
        theRule = femIntegrationCreate(3, FEM_TRIANGLE);
        theSpace = femDiscreteCreate(3, FEM_TRIANGLE);
    }
	
    // initialise le solveur
	int *dirichlet = malloc(theMotor->size * sizeof(int));
	for (int i = 0; i < theMotor->size; i++) {
        theMotor->a[i] = 0;
		dirichlet[i] = 0;
	}
	int band = femMeshComputeBand(theFemMesh);
	femSolver *theSolver = femSolverBandCreate(theMotor->size, theMotorMesh->nLocalNode, band);

	for (int i = 0; i < theEdges->nEdge; i++) {
		if (theEdges->edges[i].elem[1] < 0) {
			dirichlet[theEdges->edges[i].node[0]] = 1;
			dirichlet[theEdges->edges[i].node[1]] = 1;
		}
	}

    // initialise les structures utiles �  la résolution
	double Xloc[4], Yloc[4], Uloc[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
	int map[4], ctr[4];
    
	double **A = theSolver->local->A;
	double *Aloc = theSolver->local->A[0];
	double *Bloc = theSolver->local->B;
	for (int iElem = 0; iElem < theMotorMesh->nElem; iElem++){
        // récupère js et mu du domaine de l'élément
		double js = theMotor->js[theMotorMesh->domain[iElem]];
        double mu = theMotor->mu[theMotorMesh->domain[iElem]];

        // initialise composante �  0
		for (int i = 0; i < theSpace->n; i++){
            Bloc[i] = 0;
        }
		for (int i = 0; i < (theSpace->n) * (theSpace->n); i++){
            Aloc[i] = 0;
        }

        // récupère valeur local de l'élément
		myFemDiffusionMeshLocal(theMotor, iElem, map, ctr, Xloc, Yloc, Uloc, number, dirichlet);  
        
        // boucle d'intégration
        for (int iInteg = 0; iInteg < theRule->n; iInteg++) {
            // paramètre d'intégration
			double weight = theRule->weight[iInteg];
			double xsi = theRule->xsi[iInteg];
			double eta = theRule->eta[iInteg];

            // fonction de forme
			femDiscretePhi2(theSpace, xsi, eta, phi);
			femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            // changement de variable
			double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
			for (int i = 0; i < theSpace->n; i++) {
				dxdxsi += Xloc[i] * dphidxsi[i];
				dxdeta += Xloc[i] * dphideta[i];
				dydxsi += Yloc[i] * dphidxsi[i];
				dydeta += Yloc[i] * dphideta[i];
			}
			double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
			for (int i = 0; i < theSpace->n; i++) {
				dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
				dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
			}

            // assemblage de la matrice
			for (int i = 0; i < theSpace->n; i++) {
				for (int j = 0; j < theSpace->n; j++) {
					A[i][j] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight * 1.0 / mu;
				}
				Bloc[i] += phi[i] * jac * js * weight;
			}
		}
        // application des contraintes
		for (int i = 0; i < theSpace->n; i++)
			if (ctr[i] == 1) femFullSystemConstrain(theSolver->local, i, 0);
        // assemblage du solveur
		femSolverAssemble(theSolver, Aloc, Bloc, Uloc, map, theSpace->n);
	}
    // élimination gaussienne
	double *soluce = femSolverEliminate(theSolver);
	for (int i = 0; i < theMotor->size; i++)
		theMotor->a[i] += soluce[number[i]];

    // free le solver et le femMesh (change �  l'itération suivante)
    free(dirichlet);
    femSolverFree(theSolver);
    free(number);
    free(theFemMesh);

    return;
}


void motorFree(motor *theMotor)
{
    // free variable général de motor.c et theMotor
    femIntegrationFree(theRule);
    femDiscreteFree(theSpace);
    femEdgesFree(theEdges);
    free(theMotor->mu);
    free(theMotor->js);
    free(theMotor->movingNodes);
    free(theMotor->a);
    free(theMotor);
    return;
}

