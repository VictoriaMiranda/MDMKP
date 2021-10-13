/*************************************************************************************/
/****  0. Header files, data structures, and global varialbes  ***********************/
/*************************************************************************************/ 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <math.h>
#include <ctype.h>
using namespace std;

typedef struct{
	int *X;     // array binario que representa la solucion 
	int *S;     // objetos seleccionados
	int *NS;    // objetos no seleccionados
	int IN, ON; // numero de objetos in y out
	int *IW;    // peso en la mochila (por cada restriccion)
	long int f; // funcion objetivo
} Solution; 

typedef struct Neighbor{
	int f;      // funcion objetivo de la solucion candidata
    int type;   // 1: candidato de vecindario N1 (bit-flip) 2: candidato de vecindario N2 o N3 (swap)
    int IO;     // 0: candidato generado con add-bit 1: candidato generado con drop-bit
    int k;      // indice en S o NS de bit que cambio a 0 o 1 
    int x;      // bit que cambio a 0 o 1
    int x1;     // indice en S de bit que cambia a 0 (swap)
    int y1;     // indice en NS de bit que cambia a 1 (swap) 
} Neighbor;

string File_Name; 
double time_limit, time_to_target, AvgTime; // tiempo maximo en realizar ambas fases, time_to_target no se usa(?)
double time_one_run, starting_time;  

int f, f_best;     // mejor f.o. de un vecindario, mejor f.o. actual de la fase
int *P;            // beneficios
int *B;            // capacidades maximas y minimas
int *ES;
float *E;          // medidas de eficiencia
int **Matrix;
int totalCore;
int media;
float core;
int seed;
int **R;           // matriz de pesos
int N, M, Q, BKV;       // nº objetos, restricciones de capacidad minima, y maxima
Solution SC;       // solucion actual
Solution OS;
Solution S_BEST;   // mejor solucion actual
Solution G_BEST;   // mejor solucion de la ejecucion
int *W1, *W2, *W3; // valores de i^{gamma} en vectores hash
int *TL;
long int L;        // largo de las listas
int alpha  = 5000; // profundidad de TS
int lanbda;  // peso de la penalizacion en funcion de evaluacion
int Itt;
int resta;
float pct_L;

bool Find(int* Array, int number) {
    int i;
    for (i = 0; i < L; i++) {
        if (Array[i] == number) return true;
    }
    return false;
}

/***********************************************************************************/
/***********************          1. Initializing           ************************/
/***********************************************************************************/ 
//1.1 Inputing 
void Initializing(float min, float max) {
    int i,j, x1, x2; 
    totalCore = 0; 

	ifstream FIC1;
	ifstream FIC2; 
    FIC1.open("benchmark/250/"+File_Name);
    if (FIC1.fail()) {
        //std::cout << "### Erreur open, File_Name " << File_Name << endl;
        exit(0);
    } 
	FIC1 >> N >> M >> Q >> BKV;
 
	P = new int [N];
    E = new float [N];
    ES = new int [N];
	B = new int [M+Q]; 
	 
	R = new int *[M+Q];
	for(i=0;i<M+Q;i++) 
		R[i] = new int [N];
		 		  
    while (!FIC1.eof()) {
        for(i=0;i<N;i++) FIC1 >> P[i];
		for(i=0;i<M;i++) 
		    for(j=0;j<N;j++) FIC1 >> R[i][j];
		for(j=0;j<M;j++) FIC1 >> B[j];   
	       
		for(i=M;i<M+Q;i++) 
		    for(j=0;j<N;j++) FIC1 >> R[i][j];  
		for(j=M;j<M+Q;j++) FIC1 >> B[j];     
    }	       
    FIC1.close();
    
    FIC2.open("eficiencias/"+File_Name);
    if (FIC2.fail()) {
        //std::cout << "### Erreur open, File_Name " << "eficiencias/"+File_Name << endl;
        exit(0);
    } 
    for (j = 0; j < N; j++) {
        FIC2 >> E[j];
        if (E[j] >= min && E[j] <= max) {
            ES[totalCore] = j;
            totalCore++;
        }
    }
    //std::cout << "Elementos core: " << totalCore << endl;
    FIC2.close();
    L = pct_L*N;
}

void WriteMatrix(string instancia) {
    int i, j;
    std::ofstream outfile ("matrix_hib-"+instancia);
    for (i = 0; i < totalCore; i++) {
        for (j = 0; j < totalCore; j++) {
            outfile << Matrix[ES[i]][ES[j]] << " ";
        }
        outfile << endl;
    }
    outfile << media;
    outfile.close();
} 

void GetMedia() {
    int i, j, suma;
    suma = 0;
    for (i=0; i < totalCore; i++) {
        for (j=0; j<totalCore; j++) {
            suma += Matrix[ES[i]][ES[j]];
        }
    }
    media = suma/(totalCore*totalCore);
}

void AssignMemory() {
	SC.X      = new int [N];
	SC.S      = new int [N];
	SC.NS     = new int [N]; 
	SC.IW     = new int [M+Q];

	S_BEST.X  = new int [N];
    S_BEST.S  = new int [N];
	S_BEST.NS = new int [N];
	S_BEST.IW = new int [M+Q]; 

    OS.X      = new int [N];
	OS.S      = new int [N];
	OS.NS     = new int [N]; 
	OS.IW     = new int [M+Q]; 

	G_BEST.X  = new int [N];
    G_BEST.S  = new int [N];
    G_BEST.NS = new int [N];
    G_BEST.IW = new int [M+Q];
     
    TL = new int [L];

}

/**************************************************************************/
/************     2. Initial Solution (Random and Opposite)    ************/
/**************************************************************************/ 
void pseudoRandomSol(Solution &S) {
	int i, j, k;
	int r;     // indice random de los candidatos a elegir, objeto random elegido, numero de candidatos restantes
	int c1, c2;           // indice dinamico de S y NS (empieza en 0, avanza a medida que se agregan objetos)


	for(i=0;i<N;i++) S.S[i] = -1; 
	for(i=0;i<N;i++) S.NS[i] = -1;
    for(i=0;i<M+Q;i++) S.IW[i] = 0;

	S.IN = 0;   // no hay objetos dentro
	S.ON = 0;   // todos los objetos estan fuera
	S.f = 0;  
	c1 = 0;
	c2 = 0; 

    for (j = 0; j < N; j++) {
        if (E[j] < 1) {
            S.X[j] = 1;
            for (i = 0; i < M + Q; i++) S.IW[i] += R[i][j]; 
            S.IN++;
            S.S[c1] = j;
            S.f += P[j]; 
            c1++;
        }
        else if (E[j] > 1) {
            S.X[j] = 0;
            S.ON++;
            S.NS[c2] = j;
            c2++;
        }
        else {
            r = rand() % 2;
            S.X[j] = r;

            if (r == 1) {
                for (i = 0; i < M + Q; i++) S.IW[i] += R[i][j]; 
                S.IN ++;
                S.S[c1] = j; 
                S.f += P[j]; 
                c1++;
            }
            else {
                S.ON ++;
                S.NS[c2] = j;
                c2++;
            }
        }
    }
}

/*****************************************************************************/
/************************    5. Certification of the solution   **************/
/*****************************************************************************/ 
int proof(Solution &S) {
    int i, j; 
    int count, p; 
    int assign; 

    for(i=0; i<M+Q; i++) S.IW[i] = 0;   // peso en 0
    for(i=0;i<M+Q;i++) {
        for(j=0;j<N;j++) {
            S.IW[i] += S.X[j]*R[i][j]; // actualizar peso con la solucion encontrada 
        }
    }
   
    assign = 0; 
    for(i=0;i<M;i++) {
        if(S.IW[i]>B[i]) {
            assign = 1;
            //printf("<= errer ! \n");
        } // violacion restriccion capacidad maxima
    } 
    for(i=M;i<M+Q;i++) {
        if (S.IW[i]<B[i]) {
            assign = 1;
            //printf(">= errer ! \n");
        } // violacion restriccion capacidad minima
    }

    count = 0; // cantidad de objetos en solucion encontrada (no se usa realmente (?))
    p=0;       // beneficio de solucion encontrada
    for(i=0;i<N;i++) if(S.X[i]==1) count++;
    for(i=0;i<N;i++) {
        if(S.X[i]==1) {
            p+=P[i];
        }
    }
    return assign; // 0: factible, 1: no factible
}

/*****************************************************************************/
/*****************     4. Solution-based Tabu Search  procedures *************/
/*****************************************************************************/ 
//-----------------------------------------------------------------------------
// 4.1 Build initial information for the tabu lists
//-----------------------------------------------------------------------------
void Build_Information() {
    int i;
    for(i=0;i<L;i++) {
        TL[i] = -1;
	} 
}
//-----------------------------------------------------------------------------
// 4.2 Tabu Search with the union neighborhood of N1 and N2
//-----------------------------------------------------------------------------
void Tabu_Search(Solution &S, float min, float max) {
    int i, j, j1, k, k0, x, y, v, u, l;
    int select, swap;                     // solucion candidata seleccionada en la iteracion
    int non_improve = 0 ;                 // condicion de parada
    long int tabu_best_fc, best_fc, f_c;  // mejor f.e. tabu, f.e. no tabu, y f.o. actual
    int num_tabu_best, num_best;          // numero de soluciones candidatas tabu y no-tabu
    int n_add, n_drop;
    Neighbor best[50];                    // mejores soluciones candidatas no tabu (= f.e.)
    Neighbor tabu_best[50];               // mejores soluciones candidatas tabu (= f.e.)
	double current_time, starting_time, to_best; 
	long int sum_penalty;                 // penalizacion en f.e. P(s)
	long int delt;                        // -lambda*P(s)
	int flag;                             // si solucion es factible, flag = 1 
	 
    f_best = -999999;
    f = S.f;  
    Build_Information();
	   
    // mejor solucion = solucion actual
	S_BEST.IN = S.IN;
	S_BEST.ON = S.ON; 
	S_BEST.f= S.f; 
	for(i=0;i<N;i++)   S_BEST.S[i]  = S.S[i];
	for(i=0;i<N;i++)   S_BEST.NS[i] = S.NS[i];
	for(i=0;i<N;i++)   S_BEST.X[i]  = S.X[i]; 
	for(i=0;i<M+Q;i++) S_BEST.IW[i] = S.IW[i]; 

    int n_tl = 0;
	starting_time = clock(); 
	flag  = 0; 

    // mientras exista una mejor en las ultima 2*5000 iteraciones o, aun no se haya encontrado una solucion factible
    while( non_improve < alpha || flag == 0 ) {
        //std::cout << S.f << endl;
        tabu_best_fc = -999999999;  
        best_fc = -999999999;
        num_tabu_best = 0;  
        num_best = 0;

		//1). the add neighborhood
		for(k = 0; k < S.ON; k ++) { // agrego objeto con indice k en NS
            j = S.NS[k];
			sum_penalty = 0; 
			for(i=0;i<M;i++)   if(S.IW[i] + R[i][j] > B[i]) sum_penalty += (S.IW[i] + R[i][j] - B[i]);   
			for(i=M;i<M+Q;i++) if(S.IW[i] + R[i][j] < B[i]) sum_penalty += (B[i] - S.IW[i] - R[i][j]);  
		    delt = -1.0*lanbda*sum_penalty;   // calculo su penalizacion
		    
            // como varia f.o.
			f_c = S.f + P[j];  
		
            if (Find(TL,j) == true) {  // si es tabu
                // si la f.e. es mejor que la mejor f.e. tabu ...
                if( f_c + delt > tabu_best_fc ) { 
                    tabu_best[ 0 ].x = j; 
                    tabu_best[ 0 ].type = 1;       // la agrego a las mejores soluciones tabu ...
                    tabu_best[ 0 ].IO = 0; 
                    tabu_best[ 0 ].k = k; 
                    tabu_best[ 0 ].f = f_c; 
                    tabu_best_fc = f_c + delt;
                    num_tabu_best = 1;             // "elimino" las otras
                }
                // si la f.e. es igual que la mejor f.e. tabu ...
                else if ( f_c + delt == tabu_best_fc && num_tabu_best < 50 ) {
                    tabu_best[ num_tabu_best ].x    = j;    
					tabu_best[ num_tabu_best ].type = 1;   // la agrego a las mejores soluciones tabu ...
					tabu_best[ num_tabu_best ].IO = 0; 
					tabu_best[ num_tabu_best ].k = k; 
					tabu_best[ num_tabu_best ].f = f_c; 
                    num_tabu_best++;   
                }                                              
            } 

			else {                          // si no es tabu
                // si la f.e. es mejor que la mejor f.e. no tabu ...
                if( f_c + delt > best_fc ) {
                    best[ 0 ].x = j; 
                    best[ 0 ].type = 1;                   // la agrego a las mejores soluciones no tabu ...
                    best[ 0 ].IO = 0;
                    best[ 0 ].k = k;
                    best[ 0 ].f = f_c ;
                    best_fc  = f_c + delt; 
                    num_best = 1 ;                        // "elimino" las otras
                }
                // si la f.e. es igual que la mejor f.e. no tabu ...
                else if( f_c + delt == best_fc && num_best < 50 ) {
                    best[ num_best ].x = j;  
                    best[ num_best ].type = 1;            // la agrego a las mejores soluciones no tabu ...
                    best[ num_best ].IO = 0; 
                    best[ num_best ].k = k; 
                    best[ num_best ].f = f_c;  
                    num_best++ ;  
                }
            } 	
		} 
        
        //2) the drop neighborhood  
        for( k = 0; k < S.IN; k ++) {     // saco objeto con indice k en S
            j = S.S[k];
		    sum_penalty = 0;  
			for(i=0;i<M;i++)   if(S.IW[i] - R[i][j] > B[i]) sum_penalty += (S.IW[i] - R[i][j] - B[i]); 
			for(i=M;i<M+Q;i++) if(S.IW[i] - R[i][j] < B[i]) sum_penalty += (B[i] - S.IW[i] + R[i][j]); 
		    delt = -1.0*lanbda*sum_penalty; 
			
			f_c = S.f - P[j] ; 
			
            if (Find(TL,j) == true) { // si es tabu
                // si la f.e. es mejor que la mejor f.e. tabu ...
				if( f_c + delt > tabu_best_fc ) {
                    tabu_best[ 0 ].x = j ; 
                    tabu_best[ 0 ].type = 1 ; 
                    tabu_best[ 0 ].IO = 1; 
                    tabu_best[ 0 ].k = k; 
                    tabu_best[ 0 ].f = f_c; 
                    tabu_best_fc = f_c + delt; 
                    num_tabu_best = 1 ;
                }
                // si la f.e. es igual que la mejor f.e. tabu ...
                else if( f_c + delt == tabu_best_fc && num_tabu_best < 50 ) {
                    tabu_best[ num_tabu_best ].x   = j;    
					tabu_best[ num_tabu_best ].type = 1;   
					tabu_best[ num_tabu_best ].IO = 1; 
					tabu_best[ num_tabu_best ].k = k; 
					tabu_best[ num_tabu_best ].f = f_c; 
                    num_tabu_best++ ;   
                }                               
            } 
			else {                                                 // si no es tabu
                // si la f.e. es mejor que la mejor f.e. no tabu ...
                if( f_c + delt > best_fc ) {
                    best[ 0 ].x = j; 
                    best[ 0 ].type = 1;
                    best[ 0 ].IO = 1; 
                    best[ 0 ].k = k;
                    best[ 0 ].f = f_c ;
                    best_fc  = f_c + delt; 
                    num_best = 1 ; 
                }
                // si la f.e. es igual que la mejor f.e. no tabu ...
                else if( f_c + delt == best_fc && num_best < 50 ) {
                    best[ num_best ].x   = j ; 
                    best[ num_best ].type = 1 ;
                    best[ num_best ].IO = 1; 
                    best[ num_best ].k = k; 
                    best[ num_best ].f = f_c; 
                    num_best++ ;   
                }
            }              	
		}  
        //3) Evaluating the neighborhood N2 based candidate list  
        for(x = 0; x < S.IN; x++) {         // saco x
            for(y = 0; y < S.ON; y++) {     // agrego y
                sum_penalty = 0;   
			    for(i=0;i<M;i++) if(S.IW[i]+(R[i][S.NS[y]]-R[i][S.S[x]]) > B[i])   sum_penalty += (S.IW[i]+(R[i][S.NS[y]]-R[i][S.S[x]]) - B[i]) ; 
			    for(i=M;i<M+Q;i++) if(S.IW[i]+(R[i][S.NS[y]]-R[i][S.S[x]]) < B[i]) sum_penalty += (B[i]- S.IW[i] - (R[i][S.NS[y]]-R[i][S.S[x]])); 
			    delt = -1.0*lanbda*sum_penalty; 
			    
                // f.o. del swap
				f_c = S.f + (P[S.NS[y]]-P[S.S[x]]);
				
                if ((Find(TL,S.S[x]) == true) || (Find(TL,S.NS[y]) == true)) {
                    // f.e actual > mejor f.e. (tabu)
                    if ( f_c + delt > tabu_best_fc ) {
                        tabu_best[ 0 ].x1 = x ; 
                        tabu_best[ 0 ].y1 = y ; 
                        tabu_best[ 0 ].type = 2 ; 
                        tabu_best[ 0 ].IO = -1;
                        tabu_best[ 0 ].f = f_c; 
                        tabu_best_fc = f_c + delt; 
                        num_tabu_best = 1 ;
                    }
                    // f.e actual = mejor f.e. (tabu)
                    else if ( f_c + delt == tabu_best_fc && num_tabu_best < 50 ) {
                        tabu_best[ num_tabu_best ].x1   = x;    
                        tabu_best[ num_tabu_best ].y1   = y;  
						tabu_best[ num_tabu_best ].type = 2;   
						tabu_best[ num_tabu_best ].IO = -1; 
						tabu_best[ num_tabu_best ].f = f_c;  
                        num_tabu_best++ ;
                    }   
                    // si es un swap que aumenta mi funcion objetivo
                    if ( (f_c >= S.f) && (E[S.S[x]] >= min && E[S.S[x]] <= max) && (E[S.NS[y]] >= min && E[S.NS[y]] <= max) ) {
                        Matrix[S.S[x]][S.NS[y]] += 1;
                    }                   
                } 
			    else {                                                  // si no es tabu 
                    // f.e actual > mejor f.e. (no tabu)
                    if (f_c + delt > best_fc ) {
                        best[0].x1 = x; 
                        best[0].y1 = y;
                        best[0].type = 2;
                        best[0].IO = -1; 
                        best[0].f  = f_c; 
                        best_fc  = f_c + delt; 
                        num_best = 1; 
                    }
                    // f.e actual = mejor f.e. (no tabu)
                    else if( f_c + delt == best_fc && num_best < 50 ) {
                        best[num_best].x1   = x ; 
                        best[num_best].y1   = y ;
                        best[num_best].type = 2 ;
                        best[num_best].IO = -1 ; 
                        best[num_best].f = f_c ; 
                        num_best++ ;  
                    }
                    if ( (f_c >= S.f) && (E[S.S[x]] >= min && E[S.S[x]] <= max) && (E[S.NS[y]] >= min && E[S.NS[y]] <= max) ) {
                        Matrix[S.S[x]][S.NS[y]] += 1;
                    } 
                }    		               
			}
        }	
        //4) moves
	    if ((num_best == 0) && (num_tabu_best > 0)) { // criterio de aspiracion (si no existe solucion candidata no tabu)
            select = rand() % num_tabu_best ;   // selecciono aleatoriamente una mejor solucion tabu
            f = tabu_best[select].f ;

            if(tabu_best[select].type==2) {     // si es de N2 (swap)
                u = tabu_best[select].x1; 
                v = tabu_best[select].y1;

                // el mejor swap no aumenta en lista tabu
                if ( (f >= S.f) && (E[S.S[u]] >= min && E[S.S[u]] <= max) && (E[S.NS[v]] >= min && E[S.NS[v]] <= max) ) Matrix[S.S[u]][S.NS[v]] -=1;

                //actualizo listas
                TL[n_tl%L] = S.S[u];
                n_tl++; 
                TL[n_tl%L] = S.NS[v];
                n_tl++; 

                // actualizo peso en la mochila  
                for(i=0;i<M+Q;i++) {
                    S.IW[i] += (R[i][S.NS[v]] - R[i][S.S[u]]);   
				} 
                // actualizo solucion actual
                S.X[S.S[u]]  = 1 - S.X[S.S[u]]; 
                S.X[S.NS[v]] = 1 - S.X[S.NS[v]];

                // actualizo objetos seleccionados (S) y no seleccionados (NS)    
                swap    = S.S[u]; 
                S.S[u]  = S.NS[v];
                S.NS[v] = swap; 
                // actualizo f.o. actual  
                S.f = f; 
			}
			    	
			else if (tabu_best[select].type == 1) {    // si es de N1 (bit-flip)
                j  = tabu_best[ select ].x; 
			    k0 = tabu_best[ select ].k;  
			    	   
				if(tabu_best[ select ].IO == 1) {      // si un objeto sale
                    TL[n_tl%L] = j;
                    n_tl++;

                    // actualizo peso en la mochila   
                    for(i=0;i<M+Q;i++) {
                       	S.IW[i] -= R[i][j];  
				    } 
					      
                    // actualizo solucion actual
                    S.X[j] = 1 - S.X[j]; 

                    // actualizo objetos seleccionados
                    for(i=k0; i<S.IN; i++) S.S[i] = S.S[i+1];
                    S.IN -- ;

                    // actualizo objetos no seleccionados    
                    S.NS[S.ON] = j ;
                    S.ON ++ ; 

                    // actualizo f.o. actual
                    S.f = f;  
				}
				
                else if(tabu_best[ select ].IO == 0) {    // si un objeto entra

                    TL[n_tl%L] = l;
                    n_tl++;

                    // actualizo peso  
                    for(i=0;i<M+Q;i++) {
                        S.IW[i] += R[i][j];  
					} 

                    // actualizo solucion actual
                    S.X[j] = 1 - S.X[j]; 

                    // actualizo objetos no seleccionados
                    for(i=k0;i<S.ON;i++) S.NS[i] = S.NS[i+1];
                    S.ON--; 

                    // actualizo objetos seleccionados
                    S.S[S.IN] = j; 
                    S.IN++;

                    //actualizo f.o.
                    S.f = f;
                }	
			}
        } 
        else if (num_best > 0){                  // si no es tabu
            select = rand() % num_best;     // selecciono aleatoriamente una mejor solucion no tabu 
            f = best[select].f ; 

            if(best[select].type == 2) {   // si es de N2
                u = best[ select ].x1; 
                v = best[ select ].y1;  

                // el mejor swap no aumenta en lista tabu
                if ( (f >= S.f) && (E[S.S[u]] >= min && E[S.S[u]] <= max) && (E[S.NS[v]] >= min && E[S.NS[v]] <= max)) Matrix[S.S[u]][S.NS[v]] -=1;

                TL[n_tl%L] = S.S[u];
                n_tl++; 
                TL[n_tl%L] = S.NS[v];
                n_tl++; 
                       
                for(i=0;i<M+Q;i++) {
                    S.IW[i] += (R[i][S.NS[v]] - R[i][S.S[u]]);  
				} 
					   
                S.X[S.S[u]]  = 1 - S.X[S.S[u]]; 
                S.X[S.NS[v]] = 1 - S.X[S.NS[v]];
                       
                swap    = S.S[u] ; 
                S.S[u]  = S.NS[v];
                S.NS[v] = swap ;  
					    
                S.f = f;  
			}
			else if( best[select].type == 1 ) {  // si es de N1
			    j  = best[ select ].x; 
			    k0 = best[ select ].k; 
			    	    
			    if(best[ select ].IO == 1) {	 // si sale un objeto
                    TL[n_tl%L] = j;
                    n_tl++;

                    for(i=0;i<M+Q;i++) {
                        S.IW[i] -= R[i][j];  
					} 
					      
                    S.X[j] = 1 - S.X[j]; 
                    for(i=k0;i<S.IN;i++) S.S[i] = S.S[i+1];
                    S.IN -- ;
                          
                    S.NS[S.ON] = j ;
                    S.ON ++ ; 
                    S.f = f;   		 
				}
                
                else if(best[ select ].IO == 0) {  // si entra un objeto
                    TL[n_tl%L] = j;
                    n_tl++;

                    for(i=0;i<M+Q;i++) {
                       	S.IW[i] += R[i][j];  
					} 
                    S.X[j] = 1 - S.X[j]; 
                    for(i=k0;i<S.ON;i++) S.NS[i] = S.NS[i+1];
                    S.ON -- ; 
                    S.S[S.IN] = j; 
                    S.IN ++ ;
                    S.f = f;
                }	
			}
        } 
        //5. print
		current_time = (double) (1.0*(clock()-starting_time)/CLOCKS_PER_SEC); 

		sum_penalty = 0; 
        // calculo penalizacion de solucion actual
		for(i=0;i<M;i++)   if(S.IW[i] > B[i])  sum_penalty += (S.IW[i] - B[i]);
		for(i=M;i<M+Q;i++) if(S.IW[i] < B[i])  sum_penalty += (B[i] - S.IW[i]);

        // si solucion es factible 
		if(sum_penalty == 0) {
            flag = 1; 
        }

        // si f.o. de solucion actual es >= que mejor f.o. hasta el momento Y es factible 
        if( S.f >= f_best && sum_penalty == 0 ) { 
            // si f.o. de solucion actual es MEJOR que mejor f.o. hasta el momento 
            if( S.f > f_best  ) {
                // mejor f.o. = f.o. actual
                f_best = S.f;
                // mejor solucion = solucion actual
                for( i = 0 ; i < N ; i ++ )
                    S_BEST.X[ i ] = S.X[ i ];
                for( i = 0 ; i < N ; i ++ )
                    S_BEST.S[ i ] = S.S[ i ];
                for( i = 0 ; i < N ; i ++ )
                    S_BEST.NS[ i ] = S.NS[ i ]; 
                for( i = 0 ; i < M+Q ; i ++ )
                    S_BEST.IW[ i ] = S.IW[ i ];  
                S_BEST.IN = S.IN;
                S_BEST.ON = S.ON; 
		    	S_BEST.f = S.f; 

                non_improve = 0; // hay mejora :)
                to_best = (double) (1.0*(clock()-starting_time)/CLOCKS_PER_SEC); 
            }  
            // si f.o. actual es igual a la mejor: no hay mejora
            else if ( S.f == f_best )  non_improve ++ ;
        }  
        // si f.o. actual es peor que la mejor hasta el momento
        else non_improve ++ ;  
    }
       
    // termina fase1, actualizo solucion actual con la mejor encontrada en la fase
    for(i=0;i<N;i++)  S.X[i]  = S_BEST.X[i]; 
	for(i=0;i<N;i++)  S.S[i]  = S_BEST.S[i];
	for(i=0;i<N;i++)  S.NS[i] = S_BEST.NS[i]; 
	for(i=0;i<M+Q;i++)  S.IW[i] = S_BEST.IW[i]; 
	S.f  = S_BEST.f;
    S.IN = S_BEST.IN;
    S.ON = S_BEST.ON; 
}

//-------------------------------------------------------------------------------------------------
// 4.3 Tabu Search with the constrained swap neighborhood N3 and the extended evaluation function
//-------------------------------------------------------------------------------------------------
void Swap_tabu_search(Solution &S, double X, float min, float max) {
    int i, j, k, k0, x, y, v, u, iter ;
    int K, flag;
    int count; 
    int assign, Ncount, Ncount1, Ncount2;
    double pho, eta;
	int num, num1; 
    int n_drop, n_add;
	long int sum_penalty; 
	long int delt, tabu_delt, non_delt; 
    int select, swap; 
    int non_improve = 0 ;  //the stop condition of TS
	long int tabu_best_fc, best_fc, f_c ;
	int num_tabu_best, num_best;  //the number of tabu neighbors and non-tabu neighbors
    Neighbor best[ 50 ];
    Neighbor tabu_best[ 50 ]; 
    float factor = 1.0;

    Neighbor worst;
    long int worst_fc;
	double current_time, starting_time; 
	 
    Build_Information(); 
	 
	f = S.f; 
	f_best = S.f; 
	S_BEST.IN = S.IN;
	S_BEST.ON = S.ON; 
	S_BEST.f = S.f; 
	for(i=0;i<N;i++) S_BEST.S[i]  = S.S[i];
	for(i=0;i<N;i++) S_BEST.NS[i] = S.NS[i];
	for(i=0;i<N;i++) S_BEST.X[i]  = S.X[i]; 
	for(i=0;i<M+Q;i++) S_BEST.IW[i] = S.IW[i];  
	 
    int n_tl = 0;
	starting_time = clock(); 
    iter = 0; 
	current_time = (double) (1.0*(clock()-starting_time)/CLOCKS_PER_SEC); 
    while( current_time < X ) {
        flag = 0;
        tabu_best_fc = -99999999 ; 
        best_fc = -99999999 ;
        worst_fc = 99999999;
        num_tabu_best = 0 ; 
        num_best = 0 ;
		num = 0; 
		num1 = 0; 
		//2) Evaluating the neighborhood N3   
		for(x = 0; x < S.IN; x++) {
            for(y = 0; y < S.ON; y++) {
                f_c = f + (P[S.NS[y]]-P[S.S[x]]); 
			    if (f_c < factor*f_best) continue;

                else if (Matrix[S.S[x]][S.NS[y]] > media) {
                    Matrix[S.S[x]][S.NS[y]] -= resta;
                    continue;
                }
                
				sum_penalty = 0; 
			    for(i=0;i<M;i++)   if(S.IW[i]+(R[i][S.NS[y]]-R[i][S.S[x]]) > B[i])  sum_penalty += (S.IW[i] + R[i][S.NS[y]] - R[i][S.S[x]] - B[i]);
			    for(i=M;i<M+Q;i++) if(S.IW[i]+(R[i][S.NS[y]]-R[i][S.S[x]]) < B[i])  sum_penalty += (B[i] - S.IW[i] - R[i][S.NS[y]] + R[i][S.S[x]]);
			    delt = -1.0*lanbda*sum_penalty; 
			    
                if ( f_c + delt < worst_fc) {
                    worst.x1 = x;
                    worst.y1 = y ;  
                    worst.f = f_c ; 
                    worst_fc = f_c + delt ;
                }
				
				if( (Find(TL,S.S[x]) == true) || (Find(TL,S.NS[y]) == true) ) {
                    num ++; 
				    if( f_c + delt > tabu_best_fc ) {
                        tabu_best[ 0 ].x1 = x ; 
                        tabu_best[ 0 ].y1 = y ; 
                        tabu_best[ 0 ].type = 2 ; 
                        tabu_best[ 0 ].f = f_c ; 
                        tabu_best_fc = f_c + delt ; 
                        tabu_delt = delt; 
                        num_tabu_best = 1 ;
                    }
                    else if( f_c + delt == tabu_best_fc && num_tabu_best < 50 ) {
                        tabu_best[ num_tabu_best ].x1   = x;    
                        tabu_best[ num_tabu_best ].y1   = y;  
						tabu_best[ num_tabu_best ].type = 2;   
						tabu_best[ num_tabu_best ].f = f_c; 
                        num_tabu_best++ ;   
                    }  
				} 
			    else {   
                    num1++;  
					if( f_c + delt > best_fc ) {
                        best[ 0 ].x1 = x; 
                        best[ 0 ].y1 = y;
                        best[ 0 ].type = 2;
                        best[ 0 ].f = f_c; 
                        best_fc  = f_c + delt; 
                        non_delt = delt; 
                        num_best = 1 ; 
                    }
                    else if( f_c + delt == best_fc && num_best < 50 ) {
                        best[ num_best ].x1   = x ; 
                        best[ num_best ].y1   = y ;
                        best[ num_best ].type = 2 ;
                        best[ num_best ].f = f_c; 
                        num_best++ ;  
                    }
                }          
			}
        }

        if (non_improve == Itt) {
            for (j=0; j<N; j++) OS.X[j] = S.X[j];
            for (j=0; j<N; j++) OS.S[j] = S.S[j];
            for (j=0; j<N; j++) OS.NS[j] = S.NS[j];
            OS.f = worst.f;

            u = worst.x1; 
            v = worst.y1;  
                       
            for(i=0;i<M+Q;i++) OS.IW[i] = S.IW[i] + (R[i][S.NS[v]] - R[i][S.S[u]]);  
			
            OS.X[S.S[u]]  = 1 - S.X[S.S[u]]; 
            OS.X[S.NS[v]] = 1 - S.X[S.NS[v]];
                       
            OS.S[u]  = S.NS[v];
            OS.NS[v] = S.S[u] ; 
        }	 
        //4) moves
	    if((num_best == 0) && (num_tabu_best > 0) ) {
            select = rand() % num_tabu_best ; 
			f = tabu_best[select].f ;
            u = tabu_best[ select ].x1; 
            v = tabu_best[ select ].y1;  

            TL[n_tl%L] = S.S[u];
            n_tl++; 
            TL[n_tl%L] = S.NS[v];
            n_tl++; 
                       
           for(i=0;i<M+Q;i++) {
               S.IW[i] += (R[i][S.NS[v]] - R[i][S.S[u]]);   
			} 
            S.X[S.S[u]] = 1 - S.X[S.S[u]]; 
            S.X[S.NS[v]] = 1 - S.X[S.NS[v]];
                       
            swap    = S.S[u] ; 
            S.S[u]  = S.NS[v];
            S.NS[v] = swap   ;   
            S.f = f;   
        } 
           
        else if (num_best > 0) {
            select = rand() % num_best ;  
            f = best[ select ].f ; 
            u = best[ select ].x1; 
            v = best[ select ].y1;  

            TL[n_tl%L] = S.S[u];
            n_tl++; 
            TL[n_tl%L] = S.NS[v];
            n_tl++; 
                       
            for(i=0;i<M+Q;i++) S.IW[i] += (R[i][S.NS[v]] - R[i][S.S[u]]);  
					   
            S.X[S.S[u]]  = 1 - S.X[S.S[u]]; 
            S.X[S.NS[v]] = 1 - S.X[S.NS[v]];
                       
            swap    = S.S[u] ; 
            S.S[u]  = S.NS[v];
            S.NS[v] = swap   ;  
					    
            S.f = f;  
        } 
        else {
            factor -= 0.01;
        }
        //5. print
        iter ++ ;
		current_time = (double) (1.0*(clock()-starting_time)/CLOCKS_PER_SEC); 
        //std::cout << best_fc << endl;
        if( f >= f_best  ) {
           	sum_penalty = 0; 
			for(i=0;i<M;i++)   if(S.IW[i] > B[i])  sum_penalty += (S.IW[i] - B[i]);
			for(i=M;i<M+Q;i++) if(S.IW[i] < B[i])  sum_penalty += (B[i] - S.IW[i]);
			//printf("sum =%d\n",sum_penalty); 
            if( f > f_best && sum_penalty == 0 ) {
                f_best = f;
                  //std::cout << "Mejora: " << iter << endl; 
                for( i = 0 ; i < N ; i ++ )
                    S_BEST.X[ i ] = S.X[ i ] ;
                for( i = 0 ; i < N ; i ++ )
                    S_BEST.S[ i ] = S.S[ i ] ;
                for( i = 0 ; i < N ; i ++ )
                    S_BEST.NS[ i ] = S.NS[ i ] ; 
                for( i = 0 ; i < M+Q ; i ++ )
                    S_BEST.IW[ i ] = S.IW[ i ] ;  
                S_BEST.IN = S.IN;
                S_BEST.ON = S.ON; 
				S_BEST.f = f; 
                 // printf("\n swap_move :  %8d     %d     %d     %d     %lf      %lf", iter, f, f_best, non_improve, pho, current_time);
                  //non_improve = 0 ;
				time_one_run = current_time; 
            }  
            else if ( f == f_best )  non_improve++;

            else if ( sum_penalty > 0) non_improve++;
        }  
        else non_improve ++ ;

        if (non_improve > Itt) {
            for (i=0; i<M+Q; i++) S.IW[i] = OS.IW[i];
            for (j=0; j<N; j++) S.X[j] = OS.X[j];
            for (j=0; j<N; j++) S.S[j] = OS.S[j];
            for (j=0; j<N; j++) S.NS[j] = OS.NS[j];
            S.f = OS.f;
            //std::cout << "non_improve " << non_improve << endl;
            int c1 = 0;
            int c2 = 0;

            for (u = 0 ; u < S.IN ; u++ ) {
                if (E[S.S[u]] >= min && E[S.S[u]] <= max) {
                    for (v = c1; v < S.ON; v++) {
                        if (E[S.NS[v]] >= min && E[S.NS[v]] <= max) {
                            c2 = 1;
                            break;
                        }
                    }
                }
                if (c2 == 1) {
                    S.f = S.f + (P[S.NS[v]]-P[S.S[u]]); 
                    for(i=0;i<M+Q;i++) 
                        S.IW[i] += (R[i][S.NS[v]] - R[i][S.S[u]]);  
                    
                    S.X[S.S[u]]  = 1 - S.X[S.S[u]]; 
                    S.X[S.NS[v]] = 1 - S.X[S.NS[v]];
                            
                    swap    = S.S[u] ; 
                    S.S[u]  = S.NS[v];
                    S.NS[v] = swap   ;  

                    c2 = 0;
                    c1 = v + 1;
                    f = S.f;
                }
            }

            sum_penalty = 0; 
            // calculo penalizacion de solucion actual
            for(i=0;i<M;i++)   if(S.IW[i] > B[i])  sum_penalty += (S.IW[i] - B[i]);
            for(i=M;i<M+Q;i++) if(S.IW[i] < B[i])  sum_penalty += (B[i] - S.IW[i]);

            if( S.f >= f_best && sum_penalty == 0 ) { 
                // si f.o. de solucion actual es MEJOR que mejor f.o. hasta el momento 
                if( S.f > f_best  ) {
                    // mejor f.o. = f.o. actual
                    //std::cout << "mejora: " << iter << endl;
                    f_best = S.f;
                    // mejor solucion = solucion actual
                    for( i = 0 ; i < N ; i ++ )
                        S_BEST.X[ i ] = S.X[ i ];
                    for( i = 0 ; i < N ; i ++ )
                        S_BEST.S[ i ] = S.S[ i ];
                    for( i = 0 ; i < N ; i ++ )
                        S_BEST.NS[ i ] = S.NS[ i ]; 
                    for( i = 0 ; i < M+Q ; i ++ )
                        S_BEST.IW[ i ] = S.IW[ i ];  
                    S_BEST.IN = S.IN;
                    S_BEST.ON = S.ON; 
                    S_BEST.f = S.f; 
                    time_to_target = current_time;
                }  
            }
            non_improve = 0;
        }   
    }
       
    for(i=0;i<N;i++)  S.X[i]  = S_BEST.X[i]; 
	for(i=0;i<N;i++)  S.S[i]  = S_BEST.S[i];
	for(i=0;i<N;i++)  S.NS[i] = S_BEST.NS[i]; 
	for(i=0;i<M+Q;i++)  S.IW[i] = S_BEST.IW[i]; 
	S.f  = S_BEST.f;
    S.IN = S_BEST.IN;
    S.ON = S_BEST.ON; 
    //std::cout << "Total iteraciones " << iter << endl;
} // With a fixed value of K.  


/**************************     9.  Two-phase Search  ************************/
/*****************************************************************************/
void TwoPhaseSearch(float min, float max) {
    int i,j;
	int k, k0; 
	int select1, select2; 
	double start_time, s0_time, first_time, second_time; 

    Matrix = new int *[N];
	for(i=0;i<N;i++) 
		Matrix[i] = new int [N]();
	
	G_BEST.f = -99999; 
	start_time = clock();
    pseudoRandomSol(SC); 
    Tabu_Search(SC,min,max); 
    proof(SC); 
    first_time = (double)((clock()- start_time)/CLOCKS_PER_SEC); 
    GetMedia(); 
	Swap_tabu_search(SC,time_limit-first_time,min,max);  // the search of the second phase.
    proof(SC);
	time_one_run += first_time; 
	if(SC.f > G_BEST.f) {
		for(j=0;j<N;j++)   G_BEST.S[j]  = SC.S[j];
	    for(j=0;j<N;j++)   G_BEST.NS[j] = SC.NS[j]; 
	    for(j=0;j<N;j++)   G_BEST.X[j]  = SC.X[j];
	    for(j=0;j<M+Q;j++) G_BEST.IW[j] = SC.IW[j]; 
        G_BEST.f  = SC.f;
        G_BEST.IN = SC.IN;
        G_BEST.ON = SC.ON; 
	}
}
/*****************************************************************************/
/**************************   11. Main  Scheme     ***************************/
/*****************************************************************************/ 
int main(int argc, char **argv) {
    int i,j ; 
    File_Name  = argv[1];
    time_limit =  1.0*atoi(argv[2]); 
    core = 0.01*atoi(argv[3]);
    Itt = 1000*atoi(argv[4]);
    resta = atoi(argv[5]);
    pct_L = 0.01*atoi(argv[6]);
    lanbda = atoi(argv[7]);
    seed = atoi(argv[8]);
    float min = 1 - core;
    float max = 1 + core;
     
    Initializing(min,max); 
    AssignMemory(); 

    srand(seed);
    TwoPhaseSearch(min,max);
    sstd::cout << 100.0*(BKV-G_BEST.f)/BKV << endl;
    
    return 1;	
}
