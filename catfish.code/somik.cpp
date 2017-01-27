/* 
 * Program numerically solves PDE of Adsorption-modified Gierer-Meinhardt model
 * with Neuman boundary conditions. 
 *                                                                                   
 * Reaction and diffusion terms are computed sequentially. We use Runge-Kutta
 * scheme for the reaction part. ADI technique coupled with flux sweeping algorithm 
 * is used for solving two-dimensional diffusion problem.
 *                                                                                   
 * For the details see (Nesterenko et. al, Plos Comp. Biol., 2016). Please cite
 * this paper if you use the code.
 * */

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>

using namespace std;

// model parameters 
const double k1   = 0.97,
             k_1  = 0.1,
             ro   = 0.6,
             mu_u = 0.03,
             mu_v = 0.08,
             D    = 1.0;
    
// computational scheme parameters

const double tau  = 0.005,  // time step 
           stepr  = 0.05,   // grid step
           t_res  = 10.,    // output interval
               T  = 8.,  // integration time 
               L  = 40.;     // reactor length
// grid dimensions
const int      X  = (int)(L/stepr),
               Y  = X;
// settings of adsorption gradient    
const int som_m1 = (int)((X+1)*0.2),
          som_m2 = (int)((X+1)*0.7);
const double som_wup = 50.,
             som_wdown = 20.; 
   
// functions predefinitions
void RungeKutta(double** uu,double** vv,double** ww,double** ur,double** vr,double** wr);

void FluxSweepingSolver(double* array_lin,double Dd,double* prog,int npro);    
void Diffusion_DIA (double** array_rk,double** array_1,double** array_2,double Dd);

void output_3D (double** uu, double** vv, double** ww, char* File);
void output_diag (double** uu, double** vv, double** ww, char* File);
   
// Reaction parts of PDE 
double Growth_u (double u, double v, double w, double w_0);
double Growth_v (double u, double v, double w, double w_0);
double Growth_w (double u, double v, double w, double w_0);
        

    
//********************************* ПРОГРАММА **************************************************************************************************************
double uK1[X+1][Y+1], uK2[X+1][Y+1], uK3[X+1][Y+1], uK4[X+1][Y+1];
double vK1[X+1][Y+1], vK2[X+1][Y+1], vK3[X+1][Y+1], vK4[X+1][Y+1];
double wK1[X+1][Y+1], wK2[X+1][Y+1], wK3[X+1][Y+1], wK4[X+1][Y+1];
double prog_x[X+1];
double prog_y[Y+1];
double w_0_som[Y+1];
     	
int main() {    
      
      char File_3D[20], File_diag[20], chislo[4];

      double **u0 = new double*[X+1];
      double **u1 = new double*[X+1];
      double **u2 = new double*[X+1];
      double **ur = new double*[X+1];
      
      double **v0 = new double*[X+1];
      double **v1 = new double*[X+1];
      double **v2 = new double*[X+1];
      double **vr = new double*[X+1];
      
      double **w0 = new double*[X+1];
      double **wr = new double*[X+1];
            

      for (int i=0; i<=X; i++) {
          u0[i]= new double [Y+1];
              u1[i]= new double [Y+1];
              u2[i]= new double [Y+1];
              ur[i]= new double [Y+1];
              
              v0[i]= new double [Y+1];
              v1[i]= new double [Y+1];
              v2[i]= new double [Y+1];
              vr[i]= new double [Y+1]; 	
      
              w0[i]= new double [Y+1];
              wr[i]= new double [Y+1];
      }
                      
      // Initial conditions   
      for (int j = 0; j < som_m1; j++)
          w_0_som[j] = som_wup;
      for (int j = som_m2; j <= X; j++)
          w_0_som[j] = som_wdown; 
      for (int j = som_m1; j < som_m2; j++) {
          w_0_som[j]=som_wup + (som_wdown-som_wup)*1.08551*
              (tanh((j-0.5*(som_m1+som_m2))*0.008)/2.+0.461);
      }
      printf("Binding site gradient:\n");
      for (int j = 0; j <= X; j++) {
          printf("%8.3f ", w_0_som[j]);
          if (j % 10 == 0) 
            printf("\n");
      }

                      
      for (int i=0; i<=X; i++) {
          for (int j=0; j<=Y; j++) {
              u0[i][j]=0.001*(rand() % 101  + 1.);
              v0[i][j]=0.001*(rand() % 1001 + 1.);
              w0[i][j]=w_0_som[j];
          }
      }        
      
      // Initial conditions print
      
      strcpy(File_3D,"res_3D-d0000.txt");
      output_3D(u0,v0,w0,File_3D);
      
      strcpy(File_diag,"res_diag-d0000.txt");
      output_diag(u0,v0,w0,File_diag);
      
      // MAIN SOLVER
      printf("Starting calculations.\n");  

      double t = 0.; // time counter
      int n = 1;     // file counter

      while ( t <= T ) {
          printf("."); fflush(stdout);
          t += tau;
              
          RungeKutta(u0,v0,w0,ur,vr,wr);
          Diffusion_DIA (ur,u1,u2,D);
          Diffusion_DIA (vr,v1,v2,D);
          
          // Output          
          if ( t > n*t_res - 1.e-7 ) {
              printf("\n%i\t%.3f\n",n,t);
              gcvt(n,4,chislo);
                                              
              if (n<10) {
                strcpy(File_3D,"res_3D-d000");
                strcpy(File_diag,"res_diag-d000");
              }
              else if (n<100) {
                strcpy(File_3D,"res_3D-d00");
                strcpy(File_diag,"res_diag-d00");
              }
              else if (n<1000) {
                strcpy(File_3D,"res_3D-d0");
                strcpy(File_diag,"res_diag-d0");
              }
              else {
                strcpy(File_3D,"res_3D-d");	
                strcpy(File_diag,"res_diag-d");	
              }                         		 
              strcat(File_3D,chislo);	 
              strcat(File_diag,chislo);	 
              strcat(File_3D,"txt");                                                                                                             
              strcat(File_diag,"txt");                                                                                                             
                      
              output_3D(u2,v2,wr,File_3D);
              output_diag(u2,v2,wr,File_diag); 
                      
              n++;               
          }
      
          // Array transfer        
          for (int i=0; i<=X; i++) {
              for (int j=0; j<=Y; j++) {
                u0[i][j]=u2[i][j];
                v0[i][j]=v2[i][j];
                w0[i][j]=wr[i][j];      			
              }
          }    				
      } // End main cycle
      
      int _cl = clock()/CLOCKS_PER_SEC;
      printf("\nThe end. Total time: %d:%d:%d", _cl / 3600, (_cl % 3600) / 60, (_cl % 60) );
      
      for(int i=0; i<=X; i++) {
        delete []u0[i];
        delete []u1[i];
        delete []u2[i];
        delete []ur[i];

        delete []v0[i];
        delete []v1[i];
        delete []v2[i];
        delete []vr[i];

        delete []w0[i];
        delete []wr[i];        	   
      }
              
      delete []u0;
      delete []u1;
      delete []u2;
      delete []ur;
              
      delete []v0;
      delete []v1;
      delete []v2;
      delete []vr;
      
      delete []w0;
      delete []wr;
} // end main
    
/* 
 * Reaction function implementations.
 *
 * */

double Growth_u (double u, double v, double w, double w_0) {
        return((ro*(u+w_0-w)*(u+w_0-w)/v)-mu_u*u-k1*w*u+k_1*(w_0-w));
}

double Growth_v (double u, double v, double w, double w_0) {
    return(ro*(u+w_0-w)*(u+w_0-w)-mu_v*v);
}

double Growth_w (double u, double v, double w, double w_0) {
    return((-1.)*k1*w*u+(k_1+mu_u)*(w_0-w));
}     

/* 
 * Computational function implementations.
 *
 * */


void RungeKutta(double** uu,double** vv,double** ww,double** ur,double** vr,double** wr) {
  for (int i=0; i<=X; i++) {
    for (int j=0; j<=Y; j++) {
      uK1[i][j]=Growth_u(uu[i][j],vv[i][j],ww[i][j],w_0_som[j]);
      vK1[i][j]=Growth_v(uu[i][j],vv[i][j],ww[i][j],w_0_som[j]);
      wK1[i][j]=Growth_w(uu[i][j],vv[i][j],ww[i][j],w_0_som[j]);

      uK2[i][j]=Growth_u(uu[i][j]+uK1[i][j]*tau/2.,vv[i][j]+vK1[i][j]*tau/2.,ww[i][j]+wK1[i][j]*tau/2.,w_0_som[j]);
      vK2[i][j]=Growth_v(uu[i][j]+uK1[i][j]*tau/2.,vv[i][j]+vK1[i][j]*tau/2.,ww[i][j]+wK1[i][j]*tau/2.,w_0_som[j]);
      wK2[i][j]=Growth_w(uu[i][j]+uK1[i][j]*tau/2.,vv[i][j]+vK1[i][j]*tau/2.,ww[i][j]+wK1[i][j]*tau/2.,w_0_som[j]);

      uK3[i][j]=Growth_u(uu[i][j]+uK2[i][j]*tau/2.,vv[i][j]+vK2[i][j]*tau/2.,ww[i][j]+wK2[i][j]*tau/2.,w_0_som[j]);
      vK3[i][j]=Growth_v(uu[i][j]+uK2[i][j]*tau/2.,vv[i][j]+vK2[i][j]*tau/2.,ww[i][j]+wK2[i][j]*tau/2.,w_0_som[j]);
      wK3[i][j]=Growth_w(uu[i][j]+uK2[i][j]*tau/2.,vv[i][j]+vK2[i][j]*tau/2.,ww[i][j]+wK2[i][j]*tau/2.,w_0_som[j]);

      uK4[i][j]=Growth_u(uu[i][j]+uK3[i][j]*tau,vv[i][j]+vK3[i][j]*tau,ww[i][j]+wK3[i][j]*tau,w_0_som[j]);
      vK4[i][j]=Growth_v(uu[i][j]+uK3[i][j]*tau,vv[i][j]+vK3[i][j]*tau,ww[i][j]+wK3[i][j]*tau,w_0_som[j]);
      wK4[i][j]=Growth_w(uu[i][j]+uK3[i][j]*tau,vv[i][j]+vK3[i][j]*tau,ww[i][j]+wK3[i][j]*tau,w_0_som[j]);

      ur[i][j] = uu[i][j] + tau*(uK1[i][j]+uK4[i][j]+2.*(uK2[i][j]+uK3[i][j]))/6.;
      vr[i][j] = vv[i][j] + tau*(vK1[i][j]+vK4[i][j]+2.*(vK2[i][j]+vK3[i][j]))/6.;
      wr[i][j] = ww[i][j] + tau*(wK1[i][j]+wK4[i][j]+2.*(wK2[i][j]+wK3[i][j]))/6.;
    }
  }
} // end RungeKutta

    
void Diffusion_DIA (double** array_rk,double** array_1,double** array_2,double Dd) {
    double array_rk_x[X+1];
    double array_rk_y[Y+1];
    
    // over x
    for(int j=0; j<=Y; j++) {
      for(int i=0; i<=X; i++)
        array_rk_x[i]=array_rk[i][j];      

      FluxSweepingSolver(array_rk_x,D,prog_x,X);

      for(int i=0; i<=X; i++)
        array_1[i][j]=prog_x[i];                    
    }
      
    // over y
    for(int i=0; i<=X; i++) {
      for(int j=0; j<=Y; j++)
        array_rk_y[j]=array_1[i][j];   

      FluxSweepingSolver(array_rk_y,D,prog_y,Y);

      for(int j=0; j<=Y; j++)
        array_2[i][j]=prog_y[j];                        
    }  
        
}


void FluxSweepingSolver(double* array_lin,double Dd,double* prog,int npro) {
    
    double *P_prog = new double[npro];
    double *R_prog = new double[npro];
    double *A_prog = new double[npro];
    double *potok = new double[npro];

    P_prog[0] = 2*tau/stepr;
    R_prog[0] = array_lin[0];
    
    for (int i=1; i<=npro-1; i++) {
      A_prog[i] = 1.+(stepr/tau)*(P_prog[i-1]+(stepr/D));
      P_prog[i] = (1./A_prog[i])*(P_prog[i-1]+(stepr/D));
      R_prog[i] = (1./A_prog[i])*((P_prog[i-1]+(stepr/D))*((stepr/tau)*array_lin[i])+R_prog[i-1]);
    }
    
    prog[npro] = (2.*tau*D*R_prog[npro-1]+(stepr+P_prog[npro-1]*D)*stepr*array_lin[npro])/(stepr*stepr+stepr*P_prog[npro-1]*D+2.*tau*D);
    potok[npro-1] = (stepr*(array_lin[npro]-prog[npro]))/(2.*tau);
    prog[npro-1] = P_prog[npro-1]*potok[npro-1]+R_prog[npro-1];
    
    for(int i=npro-2; i>=0; i--) {
      potok[i] = potok[i+1]-(stepr/tau)*prog[i+1]+(stepr/tau)*array_lin[i+1];
      prog[i] = P_prog[i]*potok[i]+R_prog[i];
    }
    
    delete []P_prog;
    delete []R_prog;
    delete []A_prog;
    delete []potok;
}
  
/*
 * Output function implementation.
 *
 */
void output_3D (double** uu, double** vv, double** ww, char* File) {        
  FILE * out;
  if ((out = fopen(File, "wt")) == NULL)
    fprintf(stderr, "Cannot open file\n");

  fprintf(out, "x       y       u               v               w\n\n");

  for(int i=0; i<=X; i++)
    for(int j=0; j<=Y; j++)
      fprintf(out,"%.3f\t%.3f\t%.7e\t%.7e\t%.7e\n", i*stepr, j*stepr, uu[i][j], vv[i][j], ww[i][j]);

  fclose(out);
}

void output_diag (double** uu, double** vv, double** ww, char* File) {        
  FILE * out;
  if ((out = fopen(File, "wt")) == NULL)
    fprintf(stderr, "Cannot open file\n");

  fprintf(out,"x       y       u               v               w\n\n");

  for(int i=0; i<=X; i++)
    fprintf(out, "%.3f\t%.3f\t%.7e\t%.7e\t%.7e\n", i*stepr, i*stepr, uu[i][i], vv[i][i], ww[i][i]);

  fclose(out);
}
