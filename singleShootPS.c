#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct model{
    double *r;      //radial coordinate
    double *m;      //mass function
    double *f;      // f and g are real functions for potential one-form   A = f(r) + ig(r)
    double *g;      
    double *sigma;
    double *z;      //auxiliary variable: z = - sigma*N*g*r^2 / omega with z(0) = 0
    double *T;      //energy density
} model;

typedef struct parameters{
    double f0;
    double omega0;
    double c0;
    double sigma0;
    double rmax;        //maximum radial extent of the calculation
    int Nint;           //number of intervals the grid is divided into
    double h;           //step in the Runge-Kutta method: h = rmax / N
    int nzerotarget;    //wanted number of zero crossings (corresponds to the order of excited state with n=1 being ground state)
    double thresh;      //accuracy threshold for finding eigenfrequency of the proca star model
} parameters;

/****Functions****/
void printData(double*, double*, int, int, const char*);
void initialise(double, int, double, double, double, int);
void solveODE(double);
void shoot();
void rhsPSint(double*, double*, double*, double*, double,
              double, double, double, double, double);
int findFreqMinMax(double*, double*);
int countZeroCrossings(int*);
void rescale();
double VirialRelation();
void matchExterior();
void calcEnergyDensity();
double calcCharge();

// Global variables
parameters par;
model star;

/*=================================================================================================*/
int main(){
    double maxR = 200.0;
    int gridSize = 100000;
    double omega0 = 1.0;
    double f0 = 0.394;
    double sigma0 = 0.398;
    int sign_f;
    
    par.thresh = 1e-15;
    
    int nzerotarget = 0;        //0 = ground state; excited states > 1
    
    //initialise the grid
    initialise(maxR, gridSize, f0, sigma0, omega0, nzerotarget);
    
    shoot();
    
    matchExterior();
    
    rescale();
    
    //print data to file
    printData(star.r, star.m, 0, par.Nint - 1, "m.dat");
    printData(star.r, star.f, 0, par.Nint - 1, "f.dat");
    printData(star.r, star.g, 0, par.Nint - 1, "g.dat");
    printData(star.r, star.sigma, 0, par.Nint - 1, "sigma.dat");
    printData(star.r, star.z, 0, par.Nint - 1, "z.dat");
    
    int nzeros = countZeroCrossings(&sign_f);
    printf("\nNumber of zero crossings:      %d\n", nzeros);
    printf("(after rescaling) Frequency =      %f\n", par.omega0);
    
    double integral_I = VirialRelation();
    printf("Virial relation = %f\n", integral_I);
    
    calcEnergyDensity();
    printData(star.r, star.T, 0, par.Nint - 1, "T.dat");
    
    double Q = calcCharge();
    printf("Noether Charge      Q =  %f\n", Q);
    
    double ADM = star.m[par.Nint - 1];
    printf("ADM mass    M = %f\n", ADM);
    
    return 0;
}
/*===========================================================================================================*/
void shoot(){
    int cnt, success, cntmax;
    double om, ommin, ommax;
    int nzeros, sign_f;
    double rstop;
    
    printf("omega               Zero crossings                      \n");
    printf("========================================================\n");
    
    // (1) Find an upper and a lower limit for the frequency. This is done
    //   in findFreqMinMax, starting with omega=1 and then doubling or
    //   halving the value; see comments in that function. Should this
    //   search fail, we quit.
    success = findFreqMinMax(&ommin, &ommax);
    if(success == 0){
        printf("Could not find ommin and ommax in shoot\n"); 
        exit(0); 
        
    }
    else{
        printf("Using   omega in   [%22.16g,   %22.16g]\n", ommin, ommax);
    }
    
  // (2) Now we finetune the frequency such that it sits right on the
  //     threshold between nzerotarget and nzerotarget+1. This is the
  //     model we seek with nzerotarget zero crossings.
  printf("\n");
  printf("omega            Zero crossings        Sign_f\n");
  printf("=============================================\n");
  
    cnt = 0; cntmax = 100;
    nzeros = countZeroCrossings(&sign_f);
    
    while ((ommax - ommin) > par.thresh){
        om = (ommax + ommin) * 0.5;
        //int old_sign_f = sign_f;
        solveODE(om);
        nzeros = countZeroCrossings(&sign_f);
        
        printf("%10.10f               %d              %d \n", om, nzeros, sign_f);
    
        if (nzeros > par.nzerotarget) ommax = om;
        else ommin = om;
        
        cnt++;
        
        
    // If we fail to reach the threshold over many iterations, we likely
    // hit round-off error and reduce the accuracy requirement.        
        if (cnt > cntmax){
            par.thresh *= 10;
            cnt = 0;
            printf("Changed threshold to    %f\n", par.thresh);
        }
    }
    // ommin gives our best estimate for the model, since it has
    // nzero zero crossings whereas ommax has one zero crossing more.
    // We store this frequency in the parameter omega0 for further use.
    
    par.omega0 = ommin;
    solveODE(par.omega0);
    printf("%26.20g               %d           \n", om, nzeros);
}
/*============================================================================================================*/
int findFreqMinMax(double* ommin, double* ommax){
    
    int nzeros;
    int cnt = 0;            //cnt is a sanity measure to prevent being stuck in a loop forever,
    int cntmax = 100;       //it should never be triggered unless the frequency has values 2^1000 2^(-1000)...
    int sign_f;

        
    solveODE(par.omega0);
    nzeros = countZeroCrossings(&sign_f);    
    printf("%15.10g             %d                \n", par.omega0, nzeros);
    
  // Our strategy is as follows. We use the fact that the number of
  // zero crossings is a non-decreasing function of the frequency;
  // increasing omega0 gives you at least as many zero crossings as
  // before. We have a target, ntarget. We compute the number of zero
  // crossings for omega0 = 1.
  // If the resulting n>ntarget, we half omega0 until n <= ntarget
  // (note that the stable PS is always the limiting case to have
  // one more zero crossing -- the new crossing is at infinity in
  // the limiting case). Once we have n <= ntarget, the corresponding
  // two omega0 values are the brackets.
  // If the resulting n<=ntarget, we double omega0 until we have
  // n > ntarget and the ensuing omega0 values are our brackets.
    
    solveODE(par.omega0);
    nzeros = countZeroCrossings(&sign_f);
    
    printf("%15.10g             %d                \n", par.omega0, nzeros);
    
    if (nzeros > par.nzerotarget){
        *ommin = par.omega0;
        while(nzeros > par.nzerotarget && cnt < cntmax){
            *ommax = *ommin;
            *ommin /= 2;
            solveODE(*ommin);
            nzeros = countZeroCrossings(&sign_f);
            printf("%15.10g          %d           \n", *ommin, nzeros);
            cnt++;
        }
    }
    else{
        *ommax = par.omega0;
        while(nzeros <= par.nzerotarget && cnt < cntmax){
            *ommin = *ommax;
            *ommax *= 2;
            solveODE(*ommax);
            nzeros = countZeroCrossings(&sign_f);
            printf("%15.10g            %d          \n", *ommax, nzeros);
            cnt++;
        }
    }
    if (cnt == cntmax) return 0;        // Either upper or lower frequency limit has not been found
    else return 1;
    
    

}
/*============================================================================================================*/
int countZeroCrossings(int* sign_f){
    int nzeros = 0;
    int N = par.Nint;
    
    for (int i = 1; i < N; i++){
        /*if (fabs(star.f[i]) > 2 * star.f[0] || isnan(star.f[i])){
            *sign_f = ( (star.f[i-1] > 0) - (star.f[i-1] < 0) );
            // There is one problem here: We regard the present data point and
            // also the previous one as unrealiable. The previous point,
            // however, might have been counted as a zero crossing. If that
            // happened, we wish to undo this counting.
            if (star.f[i-1] * star.f[i-2] < 0) nzeros--;
            break;
        }*/
        if (isnan(star.g[i])){
            if (star.g[i-1] * star.g[i-2] < 0) nzeros--;
            break;
        }
        if (star.g[i] * star.g[i-1] < 0) nzeros++;
    }
    return nzeros;
}
/*============================================================================================================*/
void initialise(double maxR, int gridSize, double f0, double sigma0, double omega, int nzerotarget){
   
    /*This function sets the values of f0, omega0, c0, sigma0, rmax and Nint, 
    initialises the radial grid and allocates memory with boundary conditions*/
    
    par.f0 = f0;               
    par.omega0 = omega;           
    par.c0 = 0.0;                 
    par.sigma0 = sigma0;
    par.nzerotarget = nzerotarget;
    
    //rmax and Nint are specified by the user
    par.rmax = maxR;
    par.Nint = gridSize;
    double nn = (double) par.Nint;
    par.h = par.rmax/nn;
    
    //Allocate memory
    star.r     = (double*) malloc((size_t) par.Nint * sizeof(double));
    star.m     = (double*) malloc((size_t) par.Nint * sizeof(double));
    star.f     = (double*) malloc((size_t) par.Nint * sizeof(double));
    star.g     = (double*) malloc((size_t) par.Nint * sizeof(double));
    star.sigma = (double*) malloc((size_t) par.Nint * sizeof(double));
    star.z     = (double*) malloc((size_t) par.Nint * sizeof(double));
    star.T     = (double*) malloc((size_t) par.Nint * sizeof(double));

    
    //Initialising grid
    int n = par.Nint;
    double R = par.rmax;
    for (int i = 0; i < n; i++){    
        star.r[i] = R * i / (n - 1.0);
    }
    
    //Set boundary conditions
    star.f[0] = par.f0;
    star.g[0] = 0.0;
    star.m[0] = 0.0;
    star.sigma[0] = par.sigma0;
    star.z[0] = 0.0;
    
    printf("\n===================================\n");
    printf("Nint         = %d\n", par.Nint);
    printf("f0           = %f\n", par.f0);
    printf("sigma0       = %f\n", par.sigma0);
    printf("NzeroTarget  = %d\n", par.nzerotarget);
    printf("thresh        = %g\n", par.thresh);
    printf("====================================\n");
}
/*========================================================================================================*/
void solveODE(double omega){
    
    double r, f, z, m ,sigma;
    
    double k1_s, k2_s, k3_s, k4_s;
    double k1_f, k2_f, k3_f, k4_f;
    double k1_m, k2_m, k3_m, k4_m;
    double k1_z, k2_z, k3_z, k4_z;
    
    int n = par.Nint;
    double h = par.h;
    double om = omega;
    
    double rhs_m, rhs_f, rhs_z, rhs_sigma;
    
    for (int i = 1; i < n; i++){
        
/*......1st R-K step...............................*/
        r = star.r[i-1];
        f = star.f[i-1];
        z = star.z[i-1];
        sigma = star.sigma[i-1];
        m = star.m[i-1];
        
        rhsPSint(&rhs_m, &rhs_f, &rhs_z, &rhs_sigma,
                 r, m, f, z, sigma, om);
        
        k1_s = h * rhs_sigma;
        k1_f = h * rhs_f;
        k1_m = h * rhs_m;
        k1_z = h * rhs_z;
    
/*......2nd R-K step...............................*/
        
        r = star.r[i-1] + 0.5 * h;
        f = star.f[i-1] + 0.5 * k1_f;
        z = star.z[i-1] + 0.5 * k1_z;
        sigma = star.sigma[i-1] + 0.5 * k1_s;
        m = star.m[i-1] + 0.5 * k1_m;

        rhsPSint(&rhs_m, &rhs_f, &rhs_z, &rhs_sigma,
                 r, m, f, z, sigma, om);
        
        k2_s = h * rhs_sigma;
        k2_f = h * rhs_f;
        k2_m = h * rhs_m;
        k2_z = h * rhs_z;
        
/*......3rd R-K step...............................*/
        
        r = star.r[i-1] + 0.5 * h;
        f = star.f[i-1] + 0.5 * k2_f;
        z = star.z[i-1] + 0.5 * k2_z;
        sigma = star.sigma[i-1] + 0.5 * k2_s;
        m = star.m[i-1] + 0.5 * k2_m;
        
        rhsPSint(&rhs_m, &rhs_f, &rhs_z, &rhs_sigma,
                 r, m, f, z, sigma, om);
        
        k3_s = h * rhs_sigma;
        k3_f = h * rhs_f;
        k3_m = h * rhs_m;
        k3_z = h * rhs_z;
        
/*......4th R-K step...............................*/
        
        r = star.r[i-1] + h;
        f = star.f[i-1] + k3_f;
        z = star.z[i-1] + k3_z;
        sigma = star.sigma[i-1] + k3_s;
        m = star.m[i-1] + k3_m;
        
        rhsPSint(&rhs_m, &rhs_f, &rhs_z, &rhs_sigma,
                 r, m, f, z, sigma, om);
        
        k4_s = h * rhs_sigma;
        k4_f = h * rhs_f;
        k4_m = h * rhs_m;
        k4_z = h * rhs_z;

/*......updating variables.........................*/
        
        r = star.r[i];
        star.m[i] = star.m[i-1] + (k1_m + 2.0*k2_m + 2.0*k3_m + k4_m)/6.0;
        star.f[i] = star.f[i-1] + (k1_f + 2.0*k2_f + 2.0*k3_f + k4_f)/6.0;
        star.z[i] = star.z[i-1] + (k1_z + 2.0*k2_z + 2.0*k3_z + k4_z)/6.0;
        star.sigma[i] = star.sigma[i-1] + (k1_s + 2.0*k2_s + 2.0*k3_s + k4_s)/6.0;
        double N = 1 - 2 * star.m[i] / star.r[i];
        star.g[i] = (-1 * par.omega0 * star.z[i]) / (star.r[i] * star.r[i] * star.sigma[i] * N);
    }
    
}
/*==========================================================================================================*/
void rhsPSint(double* rhs_m, double* rhs_f, double* rhs_z, double* rhs_sigma, 
                double r, double m, double f, double z, double sigma, double om){
    
    double f_0 = par.f0;
    double sigma_0 = par.sigma0;
    double N = 1.0;
    
    //asymptotic behaviour as r --> 0
    if (r < 1.0e-15){
        //as r->0:
        /*   f(r) = f_0 + f_0*r^2*(1-om^2/sigma_0^2)/6
         *   m(r) = f_0^2*r^3 / (6sigma_0^2)
         *   sigma(r) = sigma_0 + f_0^2*r^2 / (2sigma_0)
         *   N(r) = 1-f_0^2*r^2 / (3sigma_0^2)
         *   z(r) = r^3*f_0*N / (3sigma_0)
         * */
        
        *rhs_f = f_0 + f_0 * (1.0 - (om * om) / (sigma_0 * sigma_0)) * r * r / 6.0;
        
        *rhs_m = f_0 * f_0 * r * r * r / (sigma_0 * sigma_0 * 6.0);
        
        *rhs_sigma = sigma_0 + (f_0 * f_0 * r * r) / (2.0 * sigma_0);
        
        N = 1.0 - (f_0 * f_0 * r * r) / (3.0 * sigma_0 * sigma_0);
        
        *rhs_z = (r * r * r * f_0 * N) / (3.0 * sigma_0);
    }
    else{
        
        N = 1.0 - 2.0 * m / r;
        
        *rhs_sigma = ((om * om * z * z) / (r * r * r * sigma) + (f * f * r)) / (N * N);
        
        *rhs_f = (z * (sigma - om * om / (sigma * N))) / (r * r);
        
        *rhs_m = (z * z) / (2.0 * r * r) + (om * om * z * z) / (sigma * sigma * N * 2.0 * r * r) + (f * f * r * r) / (N * sigma * sigma * 2.0);
        
        *rhs_z = r * r * f / (sigma * N);
    }
                    
}
/*==========================================================================================================*/
void matchExterior(){
    double c1, c2, c3c3;
    int i_stop = -1;
    int i_match;

    double N;
    
    double r, f, z, m ,sigma;
    
    double k1_s, k2_s, k3_s, k4_s;
    double k1_f, k2_f, k3_f, k4_f;
    double k1_m, k2_m, k3_m, k4_m;
    double k1_z, k2_z, k3_z, k4_z;
    
    int n = par.Nint;
    double h = par.h;
    double om = par.omega0;
    
    double rhs_m, rhs_f, rhs_z, rhs_sigma;
    
  // At this stage we have the correct frequency from the shooting
  // algorithm and the computed model. Here we will remove
  // the diverging part in the profile and replace it with a smooth
  // exterior. Note that we are not rescaling time and frequency
  // yet, since we will need the complete profile with exterior to do
  // that.
    for (int i = 1; i < n; i++){
        if (fabs(star.f[i]) > 2 * star.f[0] || isnan(star.f[i])){
            i_stop = i;
            break;
        }
    }
  // We now have a model that either drops in scalar amplitude all
  // the way to the edge of our grid (unlikely) or will diverge
  // somewhere along our grid towards either positive or negative
  // values. We regularize this divergence and replace it with a
  // smooth exterior solution as follows:
  // (1) Search backwards from the truncation point (i_stop -- for
  //     which we need to find the index) until we find a minimum
  //     in abs(f). That is the point from which we construct the
  //     exterior.
    for (int i = i_stop - 1; i > 0; i--){
        if (fabs(star.f[i]) < fabs(star.f[i+1]) && fabs(star.f[i]) < fabs(star.f[i-1])){
            i_match = i;
            break;
        }
    }
    if(i_match == 1)
        { printf("Search for matching point yielded i = 1\n\n"); exit(0); }

    printf("Matching to exterior at   r[%d] = %g\n\n", i_match, star.r[i_match]);

  // (2) We now match the scalar amplitude to an exterior function
  //
  //       f   = c1 * exp(-r * sqrt(1-om^2)) / r
  //       z   = - c2 * r * sigma * N * exp(-r * sqrt(1-om^2)) / r 
  //       log(sigma) = - c3^2 exp(-2*r*sqrt(1-om^2))/(2*r*(sqrt(1-om^2))^3)
  //
  //     and determine c1 and c2 through matching to f[i_match] and
  //     z[i_match].
  //     We could formally integrate z from f, but that deviates
  //     more strongly from the asymptotic limit z -> 0 at infinity.
    r = star.r[i_match];
  
    N = 1.0 - 2.0 * star.m[i_match] / star.r[i_match];
    c1  =   star.f[i_match] * r * exp(r * sqrt(1 - om*om));
    c2  =   (-1.0 * star.z[i_match] * sqrt(1-om*om) * exp(r * sqrt(1-om*om))) / (r * star.sigma[i_match] * N);
    c3c3 = -2*log(star.sigma[i_match])*sqrt(1-om*om)*sqrt(1-om*om)*sqrt(1-om*om)*r*exp(2*r*sqrt(1-om*om));
    
    double ADM = star.m[i_match];
  
  
    for (int i = i_match + 1; i < n; i++){
        double r = star.r[i];
    
      
      
        
/*......1st R-K step...............................*/
        r = star.r[i-1];
        sigma = star.sigma[i-1];
        //sigma = exp((-1*c3c3*exp(-2*r*sqrt(1-om*om)))/(2*sqrt(1-om*om)*sqrt(1-om*om)*sqrt(1-om*om)*r));
        m = star.m[i-1];
        //m = ADM;
        //f = star.f[i-1];
        f = c1 * exp(-1 * r * sqrt(1-om*om)) / r;
        N = 1 - 2*m/r;
        z = star.z[i-1];
        //z = (-1 * r * sigma * N * c2 * exp(-1 * r * sqrt(1-om*om))) / (sqrt(1-om*om));
        
        rhsPSint(&rhs_m, &rhs_f, &rhs_z, &rhs_sigma,
                 r, m, f, z, sigma, om);
        
        k1_s = h * rhs_sigma;
        k1_f = h * rhs_f;
        k1_m = h * rhs_m;
        k1_z = h * rhs_z;
    
/*......2nd R-K step...............................*/
        
        r = star.r[i-1] + 0.5 * h;
        sigma = star.sigma[i-1] + 0.5 * k1_s;
        //sigma = exp((-1*c3c3*exp(-2*r*sqrt(1-om*om)))/(2*sqrt(1-om*om)*sqrt(1-om*om)*sqrt(1-om*om)*r));
        m = star.m[i-1] + 0.5 * k1_m;
        //m = ADM;
        //f = star.f[i-1] + 0.5*k1_f;
        f = c1 * exp(-1 * r * sqrt(1-om*om)) / r;
        N = 1 - 2*m/r;
        z = star.z[i-1] + 0.5*k1_z;
        //z = (-1 * r * sigma * N * c2 * exp(-1 * r * sqrt(1-om*om))) / (sqrt(1-om*om));
        
        rhsPSint(&rhs_m, &rhs_f, &rhs_z, &rhs_sigma,
                 r, m, f, z, sigma, om);
        
        k2_s = h * rhs_sigma;
        k2_f = h * rhs_f;
        k2_m = h * rhs_m;
        k2_z = h * rhs_z;
        
/*......3rd R-K step...............................*/
        
        r = star.r[i-1] + 0.5 * h;
        sigma = star.sigma[i-1] + 0.5 * k2_s;
        //sigma = exp((-1*c3c3*exp(-2*r*sqrt(1-om*om)))/(2*sqrt(1-om*om)*sqrt(1-om*om)*sqrt(1-om*om)*r));
        m = star.m[i-1] + 0.5 * k2_m;
        //m = ADM;
        //f = star.f[i-1] + 0.5 * k2_f;
        f = c1 * exp(-1 * r * sqrt(1-om*om)) / r;
        N = 1 - 2*m/r;
        z = star.z[i-1] + 0.5 * k2_z;
        //z = (-1 * r * sigma * N * c2 * exp(-1 * r * sqrt(1-om*om))) / (sqrt(1-om*om));
        
        rhsPSint(&rhs_m, &rhs_f, &rhs_z, &rhs_sigma,
                 r, m, f, z, sigma, om);
        
        k3_s = h * rhs_sigma;
        k3_f = h * rhs_f;
        k3_m = h * rhs_m;
        k3_z = h * rhs_z;
        
/*......4th R-K step...............................*/
        
        r = star.r[i-1] + h;
        sigma = star.sigma[i-1] + k3_s;
        //sigma = exp((-1*c3c3*exp(-2*r*sqrt(1-om*om)))/(2*sqrt(1-om*om)*sqrt(1-om*om)*sqrt(1-om*om)*r));
        m = star.m[i-1] + k3_m;
        //m = ADM;
        //f = star.f[i-1] + k3_f;
        f = c1 * exp(-1 * r * sqrt(1-om*om)) / r;
        N = 1 - 2*m/r;
        z = star.z[i-1] + k3_z;
        //z = (-1 * r * sigma * N * c2 * exp(-1 * r * sqrt(1-om*om))) / (sqrt(1-om*om));
        
        rhsPSint(&rhs_m, &rhs_f, &rhs_z, &rhs_sigma,
                 r, m, f, z, sigma, om);
        
        k4_s = h * rhs_sigma;
        k4_f = h * rhs_f;
        k4_m = h * rhs_m;
        k4_z = h * rhs_z;

/*......updating variables.........................*/
        
        r = star.r[i];
        star.m[i] = star.m[i-1] + (k1_m + 2.0*k2_m + 2.0*k3_m + k4_m)/6.0;
        //star.m[i] = ADM;
        star.sigma[i] = star.sigma[i-1] + (k1_s + 2.0*k2_s + 2.0*k3_s + k4_s)/6.0;
        //star.sigma[i] = exp((-1*c3c3*exp(-2*r*sqrt(1-om*om)))/(2*sqrt(1-om*om)*sqrt(1-om*om)*sqrt(1-om*om)*r));
        star.f[i] = c1 * exp(-1 * r * sqrt(1-om*om)) / r;
        //star.f[i] = star.f[i-1] + (k1_f + 2.0*k2_f + 2.0*k3_f + k4_f)/6.0;
        N = 1 - 2*star.m[i]/star.r[i];
        //star.z[i] = (-1 * r * star.sigma[i] * N * c2 * exp(-1 * r * sqrt(1-om*om))) / (sqrt(1-om*om));
        star.z[i] = star.z[i-1] + (k1_z + 2.0*k2_z + 2.0*k3_z + k4_z)/6.0;
        
        star.g[i] = (-1 * par.omega0 * star.z[i]) / (star.r[i] * star.r[i] * star.sigma[i] * N);
  }
}

/*==========================================================================================================*/
void rescale(){
    
   // We want the lapse = 1 at infinity. Since we do not go to infinity,
   // our condition on the lapse is sigma = 1 at the outer radius. 
   //
   // Note that we also need to change omega such that
   //
   // omega_old / (sigma_old * sqrt(N_old)) = omega_new / (sigma_new * sqrt(N_new))

   int n = par.Nint;

   double sigma_target = 1.0;
   double sigma_org = star.sigma[n-1];
   
   printf("Sigma_org = %g     Sigma_target = %g\n", sigma_org, sigma_target);
   
   for (int i = 0; i < n; i++){
       star.sigma[i] *= 1/sigma_org;
    }
   
   double N = 1.0 - 2*star.m[n-1]/star.r[n-1];
   
   par.omega0 = par.omega0 / sigma_org;

    
   
}
/*==========================================================================================================*/
double VirialRelation(){
    double I = 0.0;
    double h = par.h;
    double N,f,g,sigma,r,m;
    double om = par.omega0;
    
    for (int i = 1; i < par.Nint - 1; i++){
        
        r = star.r[i];
        m = star.m[i];
        sigma = star.sigma[i];
        f = star.f[i];
        g = star.g[i];
        N = 1.0 - 2.0*m/r;
        
        I += h * sigma * r * r *(g*g - (f*f*(4*N - 1.0))/(sigma*sigma*N*N) - (sigma*sigma*N*N*g*g)/(om*om));
        //printf("%d ________________________    %f\n", i, I);
    }
    
    
    return I;
}
/*==========================================================================================================*/
void calcEnergyDensity(){
    
    double om = par.omega0;
    double N,f,g,sigma,r,m;
    
    for (int i = 0; i < par.Nint; i++){
        
        f = star.f[i];
        g = star.g[i];
        sigma = star.sigma[i];
        m = star.m[i];
        r = star.r[i];
        
        if (i != 0) N = 1.0 - 2*m/r;
        else N = 1.0;
        
        star.T[i] = (sigma*sigma*N*N*g*g)/(2*om*om) + 0.5 * (g*g*N + (f*f)/(N*sigma*sigma));
        //printf("T[%d]  =  %f\n", i, star.T[i]);
    }

}
/*==========================================================================================================*/
double calcCharge(){
    double Q = 0.0;
    double h = par.h;
    double N,f,g,sigma,r,m;
    double om = par.omega0;
    
    for (int i = 1; i < par.Nint - 1; i++){
        
        r = star.r[i];
        m = star.m[i];
        sigma = star.sigma[i];
        f = star.f[i];
        g = star.g[i];
        N = 1.0 - 2.0*m/r;
        
        Q += h*r*r*g*g*sigma*N;
    }
    
    //Q *= 4.0 * 3.14159;
    Q /= om;
    
    return Q;
}
/*==========================================================================================================*/
void printData(double* x, double* y, int n1, int n2, const char* ofil){
    
    FILE* ofp;
    
    ofp = fopen(ofil, "w");
    
    if (ofp == NULL) {
        printf("Error. Could not open the file %s.\n", ofil);
        exit(-1);
    }
    fprintf(ofp, "# %s\n", ofil);
    for (int i = n1; i <= n2; i++){
        //if (fabs(y[i]) > 1 || isnan(y[i])) break;
        fprintf(ofp, "%10.10f   %10.10f\n", x[i], y[i]);
    }
    fclose(ofp);
}
