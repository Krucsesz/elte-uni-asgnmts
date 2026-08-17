#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

int N = 1024;
double rho = 0.02; // nagyon ritka gáz, hogy jobban tetszen az ideális gáztörvénynek
double T = 2.0;
double L;                 

double **r, **v, **a;     

double momentum_transfer = 0.0; // falnak adott összimpulzus
int measure_steps = 100;        // hány lépésenként mérjük a nyomást

void initialize();
void initPositions();
void initVelocities();
void rescaleVelocities();
double instantaneousTemperature();
double gasdev();

void initialize() {
    r = new double* [N];
    v = new double* [N];
    a = new double* [N];
    for (int i = 0; i < N; i++) {
        r[i] = new double [3];
        v[i] = new double [3];
        a[i] = new double [3];
    }
    initPositions();
    initVelocities();
}

void computeAccelerations() {
    for (int i = 0; i < N; i++)          
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;

    for (int i = 0; i < N-1; i++)        
        for (int j = i+1; j < N; j++) { 
            double rij[3];               
            double rSqd = 0;
            for (int k = 0; k < 3; k++) {
                rij[k] = r[i][k] - r[j][k];                 
                rSqd += rij[k] * rij[k];
            }
            
            // levágás, hogy ne szálljon el a cuccos
            if(rSqd < 9.0) {
                double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
                for (int k = 0; k < 3; k++) {
                     a[i][k] += rij[k] * f;
                     a[j][k] -= rij[k] * f;
                }
            }
        }
}

void velocityVerlet(double dt) {
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++) {
            r[i][k] += v[i][k] * dt + 0.5 * a[i][k] * dt * dt;

            if (r[i][k] <= 0) {
                r[i][k] = -r[i][k];
                v[i][k] = -v[i][k];
                momentum_transfer += 2.0 * abs(v[i][k]);
            }
            if (r[i][k] >= L) {
                r[i][k] = 2.0 * L - r[i][k]; 
                v[i][k] = -v[i][k];
                momentum_transfer += 2.0 * abs(v[i][k]);
            }
            
            v[i][k] += 0.5 * a[i][k] * dt;
        }
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            v[i][k] += 0.5 * a[i][k] * dt;
}

int main() {
    initialize();
    double dt = 0.005;
    ofstream file("datas/pressure.data");
    file << "Time P_measured P_ideal\n";

    double time = 0.0;
    
    for (int i = 0; i < 20000; i++) {
        velocityVerlet(dt);
        time += dt;

        if ((i + 1) % measure_steps == 0) {
            double T_inst = instantaneousTemperature();
            
            // 1. P = F / A
            double delta_t_total = measure_steps * dt;
            double force = momentum_transfer / delta_t_total;
            double area = 6.0 * L * L; // A kocka 6 falának összfelülete
            double P_measured = force / area;

            // 2. pV = NkT -> p = N*T / V
            double volume = L * L * L;
            double P_ideal = (N * T_inst) / volume;

            file << time << " " << P_measured << " " << P_ideal << "\n";
           
            momentum_transfer = 0.0;
        }
        
        if (i % 200 == 0) rescaleVelocities();
    }
    file.close();
    cout << "Simulation complete. Pressure data saved to 'datas/pressure.data'." << endl;
}

void initPositions() {

    // compute side of cube from number of particles and number density
    L = pow(N / rho, 1.0/3);

    // find M large enough to fit N atoms on an fcc lattice
    int M = 1;
    while (4 * M * M * M < N)
        ++M;
    double a = L / M;           // lattice constant of conventional cell

    // 4 atomic positions in fcc unit cell 
    double xCell[4] = {0.25, 0.75, 0.75, 0.25};
    double yCell[4] = {0.25, 0.75, 0.25, 0.75};
    double zCell[4] = {0.25, 0.25, 0.75, 0.75};

    int n = 0;                  // atoms placed so far
    for (int x = 0; x < M; x++)
        for (int y = 0; y < M; y++)
            for (int z = 0; z < M; z++)
                for (int k = 0; k < 4; k++)
                    if (n < N) {
                        r[n][0] = (x + xCell[k]) * a;
                        r[n][1] = (y + yCell[k]) * a;
                        r[n][2] = (z + zCell[k]) * a;
                        ++n;
                    }
}

double gasdev () {
     static bool available = false;
     static double gset;
     double fac, rsq, v1, v2;
     if (!available) {
          do {
               v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
               v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
               rsq = v1 * v1 + v2 * v2;
          } while (rsq >= 1.0 || rsq == 0.0);
          fac = sqrt(-2.0 * log(rsq) / rsq);
          gset = v1 * fac;
          available = true;
          return v2*fac;
     } else {
          available = false;
          return gset;
     }
}

void initVelocities() {

    // Gaussian with unit variance
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] = gasdev();

    // Adjust velocities so center-of-mass velocity is zero
    double vCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vCM[i] += v[n][i];
    for (int i = 0; i < 3; i++)
        vCM[i] /= N;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] -= vCM[i];

    // Rescale velocities to get the desired instantaneous temperature
    rescaleVelocities();
}

void rescaleVelocities() {
    double vSqdSum = 0;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vSqdSum += v[n][i] * v[n][i];
    double lambda = sqrt( 3 * (N-1) * T / vSqdSum );
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] *= lambda;
}

double instantaneousTemperature() {
    double sum = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            sum += v[i][k] * v[i][k];
    return sum / (3 * (N - 1));
}