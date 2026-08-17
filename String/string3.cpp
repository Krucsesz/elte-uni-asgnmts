#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

#include "vector.hpp"
using namespace cpl;

int main() {
    double L = 1.0, rho = 0.01, T = 10.0;
    double c = sqrt(T / rho);
    
    double gamma = 3.0; 

    int N = 400;
    double dx = L / N;
    double dt = 1.0 * (dx / c);
    double C2 = (c * dt / dx) * (c * dt / dx);
    
    // a gamma súrlódási tényezőből előforduló értékek
    double K = (gamma * dt) / 2.0;
    double denom = 1.0 + K;
    double coeff_old = 1.0 - K;

    double tMax = 2.0;
    int steps = tMax / dt;
    int skip_frames = 20;

    Vector y_old(N + 1), y_curr(N + 1), y_new(N + 1);
    ofstream dataFile("string_friction.data");

    const double pi = 4 * atan(1.0);

    // kf.
    for (int i = 0; i <= N; i++) {
        double x = i * dx;
        y_curr[i] = sin(pi * x / L);
        y_old[i] = y_curr[i];
    }
    y_curr[0] = 0.0; y_curr[N] = 0.0;

    for (int i = 0; i <= N; i++) dataFile << y_curr[i] << (i == N ? "" : "\t");
    dataFile << "\n";

    // 1. időlépés
    for (int i = 1; i < N; i++) {
        y_new[i] = y_curr[i] + 0.5 * C2 * (y_curr[i+1] - 2*y_curr[i] + y_curr[i-1]);
    }
    y_new[0] = 0.0; y_new[N] = 0.0;
    y_old = y_curr; y_curr = y_new;

    // fő ciklus
    for (int n = 1; n < steps; n++) {
        for (int i = 1; i < N; i++) {
            y_new[i] = (2.0 * y_curr[i] - coeff_old * y_old[i] + C2 * (y_curr[i+1] - 2.0 * y_curr[i] + y_curr[i-1])) / denom;
        }
        y_new[0] = 0.0; y_new[N] = 0.0;

        if (n % skip_frames == 0) {
            for (int i = 0; i <= N; i++) dataFile << y_new[i] << (i == N ? "" : "\t");
            dataFile << "\n";
        }
        y_old = y_curr; y_curr = y_new;
    }

    dataFile.close();
    cout << " Friction simulation completed and data saved to string_friction.data" << endl;
    return 0;
}