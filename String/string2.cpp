#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
using namespace std;

#include "vector.hpp"
using namespace cpl;

void solveWave(const Vector& y_init, const string& filename, int N, double C2, int steps, int skip_frames) {
    ofstream dataFile(filename);
    Vector y_old = y_init;
    Vector y_curr = y_init;
    Vector y_new(N + 1);

    // kezdeti állapot mentése
    for (int i = 0; i <= N; i++) dataFile << y_curr[i] << (i == N ? "" : "\t");
    dataFile << "\n";

    // 1. lépés
    for (int i = 1; i < N; i++) y_new[i] = y_curr[i] + 0.5 * C2 * (y_curr[i+1] - 2*y_curr[i] + y_curr[i-1]);
    y_new[0] = 0.0; y_new[N] = 0.0;
    y_old = y_curr; y_curr = y_new;

    // fő ciklus
    for (int n = 1; n < steps; n++) {
        for (int i = 1; i < N; i++) {
            y_new[i] = 2.0 * y_curr[i] - y_old[i] + C2 * (y_curr[i+1] - 2.0 * y_curr[i] + y_curr[i-1]);
        }
        y_new[0] = 0.0; y_new[N] = 0.0;

        if (n % skip_frames == 0) {
            for (int i = 0; i <= N; i++) dataFile << y_new[i] << (i == N ? "" : "\t");
            dataFile << "\n";
        }
        y_old = y_curr; y_curr = y_new;
    }
    dataFile.close();
}

int main() {
    double L = 1.0, rho = 0.01, T = 10.0;
    double c = sqrt(T / rho);
    int N = 400;
    double dx = L / N;
    double dt = 1.0 * (dx / c);
    double C2 = (c * dt / dx) * (c * dt / dx);
    double tMax = 0.2;
    int steps = tMax / dt;
    int skip_frames = 20;

    Vector y_ic1(N + 1), y_ic2(N + 1), y_ic_sum(N + 1), y_mixed(N + 1);
    const double pi = 4 * atan(1.0);

    for (int i = 0; i <= N; i++) {
        double x = i * dx;
        
        // 5/a feladat: linearitás tesztelése különböző kezdeti állapotokkal
        // 1. pengetés: bal oldalon felfelé
        double x1 = 0.2 * L;
        if (x <= x1) y_ic1[i] = x / x1;
        else         y_ic1[i] = (L - x) / (L - x1);

        // 2. pengetés: jobb oldalon lefelé (csak hogy látványosabb legyen :3)
        double x2 = 0.75 * L;
        if (x <= x2) y_ic2[i] = -1.0 * (x / x2);
        else         y_ic2[i] = -1.0 * (L - x) / (L - x2);

        // 3. kombinált pengetés: a kettő összege
        y_ic_sum[i] = y_ic1[i] + y_ic2[i];

        // 5/b feladat: kevert normálmódusok 
        // legyenek a 2. és 3. normálmódusok keveréke, de a 3. legyen valamivel erősebb
        y_mixed[i] = sin(2 * pi * x / L) + 0.6 * sin(3 * pi * x / L);
    }

    cout << " 5/a: testing linearity \n";
    solveWave(y_ic1, "datas/lin_ic1.data", N, C2, steps, skip_frames);
    solveWave(y_ic2, "datas/lin_ic2.data", N, C2, steps, skip_frames);
    solveWave(y_ic_sum, "datas/lin_sum.data", N, C2, steps, skip_frames);
    
    cout << " 5/b: simulating mixed normal modes \n";
    solveWave(y_mixed, "datas/mixed_modes.data", N, C2, steps, skip_frames);

    cout << " Done, data saved to datas/lin_ic1.data, datas/lin_ic2.data, datas/lin_sum.data, and datas/mixed_modes.data\n";
    return 0;
}