#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
using namespace std;

#include "vector.hpp"
using namespace cpl;

void runStabilityMap() {
    cout << " Stability map generation..." << endl;
    double L = 1.0, rho = 0.01, T = 10.0;
    double c = sqrt(T / rho);
    ofstream file("stability_map.data");

    int resolution = 50;
    // dx és dt rács végigpásztázása
    for (int i = 0; i < resolution; i++) {
        double dx = 0.005 + i * (0.05 - 0.005) / (resolution - 1);
        int N = round(L / dx);
        dx = L / N; // újraszámoljuk, hogy pontosan illeszkedjen a hosszra

        for (int j = 0; j < resolution; j++) {
            double dt = 0.0001 + j * (0.002 - 0.0001) / (resolution - 1);
            
            double C2 = (c * dt / dx) * (c * dt / dx);
            
            Vector y_old(N + 1), y_curr(N + 1), y_new(N + 1);
            
            // kf.
            const double pi = 4 * atan(1.0);
            for (int k = 0; k <= N; k++) {
                y_curr[k] = sin(pi * k * dx / L);
                y_old[k] = y_curr[k];
            }
            
            bool stable = true;
            // időlépések szimulációja
            for (int step = 0; step < 200; step++) {
                for (int k = 1; k < N; k++) {
                    y_new[k] = 2.0*y_curr[k] - y_old[k] + C2*(y_curr[k+1] - 2.0*y_curr[k] + y_curr[k-1]);
                }
                y_new[0] = 0; y_new[N] = 0;
                y_old = y_curr;
                y_curr = y_new;
                
                // ha a középső pont amplitúdója > 10, akkor a rendszer felrobbant
                if (abs(y_curr[N/2]) > 10.0 || std::isnan(y_curr[N/2])) {
                    stable = false;
                    break;
                }
            }
            // dx, dt, és hogy stabil volt-e (1 vagy 0)
            file << dx << "\t" << dt << "\t" << (stable ? 1 : 0) << "\n";
        }
    }
    file.close();
    cout << " Done, stability map saved to stability_map.data\n";
}

int main() {
    int choice;
    cout << " STRING SIMULATION\n"
         << " 1: Plucked string\n"
         << " 2: Normal mode \n"
         << " 3: Stability map generation\n"
         << " Choose initial condition: ";
    cin >> choice;

    if (choice == 3) {
        runStabilityMap();
        return 0;
    }
    
    double L = 1.0;
    double rho = 0.01;
    double T = 10.0;
    double c = sqrt(T / rho);
    
    int N = 100;
    double dx = L / N;
    double dt = 1.0 * (dx / c);
    double C2 = (c * dt / dx) * (c * dt / dx); 

    double tMax = 0.2;         
    int steps = tMax / dt;
    int skip_frames = 5;

    Vector y_old(N + 1);
    Vector y_curr(N + 1);
    Vector y_new(N + 1);

    string filename = (choice == 1) ? "string_plucked.data" : "string_normal.data";
    ofstream dataFile(filename);

    const double pi = 4 * atan(1.0);

    // 1. kezdeti állapot beállítása
    for (int i = 0; i <= N; i++) {
        double x = i * dx;
        if (choice == 1) {
            // 1. feladat
            double x0 = 0.2 * L;
            if (x <= x0) y_curr[i] = x / x0;
            else         y_curr[i] = (L - x) / (L - x0);
        } else {
            // 2. feladat
            int n_mode = 2;
            y_curr[i] = sin(n_mode * pi * x / L);
        }
        y_old[i] = y_curr[i]; 
    }

    y_curr[0] = 0.0; y_curr[N] = 0.0;

    // kezdeti állapot mentése fájlba
    for (int i = 0; i <= N; i++) dataFile << y_curr[i] << (i == N ? "" : "\t");
    dataFile << "\n";

    // 2. első időlépés
    for (int i = 1; i < N; i++) {
        y_new[i] = y_curr[i] + 0.5 * C2 * (y_curr[i+1] - 2*y_curr[i] + y_curr[i-1]);
    }
    y_new[0] = 0.0; y_new[N] = 0.0;

    y_old = y_curr;
    y_curr = y_new;

    // 3. fő időlépések
    for (int n = 1; n < steps; n++) {
        for (int i = 1; i < N; i++) {
            y_new[i] = 2.0 * y_curr[i] - y_old[i] + C2 * (y_curr[i+1] - 2.0 * y_curr[i] + y_curr[i-1]);
        }
        
        y_new[0] = 0.0; y_new[N] = 0.0;

        if (n % skip_frames == 0) {
            for (int i = 0; i <= N; i++) dataFile << y_new[i] << (i == N ? "" : "\t");
            dataFile << "\n";
        }

        y_old = y_curr;
        y_curr = y_new;
    }

    dataFile.close();
    cout << " Ready, data saved to " << filename << "\n";
    return 0;
}