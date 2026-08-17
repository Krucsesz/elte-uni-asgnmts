#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct State {
    double n1, n2;
};

// diffegyenlet rendszer a két fajra
State f(State s, double r1, double r2, double k1, double k2, double alpha, double beta) {
    State ds;
    ds.n1 = r1 * s.n1 * (1.0 - (s.n1 + alpha * s.n2) / k1);
    ds.n2 = r2 * s.n2 * (1.0 - (s.n2 + beta * s.n1) / k2);
    return ds;
}

// RK4
State rk4_step(State s, double dt, double r1, double r2, double k1, double k2, double alpha, double beta) {
    State k1_s = f(s, r1, r2, k1, k2, alpha, beta);
    
    State s2 = {s.n1 + 0.5 * dt * k1_s.n1, s.n2 + 0.5 * dt * k1_s.n2};
    State k2_s = f(s2, r1, r2, k1, k2, alpha, beta);
    
    State s3 = {s.n1 + 0.5 * dt * k2_s.n1, s.n2 + 0.5 * dt * k2_s.n2};
    State k3_s = f(s3, r1, r2, k1, k2, alpha, beta);
    
    State s4 = {s.n1 + dt * k3_s.n1, s.n2 + dt * k3_s.n2};
    State k4_s = f(s4, r1, r2, k1, k2, alpha, beta);
    
    State next_s;
    next_s.n1 = s.n1 + (dt / 6.0) * (k1_s.n1 + 2.0 * k2_s.n1 + 2.0 * k3_s.n1 + k4_s.n1);
    next_s.n2 = s.n2 + (dt / 6.0) * (k1_s.n2 + 2.0 * k2_s.n2 + 2.0 * k3_s.n2 + k4_s.n2);
    return next_s;
}

// 1. kompetitív kizárás
void run_exclusion() {
    ofstream file("datas/exclusion.data");
    file << "Time n1 n2\n";

    double r1 = 1.0, r2 = 1.0;
    double k1 = 100.0, k2 = 150.0;
    double alpha = 1.0, beta = 1.0;
    
    State s = {10.0, 10.0};
    double dt = 0.1;
    double t_max = 50.0;

    for (double t = 0; t <= t_max; t += dt) {
        file << t << " " << s.n1 << " " << s.n2 << "\n";
        s = rk4_step(s, dt, r1, r2, k1, k2, alpha, beta);
    }
    file.close();
    cout << "Competitive exclusion data saved to datas/exclusion.data." << endl;
}

// 2. fázistér vizsgálat
void run_phase_space() {
    ofstream file("datas/coexistence.data");
    file << "X Y Survivors\n";

    double r1 = 1.0, r2 = 1.0;
    double k1 = 100.0, k2 = 100.0; // kapacitások

    for (double X = -50.0; X <= 50.0; X += 2.0) {
        for (double Y = -50.0; Y <= 50.0; Y += 2.0) {
            double alpha = (X + k1) / k2;
            double beta = (Y + k2) / k1;

            State s = {50.0, 50.0}; // esélyek
            double dt = 0.5;
            double t_max = 200.0;

            for (double t = 0; t <= t_max; t += dt) {
                s = rk4_step(s, dt, r1, r2, k1, k2, alpha, beta);
            }

            int survivors = 0;
            if (s.n1 > 1.0) survivors++;
            if (s.n2 > 1.0) survivors++;

            file << X << " " << Y << " " << survivors << "\n";
        }
    }
    file.close();
    cout << "Coexistence phase space data saved to datas/coexistence.data." << endl;
}

int main() {
    run_exclusion();
    run_phase_space();
    cout << "The simulations for the 2nd problem have completed." << endl;
    return 0;
}