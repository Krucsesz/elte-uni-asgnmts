#include <iostream>
#include <fstream>

using namespace std;

struct State {
    double nR, nF;
};

State f(State s, double a, double b, double c, double d, double K, double S) {
    State ds;
    double predation = (s.nR * s.nF) / (1.0 + s.nR / S); // Telítődő ragadozás
    
    ds.nR = a * s.nR * (1.0 - s.nR / K) - b * predation;
    ds.nF = c * predation - d * s.nF;
    return ds;
}

// RK4
State rk4_step(State s, double dt, double a, double b, double c, double d, double K, double S) {
    State k1 = f(s, a, b, c, d, K, S);
    
    State s2 = {s.nR + 0.5 * dt * k1.nR, s.nF + 0.5 * dt * k1.nF};
    State k2 = f(s2, a, b, c, d, K, S);
    
    State s3 = {s.nR + 0.5 * dt * k2.nR, s.nF + 0.5 * dt * k2.nF};
    State k3 = f(s3, a, b, c, d, K, S);
    
    State s4 = {s.nR + dt * k3.nR, s.nF + dt * k3.nF};
    State k4 = f(s4, a, b, c, d, K, S);
    
    State next_s;
    next_s.nR = s.nR + (dt / 6.0) * (k1.nR + 2.0 * k2.nR + 2.0 * k3.nR + k4.nR);
    next_s.nF = s.nF + (dt / 6.0) * (k1.nF + 2.0 * k2.nF + 2.0 * k3.nF + k4.nF);
    return next_s;
}

void run_mod_lv_simulation(string filename, double a, double b, double c, double d, double K, double S, double nR0, double nF0, double t_max) {
    ofstream file(filename);
    file << "Time nR nF\n";

    State s = {nR0, nF0};
    double dt = 0.01; 

    for (double t = 0; t <= t_max; t += dt) {
        file << t << " " << s.nR << " " << s.nF << "\n";
        s = rk4_step(s, dt, a, b, c, d, K, S);
    }
    file.close();
}

int main() {
    // 1. eset: csak a nyulak kapacitása véges (K=800), telítődés nincs
    run_mod_lv_simulation("datas/mod_lv_K.data", 0.4, 0.004, 0.004, 0.9, 800.0, 1e9, 500.0, 500.0, 100.0);

    // 2. eset: csak a rókák fogyasztása telítődik (S=2500), kapacitás nincs
    run_mod_lv_simulation("datas/mod_lv_S.data", 0.4, 0.004, 0.004, 0.9, 1e9, 2500.0, 500.0, 500.0, 100.0);

    cout << "Simulations completed, data saved to datas/mod_lv_K.data and datas/mod_lv_S.data" << endl;
    return 0;
}