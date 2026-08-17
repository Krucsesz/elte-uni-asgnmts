#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// nyulak (nR) és rókák (nF)
struct State {
    double nR, nF;
};

State f(State s, double a, double b, double c, double d) {
    State ds;
    ds.nR = a * s.nR - b * s.nR * s.nF;
    ds.nF = c * s.nR * s.nF - d * s.nF;
    return ds;
}

// RK4
State rk4_step(State s, double dt, double a, double b, double c, double d) {
    State k1 = f(s, a, b, c, d);
    
    State s2 = {s.nR + 0.5 * dt * k1.nR, s.nF + 0.5 * dt * k1.nF};
    State k2 = f(s2, a, b, c, d);
    
    State s3 = {s.nR + 0.5 * dt * k2.nR, s.nF + 0.5 * dt * k2.nF};
    State k3 = f(s3, a, b, c, d);
    
    State s4 = {s.nR + dt * k3.nR, s.nF + dt * k3.nF};
    State k4 = f(s4, a, b, c, d);
    
    State next_s;
    next_s.nR = s.nR + (dt / 6.0) * (k1.nR + 2.0 * k2.nR + 2.0 * k3.nR + k4.nR);
    next_s.nF = s.nF + (dt / 6.0) * (k1.nF + 2.0 * k2.nF + 2.0 * k3.nF + k4.nF);
    return next_s;
}

void run_lv_simulation(string filename, double a, double b, double c, double d, double nR0, double nF0) {
    ofstream file(filename);
    file << "Time nR nF\n";

    State s = {nR0, nF0};
    double dt = 0.01;
    double t_max = 100.0;

    for (double t = 0; t <= t_max; t += dt) {
        file << t << " " << s.nR << " " << s.nF << "\n";
        s = rk4_step(s, dt, a, b, c, d);
    }
    file.close();
}

int main() {
    // a=0.4, b=0.001, c=0.001, d=0.9
    // kf.: Nyúl=500, Róka=500
    run_lv_simulation("datas/lv_set1.data", 0.4, 0.001, 0.001, 0.9, 500.0, 500.0);

    // a=0.4, b=0.004, c=0.004, d=0.9
    // kf.: Nyúl=400, Róka=400
    run_lv_simulation("datas/lv_set2.data", 0.4, 0.004, 0.004, 0.9, 400.0, 400.0);

    cout << "Lotka-Volterra simulations completed and saved to datas/lv_set1.data and datas/lv_set2.data." << endl;
    return 0;
}