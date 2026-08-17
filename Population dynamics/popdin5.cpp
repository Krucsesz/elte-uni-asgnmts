#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// növény (N1), növényevő (N2), ragadozó (N3)
struct State {
    double n1, n2, n3;
};

State f(State s) {
    State ds;
    double r = 1.0, K = 100.0; // növény szaporodása és kapacitása
    double a = 0.1;            // növényevő fogyasztási rátája
    double b = 0.05, c = 0.2;  // növényevő szaporodási (b) és pusztulási (c) rátája
    double d = 0.1;            // ragadozó fogyasztási rátája
    double e = 0.02, f_p = 0.1;// ragadozó szaporodási (e) és pusztulási (f) rátája

    ds.n1 = r * s.n1 * (1.0 - s.n1 / K) - a * s.n1 * s.n2;
    ds.n2 = b * s.n1 * s.n2 - c * s.n2 - d * s.n2 * s.n3;
    ds.n3 = e * s.n2 * s.n3 - f_p * s.n3;
    
    return ds;
}

// RK4
State rk4_step(State s, double dt) {
    State k1 = f(s);
    
    State s2 = {s.n1 + 0.5 * dt * k1.n1, s.n2 + 0.5 * dt * k1.n2, s.n3 + 0.5 * dt * k1.n3};
    State k2 = f(s2);
    
    State s3 = {s.n1 + 0.5 * dt * k2.n1, s.n2 + 0.5 * dt * k2.n2, s.n3 + 0.5 * dt * k2.n3};
    State k3 = f(s3);
    
    State s4 = {s.n1 + dt * k3.n1, s.n2 + dt * k3.n2, s.n3 + dt * k3.n3};
    State k4 = f(s4);
    
    State next_s;
    next_s.n1 = s.n1 + (dt / 6.0) * (k1.n1 + 2.0 * k2.n1 + 2.0 * k3.n1 + k4.n1);
    next_s.n2 = s.n2 + (dt / 6.0) * (k1.n2 + 2.0 * k2.n2 + 2.0 * k3.n2 + k4.n2);
    next_s.n3 = s.n3 + (dt / 6.0) * (k1.n3 + 2.0 * k2.n3 + 2.0 * k3.n3 + k4.n3);
    return next_s;
}

int main() {
    ofstream file("datas/food_chain.data");
    file << "Time N1 N2 N3\n";

    State s = {150.0, 100.0, 20.0}; // kezdeti populációk
    double dt = 0.05;
    double t_max = 500.0;

    for (double t = 0; t <= t_max; t += dt) {
        file << t << " " << s.n1 << " " << s.n2 << " " << s.n3 << "\n";
        s = rk4_step(s, dt);
    }
    
    file.close();
    cout << "3 species population data saved to datas/food_chain.data\n";
    return 0;
}