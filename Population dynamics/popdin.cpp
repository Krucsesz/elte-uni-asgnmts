#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

// 1. f(x) = dx/dt = r*x*(1-x)
double f(double r, double x) {
    return r * x * (1.0 - x);
}

// 2. analitikus megoldás
double exact_solution(double r, double x0, double t) {
    if (x0 == 0.0) return 0.0;
    if (x0 == 1.0) return 1.0;
    if (r == 0.0) return x0;
    return 1.0 / (1.0 + (1.0 / x0 - 1.0) * exp(-r * t));
}

// 3. sima Euler-lépés
double euler_step(double r, double x, double dt) {
    return x + dt * f(r, x);
}

// 4. RK4 lépés
double rk4_step(double r, double x, double dt) {
    double k1 = f(r, x);
    double k2 = f(r, x + 0.5 * dt * k1);
    double k3 = f(r, x + 0.5 * dt * k2);
    double k4 = f(r, x + dt * k3);
    return x + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

// 5. ARK
double adaptive_rk4(double r, double x, double t_start, double t_end) {
    double t = t_start;
    double dt = 0.1;
    double tol = 1e-6;

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        double x1 = rk4_step(r, x, dt);
        
        double x_half = rk4_step(r, x, dt / 2.0);
        double x2 = rk4_step(r, x_half, dt / 2.0);

        double error = abs(x1 - x2);

        if (error < tol) {
            x = x2;
            t += dt;
            if (error < tol / 10.0) dt *= 1.5; 
        } else {
            dt *= 0.5;
        }
    }
    return x;
}

void run_comparison() {
    ofstream file("comparison.data");
    file << "Time Analytic Euler Adaptive_RK4\n";

    double r = 1.0;
    double x0 = 0.1;
    double t_max = 10.0;
    double dt_out = 0.5;

    double x_euler = x0;
    double x_rk = x0;

    for (double t = 0; t <= t_max; t += dt_out) {
        double x_exact = exact_solution(r, x0, t);
        
        file << t << " " << x_exact << " " << x_euler << " " << x_rk << "\n";

        x_euler = euler_step(r, x_euler, dt_out);
        x_rk = adaptive_rk4(r, x_rk, t, t + dt_out);
    }
    file.close();
    cout << "Data saved to comparison.data." << endl;
}

void run_page7_plots() {
    ofstream file("page7.data");
    file << "r x0 Time x\n";

    vector<double> r_values = {-1.0, -0.5, 0.0, 0.5, 1.0};
    double t_max = 10.0;
    double dt_out = 0.1;

    for (double r : r_values) {
        for (int i = 0; i <= 10; i++) {
            double x0 = i / 10.0;
            double x = x0;
            for (double t = 0; t <= t_max + 1e-5; t += dt_out) {
                file << r << " " << x0 << " " << t << " " << x << "\n";
                x = adaptive_rk4(r, x, t, t + dt_out);
            }
        }
    }
    file.close();
    cout << "Data saved to page7.data." << endl;
}

int main() {
    run_comparison();
    run_page7_plots();
    cout << "All simulations completed!" << endl;
    return 0;
}