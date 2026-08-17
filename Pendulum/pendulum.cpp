#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "vector.hpp"        // vectors with components of type double
#include "odeint.hpp"        // ODE integration routines, Runge-Kutta ...
using namespace cpl;

const double pi = 4 * atan(1.0);

const double g = 9.8;        // acceleration of gravity

double L = 1.0;              // length of pendulum
double q = 0.5;              // damping coefficient
double Omega_D = 2.0/3.0;    // frequency of driving force
double F_D = 0.9;            // amplitude of driving force
bool nonlinear;              // linear if false

Vector f(const Vector& x) {  // extended derivative vector
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    Vector f(3);             // Vector with 3 components
    f[0] = 1;

    f[1] = omega;
    if (nonlinear)
        f[2] = - (g/L) * sin(theta) - q * omega + F_D * sin(Omega_D * t);
    else
        f[2] = - (g/L) * theta - q * omega + F_D * sin(Omega_D * t);
    return f;
}

// sima Euler-módszer
void EulerStep(Vector& x, double tau, Vector derivs(const Vector&)) {
    x += tau * derivs(x);
}

// Euler-Cromer-módszer
void EulerCromerStep(Vector& x, double tau, Vector derivs(const Vector&)) {
    Vector d = derivs(x);
    x[0] += tau;
    x[2] += tau * d[2];
    x[1] += tau * x[2];
}

void runSinglePendulum() {
        cout << " Nonlinear damped driven pendulum\n"
         << " --------------------------------\n"
         << " Enter linear or nonlinear: ";
    string response;
    cin >> response;
    nonlinear = (response[0] == 'n');
    cout<< " Length of pendulum L: ";
    cin >> L;
    cout<< " Enter damping coefficient q: ";
    cin >> q;
    cout << " Enter driving frequencey Omega_D: ";
    cin >> Omega_D;
    cout << " Enter driving amplitude F_D: ";
    cin >> F_D;
    cout << " Enter theta(0) and omega(0): ";
    double theta, omega, tMax;
    cin >> theta >> omega;
    cout << " Enter integration time t_max: ";
    cin >> tMax;

    // kiválaszthatóvá tesszük a numerikus integrációs módszert
    int solver;
    cout << "\n Select Integration Method:\n"
         << " 1: Euler\n"
         << " 2: Euler-Cromer\n"
         << " 3: Runge-Kutta (RK4)\n"
         << " 4: Adaptive Runge-Kutta\n"
         << " Enter choice (1-4): ";
    cin >> solver;

    double dt = 0.05;
    double accuracy = 1e-6;
    ofstream dataFile("datas/pendulum.data");

    double t = 0;
    Vector x(3);
    x[0] = t;
    x[1] = theta;
    x[2] = omega;

    while (t < tMax) {
        if (solver == 1) {
            EulerStep(x, dt, f);
        } else if (solver == 2) {
            EulerCromerStep(x, dt, f);
        } else if (solver == 3) {
            RK4Step(x, dt, f);
        } else {
            // ARK4 lépésméret választása
            adaptiveRK4Step(x, dt, accuracy, f);
        }

        t = x[0], theta = x[1], omega = x[2];
        
        // a szöget [-pi, pi) intervallumba hozzuk, hogy a grafikon ne legyen túl kusza a nemlineáris esetben
        if (nonlinear) {
            while (theta >= pi) theta -= 2 * pi;
            while (theta < -pi) theta += 2 * pi;
        }

        double E;
        if (nonlinear) {
            E = 0.5 * L * L * omega * omega + g * L * (1.0 - cos(theta));
        } else {
            E = 0.5 * L * L * omega * omega + 0.5 * g * L * theta * theta;
        }

        // idő, szög, szögsebesség, energia
        dataFile << t << '\t' << theta << '\t' << omega << '\t' << E << '\n';
    }

    cout << " Output data to file datas/pendulum.data" << endl;
    dataFile.close();
}



double m1 = 1.0, m2 = 1.0;
double L1 = 1.0, L2 = 1.0;

Vector f_double(const Vector& x) {  
    double t  = x[0];
    double th1 = x[1];
    double th2 = x[2];
    double w1  = x[3];
    double w2  = x[4];

    Vector f(5);
    f[0] = 1.0;
    f[1] = w1;
    f[2] = w2;

    double delta = th1 - th2;
    double den1 = 2*m1 + m2 - m2*cos(2*delta);
    
    f[3] = (-g*(2*m1+m2)*sin(th1) - m2*g*sin(th1-2*th2) - 2*sin(delta)*m2*(w2*w2*L2 + w1*w1*L1*cos(delta))) / (L1 * den1);
    f[4] = (2*sin(delta)*(w1*w1*L1*(m1+m2) + g*(m1+m2)*cos(th1) + w2*w2*L2*m2*cos(delta))) / (L2 * den1);

    return f;
}

void runDoublePendulum() {
    double th1, th2, w1, w2, tMax;
    cout << " Enter initial th1 and th2 (rad): ";
    cin >> th1 >> th2;
    cout << " Enter initial w1 and w2 (rad/s): ";
    cin >> w1 >> w2;
    cout << " Enter integration time t_max: ";
    cin >> tMax;

    double dt = 0.05;
    double accuracy = 1e-6;
    ofstream dataFile("datas/double_pendulum.data");

    double t = 0;
    Vector x(5);
    x[0]=t; x[1]=th1; x[2]=th2; x[3]=w1; x[4]=w2;

    while (t < tMax) {
        adaptiveRK4Step(x, dt, accuracy, f_double);
        t = x[0];

        double T = 0.5*m1*L1*L1*x[3]*x[3] + 0.5*m2*(L1*L1*x[3]*x[3] + L2*L2*x[4]*x[4] + 2*L1*L2*x[3]*x[4]*cos(x[1]-x[2]));
        double V = -(m1+m2)*g*L1*cos(x[1]) - m2*g*L2*cos(x[2]);
        double E = T + V;

        dataFile << t << '\t' << x[1] << '\t' << x[2] << '\t' << x[3] << '\t' << x[4] << '\t' << E << '\n';
    }
    cout << " Output data to file datas/double_pendulum.data\n";
    dataFile.close();

}

int main() {
    int choice;
    cout << " Select system (1: Single Pendulum, 2: Double Pendulum): ";
    cin >> choice;

    if (choice == 1) {
        runSinglePendulum();
    } else if (choice == 2) {
        runDoublePendulum();
    } else {
        cout << " Invalid choice.\n";
    }
    
    return 0;
}




