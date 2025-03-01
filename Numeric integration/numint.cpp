#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Integrálandó példafüggvény
double f(double x) {
    return cos(x) * exp(-x * x);
}

// Kompozit trapéz módszer
double integrate(int n, double x0, double x1) { // n: a részek száma, x0 és x1: az integrálási intervallum
    double h = (x1 - x0) / n; // lépésköz
    double sum = 0.5 * (f(x0) + f(x1)); // kezdő és záró értékek fele

    // A közbenső értékek összege
    for (int i = 1; i < n; ++i) {
        double x = x0 + i * h;
        sum += f(x); // összegzés
    }

    return sum * h; // végső integrál érték
}

int main() {
    double x0 = -1.0;
    double x1 = 3.0;

    // 16 tizedesjegy pontosság
    cout << fixed << setprecision(16);

    // Különböző n értékek a konvergencia ellenőrzéséhez
    int n_values[] = {10, 100, 1000, 10000, 100000};

    // Az integrál eredményének kiszámítása és kiírása a különböző n értékeinkre
    for (int n : n_values) {
        double result = integrate(n, x0, x1); // integrál érték
        cout << "n = " << n << ", az integrál eredménye: " << result << endl;
    }

    return 0;
}