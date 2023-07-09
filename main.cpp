#include <iostream>
#include <cmath>
#include <fstream>
#define PI 3.14
using namespace std;

double Euler_with_const_B (double x0, double y0, double Vy, double Vx, double dt, double t1, double B, double m, double q, double z0, double Vz) {
    ofstream out("Euler_with_const_B.dat");
    for(double t0 = 0; t0 < t1; t0+=dt){
        out << x0 << "\t" << y0 << "\t" << z0 << endl;
        x0 = x0 + Vx * dt;
        y0 = y0 + Vy * dt;
        z0 = z0 + Vz * dt;
        double Vynew = Vy - q * B/m * Vx * dt;
        Vx = Vx + q * B/m * Vy * dt;
        Vy = Vynew;
    }
    out.close();
    return 0;
}

double Euler_without_const_B (double x0, double y0, double Vy, double Vx, double dt, double t1, double m, double q, double z0, double Vz, double I) {
    ofstream out("Euler_withOUT_const_B.dat");
    for(double t0 = 0; t0 < t1; t0+=dt){
        out << x0 << "\t" << y0 << "\t" << z0 << endl;
        x0 = x0 + Vx * dt;
        y0 = y0 + Vy * dt;
        z0 = z0 + Vz * dt;
        double Vynew = Vy - 2 * q * Vx * I/(z0 * m) * dt;
        double Vxnew = Vx + 2 * Vy * q * I/(m * z0) * dt;
        Vx = Vxnew;
        Vy = Vynew;
    }
    out.close();
    return 0;
}

int main() {
    double x0 = 0, y0 = 0, z0 = 0, Vz = 10e3, V0, alpha, t1=10e-11, dt = 10e-16, B, m = 9.11e-31, Vy, Vx, q = -1.6e-19, I;
    //double x0 = 0, y0 = 0, z0 = 0, Vz = 10e3, V0 = 2e4, alpha = 47.2 * 180 / PI, t1=10e-12, dt = 10e-18, B = 2e4, m = 3.7e-27, q = 1.6e-19, I = 1e-3; //proton
    //double x0 = 0, y0 = 0, z0 = 0, Vz = 10e3, V0 = 21e4, alpha = 30 * 180 / PI, t1=10e-11, dt = 10e-16, B = 2, m = 9.11e-31, q = -1.6e-19, I = 1e-6; //electron
    cout << "Please, enter start speed" << endl;
    cin >> V0;
    cout << "Please, enter value of the angle between magnetic field and speed" << endl;
    cin >> alpha;
    alpha = alpha * 180 / PI;
    Vy = V0*sin(alpha); Vx = V0*cos(alpha);
    cout << "Please, enter B" << endl;
    cin >> B;
    cout << "Please, enter I" << endl;
    cin >> I;
    Euler_with_const_B (x0, y0, Vy, Vx, dt, t1, B, m, q, z0, Vz);
    Euler_without_const_B (x0, y0, Vy, Vx, dt, t1, m, q, z0, Vz, I);
}




