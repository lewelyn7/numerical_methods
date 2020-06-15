#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#define DERIVATIVE_STEP 0.00001
using namespace std::placeholders;
using namespace std;
double derivative(std::function<double(double)> f, double x){
    return (f(x + DERIVATIVE_STEP) - f(x))/2;
}

double ftest(double x){
    return sin(x);
}

vector<double> generate_mesh(double a, double b, int n){
    double step = (b - a)/n;
    vector<double> result(n);
    for(int i = 0; i < n; i++){
        result[i] = a + i*step;
    }
    return result;
}

vector<double> Euler_method(vector<double> mesh, function<double(double, double)> ftu, double u0){
    int n = mesh.size();
    double step = mesh[1] - mesh[0];
    vector<double> result(n);
    result[0] = u0;
    for(int i = 1; i < n; i++){
        result[i] = result[i-1] + ftu(mesh[i-1], result[i-1])*step;
    }
    return result;
}
vector<double> backward_Euler_method(vector<double> mesh, function<double(double, double)> ftu, double u0){
    int n = mesh.size();
    double step = mesh[1] - mesh[0];
    vector<double> result(n);
    result[0] = u0;
    for(int i = 1; i < n; i++){
        result[i] = result[i-1];
        result[i] = result[i-1] + ftu(mesh[i], result[i])*step;
    }
    return result;
}

vector<double> Rung_Kuta_2(vector<double> mesh, function<double(double, double)> ftu, double u0){
    int n = mesh.size();
    double step = mesh[1] - mesh[0];
    vector<double> result(n);
    double c1 = 0;
    double c2 = 1;
    auto k1 = [ftu](double t, double u, double h){ return ftu(t, u);};
    auto k2 = [ftu](double t, double u, double h){ return ftu(t + 0.5*h, u + 0.5*h*ftu(t, u));};

    result[0] = u0;
    for(int i = 1; i < n; i++){
        result[i] = result[i-1] + step*(c1*k1(mesh[i-1], result[i-1], step) + c2*k2(mesh[i-1], result[i-1], step));
    }
    return result;
}
vector<double> Rung_Kuta_4(vector<double> mesh, function<double(double, double)> ftu, double u0){
    int n = mesh.size();
    double step = mesh[1] - mesh[0];
    vector<double> result(n);
    double c1 = 1;
    double c2 = 2;
    double c3 = 2;
    double c4 = 1;
    auto k1 = [ftu](double t, double u, double h){ return ftu(t, u);};
    auto k2 = [ftu, k1](double t, double u, double h){ return ftu(t + 0.5*h, u + 0.5*h*k1(t, u, h));};
    auto k3 = [ftu, k2](double t, double u, double h){ return ftu(t + 0.5*h, u + 0.5*h*k2(t, u, h));};
    auto k4 = [ftu, k3](double t, double u, double h){ return ftu(t + h, u + h*k3(t, u, h));};

    result[0] = u0;
    for(int i = 1; i < n; i++){
        result[i] = result[i - 1] + step / 6.0 * (c1 * k1(mesh[i - 1], result[i - 1], step) + c2 * k2(mesh[i - 1], result[i - 1], step) + c3 * k3(mesh[i - 1], result[i - 1], step) + c4 * k4(mesh[i - 1], result[i - 1], step));
    }
    return result;
}
vector<double> calc_values(vector<double> mesh, function<double(double)> f){
    int n = mesh.size();
    vector<double> result(n);
    for(int i = 0; i < n; i++){
        result[i] = f(mesh[i]);
    }
    return result;
}
void print_vec(vector<double> vec){
    for(int i = 0; i < vec.size(); i++){
        cout << vec[i] << endl;
    }
}
double ftu_test(double t, double u){
    return cos(t);
}

double ftukm(double t, double u, double k, double m){
    return k*k*m*sin(m*t)*cos(m*t) + k*m*u*sin(m*t);
}

double exact_ftukm(double k, double m , double x){
    return exp(-1*k*cos(m*x)) - k*cos(m*x) + 1;
}
void save_result(vector<double> mesh, vector<double> res, string filename){
    ofstream file;
    file.open(filename, ios::out);
    int n = res.size();
    for(int i = 0; i < n; i++){
        file << mesh[i] << " " << res[i] << endl;
    }
    file.close();

}

void test(){
    double k = 1.0;
    double m = 1.0;
    double x0 = 0;
    double xk = 18.0;

    auto ftu = bind(ftukm, _1, _2, k ,m);
    auto exact_sol = bind(exact_ftukm, k, m, _1);

    double a = exact_sol(x0);


    vector<int> n = {10, 20, 40, 50, 100, 300, 500, 1000};
    for(int i = 0; i < n.size(); i++){
        vector<double> mesh = generate_mesh(x0, xk, n[i]);
        vector<double> resEuler = Euler_method(mesh, ftu, a);
        vector<double> resRG = Rung_Kuta_4(mesh, ftu, a);
        save_result(mesh, resEuler, "euler_" + to_string(n[i]) + ".txt");
        save_result(mesh, resRG, "rung_kuta_" + to_string(n[i]) + ".txt");
    }

    vector<double> mesh = generate_mesh(x0, xk, 200);
    vector<double> res_exact(mesh.size());
    for(int i = 0; i < mesh.size(); i++) res_exact[i] = exact_sol(mesh[i]);
    save_result(mesh, res_exact, "exact_values.txt");
    
}
int main(){
    // vector<double> mesh = generate_mesh(0,6.28, 100);
    // vector<double> res = Rung_Kuta_4(mesh, ftu_test, 0);
    // print_vec(res);
    // save_result(mesh, res, "wynik.txt");
    test();
}