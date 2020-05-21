#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include "solve.cpp"
#include "trigo.cpp"
#define K 2
#define M 3






double f1(double x){
    double k = K;
    double m = M;
    return sin(k*x/M_PI)*exp(-m*x/M_PI);
}

double f2(double x){
    double k = K;
    double m = M;
    return exp(k*cos(m*x));
}

double approx_func(double x, int n){
    return pow(x, (double) n);
}

void func_to_file(double (*f)(double), double a, double b, double n, std::string name){

    double step = my_abs(b-a)/n;
    std::ofstream file;
    file.open(name);
    for(double i = a; i < b; i +=step){

        file << i << " " << f(i) << std::endl;

    }
    file.close();
}



void points_to_file(arr2d points , std::string name){

    std::ofstream file;
    file.open(name  + ".txt");
    int n = points.size();
    for(int i = 0; i < n; i++){

        file << points[i][0] << " " << points[i][1] << std::endl;

    }
    file.close();

}

arr2d create_A(arr2d &points, double(*appr)(double x, int n)){
    
    arr2d res;
    int n = points.size();
    res.resize(n);

    for(int i = 0; i < n; i++){

        res[i].resize(n);

        for(int j = 0; j < n; j++){

            double sum = 0;
            for(int s = 0; s < n; s++) sum += appr(points[s][0], i+j);

            res[i][j] = sum;

        }

    }
    return res;
}

std::vector<double> create_B(arr2d &points, double(*appr)(double x, int n)){
    std::vector<double> res;
    
    int n = points.size();
    res.resize(n);

    for(int i = 0; i < n; i++){

        double sum = 0;
        for(int s = 0; s < n; s++) sum += points[s][1]*appr(points[s][0], i);
        res[i] = sum;

    }

    return res;

}

arr2d generatePoints(double (*f)(double x), double a, double b, double n){

    arr2d res;
    double step = my_abs(b-a)/n;


    int it = 0;
    for(double i = a; i < b; i +=step){
        std::vector<double> a;
        a.resize(2);
        a[0] = i;
        a[1] = f(i);
        res.push_back(a);
        it++;
    }
    return res;
}

double result_func(std::vector<double> A, double (*appr)(double x, int n), double x){
    double sum = 0;
    for(int i = 0; i < (int) A.size(); i++){
        sum += A[i]*appr(x,i);
    }
    return sum;
}


arr2d generatePoints2(std::vector<double> &A, double a, double b, double n){

    arr2d res;
    double step = my_abs(b-a)/n;


    int it = 0;
    for(double i = a; i < b; i +=step){
        std::vector<double> a;
        a.resize(2);
        a[0] = i;
        a[1] = result_func(A, approx_func, i);
        res.push_back(a);
        it++;
    }
    return res;
}

void print(std::vector<double> x){
        for(int i = 0; i < (int)x.size(); i++){
        std::cout << x[i] << "  ";
    }
    std::cout << std::endl;
}


std::vector<double> solve_coeff(arr2d &points){

    
    arr2d A = create_A(points, approx_func);
    std::vector<double> B = create_B(points, approx_func);

    std::vector<int> P(A.size()+1,0);
    std::vector<double> X(A.size(),0);
    
    LUPDecompose(A,A.size(), 1.0e-7, P);

    LUPSolve(A, P, B, A.size(), X);
    return X;
}

double approx_error_mean(arr2d &points1, arr2d &points2){
    double err = 0;
    
    for(int i = 0; i < points1.size(); i++){
        err += my_abs(points1[i][1] - points2[i][1]);
    }
    return err/(double)points1.size();
}
double approx_error_max(arr2d &points1, arr2d &points2){
    double err = 0;
    double max = 0;
    for(int i = 0; i < points1.size(); i++){
        err = my_abs(points1[i][1] - points2[i][1]);
        if( err > max) max = err;
    }
    return max;
}
int main(){
    
    // double n = 100;
    // double a = -10;
    // double b = 10;
    
    // double step = my_abs(b-a)/n;



    // arr2d A({{5,3,2}, {1, 2, 0}, {3, 0, 4}});
    // std::vector<double> B({10, 5, -2});

    std::ofstream raportf;
    raportf.open("raport_f1.txt");

    for(int i = 2; i < 10; i+=1){

        arr2d discrete_points = generatePoints(f1, 0, 5, i);
        std::vector<double> A = solve_coeff(discrete_points);

        
        arr2d approx_points = generatePoints2(A, 0, 6, i);
        

        arr2d graphing_points = generatePoints2(A, 0, 5, 200);
        arr2d error_points = generatePoints(f1, 0, 5, 200);
        raportf << "Aproksymacja funkcji f1, liczba punktów dyskretyzacji: " << i << " błąd średni: " 
        << approx_error_mean(graphing_points, error_points) << " maksymalny błąd: "
         << approx_error_max(graphing_points, error_points) << std::endl;

        points_to_file(graphing_points, "f1_approx" + std::to_string(i));

    }
    //  for(int i = 2; i < 30; i+=4){

    //     arr2d discrete_points = generatePoints(f1, 0, 5, i);

    //     // arr2d approx_points = generatePoints3(discrete_points, 0, 6, i);
        
    //     // raportf << "Aproksymacja funkcji f2, liczba punktów dyskretyzacji: " << i << " błąd: " << approx_error(discrete_points, approx_points) << std::endl;

    //     arr2d graphing_points = generatePoints3(discrete_points, 0, 5, 200, i);
    //     // points_to_file(graphing_points, "f2_approx_trig" + std::to_string(i));
    //     arr2d error_points = generatePoints(f1, 0, 5, 200);
    //     raportf << "Aproksymacja funkcji f1, liczba punktów dyskretyzacji: " << i << " błąd średni: " 
    //     << approx_error_mean(graphing_points, error_points) << " maksymalny błąd: "
    //      << approx_error_max(graphing_points, error_points) << std::endl;

    // }   
    raportf.close();



    // func_to_file(f1, 0, 5, 500, "f1.txt");
    // func_to_file(f2, 0, 5, 500, "f2.txt");

}