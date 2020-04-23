#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>

double test_f1(double x){
    // a = 3/2pi b = 2pi
    return cos(x) * cosh(x) - 1;
}

double test_f2(double x ){

    // a = 0, b = pi/2
    return (1/x) - tan(x);

}

double test_f3(double x){

    //a = 1, b = 3
    return pow(2, -x) + pow(M_E, x) + 2*cos(x) - 6;
}

double my_abs(double x){
    if(x > 0){
        return x;
    }else{
        return -x;
    }
}

double derivative(double (*func)(double x), double x0){
    double h = 1.0e-10;
    return (func(x0+h) - func(x0))/h;
}

double zad1(double a, double b, double eps, double (func)(double x), int &iter_num){    
    
    iter_num = 0;

    while( my_abs(a-b) > eps){
        iter_num++;
        // std::cout << "iteration: " << iter_num << std::endl;
        double x1 = (a+b) / 2;

        if( my_abs(func(x1)) <= eps){
            break;
        }else if( func(x1) * func(a) < 0){
            b = x1;
        }else{
            a = x1;
        }
    }
    return ((a+b) / 2);
}

double zad2(double a, double b, double eps, double(*func)(double x), int &iter_done, int iter_todo){
    double x1;
    double x2;

    x1 = a;
    x2 = b;
    iter_done = 0;
    while(my_abs(func(x2)) > eps && my_abs(x2 - x1) >= eps && iter_done < iter_todo){
        x1 = x2;
        x2 = x1 - func(x1)/derivative(func, x1);
        iter_done++;
    }
    return x2;

}

double zad3(double a, double b, double eps, double(*func)(double x), int &iter_done, int iter_todo){
    double x1;
    double x2;
    double x3;

    x1 = a;
    x2 = b;
    x3 = a;
    iter_done = 0;
    while(my_abs(func(x3)) >= eps && my_abs(x3 - x2) >= eps && iter_done < iter_todo){
        x1 = x2;
        x2 = x3;
        x3 = (func(x2)*x1 - func(x1)*x2) / (func(x2) - func(x1));
        iter_done++;
    }
    return x3;

}
int main(){

    int mem;
    std::string headers[3] = {"10^-7", "10^-15", "10^-33"};
    double eps[3] = { 1.0e-7, 1.0e-15, 1.0e-33};
    double (*f_pointers[3])(double x);
    f_pointers[0] = test_f1;
    f_pointers[1] = test_f2;
    f_pointers[2] = test_f3;
    double a_b[3][2] = {{(3.0/2.0)*M_PI,  2*M_PI}, 
                        { 0, M_PI/2.0},
                        { 1.0, 3.0}};

    std::cout << std::setprecision(10);

    for(int i = 0; i < 3; i++){

        std::cout << "FUNCTON: " << i << std::endl << std::endl;

        for(int j = 0; j < 2; j++){
            int iter_num;
            std::cout << "epsilon: " << headers[j] << std::endl;
            std::cout << "result: " << zad1(a_b[i][0], a_b[i][1], eps[j], f_pointers[i], iter_num) << " steps: " << iter_num << std::endl;
            std::cout << std::endl;
        }
    }

    std::cout << "-----------NEWTON METHOD--------------" << std::endl;
    for(int i = 0; i < 3; i++){

        std::cout << "FUNCTION: " << i << std::endl << std::endl;

        for(int j = 0; j < 2; j++){
            int iter_num;
            std::cout << "epsilon: " << headers[j] << std::endl;
            std::cout << "result: " << zad2(a_b[i][0], a_b[i][1], eps[j], f_pointers[i], iter_num, 1000) << " steps: " <<  iter_num << std::endl;
            std::cout << std::endl;
        }
    }
    
    std::cout << "-----------SECANT METHOD--------------" << std::endl;
    for(int i = 0; i < 3; i++){

        std::cout << "FUNCTION: " << i << std::endl << std::endl;

        for(int j = 0; j < 2; j++){
            int iter_num;
            std::cout << "epsilon: " << headers[j] << std::endl;
            std::cout << "result: " << zad3(a_b[i][0], a_b[i][1], eps[j], f_pointers[i], iter_num, 1000) << " steps: " <<  iter_num << std::endl;
            std::cout << std::endl;
        }
    }
}