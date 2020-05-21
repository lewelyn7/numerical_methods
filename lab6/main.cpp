#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include "solve.cpp"

double my_abs(double x){
    if(x > 0){
        return x;
    }else{
        return x*-1.0;
    }
}
void LUD(arr2d &A, arr2d &L, arr2d &U, arr2d &D){
    int n = A.size();

    for(int i =0; i < n; i++){
        D[i][i] = A[i][i];
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < i; j++){
            L[i][j] = A[i][j];
        }
    }

    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; j++){
            U[i][j] = A[i][j];
        }
    }
}

void zero_arr2d( arr2d &A){
    int n = A.size();

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            A[i].resize(n);
            A[i][j] = 0;


        }
    }
}

void print_matrix(arr2d &A){
    int n = A.size();

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void print_matrix1d(arr1d &A){
    int n = A.size();

    for(int j = 0; j < n; j++){
        std::cout << A[j] << std::endl;
    }
}

arr1d Jacobi_solve(arr2d &A, arr1d &B, int k){

    int n = A.size();

    arr1d X1(n, 0);
    arr1d X2(n, 0);

    for(int i = 0; i < k; i++){
        for(int j = 0; j < n; j++){
            double sum = 0;
            for(int s = 0; s < n; s++){
                if( s != j){
                    sum += A[j][s]*X1[s];
                }
            }

            X2[j] = (1.0/A[j][j])*(B[j] - sum); 
        }
        for(int j = 0; j < n; j++){
            X1[j] = X2[j];
        }
    }


    return X2;
}

arr1d Gauss_Seidel_solve(arr2d &A, arr1d &B, int k){

    int n = A.size();


    arr1d X1(n, 0);
    arr1d X2(n, 0);

    for(int i = 0; i < k; i++){
        for(int j = 0; j < n; j++){
            double sum1 = 0;
            for(int s = 0; s < j; s++){
                sum1 += A[j][s]*X2[s];
            }
            double sum2 = 0;
            for(int s = j+1; s < n; s++){
                sum2 += A[j][s]*X1[s];
            }
            X2[j] = (1.0/A[j][j])*(B[j] - sum1 - sum2); 
        }
        for(int j = 0; j < n; j++){
            X1[j] = X2[j];
        }
    }

    return X2;
}

arr1d SOR_solve(arr2d &A, arr1d &B, int k){

    int n = A.size();
    double w = 1.5;

    arr1d X1(n, 0);
    arr1d X2(n, 0);

    for(int i = 0; i < k; i++){
        for(int j = 0; j < n; j++){
            double sum1 = 0;
            for(int s = 0; s < j; s++){
                sum1 += A[j][s]*X2[s];
            }
            double sum2 = 0;
            for(int s = j+1; s < n; s++){
                sum2 += A[j][s]*X1[s];
            }
            X2[j] = (1.0 - w)*X1[j] + (w/A[j][j])*(B[j] - sum1 - sum2); 
        }
        for(int j = 0; j < n; j++){
            X1[j] = X2[j];
        }
    }


    return X2;
}

arr1d solve_normally(arr2d &A, arr1d &B){

    std::vector<int> P(A.size()+1,0);
    std::vector<double> X(A.size(),0);
    
    LUPDecompose(A,A.size(), 1.0e-7, P);

    LUPSolve(A, P, B, A.size(), X);
    return X;
}

double calc_error_sum(arr1d &X1, arr1d &X2){

    double error_sum = 0;
    int n = X1.size();
    for(int i = 0; i < n; i++){
        error_sum += my_abs(X1[i] - X2[i]);
    }

    return error_sum;
}

void run_tests(arr1d (*test_func)(arr2d &A, arr1d &B, int k), int iter){
    arr2d A1({{5 ,2 ,1 ,1}, {2 ,6 ,2 ,1}, {1, 2, 7, 1}, {1, 1, 2, 8}});
    arr2d A2({{2,1}, {5,7}});
    arr2d A3({{3.23,2.22,1.11}, {3.34, 4.78, -1.2}, {-1, -1, 1.96}});
    arr2d A4({{333, 90, -20}, {70,120,40}, {20,40,85}});
    arr2d A5({{45, 12, 1, 6, 10},
              {7, 34, 0, 2, 6},
               {6, 3, 62, 2 ,1},
               {5, 4, 1, 34, 3},
               {15, 13, 8, 9, 53}});

    arr2d A6({{29, 2, 12, 7},
              {17, 54, 18, 21},
              {7, 9, 39, 9},
              {2, 4, 1, 19}});


    arr1d B1({29, 31, 26, 19});
    arr1d B2({11,13});
    arr1d B3({7,14,2});
    arr1d B4({90,70,80});
    arr1d B5({33, 17, 5 ,1, 22});
    arr1d B6({1, 22, 1, 0});


    arr1d X;
    arr1d Xn;

    std::cout << "1. blad: ";
    X = test_func(A1, B1, iter);
    Xn = solve_normally(A1, B1);
    std::cout << calc_error_sum(X, Xn) << std::endl;

    std::cout << "2. blad: ";
    X = test_func(A2, B2, iter);
    Xn = solve_normally(A2, B2);
    std::cout << calc_error_sum(X, Xn) << std::endl;
    
    std::cout << "3. blad: ";
    X = test_func(A3, B3, iter);
    Xn = solve_normally(A3, B3);
    std::cout << calc_error_sum(X, Xn) << std::endl;
     
    std::cout << "4. blad: ";
    X = test_func(A4, B4, iter);
    Xn = solve_normally(A4, B4);
    std::cout << calc_error_sum(X, Xn) << std::endl;
    
    std::cout << "5. blad: ";
    X = test_func(A5, B5, iter);
    Xn = solve_normally(A5, B5);
    std::cout << calc_error_sum(X, Xn) << std::endl;
    
    std::cout << "6. blad: ";
    X = test_func(A6, B6, iter);
    Xn = solve_normally(A6, B6);
    std::cout << calc_error_sum(X, Xn) << std::endl;

}

double conv_test(arr1d (*test_func)(arr2d &A, arr1d &B, int k), int iter){
    arr2d A3({{3.23,2.22,1.11}, {3.34, 4.78, -1.2}, {-1, -1, 1.96}});


    arr1d B3({7,14,2});

    arr1d X;
    arr1d Xn;

    // std::cout << "1. blad: ";
    X = test_func(A3, B3, iter);
    Xn = solve_normally(A3, B3);
    // std::cout << calc_error_sum(X, Xn) << std::endl;
    return calc_error_sum(X, Xn);

}
int main(){

    run_tests(SOR_solve, 10);

    std::ofstream sor_file;
    std::ofstream Gauss_Seidel_file;
    std::ofstream Jacobi_file;
    sor_file.open("sor_file.txt", std::ios::out);
    Gauss_Seidel_file.open("Gauss_Seidel.txt", std::ios::out);
    Jacobi_file.open("jacobi_file.txt", std::ios::out);


    for(int i = 0; i < 21; i++){
        Jacobi_file << i << " " << conv_test(Jacobi_solve, i) << std::endl;
        sor_file << i << " " << conv_test(SOR_solve, i) << std::endl;
        Gauss_Seidel_file << i << " " << conv_test(Gauss_Seidel_solve, i) << std::endl;
    }

    sor_file.close();
    Gauss_Seidel_file.close();
    Jacobi_file.close();
                                       
}