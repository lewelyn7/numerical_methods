#include "aghMatrix.h"
#include <iostream>

int main() 
{
    // initialize matrices using init value
    AGHMatrix<double> mat1(5, 5, 1.2);
    AGHMatrix<double> mat2(5, 5, 2.8);

    // Uncomment when implemented
    AGHMatrix<double> mat3 = mat1 + mat2;
    std::cout << mat3;

    // initialize matrix using specified values
    std::vector<std::vector<double>> init { { 1.0, 2.0, 3.0 }, 
                                            { 4.0, 5.0, 6.0 }, 
                                            { 7.0, 8.0, 9.0 } }; 

    // initialize matrix using specified values
    std::vector<std::vector<double>> init_LU {{ 5.0, 3.0, 2.0 }, 
                                            { 1.0, 2.0, 0.0 }, 
                                            { 3.0, 0.0, 4.0 }};

    // AGHMatrix<double> mat4(init_LU);

    // std::cout << LUdet(mat4);


    // initialize matrix using specified values
    std::vector<std::vector<double>> init_cholesky {{ 4.0, 12.0, -16.0 }, 
                                                    { 12.0, 37.0, -43.0 }, 
                                                    { -16.0, -43.0, 98.0 }};

    std::vector<std::vector<double>> gauss_check{{0.0001, -5.0300, 5.8090, 7.8320, 9.5740},
               {2.2660, 1.9950,  1.2120, 8.0080, 7.2190},
               {8.8500, 5.6810,  4.5520, 1.3020, 5.7300},
               {6.7750, -2.253,  2.9080, 3.9700, 6.2910}};


    std::vector<std::vector<double>> jacoby_check{{2,1,11},{5,7,13}};

    // Jeśli się korzysta z implementacji laboratoryjnej
    AGHMatrix<double> mat5(init_cholesky);
    AGHMatrix<double> mat6(gauss_check);
    AGHMatrix<double> mat7(jacoby_check);
    std::cout << std::endl;
    // std::cout << Gauss_alg(mat6);
    std::cout << std::endl << gauss_solve(mat6);

    std::cout << std::endl << Jacoby(mat7, 100);
    return 0;

}