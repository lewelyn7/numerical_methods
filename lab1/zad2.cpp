#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cfloat>

std::fstream freport;


void zad2_1(int n){
    float sum = 0;
    for(int i = 1; i <= n; i++){
        sum += 1.0f/pow(2.0f,i + 1.0f);
    }
    freport << "n: " << n << " sum: " << sum << std::endl;
}
void zad2_2(int n){
    float sum = 0;
    for(int i = n; i >= 1; i--){
        sum += 1.0f/pow(2.0f,i + 1.0f);
        // std::cout << sum << std::endl;
    }
    freport << "n: " << n << " sum: " << sum << std::endl;
}

void zad2_31(int n){
    double sum = 0;
    for(int i = 1; i <= n; i++){
        sum += 1.0/pow(2.0,i + 1.0);
    }
    freport << "n: " << n << " sum: " << sum << std::endl;
}
void zad2_32(int n){
    double sum = 0;
    for(int i = n; i >= 1; i--){
        sum += 1.0/pow(2.0,i + 1.0);
        // std::cout << sum << std::endl;
    }
    freport << "n: " << n << " sum: " << sum << std::endl;
}

float KahanSum(std::vector<float> input){
    float sum_good_final = 2000000;

    float sum = 0.0;
    float err = 0.0;
    for(int i = 1; i < input.size(); i++){
        float y = input[i] - err;
        float temp = sum + y;
        err = (temp - sum) - y;
        sum = temp;
    }

    return sum;

}
float KahanSumD(std::vector<double> input){
    double sum_good_final = 2000000;

    double sum = 0.0;
    double err = 0.0;
    for(int i = 0; i < input.size(); i++){
        double y = input[i] - err;
        double temp = sum + y;
        err = (temp - sum) - y;
        sum = temp;
    }
    std::cout << "err: " << err << std::endl;
    return sum;

}
void KahanWrapperF(int n){
    std::vector<float> vec;
    for(int i = 0; i <= n; i++){
        vec.push_back(1.0f/pow(2.0f,i + 1.0f));
    }

    freport << "n: " << n << " sum: " << KahanSum(vec) << std::endl;

}
void KahanWrapperD(int n){
    std::vector<double> vec;
    for(int i = 1; i <= n; i++){
        vec.push_back(1.0/pow(2.0,i + 1.0));
    }

    freport << "n: " << n << " sum: " << KahanSumD(vec) << std::endl;

}
void zad4(void){

    float EPS = 0.5f;
	float prev_epsilon; 
	while ((1.0f+EPS) != 1.0f) 
	{ 
		prev_epsilon = EPS; 
		EPS /=2.0f; 
	} 
    std::cout << "epsilon for float: " << prev_epsilon << std::endl;

    double EPSd = 0.5;
	double prev_epsilond; 
	while ((1.0+EPSd) != 1.0) 
	{ 
		prev_epsilond = EPSd; 
		EPSd /=2.0; 
	} 
    std::cout << "epsilon for double: " << prev_epsilond;
}
int main(){
    freport.open("zad2_report.txt", std::ios::out);
    
    std::vector<int> vec{50, 100, 200, 500, 800};

    for(int i = 0; i < vec.size() ; i++){
        freport << "zad2.1 " << std::endl;
        zad2_1(vec[i]);
        freport << "zad2.2 " << std::endl;
        zad2_2(vec[i]);
        freport << "zad2.4 " << std::endl;
        KahanWrapperF(vec[i]);
        freport << std::endl;
    }

    freport << std::endl << "-----------DOUBLE-------------" << std::endl << std::endl;

    for(int i = 0; i < vec.size() ; i++){
        freport << "zad2.3.1 " << std::endl;
        zad2_31(vec[i]);
        freport << "zad2.3.2 " << std::endl;
        zad2_32(vec[i]);
        freport << "zad2.4 " << std::endl;
        KahanWrapperD(vec[i]);
        freport << std::endl;
    }
    zad4();
    freport.close();

}