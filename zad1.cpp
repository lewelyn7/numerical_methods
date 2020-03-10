#include <iostream>
#include <vector>
#include <fstream>

double abs(double a){
    if (a < 0){
        return a*(-1);
    }
    return a;
}


int main(){
    std::fstream freport;

    freport.open("zad1_report.txt", std::ios::out);

    std::vector<float> vec;
    float x = 2000000;
    float sum = 0;

    long int sum_good = 0;

    for(int i = 0; i < (int)1e7;i++){
        vec.push_back(0.2f);
    }
    for(int i = 0; i  < vec.size(); i++){
        sum += vec[i];
        sum_good += 2;
        if(i%25000 == 0){
             //std::cout << "bezwzgledny: " << sum - ((double)sum_good)/10 << std::endl;
             std::cout << "wzgledny: " << abs((sum -  ((double)sum_good)/10)/(((double)sum_good)/10) * 100) << std::endl;    
             freport << abs((sum -  ((double)sum_good)/10)/(((double)sum_good)/10) * 100) << std::endl;    

        }
    }

    std::cout << "na koncu" << std::endl;
    std::cout << "bezwzgledny: " << sum - x << std::endl;
    std::cout << "wzgledny: " << (sum - x)/x * 100 << std::endl;

    return 0;

}