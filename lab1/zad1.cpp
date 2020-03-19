#include <iostream>
#include <vector>
#include <fstream>

double abs(double a){
    if (a < 0){
        return a*(-1);
    }
    return a;
}

void zad1_23(std::vector<float> vec){
    std::fstream freport;
    freport.open("zad1_3_report.txt", std::ios::out);
    long int sum_good = 0;

    float sum_good_final = 2000000;
    float sum = 0;

    for(int i = 0; i  < vec.size(); i++){
        sum += vec[i];
        sum_good += 2;
        if(i%25000 == 0){
             //std::cout << "bezwzgledny: " << sum - ((double)sum_good)/10 << std::endl;
            //  std::cout << "wzgledny: " << abs((sum -  ((double)sum_good)/10)/(((double)sum_good)/10) * 100) << std::endl;    
             freport << abs((sum -  ((double)sum_good)/10)/(((double)sum_good)/10) * 100) << std::endl;    

        }
    }
    std::cout << "sumowanie normalne: " << std::endl;
    std::cout << "na koncu " << sum << std::endl;
    std::cout << "bezwzgledny: " << abs(sum - sum_good_final) << std::endl;
    std::cout << "wzgledny: " << abs((sum - sum_good_final)/sum_good_final * 100) << "%" <<  std::endl << std::endl;
}

void zad1_45(std::vector<float> vec){
    std::fstream freport;
    freport.open("zad1_3_report.txt", std::ios::out);
    float sum_good_final = 2000000;
    float sum = 0;

    for(int j = 2; j < vec.size() * 2; j*=2){
        for(int i = 0; i < vec.size() - j/2; i+=j){
            
            vec[i] = vec[i] + vec[i+j/2];

        }    
    }

    sum = vec[0];
    std::cout << "sumowanie rekurencyjne: " << std::endl;
    std::cout << "na koncu " << sum << std::endl;
    std::cout << "bezwzgledny: " << abs(sum - sum_good_final) << std::endl;
    std::cout << "wzgledny: " << abs((sum - sum_good_final)/sum_good_final * 100) << "%" <<  std::endl <<  std::endl;
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

    std::cout << "sumowanie Kahana: " << std::endl;
    std::cout << "na koncu " << sum << std::endl;
    std::cout << "bezwzgledny: " << abs(sum - sum_good_final) << std::endl;
    std::cout << "wzgledny: " << abs((sum - sum_good_final)/sum_good_final * 100) << "%" <<  std::endl <<  std::endl;
    return sum;

}

int main(){
    std::vector<float> vec;
    std::cout << std::fixed;
    for(int i = 0; i < (int)1e7;i++){
        vec.push_back(0.2f);
    }
    // for(int i = 0; i < (int) 123; i++){
    //     vec.push_back(1.0f);
    // }
    zad1_23(vec);
    zad1_45(vec);
    KahanSum(vec);

    return 0;

}