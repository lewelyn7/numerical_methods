#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <chrono>

class DFT{
    public:
        static std::vector<std::vector<double>> dft(std::vector<double> data){
            std::vector<std::vector<double>>  result;
            result.resize(data.size());
            int n = data.size();
            for(int i = 0; i < n; i++){
                result[i].resize(2);
                double sum = 0;
                double sumi = 0;
                for(int j = 0; j < n; j++){
                    sum += data[j]*(cos(2.0*M_PI*(double)i*(double)j/(double)n));
                    sumi += -1.0*data[j]*(sin(2.0*M_PI*(double)i*(double)j/(double)n));
                }
                result[i][0] = sum;
                result[i][1] = sumi;

            }
            return result;
        }
};

class FFT{
    public:
        static std::vector<std::vector<double>> fft(std::vector<double> data){
            std::vector<std::vector<double>>  result;
            int n = data.size();
            std::vector<double> data_even;
            std::vector<double> data_odd;
            if(n <= 32){
                return DFT::dft(data);
            }

            for(int i = 0;i < n; i+=2) data_even.push_back(data[i]);
            std::vector<std::vector<double>> x_even = FFT::fft(data_even);

            for(int i = 1;i < n; i+=2) data_odd.push_back(data[i]);
            std::vector<std::vector<double>> x_odd = FFT::fft(data_odd);


            std::vector<std::vector<double>> factor;
            for(int i = 0; i < n; i++){
                std::vector<double> item;
                item.push_back(cos(-2.0*M_PI*(double)i/(double)n));
                item.push_back(sin(-2.0*M_PI*(double)i/(double)n));
                factor.push_back(item);
            }

            result.resize(n);
            for(int i = 0; i < n/2; i++){
                result[i].resize(2);
                //ac - bd + (bc +ad)i
                result[i][0] = x_odd[i][0]*factor[i][0] - x_odd[i][1]*factor[i][1] + x_even[i][0];
                result[i][1] = x_odd[i][1]*factor[i][0] + x_odd[i][0]*factor[i][1] + x_even[i][1];

            }
            for(int i = 0; i < n/2; i++){
                result[i+n/2].resize(2);
                //ac - bd + (bc +ad)i
                result[i+n/2][0] = x_odd[i][0]*factor[i+n/2][0] - x_odd[i][1]*factor[i+n/2][1] + x_even[i][0];
                result[i+n/2][1] = x_odd[i][1]*factor[i+n/2][0] + x_odd[i][0]*factor[i+n/2][1] + x_even[i][1];

            }
            return result;
        }
};
void print_vector(std::vector<double> vec){
    int n  = vec.size();
    for(int i = 0; i < n; i++){
        std::cout << vec[i] << std::endl;
    }
}

void save_complex(std::string filename, std::vector<std::vector<double>> arr){
    int n = arr.size();
    std::ofstream file;
    file.open(filename, std::ios::out);
    for(int i = 0; i < n; i++){
        file << arr[i][0] << " " << arr[i][1] << std::endl;
    }
    file.close();
}

std::vector<double> read_from_file(std::string filename){
    std::ifstream file;
    file.open(filename, std::ios::in);
    std::vector<double> input;

    while(file){
        double val;
        file >> val;
        input.push_back(val);
 
    }
    input.pop_back();
    return input;
}

void test(void){
    std::ofstream fstats;
    std::ofstream dstats;
    fstats.open("stats_f.txt", std::ios::out);
    dstats.open("stats_d.txt", std::ios::out);

    for(int i = 0; i < 10; i++){
        std::vector<double> dataset = read_from_file("data" + std::to_string(i) + ".txt");
        std::chrono::steady_clock::time_point start_fft = std::chrono::steady_clock::now();
        std::vector<std::vector<double>> fresult = FFT::fft(dataset);
        std::chrono::steady_clock::time_point stop_fft = std::chrono::steady_clock::now();

        std::chrono::steady_clock::time_point start_dft = std::chrono::steady_clock::now();
        std::vector<std::vector<double>> dresult = DFT::dft(dataset);
        std::chrono::steady_clock::time_point stop_dft = std::chrono::steady_clock::now();

        fstats << dataset.size() << " ";
        fstats << std::chrono::duration_cast<std::chrono::microseconds>(stop_fft - start_fft).count() << std::endl;

        dstats << dataset.size() << " ";
        dstats << std::chrono::duration_cast<std::chrono::microseconds>(stop_dft - start_dft).count() << std::endl;
    }
    fstats.close();
    dstats.close();
}

std::vector<std::vector<double>> amplitude_and_phase(std::vector<std::vector<double>> data){
    std::vector<std::vector<double>> result;
    int n = data.size();
    double T = 24.0;
    result.resize(n);
    for(int i = 0; i < n; i++){
        result[i].resize(2);
        result[i][0] = T*i/(double)n;
        result[i][1] = sqrt(data[i][0]*data[i][0] + data[i][1]*data[i][1]) * 1.0/(double)n;
    }
    return result;
}

int main(){

    // std::vector<double> dataset = {1, 2, 3, 4, 5, 1, 2 ,3 ,4 ,5};
    std::vector<double> dataset = read_from_file("weather_dataset.txt");
    // print_vector(dataset);
    std::vector<std::vector<double>> result = DFT::dft(dataset);
    save_complex("weather_result.txt", amplitude_and_phase(result));
    // for(int i = 0; i < result.size(); i++){
    //     std::cout << result[i][0] <<  "  " << result[i][1] << std::endl;
    // }
    // test();
}