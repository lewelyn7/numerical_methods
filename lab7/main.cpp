#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <fstream>
#include <string>
#include <random>

double my_abs(double x){
    if(x > 0){
        return x;
    }else{
        return -x;
    }
}

class IIntegrationEngine
{
    public:
        virtual double integrate(double a, double b, std::function<double(double)> f) = 0;
        virtual void set_precision(int n) {};
};

class RectangleIntegration: public IIntegrationEngine
{
    private:
        int n = 1;
    public:
        virtual double  integrate(double a, double b, std::function<double(double)> f)
        {
            std::vector<double> points;
            points.resize(n);
            for(int i = 0; i < n; i++){
                points[i] = (a + ((double)i)/((double)n) * (b -a));
            }

            double dx = (b - a)/((double)n);
            std::vector<double> values;
            values.resize(n);
            for(int i = 0; i < n; i++){
                values[i] = (f(points[i]));
            }
            double sum = 0;
            for(int i = 0; i < n; i++){
                sum += values[i];
            }
            sum *= dx;
            return sum;
        }
        virtual void set_precision(int x)
        {
            this->n = x;
        }
};

class TrapeziumIntegration: public IIntegrationEngine
{
    private:
        int n = 1;
    public:
        virtual double  integrate(double a, double b, std::function<double(double)> f)
        {
            std::vector<double> points;
            points.resize(n+1);
            for(int i = 0; i <= n; i++){
                points[i] = (a + ((double)i)/((double)n) * (b -a));
            }

            double dx = (b - a)/((double)n);
            std::vector<double> values;
            values.resize(n+1);
            for(int i = 0; i <= n; i++){
                values[i] = (f(points[i]));
            }
            double sum = 0;
            for(int i = 0; i < n; i++){
                sum += values[i] + values[i+1];
            }
            sum *= dx/2;
            return sum;
        }
        virtual void set_precision(int x)
        {
            this->n = x;
        }
};
class SimpsonIntegration: public IIntegrationEngine
{
    private:
        int n = 1;
    public:
        virtual double  integrate(double a, double b, std::function<double(double)> f)
        {
            std::vector<double> points;
            points.resize(n+1);
            for(int i = 0; i <= n; i++){
                points[i] = (a + ((double)i)/((double)n) * (b -a));
            }

            std::vector<double> t_points;
            t_points.resize(n);
            for(int i = 0; i < n; i++){
                t_points[i] = (points[i+1] + points[i])/2;
            }
            double dx = (b - a)/((double)n);
            std::vector<double> values;
            values.resize(n+1);
            for(int i = 0; i <= n; i++){
                values[i] = (f(points[i]));
            }
            std::vector<double> t_values;
            t_values.resize(n);
            for(int i = 0; i < n; i++){
                t_values[i] = (f(t_points[i]));
            }
            double sum = values[0] + values[n];
            double partial_sum1 = 0;
            for(int i = 1; i < n; i++){
                partial_sum1 += values[i];
            }
            double partial_sum2 = 0;
            for(int i = 0; i < n; i++){
                partial_sum2 += t_values[i];
            }

            sum = (dx/6.0)*(sum + 2.0*partial_sum1 + 4.0*partial_sum2);
            return sum;
        }
        virtual void set_precision(int x)
        {
            this->n = x;
        }
};

class MCIntegration: public IIntegrationEngine
{
    private:
        int n = 1;
    public:
        virtual double  integrate(double a, double b, std::function<double(double)> f)
        {
            std::default_random_engine gen; 
            std::uniform_real_distribution<double> distr(a,b); 
            double dx = b-a;
            double sum = 0;
            for(int i = 0; i < n; i++){
                double x = distr(gen);
                sum += f(x);
            }
            return sum*dx/(double)n;
        }
        virtual void set_precision(int x)
        {
            this->n = x;
        }
};

double pi_MC(int N){ 
    std::default_random_engine gen; 
    std::uniform_real_distribution<double> distr(0,1.0); 
    int hits = 0;
    for (int i = 0; i < N; i++) {  
        double x = distr(gen); 
        double y = distr(gen); 
        double l = sqrt(x * x + y * y); 
        if (l <= 1) hits++; 
    } 
    double pi = double(hits) / (double)N * 4.0;
    return pi; 
} 

void make_tests(){
    std::function<double(double)> fs[3];
    fs[0] = [](double x){ return exp(x);};
    fs[1] = [](double x){ return sin(x*x*x);};
    fs[2] = [](double x){ return log(x)*x*x;};
    double ab[3][2] = {
        {1, 3},
        {1, 3},
        {1, 3}
    };
    double exact_values[3] = {17.367, 0.222579, 6.9986};



    RectangleIntegration rect;
    TrapeziumIntegration trap;
    SimpsonIntegration simp;


    for(int i = 0; i < 3; i++){
        std::ofstream frect;
        std::ofstream ftrap;
        std::ofstream fsimpson;
        frect.open("rect" + std::to_string(i) + ".txt", std::ios::out);
        ftrap.open("trap" + std::to_string(i) + ".txt", std::ios::out);
        fsimpson.open("simpson" + std::to_string(i) + ".txt", std::ios::out);
        for(int j = 5; j < 40; j+=5){
            rect.set_precision(j);
            trap.set_precision(j);
            simp.set_precision(j);

            frect << j << " " << my_abs(rect.integrate(ab[i][0], ab[i][1], fs[i]) - exact_values[i]) << std::endl;
            ftrap << j << " " << my_abs(trap.integrate(ab[i][0], ab[i][1], fs[i]) - exact_values[i]) << std::endl;
            fsimpson << j << " " << my_abs(simp.integrate(ab[i][0], ab[i][1], fs[i]) - exact_values[i]) << std::endl;
        }
        frect.close();
        ftrap.close();
        fsimpson.close();
    }

    //Monte Carlo test
    std::vector<double> monte_values({10, 100, 1000, 10000, 100000, 1000000, 10000000});
    MCIntegration mcinter;
    for(int i = 0; i < 3; i++){
        std::ofstream fmonte;
        fmonte.open("monte" + std::to_string(i) + ".txt", std::ios::out);
        
        for(int j = 0; j < monte_values.size(); j++){
            mcinter.set_precision(monte_values[j]);
            fmonte << j*5+5 << " " << my_abs(mcinter.integrate(ab[i][0], ab[i][1], fs[i]) - exact_values[i]) << std::endl;
        }

        fmonte.close();
    }

}


void pi_test(){
    std::vector<double> values({10000, 100000, 1000000, 10000000, 100000000});
    double exact_value = M_PI;
    std::ofstream pif;
    pif.open("pi.txt", std::ios::out);
    for(int i = 0; i < values.size(); i++){
        double calc_pi = pi_MC(values[i]);
        std::cout << values[i] << "  " << calc_pi << "  " << my_abs(calc_pi - exact_value) << std::endl; 
        pif << std::fixed << i << "  " << my_abs(calc_pi - exact_value) << std::endl; 
    }
    pif.close();
}   
int main(){
    make_tests();
    // pi_test();
    
    // std::function<double(double)> func1 = [](double x){ return exp(x);};
    // RectangleIntegration integration;
    // TrapeziumIntegration integration2;
    // SimpsonIntegration integration3;
    // integration.set_precision(10);
    // integration2.set_precision(10);
    // integration3.set_precision(10);
    // std::cout << func1(1) << std::endl << integration.integrate(0,3,func1) << std::endl  << integration2.integrate(0,3,func1) << std::endl  << integration3.integrate(0,3,func1);


}