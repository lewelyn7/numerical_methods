#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

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
void make_tests(){
    std::function<double(double)> fs[3];
    std::function<double(double)> fs[0] = [](double x){ return exp(x);};

    for(int i = 0; i < 3; i++){
        for(int j = 5; j < 25; j+=5){
            
        }
    }

}
int main(){

    std::function<double(double)> func1 = [](double x){ return exp(x);};
    RectangleIntegration integration;
    TrapeziumIntegration integration2;
    SimpsonIntegration integration3;
    integration.set_precision(10);
    integration2.set_precision(10);
    integration3.set_precision(10);
    std::cout << func1(1) << std::endl << integration.integrate(0,3,func1) << std::endl  << integration2.integrate(0,3,func1) << std::endl  << integration3.integrate(0,3,func1);


}