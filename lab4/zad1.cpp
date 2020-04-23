#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

double func1(double x){
    double k = 2;
    double m = 1.2;
    return sin((k*x)/M_PI)*exp(-(m*x)/M_PI);
    // return exp(k*cos(m*x));
    // return sin(x);
}
double func1prim(double x){
    double h = 0.001;
    return (func1(x) - func1(x+h)) / h;
}
struct point{
    double x;
    double y;
};

int fact(int n) 
{ 
    int f = 1; 
    for (int i = 2; i <= n; i++) 
        f *= i; 
    return f; 
} 


double uCoefNewton(double u, int n) 
{ 
    double temp = u; 
    for (int i = 1; i < n; i++) 
        temp = temp * (u - i); 
    return temp; 
} 

std::vector<std::vector<double>> Newton_diff_calc(std::vector<struct point> points){
    
    std::vector<std::vector<double>> diff_matrix(points.size());
    for(int i = 0; i < points.size(); i++){
        diff_matrix[i].resize(points.size());
        diff_matrix[i][0] = points[i].y;
    }
    for (int i = 1; i < points.size(); i++) { 
        for (int j = 0; j < points.size() - i; j++) 
            diff_matrix[j][i] = diff_matrix[j + 1][i - 1] - diff_matrix[j][i - 1]; 
    } 
    return diff_matrix;   
}

double proterm(int i, double value, std::vector<struct point> x) 
{ 
    double pro = 1; 
    for (int j = 0; j < i; j++) { 
        pro = pro * (value - x[j].x); 
    } 
    return pro; 
} 

std::vector<std::vector<double>> Hermite_diff_calc(std::vector<struct point> points){
    
    std::vector<std::vector<double>> diff_matrix(points.size());
    for(int i = 0; i < points.size(); i++){
        diff_matrix[i].resize(points.size());
        diff_matrix[i][0] = points[i].y;
    }

    for (int i = 1; i <  points.size(); i++) { 
        for (int j = 0; j <  points.size() - i; j++) { 
            if(points[i+j].x == points[j].x){
                diff_matrix[j][i] = func1prim(points[j].x);
            }else{
            diff_matrix[j][i] = (diff_matrix[j][i - 1] - diff_matrix[j + 1][i - 1]) / (points[j].x - points[i + j].x); 
            }
        } 
    } 

    return diff_matrix;   
}

double Lagrange(std::vector<point> points, int x) 
{ 
    double result = 0;
  
    for (int i=0; i<points.size(); i++) 
    { 
 
        double term = points[i].y; 
        for (int j=0;j<points.size();j++) 
        { 
            if (j!=i) 
                term = term*(x - points[j].x)/double(points[i].x - points[j].x); 
        } 

        result += term; 
    } 
  
    return result; 
} 

double Newton(std::vector<std::vector<double>> diff_matrix, std::vector<struct point> points, double x){

    double sum = diff_matrix[0][0]; 
    double u = (x - points[0].x) / (points[1].x - points[0].x); 
    for (int i = 1; i < diff_matrix.size(); i++) { 
        sum = sum + (uCoefNewton(u, i) * diff_matrix[0][i]) / fact(i); 
    }  

    return sum;
}

std::vector<struct point> calc_func_Newton(std::vector<struct point> points, double step){

    std::vector<std::vector<double>> diffs = Newton_diff_calc(points);
    std::vector<struct point> result;
    for(double i = points[0].x; i < points[points.size()-1].x; i += step){
        struct point pt;
        pt.x = i;
        pt.y = Newton(diffs, points, i);
                    std::cout << "x: " << pt.x << " y: " << pt.y << std::endl;
        result.push_back(pt);
    }

    return result;

}
std::vector<struct point> calc_func_Lagrange(std::vector<struct point> points, double step){

    std::vector<struct point> result;
    for(double i = points[0].x; i < points[points.size()-1].x; i += step){
        struct point pt;
        pt.x = i;
        pt.y = Lagrange(points, i);
        result.push_back(pt);
    }

    return result;

}

double cheby_point(double k, double n, double a, double b){
    k += 1;
    double val = 0.5*(a+b) + 0.5*(b-a)*cos(((2*k-1)*M_PI)/(2*n));
    return val;
}

std::vector<struct point> Hermite_prep(std::vector<struct point> points){
    // std::vector<std::vector<double>> temp;
    std::vector<double> z;
    std::vector<struct point> pts;
    z.resize(points.size()*2);
    pts.resize(points.size()*2);
    for(int i = 0; i < points.size(); i++){
        z[i*2] = z[i*2+1] = points[i].x;
    }

    for(int i = 0; i < z.size(); i++){
        pts[i].x = z[i];
        pts[i].y = func1(z[i]);
    }

    return pts;
}

double Hermite(std::vector<std::vector<double>> y, std::vector<struct point> x, double value) 
{ 
    float sum = y[0][0]; 
  
    for (int i = 1; i < y.size(); i++) { 
      sum = sum + (proterm(i, value, x) * y[0][i]); 
    } 
    return sum; 
} 


std::vector<struct point> calc_func_Hermite(std::vector<struct point> points, double step){
    std::vector<struct point> hermite_prep = Hermite_prep(points);
    std::vector<std::vector<double>> diffs = Hermite_diff_calc(hermite_prep);
    std::vector<struct point> result;
    for(double i = points[0].x; i < points[points.size()-1].x; i += step){
        struct point pt;
        pt.x = i;
        pt.y = Hermite(diffs, hermite_prep, i);
        std::cout << "x: " << pt.x << " y: " << pt.y << std::endl;
        result.push_back(pt);
    }

    return result;

}

void print_vector(std::vector<struct point> points, std::string name){
    std::ofstream outfile;
    outfile.open(name);
    for(int i = 0; i < points.size(); i++){
        outfile << points[i].x << " " << points[i].y << std::endl;
    }
}

int main(){
    std::ofstream func1file;
    
    func1file.open("func1.txt");
    // double ttt[] = {0.7071, 0.7660, 0.8192, 0.8660};
    // int liczba_wezlow[] = { 2, 3, 5, 20};
    std::vector<int> liczba_wezlow = {2,3,5,20};

    for(int i = 0; i < liczba_wezlow.size(); i++){
        std::vector<struct point> points;
        double step = 20.0/(double)liczba_wezlow[i];
        for(double i  = 0 ; i <= 20; i+=step){
            struct point pt;
            pt.x = i;
            pt.y = func1(i);
            points.push_back(pt);
        }

        // std::vector<double> xs;
        // for(double j  = 0 ; j < (double)liczba_wezlow[i]; j++){
        //     xs.push_back(cheby_point(j,liczba_wezlow[i], 0, 20));
        // }
        // std::sort(xs.begin(), xs.end());
        // for(double j  = 0 ; j < (double)liczba_wezlow[i]; j++){
        //     struct point pt;            
        //     pt.x = xs[j];
        //     pt.y = func1(pt.x);

        //     points.push_back(pt);           
        // }


        // print_vector(points, "Cheby_out" + std::to_string(liczba_wezlow[i]) + ".txt");
        // print_vector(calc_func_Newton(points, 0.5), "Cheby_out_Newton"  + std::to_string(liczba_wezlow[i]) +  ".txt");
        // print_vector(calc_func_Lagrange(points, 0.5), "Cheby_out_Lagrange"  + std::to_string(liczba_wezlow[i]) + ".txt");
        print_vector(calc_func_Hermite(points, 0.5), "out_Hermite"  + std::to_string(liczba_wezlow[i]) + ".txt");

    }


    
}
