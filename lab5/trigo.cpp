
double my_abs(double x){
    if(x > 0){
        return x;
    }else{
        return x*-1.0;
    }
}
double a_coef( int n, arr2d &points){
    double sum = 0;
    for(int i = 0; i < points.size(); i++){
        sum += points[i][1]*cos(points[i][0]*n);
    }
    return sum/points.size()*2;
}
double b_coef( int n, arr2d &points){
    double sum = 0;
    for(int i = 0; i < points.size(); i++){
        sum += points[i][1]*sin(points[i][0]*n);
    }
    return sum/points.size()*2;
}

double result_func_trigo(double x, int n, arr2d &points){
    double sum = 0.5*a_coef(0,points);
    for(int i = 1; i < n; i++){
        sum += (a_coef(i, points)*cos(i*x) + b_coef(i,points)*sin(i*x));
    }
    return sum;
}
arr2d generatePoints3(arr2d &points, double a, double b, double n, int approx_lvl){

    arr2d res;
    double step = my_abs(b-a)/n;


    int it = 0;
    for(double i = a; i < b; i +=step){
        std::vector<double> a;
        a.resize(2);
        a[0] = i;
        a[1] = result_func_trigo(i, approx_lvl, points);
        res.push_back(a);
        it++;
    }
    return res;
}