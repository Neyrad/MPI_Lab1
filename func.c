double t_max = 100;
double x_max = 120;
double t_step = 1.0;
double x_step = 1.0;

double func(double t, double x){
    return x*t;
}

double fi(double x){
    return x*x*x / 12;
}

double ksi(double t){
    return t*t*t / 12;
}