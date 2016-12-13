//
//  main.cpp
//  p4
//
//  Created by 徐徽 on 16/4/30.
//  Copyright © 2016年 Hui Xu. All rights reserved.
//

#include <iostream>
#include <math.h>

#define rank 1
#define e 2.71828
#define Tol1 1e-1
#define Tol2 1e-3
#define N 2
#define R1 1e4
#define R2 1e4
#define R3 1e4
#define C1 1e-12
#define C2 1e-12
#define print_time 1
#define VDD 5
#define RG 1e4
#define RL 1e4

////////////////////////////////////////////////////////////////////////////////
//  validation
////////////////////////////////////////////////////////////////////////////////

double dfda(double x_, double t){//, double x[]) {
    //x[0]=4*pow(e,(.8*t))-x_[0]*.5;
    return 4*pow(e,(.8*t))-x_*.5;
}

double f0(double x, double x_, double t, double h) {
    return x-x_-h*dfda(x, t+1);
    /*
    double d[1];
    dfda(x, t+1, d);
    f[0] = x[0]-x_[0]-h*d[0];
    return sqrt(f[0]*f[0]);*/
}

// f1 = x_i+1 - x_i - h * (f(x_i+1,t_i+1) + f(x_i, t_i))/2
double f1(double x, double x_, double t, double h){
    return x-x_-0.5*h*(dfda(x, t+1)+dfda(x_,t));
}

// forward euler
double foward(double x_, double t, double h) {
    return x_ + dfda(x_, t) * h;
}


// backward euler
double backward(double x_, double t, double x, double h) {
    //double step[rank];
    //dfda(x_, a, step);
    double step = 0.1;
    
    while (sqrt(f0(x, x_, t, h)*f0(x, x_, t, h))>0.0001 ){
        //std::cout<<"f "<<f0(x, x_, t, h)<<" "<<std::endl;
        double j = (f0(x+0.000001*x, x_, t, h)-f0(x, x_, t, h))/0.000001*x;
        //std::cout<<"j "<<j<<" "<<std::endl;
        step = -1 * f0(x, x_, t, h)/j;
        //std::cout<<"delta "<<step<<" "<<std::endl;
        
        double min = f0(x, x_, t, h);
        bool cont = true;
        double step_ = step;
        while (cont) {
            cont = false;
            if (f0(x+step_, x_, t, h)<min) {
                min = f0(x, x_, t, h);
                cont =true;
            }
            step_ /= 2;
        }
        step = step_ * 2;
        x += step_;
         
        //x += step;
    }
    return x;
    //evolve(x_, step, h, x);
}


// trapezoidal euler
double trapezoidal(double x_, double t, double x, double h) {
    double step = 0.1;
    
    while (sqrt(f1(x, x_, t, h)*f1(x, x_, t, h))>0.0001){
        //std::cout<<"f "<<f1(x, x_, t, h)<<" "<<std::endl;
        double j = (f1(x+0.00001*x, x_, t, h)-f1(x, x_, t, h))/0.00001*x;
        //std::cout<<"j "<<j<<" "<<std::endl;
        step = -1 * f1(x, x_, t, h)/j;
        //std::cout<<"step "<<sqrt(step*step)<<" "<<std::endl;
        
        double min = f1(x, x_, t, h);
        bool cont = true;
        double step_ = step;
        while (cont) {
            cont = false;
            if (f1(x+step_, x_, t, h)<min) {
                min = f1(x, x_, t, h);
                cont =true;
            }
            step_ /= 2;
        }
        step = step_ * 2;
        x += step_;
    }
    return x;
}

// RK4

double rk4(double x, double t, double h) {
    double k1 = dfda(x, t);
    double k2 = dfda(x+k1*h/2, t+h/2);
    double k3 = dfda(x+k2*h/2, t+h/2);
    double k4 = dfda(x+k3*h, t+h);
        //std::cout << k1 << " "<<k2 << " "<< k3 << " "<<k4 << std::endl;
    double x_ = x + (k1+2*k2+2*k3+k4)/6;
    //std::cout << x_ << std::endl;
    return x_;
}

// RK34
double rk34(double x, double t, double h) {
    // rk3
    double x_4;
    while (true) {
        double k1 = dfda(x, t);
        double k2 = dfda(x+k1*h/2, t+h/2);
        double k3 = dfda(x+k2*3*h/4, t+3*h/4);
        double x_3 = x + (2*k1+3*k2+4*k3)/9;
        // rk4
        double k4 = dfda(x+k3*h, t+h);
        // double k4 = dfda(x_3, t+h);
        x_4 = x + (7*k1+6*k2+8*k3+3*k4)*h/24;
        //std::cout << k1 << " "<<k2 << " "<< k3 << " "<<k4 << std::endl;
        double E = x_3 - x_4;
        if (sqrt(E*E)/sqrt(x_4*x_4)>Tol1) {
            std::cout<< "make it smaller" << std::endl;
            h = h/2;
        } else if (sqrt(E*E)/sqrt(x_4*x_4)<Tol2) {
            h = h*2;
            std::cout<< "make it larger" << std::endl;
        } else {
            break;
        }

    }
    //std::cout << x_ << std::endl;
    return x_4;
}






void validate() {
    //std::cout<<2+3.2/1.3*pow(e,0.8)+0.7*pow(e,-0.5)<<std::endl;
    //std::cout<<4*pow(e,0.8)-3.95<<std::endl;
    //std::cout<<4/1.3*0.8*pow(e,0.8)+1.3/0.7*pow(e,-0.5)<<std::endl;
    std::cout<<"True"<<std::endl;
    std::cout<<"6.1946\t14.843\t33.677\t75.339"<<std::endl;
    
    // validate
    double x = 2;
    double t = 0;
    double h = 1;
    // test forward
    std::cout<<"Forward"<<std::endl;
    for (int i=0; i<5; i++) {
        t = i;
        x= foward(x, t, h);
        std::cout<<x<<"\t";
    }
    
    x=2;
    t=0;
    h=1;
    
    // test backward
    std::cout<<"\nBackward"<<std::endl;
    for (int i=0; i<5; i++) {
        t = i;
       x = backward(x, t, x, h);
        std::cout<<x<<"\t";
    }
    x=2;
    t=0;
    h=1;
    
    
    // test trapezoidal
    std::cout<<"\nTrapezoidal"<<std::endl;
    for (int i=0; i<5; i++) {
        t = i;
        x = trapezoidal(x, t, x, h);
        std::cout<<x<<"\t";
    }
    
    x=2;
    t=0;
    h=1;
    
    // test rk4
    std::cout<<"\nRK4"<<std::endl;
    for (int i=0; i<5; i++) {
        t = i;
        x = rk4(x, t, h);
        std::cout<<x<<"\t";
    }
    
    x=2;
    t=0;
    h=1;
    // test rk34
    std::cout<<"\nRK34-adaptive"<<std::endl;
    for (int i=0; i<5; i++) {
        t = i;
        x = rk34(x, t, h);
        std::cout<<x<<"\t";
    }
    
    std::cout<<"\n";
    
    
    
}

////////////////////////////////////////////////////////////////////////////////
//  task 3
////////////////////////////////////////////////////////////////////////////////

double i(double t) {
    t = t*1e9;
    
    //if ( (t<0) or (t>100))
        //std::cout<<t<<" T out of range"<<std::endl;
    
    while (t>20)
        t -= 20;
    
    if (t<1)
        return .0001*t;
    else if (t<10)
        return .0001;
    else if (t<11)
        return .0001-(t-10)*.0001;
    else
        return 0;
}

/*
void dfda(double x[], double t, double d[]) {
    //std::cout<<"x0 "<<x[0]<<std::endl;
    //std::cout<<"x1 "<<x[1]<<std::endl;
    d[0] = -1/C1*( x[0]/R1 + (x[0]-x[1])/R2 - i(t) );
    d[1] = -1/C2*( (x[1]-x[0])/R2 + x[1]/R3 );
}*/

double Ids(double x[]) {
    double Is = 5e-6, K = 0.7, Vth = 1.5 , Vgs = x[0], Vds = x[1], VT = 26e-3;
    return Is*(pow(log(1+exp(K*(Vgs-Vth)/2/VT)),2)-pow(log(1+exp((K*(Vgs-Vth)-Vds)/2/VT)),2));
}

void dfda(double x[], double t, double d[]) {
    d[0] = -1 * x[0]/RG/C1 + i(t);
    d[1] = -1 * Ids(x)/C2 - x[1]/RL/C2 + VDD/RL/C2;
}


// forward euler
void Foward(double x_[], double t, double h, double x[]) {
    double d[N];
    dfda(x_, t, d);
    
    for (int i=0; i<N; i++) {
        x[i] = x_[i] + d[i] * h;
    }
}

double f_tr(double x[], double x_[], double t, double h, double f[]) {
    double d[N];
    dfda(x, t+h, d);
    double d_[N];
    dfda(x_, t, d_);
    //std::cout<<"d0 "<<d[0]<<std::endl;
    //std::cout<<"d1 "<<d[1]<<std::endl;
    for (int i=0; i<N; i++) {
        //f[i] = x[i] - x_[i] - h * d[i];
        f[i] = x[i] - x_[i] - 0.5* h * (d[i]+d_[i]);
    }
    //std::cout<<"f0 "<<f[0]<<std::endl;
    //std::cout<<"f1 "<<f[1]<<std::endl;
    double s = 0;
    for (int i=0; i<N; i++) {
        s += f[i] * f[i];
    }
    return sqrt(s);
}

double f_back(double x[], double x_[], double t, double h, double f[]) {
    double d[N];
    dfda(x, t+h, d);
    //std::cout<<"d0 "<<d[0]<<std::endl;
    //std::cout<<"d1 "<<d[1]<<std::endl;
    for (int i=0; i<N; i++) {
        f[i] = x[i] - x_[i] - h * d[i];
    }
    //std::cout<<"f0 "<<f[0]<<std::endl;
    //std::cout<<"f1 "<<f[1]<<std::endl;
    double s = 0;
    for (int i=0; i<N; i++) {
        s += f[i] * f[i];
    }
    return sqrt(s);
}


double jacobian_bw(double x[], double x_[], double t, double h, double j[], double f[]) {
    double norm = f_back(x, x_, t, h, f);
    //std::cout<<norm<<std::endl;
    
    // dx0
    double x0[N] = {1.00001*x[0], x[1]};
    double f0[N];
    f_back(x0, x_, t, h, f0);
    // dx1
    double x1[N] = {x[0], 1.00001*x[1]};
    double f1[N];
    f_back(x1, x_, t, h, f1);
    // compute J
    //std::cout<<f0[0]<<" <> "<<f[0]<<std::endl;
    j[0] = (f0[0]-f[0])/0.00001/x[0];
    j[1] = (f1[0]-f[0])/0.00001/x[1];
    j[2] = (f0[1]-f[1])/0.00001/x[0];
    j[3] = (f1[1]-f[1])/0.00001/x[1];
    return norm;
}

double jacobian_tr(double x[], double x_[], double t, double h, double j[], double f[]) {
    double norm = f_tr(x, x_, t, h, f);
    //std::cout<<norm<<std::endl;
    
    // dx0
    double x0[N] = {1.00001*x[0], x[1]};
    double f0[N];
    f_tr(x0, x_, t, h, f0);
    // dx1
    double x1[N] = {x[0], 1.00001*x[1]};
    double f1[N];
    f_tr(x1, x_, t, h, f1);
    // compute J
    //std::cout<<f0[0]<<" <> "<<f[0]<<std::endl;
    j[0] = (f0[0]-f[0])/0.00001/x[0];
    j[1] = (f1[0]-f[0])/0.00001/x[1];
    j[2] = (f0[1]-f[1])/0.00001/x[0];
    j[3] = (f1[1]-f[1])/0.00001/x[1];
    return norm;
}

double f_rkad(double x[], double x_[], double t, double h, double f[], double k1[], double k2[], double k3[], double k4[]) {
    dfda(x_, t, k1);
    double xx[] = {x_[0]+k1[0]*h/2,x_[1]+k1[1]*h/2};
    //std::cout << xx[0] << " xx0 "<<xx[1] << " xx1"<<std::endl;
    dfda(xx, t+h/2, k2);
    xx[0] = x_[0]+k2[0]*3*h/4;
    xx[1] = x_[1]+k2[1]*3*h/4;
    //std::cout << xx[0] << " xx0 "<<xx[1] << " xx1"<<std::endl;
    dfda(xx, t+h/2, k3);
    //std::cout << xx[0] << " xx0 "<<xx[1] << " xx1"<<std::endl;
    
    
    dfda(x, t+h, k4);
    for (int i=0; i<N; i++) {
        //std::cout << k1 << " "<<k2 << " "<< k3 << " "<<k4 << std::endl;
        f[i] = x_[i] + h*(7*k1[i]+6*k2[i]+8*k3[i]+3*k4[i])/24 - x[i];
        //std::cout << x_ << " x_ "<<(k1+2*k2+2*k3+k4)/6<< "  phih "<<std::endl;
    }
    double s = 0;
    for (int i=0; i<N; i++) {
        s += f[i] * f[i];
    }
    return sqrt(s);
}

double jacobian_rkad(double x[], double x_[], double t, double h, double j[], double f[], double k1[], double k2[], double k3[], double k4[]) {
    double norm = f_rkad(x, x_, t, h, f, k1, k2, k3, k4);
    // dx0
    double x0[N] = {1.00001*x[0], x[1]};
    double f0[N];
    f_rkad(x0, x_, t, h, f0, k1, k2, k3, k4);
    // dx1
    double x1[N] = {x[0], 1.00001*x[1]};
    double f1[N];
    f_rkad(x1, x_, t, h, f1, k1, k2, k3, k4);
    // compute J
    //std::cout<<f0[0]<<" <> "<<f[0]<<std::endl;
    j[0] = (f0[0]-f[0])/0.00001/x[0];
    j[1] = (f1[0]-f[0])/0.00001/x[1];
    j[2] = (f0[1]-f[1])/0.00001/x[0];
    j[3] = (f1[1]-f[1])/0.00001/x[1];
    return norm;
}

void inverted(double j[]) {
    double ij[N*N];
    //std::cout<<"j0 "<<j[0]<<" "<<j[1]<<" "<<j[2]<<" "<<j[3]<<std::endl;
    double delta = j[0]*j[3]-j[1]*j[2];
    ij[0] = j[3]/delta;
    ij[1] = -1*j[1]/delta;
    ij[2] = -1*j[2]/delta;
    ij[3] = j[0]/delta;
    j[0] = ij[0];
    j[1] = ij[1];
    j[2] = ij[2];
    j[3] = ij[3];
    
    //std::cout<<"j0 "<<j[0]<<" "<<j[1]<<" "<<j[2]<<" "<<j[3]<<std::endl;
}

void step(double f[], double j[], double del[]) {
    inverted(j);
    for (int i=0; i<N; i++) {
        double tmp = 0;
        for (int k=0; k<N; k++) {
            tmp += f[k]*j[i*N+k];
            }
        del[i] = -1*tmp;
    }
}


// backward euler
void Backward(double x_[], double t, double h, double x[]) {
    double delta[N];
    double f[N];
    double norm = 0.1;
    while (norm>0.001 ){
        //std::cout<<"f "<<f0(x, x_, t, h)<<" "<<std::endl;
        
        // compute J
        double j[N*N];
        norm = jacobian_bw(x, x_, t, h, j, f);
        
        
        //std::cout<<"norm "<<norm<<" "<<std::endl;
        
        // compute step
        step(f, j, delta);
        //std::cout<<"del0 "<<delta[0]<<std::endl;
        //std::cout<<"del1 "<<delta[1]<<std::endl;
    
        //step = -1 * f0(x, x_, t, h)/j;
        //std::cout<<"delta "<<step<<" "<<std::endl;
        
        
        double min = norm;
        bool cont = true;
        //std::cout<<"linesearch"<<std::endl;
        while (cont) {
            cont = false;
            double tmp = f_back(x, x_, t, h, f);
            //std::cout<<"min "<<min<<" temp "<<tmp<<std::endl;
            if (tmp<min) {
                //std::cout<<"find!"<<std::endl;
                min = tmp;
                cont =true;
                for (int i=0; i<N; i++) {
                    delta[i] = delta[i]/2;
                }
            }
        }
        //std::cout<<"linesearch end"<<std::endl;
        //std::cout<<"x0 "<<x[0]<<std::endl;
        //std::cout<<"x1 "<<x[1]<<std::endl;

        for (int i=0; i<N; i++) {
            x[i] += delta[i];
        }
        //std::cout<<"del0 "<<delta[0]<<std::endl;
        //std::cout<<"del1 "<<delta[1]<<std::endl;
        //x += step;
    }
    //return x;
    //evolve(x_, step, h, x);
}

// Heun's
void TR(double x_[], double t, double h, double x[]) {
    double delta[N];
    double f[N];
    double norm = 0.1;
    while (norm>0.001) {
        //std::cout<<"f "<<norm<<" "<<std::endl;
        
        // compute J
        double j[N*N];
        norm = jacobian_tr(x, x_, t, h, j, f);
        
        
        //std::cout<<"norm "<<norm<<" "<<std::endl;
        
        // compute step
        step(f, j, delta);
        //std::cout<<"del0 "<<delta[0]<<std::endl;
        //std::cout<<"del1 "<<delta[1]<<std::endl;
        
        //step = -1 * f0(x, x_, t, h)/j;
        //std::cout<<"delta "<<step<<" "<<std::endl;
        
        
        double min = norm;
        bool cont = true;
        //std::cout<<"linesearch"<<std::endl;
        while (cont) {
            cont = false;
            double tmp = f_tr(x, x_, t, h, f);
            //std::cout<<"min "<<min<<" temp "<<tmp<<std::endl;
            if (tmp<min) {
                //std::cout<<"find!"<<std::endl;
                min = tmp;
                cont =true;
                for (int i=0; i<N; i++) {
                    delta[i] = delta[i]/2;
                }
            }
        }
        //std::cout<<"linesearch end"<<std::endl;
        //std::cout<<"x0 "<<x[0]<<std::endl;
        //std::cout<<"x1 "<<x[1]<<std::endl;
        
        for (int i=0; i<N; i++) {
            x[i] += delta[i];
        }
        //std::cout<<"del0 "<<delta[0]<<std::endl;
        //std::cout<<"del1 "<<delta[1]<<std::endl;
        //x += step;
    }
    //return x;
    //evolve(x_, step, h, x);
}

// RK4
void RK4(double x[], double t, double h) {
    double d1[N],d2[N],d3[N],d4[N];
    dfda(x, t, d1);
    double xx[] = {x[0]+d1[0]*h/2,x[1]+d1[1]*h/2};
    //std::cout << xx[0] << " xx0 "<<xx[1] << " xx1"<<std::endl;
    dfda(xx, t+h/2, d2);
    xx[0] = x[0]+d2[0]*h/2;
    xx[1] = x[1]+d2[1]*h/2;
    //std::cout << xx[0] << " xx0 "<<xx[1] << " xx1"<<std::endl;
    dfda(xx, t+h/2, d3);
    xx[0] = x[0]+d3[0]*h;
    xx[1] = x[1]+d3[1]*h;
    //std::cout << xx[0] << " xx0 "<<xx[1] << " xx1"<<std::endl;
    dfda(xx, t+h, d4);
    
    for (int i=0; i<N; i++) {
        double k1 = d1[i];
        //std::cout <<k1<<" "<<i<<std::endl;
        double k2 = d2[i];
        double k3 = d3[i];
        double k4 = d4[i];
        //std::cout << k1 << " "<<k2 << " "<< k3 << " "<<k4 << std::endl;
        x[i] += h*(k1+2*k2+2*k3+k4)/6;
        //std::cout << x_ << " x_ "<<(k1+2*k2+2*k3+k4)/6<< "  phih "<<std::endl;
    }

}

void ForwardEuler() {
    double h[] = {2e-10, 1e-9};
    double x[] = {2.5, 2.5};
    double t = 0;
    int j = 1;
    
    
    // forward euler
    std::cout<<"\nFW"<<std::endl;
    for (int i=0; i<1e-7/h[j]; i++) {
        //std::cout<<t<<" t ";
        Foward(x, t, h[j], x);
        for(int m=0; m<print_time; m++) {
            for (int k=0; k<N; k++) {
                std::cout<<x[k]<<"\t";
            }
            std::cout<<std::endl;
        }
        t += h[j];
    }
}

void BackwardEuler() {
    double h[] = {2e-10, 1e-9};
    double x[] = {2.5, 2.5};
    double t = 0;
    int j = 1;

    std::cout<<"BW"<<std::endl;

    // backward euler
    for (int i=0; i<1e-7/h[j]; i++) {
        //std::cout<<t<<" t ";
        double guess[N] = {x[0]+0.01, x[1]+0.01};
        Backward(x, t, h[j], guess);
        for(int m=0; m<print_time; m++) {
            for (int k=0; k<N; k++) {
                std::cout<<guess[k]<<"\t";
                x[k] = guess[k];
            }
            std::cout<<"\n";
        }
        //std::cout<<"\n"<<t<<"                *******"<<i<<std::endl;
        t += h[j];
    }
}

void TrapezoidalEuler() {
    double h[] = {2e-10, 1e-9};
    double x[] = {2.5, 2.5};
    double t = 0;
    int j = 1;
    
    std::cout<<"TR"<<std::endl;

    // trapezoidal
    for (int i=0; i<1e-7/h[j]; i++) {
        //std::cout<<t<<" t ";
        double guess[N] = {x[0]+0.01, x[1]+0.01};
        TR(x, t, h[j], guess);
        for (int m=0; m<print_time; m++) {
            for (int k=0; k<N; k++) {
                std::cout<<x[k]<<"\t";
                x[k] = guess[k];
            }
            std::cout<<"\n";
        }
        //std::cout<<"\n"<<t<<"                *******"<<i<<std::endl;
        t += h[j];
    }
}

void RungeKutta4() {
    double h[] = {2e-10, 1e-9};
    double x[] = {2.5, 2.5};
    double t = 0;
    int j = 1;
    
    std::cout<<"RK4"<<std::endl;

    // forward euler
    for (int i=0; i<1e-7/h[j]; i++) {
        //std::cout<<t<<" t ";
        //double guess[N] = {x[0]+0.0001, x[1]+0.0001};
        RK4(x, t, h[j]);
        for(int m=0; m<print_time; m++) {
            for (int k=0; k<N; k++) {
                std::cout<<x[k]<<"\t";
                //x[k] = guess[k];
            }
            std::cout<<"\n";
        }
        //std::cout<<"\n"<<t<<"                *******"<<i<<std::endl;
        t += h[j];
    }
}

// rk34 adaptive
void RK34ad(double x_[], double t, double h, double x[]) {
    double f[N], E[N], k1[N], k2[N], k3[N], k4[N];
    double lastcheck = 0;
    while (true) {
        
        double norm = 0.1;
        while (norm>0.001) {
            // get f
            double j[N*N];
            // compute j
            norm = jacobian_rkad(x, x_, t, h, j, f, k1, k2, k3, k4);
        
            double delta[N];
            step(f, j, delta);
        
            double min = norm;
            bool cont = true;
            //std::cout<<"linesearch"<<std::endl;
            while (cont) {
                cont = false;
                double tmp = f_rkad(x, x_, t, h, f, k1, k2, k3, k4);
                //std::cout<<"min "<<min<<" temp "<<tmp<<std::endl;
                if (tmp<min) {
                    //std::cout<<"find!"<<std::endl;
                    min = tmp;
                    cont =true;
                    for (int i=0; i<N; i++) {
                        delta[i] = delta[i]/2;
                    }
                }
            }
            //std::cout<<"linesearch end"<<std::endl;
            //std::cout<<"x0 "<<x[0]<<std::endl;
            //std::cout<<"x1 "<<x[1]<<std::endl;
        
            for (int i=0; i<N; i++) {
                x[i] += delta[i];
            }
        }
        //std::cout<<k1[0]<< " " << k2[0]<< " "<<k3[0]<<" "<<k4[0]<<std::endl;
        // error estimation
        for (int i=0; i<N; i++) {
            E[i] = (-5*k1[i]+6*k2[i]+8*k3[i]-9*k4[i])*h/72;
            //std::cout<<i<< " -- " << E[i] << " "<<h<<std::endl;
        }
        int check=0;
        
        for (int i=0; i<N; i++) {
            if (sqrt(E[i]*E[i])/sqrt(x[i]*x[i])>Tol1) {
                //std::cout<< "--" << std::endl;
                check--;
            } else if (sqrt(E[i]*E[i])/sqrt(x[i]*x[i])<Tol2) {
                check++;
                //std::cout<< "++" << std::endl;
            } else {
                //std::cout<< "!" << std::endl;
            }
        }
        
        // adaptive
        if (check<0) {
            h = h/2;
            //std::cout<< "make it smaller" << std::endl;
            if (check != lastcheck) break;
        } else if (check>0) {
            h = h*2;
            //std::cout<< "make it larger" << std::endl;
            if (check != lastcheck) break;
        } else {
            //std::cout<< "goood" << std::endl;
            lastcheck = check;
            break;
        }
    }
}

// rk34
void RK34(double x_[], double t, double h, double x[]) {
    double f[N], k1[N], k2[N], k3[N], k4[N];
        double norm = 0.1;
        while (norm>0.001) {
            // get f
            double j[N*N];
            // compute j
            norm = jacobian_rkad(x, x_, t, h, j, f, k1, k2, k3, k4);
            
            double delta[N];
            step(f, j, delta);
            
            double min = norm;
            bool cont = true;
            //std::cout<<"linesearch"<<std::endl;
            while (cont) {
                cont = false;
                double tmp = f_rkad(x, x_, t, h, f, k1, k2, k3, k4);
                //std::cout<<"min "<<min<<" temp "<<tmp<<std::endl;
                if (tmp<min) {
                    //std::cout<<"find!"<<std::endl;
                    min = tmp;
                    cont =true;
                    for (int i=0; i<N; i++) {
                        delta[i] = delta[i]/2;
                    }
                }
            }
            //std::cout<<"linesearch end"<<std::endl;
            //std::cout<<"x0 "<<x[0]<<std::endl;
            //std::cout<<"x1 "<<x[1]<<std::endl;
            
            for (int i=0; i<N; i++) {
                x[i] += delta[i];
            }
        }
        //std::cout<<k1[0]<< " " << k2[0]<< " "<<k3[0]<<" "<<k4[0]<<std::endl;
        // error estimation

}


void RungeKutta34() {
    double h[] = {2e-10, 1e-9};
    for (int j=1; j<2; j++) {
    double x[] = {2.5, 2.5};
    double t = 0;
        std::cout<<"RK34!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    // forward euler
    for (int i=0; i<1e-7/h[j]; i++) {
        //std::cout<<t<<" t ";
        double guess[N] = {x[0]+0.0001, x[1]+0.0001};
        RK34(x, t, h[j], guess);
        for (int m=0; m<print_time; m++) {
            for (int k=0; k<N; k++) {
                std::cout<<x[k]<<"\t";
                x[k] = guess[k];
            }
            std::cout<<"\n";
        }
        //std::cout<<"\n"<<t<<"                *******"<<i<<std::endl;
        t += h[j];
    }
    }
}

void RungeKutta34_adaptive() {
    double h[] = {2e-10, 1e-9};
    for (int j=1; j<2; j++) {
        double x[] = {2.5, 2.5};
        double t = 0;
        std::cout<<"RK34 adaptive!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
        // forward euler
        for (int i=0; i<1e-7/h[j]; i++) {
            //std::cout<<t<<" t ";
            double guess[N] = {x[0]+0.0001, x[1]+0.0001};
            RK34ad(x, t, h[j], guess);
            for(int m=0; m<print_time; m++){
                for (int k=0; k<N; k++) {
                    std::cout<<x[k]<<"\t";
                    x[k] = guess[k];
                }
                std::cout<<"\n";
            }
            //std::cout<<"\n"<<t<<"                *******"<<i<<std::endl;
            t += h[j];
        }
    }
}

void task3() {
    //ForwardEuler();
    //BackwardEuler();
    //TrapezoidalEuler();
    //RungeKutta4();
    //RungeKutta34();
    RungeKutta34_adaptive();
}

int main(int argc, const char * argv[]) {
    //validate();
    task3();
    
    

    
    return 0;
}

