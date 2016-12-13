//
//  main.cpp
//  QUASI_NEWTON
//
//  Created by 徐徽 on 16/4/14.
//  Copyright © 2016年 Hui Xu. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <iostream>
#include <fstream>

#define M 199 // size of input
#define P 1.0000001 // x+delta(x) ~ x * P
#define Q 0.0000001 // delta(x) ~ x * Q
#define VT 0.026

double last_f0 = 0;
double last_f1 = 0;
double last_f2 = 0;

double vectorNorm (double v1[], int size) {
    double sum = 0;
    for (int i=0; i<size; i++)
        sum += v1[i]*v1[i];
    return sqrt(sum);
}

void vectorMatrixMinusMultiply(double a[], double b[], double c[]) {
    for (int i=0; i<3; i++) {
        double temp = 0;
        for (int j=0; j<3; j++) {
            temp += a[i*3+j]*b[j];
        }
        c[i] = -1*temp;
    }
}

void computeInvertedMatrix3(double in[], double out[]) {
    double s11 = in[0], s12 = in[1], s13 = in[2];
    double s21 = in[3], s22 = in[4], s23 = in[5];
    double s31 = in[6], s32 = in[7], s33 = in[8];
    double delta = (s11*s22*s33-s11*s23*s32-s21*s12*s33+s21*s13*s32+s31*s12*s23-s31*s13*s22);
    /*std::cout<< "delta[]::"<<std::endl;
    std::cout<< s11*s22*s33<<" "<<s11*s23*s32<<" "<<s21*s12*s33<<std::endl;
    std::cout<< s21*s13*s32<<" "<<s31*s12*s23<<" "<<s31*s13*s22<<std::endl;
    std::cout<< "delta::"<< delta <<std::endl;*/
    out[0] = (s22*s33-s23*s32)/delta;
    out[1] = -(s12*s33-s13*s32)/delta;
    out[2] = (s12*s23-s13*s22)/delta;
    out[3] = -(s21*s33-s23*s31)/delta;
    out[4] = (s11*s33-s13*s31)/delta;
    out[5] = -(s11*s23-s13*s21)/delta;
    out[6] = (s21*s32-s22*s31)/delta;
    out[7] = -(s11*s32-s12*s31)/delta;
    out[8] = (s11*s22-s12*s21)/delta;
}

double computev (double Vds[], double Vgs[], double Ids[], int size, double Is, double K, double Vth) {
    //std::cout << "Is!!!!!!!!!!!!" <<Is<<std::endl;
    //std::cout << "k!!!!!!!!!!!!" <<K<< std::endl;
    //std::cout << "vth!!!!!!!!!!!!" <<Vth<< std::endl;
    double V = 0;
    
    for (int i=0;i<size;i++) {
        double IF = Is* pow( log( 1+exp( ( K*(Vgs[i]-Vth)/(2*VT) ) ) ) ,2) ;
        //std::cout <<  K*(Vgs[i]-Vth)/(2*VT)<<" "<< 1+exp( ( K*(Vgs[i]-Vth)/(2*VT) ) ) <<" "<<log( 1+exp( ( K*(Vgs[i]-Vth)/(2*VT) ) ) )  << std::endl;
        double IR = Is* pow( log( 1+exp( ( K*(Vgs[i]-Vth)-Vds[i])/(2*VT) ) ) ,2);
        double ID = IF - IR;
        
          // std::cout << "round "<< i<<" IF "<< IF << std::endl;
          // std::cout << "IR "<< IR << std::endl;

        V += pow((ID-Ids[i]),2);
    }
    
    return V;
}

double computev (double t[], double v[], int size, double a0, double a1, double a2) {
    double V = 0;
    
    for (int i=0;i<size;i++) {
        double V2 = a0*exp(-1*t[i]/a1)+(v[0]-a0)*exp(-1*t[i]/a2);
        // std::cout << "round "<< i<<" IF "<< IF << std::endl;
        // std::cout << "IR "<< IR << std::endl;
        
        V += pow((V2-v[i]),2);
    }
    
    return V;
}


double computef0(double Vds[], double Vgs[], double Ids[], int size, double Is, double K, double Vth) {
    double f0;
    f0 = (computev(Vds, Vgs, Ids, size, P*Is, K, Vth)-computev(Vds, Vgs, Ids, size, Is, K, Vth))/(Q*Is);
    //std::cout<< "1: "<<computev(Vds, Vgs, Ids, size, P*Is, K, Vth)<<" 2: "<<computev(Vds, Vgs, Ids, size, Is, K, Vth)<<std::endl;
    return f0;
}

double computef0(double x0[], double x1[], int size, double a0, double a1, double a2) {
    double a[] = {a0, a1, a2};
    return (computev(x0, x1, size, P*a[0], a[1], a[2])-computev(x0, x1, size, a[0], a[1], a[2]))/(Q*a[0]);
}

double computef1(double x0[], double x1[], int size, double a0, double a1, double a2) {
    double a[] = {a0, a1, a2};
    return (computev(x0, x1, size, a[0], P*a[1], a[2])-computev(x0, x1, size, a[0], a[1], a[2]))/(Q*a[1]);
}

double computef2(double x0[], double x1[], int size, double a0, double a1, double a2) {
    double a[] = {a0, a1, a2};
    return (computev(x0, x1, size, a[0], a[1], P*a[2])-computev(x0, x1, size, a[0], a[1], a[2]))/(Q*a[2]);
}

double computef1(double Vds[], double Vgs[], double Ids[], int size, double Is, double K, double Vth) {
    double f1;
    f1 = (computev(Vds, Vgs, Ids, size, Is, P*K, Vth)-computev(Vds, Vgs, Ids, size, Is, K, Vth))/(Q*K);
    return f1;
}

double computef2(double Vds[], double Vgs[], double Ids[], int size, double Is, double K, double Vth) {
    double f2;
    f2 = (computev(Vds, Vgs, Ids, size, Is, K, P*Vth)-computev(Vds, Vgs, Ids, size, Is, K, Vth))/(Q*Vth);
    return f2;
}

void computeF(double Vds[], double Vgs[], double Ids[], int size, double Is, double K, double Vth, double f[]) {
    f[0] = computef0(Vds, Vgs, Ids, size, Is, K, Vth);
    f[1] = computef1(Vds, Vgs, Ids, size, Is, K, Vth);
    f[2] = computef2(Vds, Vgs, Ids, size, Is, K, Vth);
    if (f[0]!=f[0])
        std::cout << "!!!!!!!!!!!!" << std::endl;
}

void computeF(double t[], double v2[], int size, double a0, double a1, double a2, double f[]) {
    f[0] = computef0(t, v2, size, a0, a1, a2);
    f[1] = computef1(t, v2, size, a0, a1, a2);
    f[2] = computef2(t, v2, size, a0, a1, a2);
    if (f[0]!=f[0])
        std::cout << "!!!!!!!!!!!!" << std::endl;
}

void computePace (double f[], double Ji[], double out[]) {
    vectorMatrixMinusMultiply(Ji, f, out);
}

void secantJ (double Vds[], double Vgs[], double Ids[], int size, double J[], double a[], double delta[]) {
    
    J[0]=(computef0(Vds, Vgs, Ids, size, P*a[0],a[1],a[2])-computef0(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[0]);
    std::cout << "computef0(Vds, Vgs, Ids, size, P*a[0],a[1],a[2]) " << computef0(Vds, Vgs, Ids, size, P*a[0],a[1],a[2])<<std::endl;
    std::cout << "computef0(Vds, Vgs, Ids, size, a[0],a[1],a[2]) " << computef0(Vds, Vgs, Ids, size, a[0],a[1],a[2])<<std::endl;
    J[1]=0;
    //std::cout << "j0: " << J[1]<<std::endl;
    J[2]=0;
    //std::cout << "j0: " << J[2]<<std::endl;
    J[3]=0;
    J[4]=(computef1(Vds, Vgs, Ids, size, a[0],P*a[1],a[2])-computef1(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[1]);
    J[5]=0;
    J[6]=0;
    J[7]=0;
    J[8]=(computef2(Vds, Vgs, Ids, size, a[0],a[1],P*a[2])-computef2(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[2]);
    
}

void quasiJ (double x0[], double x1[], int size, double J[], double a[]) {
    // P=1.00001 Q=0.00001
        J[0]=(computef0(x0, x1, size, P*a[0],a[1],a[2])-computef0(x0, x1, size, a[0],a[1],a[2]))/(Q*a[0]);
        J[1]=(computef0(x0, x1, size, a[0],P*a[1],a[2])-computef0(x0, x1, size, a[0],a[1],a[2]))/(Q*a[1]);
        J[2]=(computef0(x0, x1, size, a[0],a[1],P*a[2])-computef0(x0, x1, size, a[0],a[1],a[2]))/(Q*a[2]);
        J[3]=(computef1(x0, x1, size, P*a[0],a[1],a[2])-computef1(x0, x1, size, a[0],a[1],a[2]))/(Q*a[0]);
        J[4]=(computef1(x0, x1, size, a[0],P*a[1],a[2])-computef1(x0, x1, size, a[0],a[1],a[2]))/(Q*a[1]);
        J[5]=(computef1(x0, x1, size, a[0],a[1],P*a[2])-computef1(x0, x1, size, a[0],a[1],a[2]))/(Q*a[2]);
        J[6]=(computef2(x0, x1, size, P*a[0],a[1],a[2])-computef2(x0 ,x1, size, a[0],a[1],a[2]))/(Q*a[0]);
        J[7]=(computef2(x0, x1, size, a[0],P*a[1],a[2])-computef2(x0, x1, size, a[0],a[1],a[2]))/(Q*a[1]);
        J[8]=(computef2(x0, x1, size, a[0],a[1],P*a[2])-computef2(x0, x1, size, a[0],a[1],a[2]))/(Q*a[2]);

}

void quasiJ (double Vds[], double Vgs[], double Ids[], int size, double J[], double a[]) {
    J[0]=(computef0(Vds, Vgs, Ids, size, P*a[0],a[1],a[2])-computef0(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[0]);
    //std::cout << "j0: " << J[0]<<std::endl;
    J[1]=(computef0(Vds, Vgs, Ids, size, a[0],P*a[1],a[2])-computef0(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[1]);
    //std::cout << "j0: " << J[1]<<std::endl;
    J[2]=(computef0(Vds, Vgs, Ids, size, a[0],a[1],P*a[2])-computef0(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[2]);
    //std::cout << "j0: " << J[2]<<std::endl;
    J[3]=(computef1(Vds, Vgs, Ids, size, P*a[0],a[1],a[2])-computef1(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[0]);
    J[4]=(computef1(Vds, Vgs, Ids, size, a[0],P*a[1],a[2])-computef1(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[1]);
    J[5]=(computef1(Vds, Vgs, Ids, size, a[0],a[1],P*a[2])-computef1(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[2]);
    J[6]=(computef2(Vds, Vgs, Ids, size, P*a[0],a[1],a[2])-computef2(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[0]);
    J[7]=(computef2(Vds, Vgs, Ids, size, a[0],P*a[1],a[2])-computef2(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[1]);
    J[8]=(computef2(Vds, Vgs, Ids, size, a[0],a[1],P*a[2])-computef2(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[2]);

}

void vectorAdd (double v1[], double v2[], int size, int time) {
    for (int i=0; i<size; i++) {
        v1[i] = v1[i] + v2[i]*time;
    }
}

void vectorAdd (double v1[], double v2[], int size) {
    vectorAdd(v1, v2, size, 1);
}

void vectorMinus (double v1[], double v2[], int size, int time) {
    for (int i=0; i<size; i++) {
        v1[i] = v1[i] - v2[i]*time;
    }
}

void vectorCopy (double v1[], double v2[], int size) {
    for (int i=0; i<size; i++) {
        v1[i] = v2[i];
    }
}

int lineSearch (double Vds[], double Vgs[], double Ids[], int size, double abc_[], double delta_abc[], double f[]) {
    double tempf[3];
    computeF(Vds, Vgs, Ids, size, abc_[0], abc_[1], abc_[2], tempf);
    double minf = vectorNorm(tempf, 3);
    bool cont = true;
    while (cont) {
        cont = false;
        std::cout<<"delta "<<delta_abc[0]<<" abc "<<abc_[0]<<std::endl;
        vectorAdd(abc_, delta_abc, size);
        std::cout<<"abc "<<abc_[0]<<std::endl;
        computeF(Vds, Vgs, Ids, size, abc_[0], abc_[1], abc_[2], tempf);
        std::cout<<"!@#  "<<vectorNorm(tempf, 3)<<std::endl;
        if (vectorNorm(tempf, 3)<minf) {
            std::cout<<"!!!"<<std::endl;
            minf = vectorNorm(tempf, 3);
            cont =true;
        }
        vectorMinus(abc_, delta_abc, size, 1);
        delta_abc[0]/=2;
        delta_abc[1]/=2;
        delta_abc[2]/=2;
        
        //step_ /= 2;
    }
    delta_abc[0]*=2;
    delta_abc[1]*=2;
    delta_abc[2]*=2;
    return 0;
}


int lineSearch (double t[], double v2[], int size, double abc_[], double delta_abc[], double f[]) {
    double minf = 1000;
    bool cont = true;
    double tempf[3];
    while (cont) {
        cont = false;
        computeF(t, v2, size, abc_[0], abc_[1], abc_[2], tempf);
        std::cout<<"!@#  "<<vectorNorm(tempf, 3)<<std::endl;
        if (vectorNorm(tempf, 3)<minf) {
            std::cout<<"!!!"<<std::endl;
            minf = vectorNorm(tempf, 3);
            cont =true;
        }
        delta_abc[0]/=2;
        delta_abc[1]/=2;
        delta_abc[2]/=2;
        //step_ /= 2;
    }
    delta_abc[0]*=2;
    delta_abc[1]*=2;
    delta_abc[2]*=2;
    return 0;
}

/*
int lineSearch (double Vds[], double Vgs[], double Ids[], int size, double abc_[], double delta_abc[], double f[]) {
    int T=10,t=0;
    double part[3];
    for (int i=0; i<3; i++) {
        part[i]=delta_abc[i]/T;
    }
    
    double minf = 1000;
    double tempf[3];
    
    for (int i=0; i<T; i++) {
        vectorAdd(abc_, part, 3, i);
        computeF(Vds, Vgs, Ids, size, abc_[0], abc_[1], abc_[2], tempf);
        
        if (vectorNorm(tempf, 3)<minf) {
            minf = vectorNorm(tempf, 3);
            //vectorCopy(mindelta,delta_abc, 3);
            t = i;
        }
        vectorMinus(abc_, part, 3, i);
    }
    //if (t!=0) {
        //vectorCopy(delta_abc, mindelta, 3);
        //vectorCopy(f, minf, 3);
    //}
    return t;
}
*/

double relativeStepSize (double a[], double del[], int size) {
    double c = 0;
    for (int i=0; i<size; i++) {
        c += del[i]*del[i]/a[i]/a[i];
    }
    return sqrt(c);
}

void vectorMultiply (double v[], int size, int t) {
    for (int i=0; i<size; i++) {
        v[i] = v[i]*t/10;
    }
}

void quasi_newton() {
    //double Vds[M]; // load from file
    //double Vgs[M]; // load from file
    //double Ids[M]; // load from file
    double t[M];
    double v2[M];

    
    std::ifstream file_reader("/Users/xuhui/Documents/xcodecpp/QUASI_NEWTON/h1.txt" );
    
    if(!file_reader.is_open()){
        std::cout << "Could not open file!"<<'\n';
    }
    
    // skip the metadata
    std::string x;
    std::string y;
    //std::string z;
    file_reader >> x;
    file_reader >> y;
    //file_reader >> z;
    
    // read data
    int i = 0;
    while(!file_reader.eof()){
        file_reader >> x;
        file_reader >> y;
        //file_reader >> z;
        //Vds[i] = atof(x.c_str());
        //Vgs[i] = atof(y.c_str());
        //Ids[i] = atof(z.c_str());

        i++;
    }
    
    std::cout<<i<<std::endl;
    int size = i;
    file_reader.close();
    
    // initial guess
    
    /*
    double Is = 0.0000001;
    double K = 1;
    double Vth = 1;
    double a[]={Is, K, Vth};
    */
    
    //double a[] = {1, 1, 1};

    double delta[3];
    
    double J[9];
    
    double J_inv[9];
    
    double f[3] = {0.1, 0.1, 0.1};

    while ( vectorNorm(f,3) > 0.0001 ) {
        //std::cout << "v "<<computev(t, v2, size, a[0], a[0], a[2])<< std::endl;
        //std::cout << "Deviation: "<<relativeStepSize(a,delta,3)<<std::endl;
        
        // compute Jarcobian matrix
        //quasiJ(Vds, Vgs, Ids, size, J, a);
        quasiJ(t, v2, size, J, a);
        
        
        std::cout<< "J[]::"<<std::endl;
        std::cout<< J[0]<<" "<<J[1]<<" "<<J[2]<<std::endl;
        std::cout<< J[3]<<" "<<J[4]<<" "<<J[5]<<std::endl;
        std::cout<< J[6]<<" "<<J[7]<<" "<<J[8]<<std::endl;
        
        computeInvertedMatrix3(J, J_inv);
        
        
        std::cout<< "Jinv[]::"<<std::endl;
        std::cout<< J_inv[0]<<" "<<J_inv[1]<<" "<<J_inv[2]<<std::endl;
        std::cout<< J_inv[3]<<" "<<J_inv[4]<<" "<<J_inv[5]<<std::endl;
        std::cout<< J_inv[6]<<" "<<J_inv[7]<<" "<<J_inv[8]<<std::endl;
        
        
        // evaluate f(x_i)
        //computeF(Vds, Vgs, Ids, size, a[0], a[1], a[2], f);
        computeF(t, v2, size, a[0], a[1], a[2], f );
        std::cout<< "before line search: "<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<< vectorNorm(f,3)<<std::endl;
        
        // compute step size
        computePace(f, J_inv, delta);
        
        // line search
        //double t = lineSearch(Vds, Vgs, Ids, size, a, delta, f);
        double tt = lineSearch(t, v2, size, a, delta, f);
        
        //vectorMultiply(delta, 3, t);
        
        //if (t==0) break;
        
        // line search decision
        //std::cout << "percent of origin.: "<<t/10<< std::endl;
        
        // set x_i+1 = x_i + delta x_i
        a[0] += delta[0];
        a[1] += delta[1];
        a[2] += delta[2];
        
        
        std::cout<< "after line search: "<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<vectorNorm(f,3)<<std::endl;
        std::cout<< delta[0]/a[0]<<" delta "<<delta[1]/a[1]<<" "<<delta[2]/a[2]<<std::endl;
        
    }
    
    
    /*std::cout << "Is << " << a[0] << std::endl;
    std::cout << "K << " << a[1] << std::endl;
    std::cout << "Vth << " << a[2] << std::endl;
    */
    std::cout << "a[0] << " << a[0] << std::endl;
    std::cout << "a[1] << " << a[1] << std::endl;
    std::cout << "a[2] << " << a[2] << std::endl;
    
}

void secant() {
    double Vds[M]; // load from file
    double Vgs[M]; // load from file
    double Ids[M]; // load from file
    
    std::ifstream file_reader("/Users/xuhui/Documents/xcodecpp/496p2/p3/496p3/transconductance.txt" );
    
    if(!file_reader.is_open()){
        std::cout << "Could not open file!"<<'\n';
    }
    
    // skip the metadata
    std::string x;
    std::string y;
    std::string z;
    file_reader >> x;
    file_reader >> y;
    file_reader >> z;
    
    // read data
    int i = 0;
    while(!file_reader.eof()){
        file_reader >> x;
        file_reader >> y;
        file_reader >> z;
        Vds[i] = atof(x.c_str());
        Vgs[i] = atof(y.c_str());
        Ids[i] = atof(z.c_str());
        i++;
    }
    
    std::cout<<i<<std::endl;
    int size = i;
    file_reader.close();
    
    // initial guess
    double Is = 0.0000001;
    double K = 1;
    double Vth = 1;
    double a[]={Is, K, Vth};
    
    double delta[] = {0.1, 0.2, 0.3};
    
    last_f0 = 0.001;
    last_f1 = 0.001;
    last_f2 = 0.001;
    
    double J[9];
    
    double J_inv[9];
    
    double f[3] = {100, 0.001, 0.001};
    
    double last_t = 1;
    
    while ( vectorNorm(f,3) > 0.0001 ) {
        
        secantJ(Vds, Vgs, Ids, size, J, a, delta);
        
        computeInvertedMatrix3(J, J_inv);
        
        
        std::cout<< "scantJ done"<<std::endl;
        std::cout<< "J[]::"<<std::endl;
        std::cout<< J[0]<<" "<<J[1]<<" "<<J[2]<<std::endl;
        std::cout<< J[3]<<" "<<J[4]<<" "<<J[5]<<std::endl;
        std::cout<< J[6]<<" "<<J[7]<<" "<<J[8]<<std::endl;
        std::cout<< "Jinv[]::"<<std::endl;
        std::cout<< J_inv[0]<<" "<<J_inv[1]<<" "<<J_inv[2]<<std::endl;
        std::cout<< J_inv[3]<<" "<<J_inv[4]<<" "<<J_inv[5]<<std::endl;
        std::cout<< J_inv[6]<<" "<<J_inv[7]<<" "<<J_inv[8]<<std::endl;
        
        computeF(Vds, Vgs, Ids, size, a[0], a[1], a[2], f);
        
        //std::cout<< "computeF done"<<std::endl;
        
        last_f0 = f[0];
        last_f1 = f[1];
        last_f2 = f[2];
        
        //std::cout<< f[0]<<" f "<<f[1]<<" "<<f[2]<<std::endl;
        
        // compute delta x using initial guess <xi>
        computePace(f,J_inv, delta);
        
        //std::cout<< delta[0]<<" delta "<<delta[1]<<" "<<delta[2]<<std::endl;
        
        double t = lineSearch(Vds, Vgs, Ids, size, a, delta, f);
        
        if (t == 0 && last_t == 0) break;
        
        last_t = t;
        //std::cout << "t: "<<t<< std::endl;
        
        vectorMultiply(delta, 3, t);
        
        // revise <xi+delta x>
        a[0] += delta[0];
        a[1] += delta[1];
        a[2] += delta[2];
        
        std::cout<< delta[0]<<" delta "<<delta[1]<<" "<<delta[2]<<std::endl;
        std::cout << "Is <~ " << a[0] << std::endl;
        std::cout << "K <~ " << a[1] << std::endl;
        std::cout << "Vth <~ " << a[2] << std::endl;
        
    }
    
    std::cout << "Is <~ " << a[0] << std::endl;
    std::cout << "K <~ " << a[1] << std::endl;
    std::cout << "Vth <~ " << a[2] << std::endl;
    
}

int main(int argc, const char * argv[]) {
    // insert code here...
    
    quasi_newton();
    //secant();
    std::cout << "Hello, World!\n";
    
    return 0;
    
}
