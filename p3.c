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

#define M 1011 // size of input
#define P 1.000000001 // x+delta(x) ~ x * P
#define Q 0.000000001 // delta(x) ~ x * Q
#define VT 0.026

double last_f0 = 0;
double last_f1 = 0;
double last_f2 = 0;

double vectorNorm (double v[], int size) {
    double sum = 0;
    for (int i=0; i<size; i++)
        sum += v[i]*v[i];
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
    out[3] = -(s12*s33-s13*s32)/delta;
    out[6] = (s12*s23-s13*s22)/delta;
    out[1] = -(s21*s33-s23*s31)/delta;
    out[4] = (s11*s33-s13*s31)/delta;
    out[7] = -(s11*s23-s13*s21)/delta;
    out[2] = (s21*s32-s22*s31)/delta;
    out[5] = -(s11*s32-s12*s31)/delta;
    out[8] = (s11*s22-s12*s21)/delta;
}

double computes (double Vds, double Vgs, int size, double Is, double K, double Vth) {
    double IF = Is* pow( log( 1+exp( ( K*(Vgs-Vth)/(2*VT) ) ) ) ,2) ;
    //std::cout <<  K*(Vgs[i]-Vth)/(2*VT)<<" "<< 1+exp( ( K*(Vgs[i]-Vth)/(2*VT) ) ) <<" "<<log( 1+exp( ( K*(Vgs[i]-Vth)/(2*VT) ) ) )  << std::endl;
    double IR = Is* pow( log( 1+exp( ( K*(Vgs-Vth)-Vds)/(2*VT) ) ) ,2);
    double ID = IF - IR;
    
    return ID;
}

double computev (double Vds[], double Vgs[], double Ids[], int size, double Is, double K, double Vth) {
    //std::cout << "Is!!!!!!!!!!!!" <<Is<<std::endl;
    //std::cout << "k!!!!!!!!!!!!" <<K<< std::endl;
    //std::cout << "vth!!!!!!!!!!!!" <<Vth<< std::endl;
    double V = 0;
    
    for (int i=0;i<size;i++) {
        
        double ID = computes(Vds[i], Vgs[i], size, Is, K, Vth);
        
        // std::cout << "round "<< i<<" IF "<< IF << std::endl;
        // std::cout << "IR "<< IR << std::endl;
        
        V += pow((ID-Ids[i]),2);
    }
    
    return V;
}


double computef0(double Vds[], double Vgs[], double Ids[], int size, double Is, double K, double Vth) {
    double f0;
    //std::cout<< Is << " "<< Q*Is << std::endl;
    f0 = (computev(Vds, Vgs, Ids, size, P*Is, K, Vth)-computev(Vds, Vgs, Ids, size, Is, K, Vth))/(Q*Is);
    //std::cout<< "1: "<<computev(Vds, Vgs, Ids, size, P*Is, K, Vth)<<" 2: "<<computev(Vds, Vgs, Ids, size, Is, K, Vth)<<std::endl;
    return f0;
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

void computePace (double f[], double Ji[], double out[]) {
    vectorMatrixMinusMultiply(Ji, f, out);
}

/*
 void secantJ (double Vds[], double Vgs[], double Ids[], int size, double J[], double a[], double delta[]) {
 double f0 = computef0(Vds, Vgs, Ids, size, a[0], a[1], a[2]);
 
 double f1 = computef1(Vds, Vgs, Ids, size, a[0], a[1], a[2]);
 
 double f2 = computef2(Vds, Vgs, Ids, size, a[0], a[1], a[2]);
 
 J[0]=(f0-last_f0)/delta[0];
 //std::cout<<"J0"<<J[0]<<std::endl;
 J[1]=0;
 J[2]=0;
 J[3]=0;
 J[4]=(f1-last_f1)/delta[1];
 //std::cout<<"J4"<<J[4]<<std::endl;
 J[5]=0;
 J[6]=0;
 J[7]=0;
 J[8]=(f2-last_f2)/delta[2];
 //std::cout<<"J8"<<J[8]<<std::endl;
 }
 */

void secantJ (double Vds[], double Vgs[], double Ids[], int size, double J[], double a[], double delta[]) {
    
    J[0]=(computef0(Vds, Vgs, Ids, size, P*a[0],a[1],a[2])-computef0(Vds, Vgs, Ids, size, a[0],a[1],a[2]))/(Q*a[0]);
    //std::cout << "j0: " << J[0]<<std::endl;
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

void sensitivity (double Vds[], double Vgs[], double Ids[], int size, double Is, double K, double Vth, double sens[]) {
    //std::cout << Vds[0]<<Vgs[0]<<Ids[0];
    int i = 16;
    double ID = computes(Vds[i], Vgs[i], size, Is, K, Vth);
    std::cout << ID;
    double ID0 = computes(Vds[i], Vgs[i], size, P*Is, K, Vth);
    double ID1 = computes(Vds[i], Vgs[i], size, Is, P*K, Vth);
    double ID2 = computes(Vds[i], Vgs[i], size, Is, K, P*Vth);
    sens[0] = (ID0/ID)/(Q);
    sens[1] = (ID1/ID)/(Q);
    sens[2] = (ID2/ID)/(Q);
    
}

void search (double Vds[], double Vgs[], double Ids[], int size, double re[]) {
    double Is[8] = {1e-8, 3*1e-8, 1e-7, 3*1e-7, 1e-6, 3*1e-6, 1e-5, 3*1e-5};
    double K[7] = {1, 1.5, 2, 2.5, 3, 3.5, 4};
    double Vth[10] = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    re[0] = Is[0];
    re[1] = K[0];
    re[2] = Vth[0];
    double minv = computev(Vds, Vgs, Ids, size, Is[0], K[0], Vth[0]);
    
    for (int i=0; i<8; i++) {
        for (int j=0; j<7; j++) {
            for (int k=0; k<10; k++) {
                double tmp = computev(Vds, Vgs, Ids, size, Is[i], K[j], Vth[k]);
                //std::cout << i << " "<< j << " " << k << std::endl;
                if (tmp < minv) {
                    minv = tmp;
                    re[0] = Is[i];
                    re[1] = K[j];
                    re[2] = Vth[k];
                    
                }
            }
        }
    }
    
    
}

void quasi_newton() {
    double Vds[M]; // load from file
    double Vgs[M]; // load from file
    double Ids[M]; // load from file
    
    std::ifstream file_reader("/Users/xuhui/Documents/xcodecpp/496p2/p3/496p3/outputCharacteristics.txt" );
    
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
        Vgs[i] = atof(x.c_str());
        Vds[i] = atof(y.c_str());
        Ids[i] = atof(z.c_str());
        i++;
    }
    /*std::ifstream file_reader("/Users/xuhui/Documents/xcodecpp/496p2/p3/496p3/transconductance.txt" );
     
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
     }*/
    
    std::cout<<i<<std::endl;
    int size = i;
    file_reader.close();
    
    // initial guess
    double Is = 0.0000001;
    double K = 1;
    double Vth = 1;
    double a[]={Is, K, Vth};
    
    std::cout << "V: "<< computev(Vds, Vgs, Ids, size, a[0], a[1], a[2]) << std::endl;
    
    double delta[3];
    
    double J[9];
    
    double J_inv[9];
    
    double f[3] = {0.1, 0.1, 0.1};
    
    while ( vectorNorm(f,3) > 0.0001 ) {
        
        std::cout << "Absolute Deviation: "<< computev(Vds, Vgs, Ids, size, a[0], a[1], a[2]) << std::endl;
        std::cout << "Relative Step Size:                    "<< relativeStepSize(a,delta,3) << std::endl;
        
        // compute Jarcobian matrix
        quasiJ(Vds, Vgs, Ids, size, J, a);
        
        
        std::cout<< "J[]::"<<std::endl;
        std::cout<< J[0]<<" "<<J[1]<<" "<<J[2]<<std::endl;
        std::cout<< J[3]<<" "<<J[4]<<" "<<J[5]<<std::endl;
        std::cout<< J[6]<<" "<<J[7]<<" "<<J[8]<<std::endl;
        
        computeInvertedMatrix3(J, J_inv);
        
        // evaluate f(x_i)
        computeF(Vds, Vgs, Ids, size, a[0], a[1], a[2], f);
        std::cout<< "before line search: "<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<< vectorNorm(f,3)<<std::endl;
        
        // compute step size
        computePace(f, J_inv, delta);
        
        // line search
        double t = lineSearch(Vds, Vgs, Ids, size, a, delta, f);
        
        //vectorMultiply(delta, 3, t);
        
        // line search decision
        //std::cout << "percent of origin.: "<<t/10<< std::endl;
        
        // set x_i+1 = x_i + delta x_i
        a[0] += delta[0];
        a[1] += delta[1];
        a[2] += delta[2];
        
        std::cout << "Is << " << a[0] << std::endl;
        std::cout << "K << " << a[1] << std::endl;
        std::cout << "Vth << " << a[2] << std::endl;
        
        
        //std::cout<< "after line search: "<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<vectorNorm(f,3)<<std::endl;
        std::cout<< delta[0]/a[0]<<" delta/ "<<delta[1]/a[1]<<" "<<delta[2]/a[2]<<std::endl;
        
    }
    
    double sen[3];
    sensitivity(Vds, Vgs, Ids, size, a[0], a[1], a[2], sen);
    std::cout << "Parameter Sensitivity: "<< sen[0] << " " << sen[1] << " " << sen[2] << std::endl;
    
    std::cout << "Is << " << a[0] << std::endl;
    std::cout << "K << " << a[1] << std::endl;
    std::cout << "Vth << " << a[2] << std::endl;
    
    // 7 enum
    double re[3];
    search(Vds, Vgs, Ids, size, re);
    std::cout << "Search Result: "<< re[0] << " " << re[1] << " " << re[2] << std::endl;
    
    
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
    
    //std::cout<<i<<std::endl;
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
    
    while ( vectorNorm(f,3) > 0.000000001 ) {
        
        //std::cout << "Absolute Deviation: "<< computev(Vds, Vgs, Ids, size, a[0], a[1], a[2]) << std::endl;
        //std::cout << "Relative Step Size:                     "<< relativeStepSize(a,delta,3) << std::endl;
        
        secantJ(Vds, Vgs, Ids, size, J, a, delta);
        
        computeInvertedMatrix3(J, J_inv);
        
        /*
         std::cout<< "scantJ done"<<std::endl;
         std::cout<< "J[]::"<<std::endl;
         std::cout<< J[0]<<" "<<J[1]<<" "<<J[2]<<std::endl;
         std::cout<< J[3]<<" "<<J[4]<<" "<<J[5]<<std::endl;
         std::cout<< J[6]<<" "<<J[7]<<" "<<J[8]<<std::endl;
         std::cout<< "Jinv[]::"<<std::endl;
         std::cout<< J_inv[0]<<" "<<J_inv[1]<<" "<<J_inv[2]<<std::endl;
         std::cout<< J_inv[3]<<" "<<J_inv[4]<<" "<<J_inv[5]<<std::endl;
         std::cout<< J_inv[6]<<" "<<J_inv[7]<<" "<<J_inv[8]<<std::endl;*/
        
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
        
        //if (t == 0 && last_t == 0) break;
        
        //last_t = t;
        //std::cout << "t: "<<t<< std::endl;
        if (t!=0)
            vectorMultiply(delta, 3, t);
        
        // revise <xi+delta x>
        a[0] += delta[0];
        a[1] += delta[1];
        a[2] += delta[2];
        
        //std::cout<< delta[0]<<" delta "<<delta[1]<<" "<<delta[2]<<std::endl;
        
    }
    
    std::cout << "Is <~ " << a[0] << std::endl;
    std::cout << "K <~ " << a[1] << std::endl;
    std::cout << "Vth <~ " << a[2] << std::endl;
    
    
}

int main(int argc, const char * argv[]) {
    
    // insert code here...
    std::cout << "Quasi-newton\n";
    quasi_newton();
    std::cout << "Secant\n";
    //secant();
    
    return 0;
    
}
