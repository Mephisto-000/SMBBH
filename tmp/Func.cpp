#include "Func.h"


double init_circle_v(double M1, double R, double G){
    return sqrt(G*M1 / R);
}


double init_elliptical_v(double e, double init_circle_v){
    return sqrt((1 - e)*pow(init_circle_v, 2));
}


double Two_Body_Func(int eq, double t, double w[], double G, double mu){
    double M;
    double r = sqrt(pow(w[1], 2) + pow(w[2], 2) + pow(w[3], 2));

    if (eq == 1) {M = w[4];}
    else if (eq == 2) {M = w[5];}
    else if (eq == 3) {M = w[6];}

    else if (eq == 4) {M = -mu*w[1] / pow(r, 3);}
    else if (eq == 5) {M = -mu*w[2] / pow(r, 3);}
    else if (eq == 6) {M = -mu*w[3] / pow(r, 3);}
    else {return 0;}

    return M;
}


double Barycentric(int eq, double w[], double M1, double M2, double pi){
    double M;
    double ratio_M1 = M2 / (M1 + M2);
    double ratio_M2 = M1 / (M1 + M2);
    double c = (M1 + M2)*2 / pi;
    double r1 = sqrt(pow(w[1], 2) + pow(w[2], 2) + pow(w[3], 2));
    double r2 = sqrt(pow(w[4], 2) + pow(w[5], 2) + pow(w[6], 2));

    // M1_pos : 
    if (eq == 1) {M = w[1]*ratio_M1;}
    else if (eq == 2) {M = w[2]*ratio_M1;}
    else if (eq == 3) {M = w[3]*ratio_M1;}
    // M2_pos : 
    else if (eq == 4) {M = -w[1]*ratio_M2;}
    else if (eq == 5) {M = -w[2]*ratio_M2;}
    else if (eq == 6) {M = -w[3]*ratio_M2;}

    // M1_v : 
    else if (eq == 7) {M = w[4]*ratio_M1 + (c*w[1]/r1)*(1/(r1 + pow(r1, 3)) - atan(r1)/pow(r1, 2));}
    else if (eq == 8) {M = w[5]*ratio_M1 + (c*w[2]/r1)*(1/(r1 + pow(r1, 3)) - atan(r1)/pow(r1, 2));}
    else if (eq == 9) {M = w[6]*ratio_M1 + (c*w[3]/r1)*(1/(r1 + pow(r1, 3)) - atan(r1)/pow(r1, 2));}
    // M2_v : 
    else if (eq == 10) {M = -w[4]*ratio_M2 + (c*w[4]/r2)*(1/(r2 + pow(r2, 3)) - atan(r2)/pow(r2, 2));}
    else if (eq == 11) {M = -w[5]*ratio_M2 + (c*w[5]/r2)*(1/(r2 + pow(r2, 3)) - atan(r2)/pow(r2, 2));}
    else if (eq == 12) {M = -w[6]*ratio_M2 + (c*w[6]/r2)*(1/(r2 + pow(r2, 3)) - atan(r2)/pow(r2, 2));}

    else {return 0;}


    return M;
}


