#include "Func.h"

using namespace std;


int main(){

    double pi = 3.1415926;
    double G = 4*pow(pi, 2);
    double M1 = 0.5;
    double M2 = 0.5; // 1

    double mu = G*(M1 + M2);
    double R;
    double e = 0.5;

    string result_name;

    cout << "Please input the result file name : ";
    cin>> result_name;
    cout<< "Please input the R : ";
    cin>> R;

    double v_c_init = init_circle_v(M2, R, G);
    double v_e_init = init_elliptical_v(e, v_c_init);


    double t = 0.0;
    double dt = 0.0005;  // 0.0005
    int t_length = round(864*pow(pow(R, 0.5), 3));  // 864 
    int eq = 6;
    

    double w[eq+1];
    double w_baryc[eq*2+1];

    double v1[eq+1], v2[eq+1], v3[eq+1], v4[eq+1];
    double k1[eq+1], k2[eq+1], k3[eq+1], k4[eq+1];
    int i, j;

    w[1] = R;
    w[2] = 0.0;
    w[3] = 0.0;
    w[4] =  0.0;
    w[5] = v_e_init;
    w[6] =  0.0;

    w_baryc[1] = w[1];
    w_baryc[2] = w[2];
    w_baryc[3] = w[3];
    w_baryc[4] = 0.0;
    w_baryc[5] = 0.0;
    w_baryc[6] = 0.0;
    
    w_baryc[7] = w[4];
    w_baryc[8] = w[5];
    w_baryc[9] = w[6];
    w_baryc[10] = 0.0;
    w_baryc[11] = 0.0;
    w_baryc[12] = 0.0;

    ofstream fout (result_name);
    if (!fout){
        cout<< "\n error" <<endl;
    }

    fout<< scientific;
    fout.precision(8);
    cout << "r1x (wb[1])" << "," << "r1y (wb[2])" << "," << "r1z (wb[3])" << "," 
        << "r2x (wb[4])" << "," << "r2y (wb[5])" << "," << "r2z (wb[6])" << "," 
        << "v1x (wb[7])" << "," << "v1y (wb[8])" << "," << "v1z (wb[9])" << ","
        << "v2x (wb[10])" << "," << "v2y (wb[11])" << "," << "v2z (wb[12])" << "\n";

    for (i = 1; i <= t_length; i++){
        for (j = 1; j <= eq; j++){
            k1[j] = dt*Two_Body_Func(j, t, w, G, mu);
            v1[j] = w[j] + 0.5*k1[j];
        }
        for (j = 1; j <= eq; j++){
            k2[j] = dt*Two_Body_Func(j, t+0.5*dt, v1, G, mu);
            v2[j] = w[j] + 0.5*k2[j];
        }
        for (j = 1; j <= eq; j++){
            k3[j] = dt*Two_Body_Func(j, t+0.5*dt, v2, G, mu);
            v3[j] = w[j] + k3[j];
        }
        for (j = 1; j <= eq; j++){
            k4[j] = dt*Two_Body_Func(j, t+dt, v3, G, mu);
        }
        for (j = 1; j <= eq; j++){
            w[j] = w[j] + (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) / 6.0;
        }



        for (j = 1; j <= eq*2+1; j++){
            w_baryc[j] = Barycentric(j, w, M1, M2, pi);
        }

        t += dt;
        fout << w_baryc[1] << "," << w_baryc[2] << "," << w_baryc[3] << "," << w_baryc[4] << "," 
        << w_baryc[5] << "," << w_baryc[6] << "," << w_baryc[7] << "," << w_baryc[8] << ","
        << w_baryc[9] << "," << w_baryc[10] << "," << w_baryc[11] << "," << w_baryc[12] << "\n";

    }


    cout<< "w[1] = " << w[1] <<endl;
    cout<< "w[2] = " << w[2] <<endl;
    cout<< "w[3] = " << w[3] <<endl;
    cout<< "w[4] = " << w[4] <<endl;
    cout<< "w[5] = " << w[5] <<endl;
    cout<< "w[6] = " << w[6] <<endl;
    cout<< "t = " << t << "\n" <<endl;
    cout<< "w_baryc[1] = " << w_baryc[1] <<endl;
    cout<< "w_baryc[2] = " << w_baryc[2] <<endl;
    cout<< "w_baryc[3] = " << w_baryc[3] <<endl;
    cout<< "w_baryc[4] = " << w_baryc[4] <<endl;
    cout<< "w_baryc[5] = " << w_baryc[5] <<endl;
    cout<< "w_baryc[6] = " << w_baryc[6] <<endl;
    cout<< "w_baryc[7] = " << w_baryc[7] <<endl;
    cout<< "w_baryc[8] = " << w_baryc[8] <<endl;
    cout<< "w_baryc[9] = " << w_baryc[9] <<endl;
    cout<< "w_baryc[10] = " << w_baryc[10] <<endl;
    cout<< "w_baryc[11] = " << w_baryc[11] <<endl;
    cout<< "w_baryc[12] = " << w_baryc[12] <<endl;

    fout.close();

    system("pause");
    return 0;
}

