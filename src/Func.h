#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;


double init_circle_v(double M1, double R, double G);

double init_elliptical_v(double e, double init_circle_v);

double Two_Body_Func(int eq, double t, double w[], double G, double mu);

double Barycentric(int eq, double w[], double M1, double M2, double pi);
