#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

///////////////////Constants////////////////////
constexpr double G = 9.80665;
constexpr double Alpha = 0.33333333333333333;       //?????????
//constexpr double Beta = 0.3;        //?????????
constexpr double eps = 0.000000000001;    //epsilon = 10^-12

///////////////////Grid_on_Ox///////////////////
constexpr double DELTA_X = 0.1;  // step on Ox
constexpr double LENGTH = 10;    // length of Ox
constexpr double X0 = LENGTH / 2;  // location of "barrier"
constexpr int X_STEPS = static_cast<int>(LENGTH / DELTA_X);  // number of division points

vector<double> X(X_STEPS), H(X_STEPS), U(X_STEPS);

//////////////////Output_files//////////////////
const string PATH_H = "heights_2.csv";
const string PATH_U = "speeds_2.csv";

/////////////csv-oriented_functions/////////////
void init_file(const string& path) {
    ofstream f(path);
    f << "time";
    for (const double x : X) {
        f << "," << x;
    }
    f << endl;
    f.close();
}

void add_string(const string& path, const vector<double>& data, const double t) {
    ofstream f(path, ios::app);
    f << t;
    for (const double value : data) {
        f << "," << value;
    }
    f << endl;
    f.close();
}

////////////math-oriented_functions/////////////
double h_half(const int x) {
    return (H[x] + H[x + 1]) / 2;
}

double u_half(const int x) {
    return (U[x] + U[x + 1]) / 2;
}

////////////WENO_math///////////////////////////

void polynom_count(const int x,  double &h2,  double &h1, double &h0, double &u2, double &u1, double &u0){
    //x  - number of left cell/dot; we count(т.е. подаем i-1, где i - номер центральной ячейки из 3-х)
    //выдает полиномы высоты и скорости по 3м точкам
    double x0, x1, x2;
    x0 = X[x] + DELTA_X/2;
    x1 = x0 + DELTA_X;
    x2 = x1 + DELTA_X;
    
    h2 = (h_half(x) - 2.0 * h_half(x+1) + h_half(x+2))/(2*DELTA_X*DELTA_X);
    u2 = (u_half(x) - 2.0 * u_half(x+1) + u_half(x+2))/(2*DELTA_X*DELTA_X);
    
    h1 = (-1)*(h_half(x) * (x1 + x2) - 2.0 * h_half(x+1) * (x0 + x2) + h_half(x+2) * (x0 + x1)) / (2*DELTA_X*DELTA_X);
    u1 = (-1)*(u_half(x) * (x1 + x2) - 2.0 * u_half(x+1) * (x0 + x2) + u_half(x+2) * (x0 + x1)) / (2*DELTA_X*DELTA_X);
    //what's right?
    //h1 = (h_half(x) * (x1 + x2) - 2.0 * h_half(x+1) * (x0 + x2) + h_half(x+2) * (x0 + x1)) / (2*DELTA_X*DELTA_X);
    //u1 = (u_half(x) * (x1 + x2) - 2.0 * u_half(x+1) * (x0 + x2) + u_half(x+2) * (x0 + x1)) / (2*DELTA_X*DELTA_X);

    h0 = (h_half(x) * x1*x2 - 2.0 * h_half(x+1) * x0*x2 + h_half(x+2) * x0*x1)/(2*DELTA_X*DELTA_X);
    u0 = (u_half(x) * x1*x2 - 2.0 * u_half(x+1) * x0*x2 + u_half(x+2) * x0*x1)/(2*DELTA_X*DELTA_X);
}

double h_center(const int x, const double hl0, const double hl1, const double hl2, const double hc0, const double hc1, const double hc2, const double hr0, const double hr1, const double hr2){
    double x1, x2, hl, hc, hr;
    x1 = X[x] + DELTA_X/2;        
    x2 = x1 * x1;
    //weights: 3/10; 3/5; 1/10;
    //cout << "\n" << x << ": ";
    hl = (x2*hl2 + x1*hl1 + hl0)*0.3;
    //cout << " hl: " << hl;
    hc = (x2*hc2 + x1*hc1 + hc0)*0.6;
    //cout << "  hc: " << hc;
    hr = (x2*hr2 + x1*hr1 + hr0)*0.1;
    //cout << "  hr: " << hr;
    return (hl + hc + hr);      //возвращает взвешенное центральное значение для ячейки
}

double u_center(const int x, const double ul0, const double ul1, const double ul2, const double uc0, const double uc1, const double uc2, const double ur0, const double ur1, const double ur2){
    double x1, x2, ul, uc, ur;
    x1 = X[x] + DELTA_X/2;        
    x2 = x1 * x1;
    //weights: 3/10; 3/5; 1/10;
    ul = (x2*ul2 + x1*ul1 + ul0)*0.3; 
    uc = (x2*uc2 + x1*uc1 + uc0)*0.6;
    ur = (x2*ur2 + x1*ur1 + ur0)*0.1;
    
    return (ul + uc + ur);      //возвращает взвешенное центральное значение для ячейки
}
//
double hc_new(int x, double dt, double h_c)
{
    //cout << (h_c - dt*(H[x+1]*U[x+1] - H[x]*U[x])/DELTA_X);    
    return h_c - dt*(H[x+1]*U[x+1] - H[x]*U[x])/DELTA_X;  //новое центральное значение для ячейки
}

double uc_new(int x, double dt, double h_c, double u_c, double h_c_new)
{
    double  uh_c;
    uh_c = h_c*u_c - dt*(H[x+1]*U[x+1]*U[x+1] - H[x]*U[x]*U[x])/DELTA_X - dt*G*(H[x+1]*H[x+1] - H[x]*H[x])/(2*DELTA_X);  
    if (h_c_new >= eps)
        return (uh_c/h_c_new);  //новое центральное значение для ячейки
    else
        return 0.0;     //нет волны - не делим на 0 и ставим скорость равной 0
    
}
////////////////////////////////////////////////
int main() {
    // Initialize grid on Ox
    for (int i = 0; i < X_STEPS; ++i) {
        X[i] = i * DELTA_X;
    }

    // Initialize output files
    init_file(PATH_H);
    init_file(PATH_U);

    // Get initial conditions
    double h1, u1, h2, u2;
    //cout << "Enter h1: ";
    //cin >> h1;
    h1 = 2.0;
    //cout << "Enter u1: ";
    //cin >> u1;
    u1 = 0;
    //cout << "Enter h2: ";
    //cin >> h2;
    h2 = 1.0;
    //cout << "Enter u2: ";
    //cin >> u2;
    u2 = 0;
        
    // Initialize H and U
    for (int i = 0; i < static_cast<int>(X0 / DELTA_X); ++i) {
        H[i] = h1;
        U[i] = u1;
    }
    for (int i = static_cast<int>(X0 / DELTA_X); i < X_STEPS; ++i) {
        H[i] = h2;
        U[i] = u2;
    }

    // Simulation loop
    double t;
    cout << "Enter time: ";
    //cin >> t;
    t = 0.5;
    double t_cur = 0;

    char c;

    while (t > 0) {
        double dt = t + 1;
        double buf = h1;
        //
        for (int x = 0; x < X_STEPS; ++x)   {
            if (H[x] != 0) {
                dt = min(dt, Alpha * DELTA_X / ( abs(U[x]) + sqrt(G * H[x]) ));
            }
        }

        vector<double> H_new(X_STEPS), U_new(X_STEPS), H_centr(X_STEPS), U_centr(X_STEPS);
        //H.insert(H.begin(), H.front());
        //H.push_back(H.back());
        //U.insert(U.begin(), U.front());
        //U.push_back(U.back());
        double h[9], u[9];
        //отсюда значения меняются
        for (int i = 0; i<9; i++){
            h[i] = 0.0;
            u[i] = 0.0;
        }

        H_new[0] = h1;
        H_new[1] = h1;
        U_new[0] = u1;
        U_new[1] = u1;
        
        H_new[X_STEPS] = h2;
        H_new[X_STEPS-1] = h2;
        //H_new[X_STEPS-2] = h2;
        U_new[X_STEPS] = u2;
        U_new[X_STEPS-1] = u2;
        //U_new[X_STEPS-2] = u2;
        
        //continue changes from here
        H_centr[0] = h1;    //1-я ячейка из X_STEPS
        H_centr[1] = h1;    //2-я ...
        U_centr[0] = u1;    //1-я ...
        U_centr[1] = u1;    //2-я ...
        
        for (int x = 0; x <= (X_STEPS-4); ++x) {
            //if (((h1 - h_half(x)) < eps)&&(abs(u_half(x) - u1) < eps))       //стартовые левые ячейки
            //{
            //    if (((h1 - h_half(x+4)) > eps)||((abs(u_half(x+4) - u1)) > eps))  //ячейки отличаются от стартовых левых => будут изменения
            //    {
                    //to not recount again and again same polynoms:                    
                    /*if ((abs(h[0]) + abs(h[1]) + abs(h[2])) == 0)
                        polynom_count(x-1, h[2], h[1], h[0], u[2], u[1], u[0]);
                    if ((abs(h[3]) + abs(h[4]) + abs(h[5])) == 0)
                        polynom_count(x, h[5], h[4], h[3], u[5], u[4], u[3]);
                    if ((abs(h[6]) + abs(h[7]) + abs(h[8])) == 0)
                        polynom_count(x+1, h[8], h[7], h[6], u[8], u[7], u[6]);
                    */                    
                    polynom_count(x+0, h[2], h[1], h[0], u[2], u[1], u[0]);
                    polynom_count(x+1, h[5], h[4], h[3], u[5], u[4], u[3]);
                    polynom_count(x+2, h[8], h[7], h[6], u[8], u[7], u[6]);
                    //получили полиномы для набора из 5 ячеек

                    H_centr[x+2] = h_center(x+2, h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], h[8]); //пересчёт центральных значений по полиномам
                    //cout << "\n" << H_centr[x+2] << "\n";
                    //cin >> c;
                    U_centr[x+2] = u_center(x+2, u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]); //подаём номер левого значения ячейки    
                    
                    /*
                    buf = H_centr[x+2];
                    H_centr[x+2] = hc_new(x+1, dt, buf);
                    cout << "\n" << H_centr[x+2] << "\n";
                    U_centr[x+2] = uc_new(x+1, dt, buf, U_centr[x+2], H_centr[x+2]);
                    
                    H_new[x+1] = (H_centr[x+1] + H_centr[x+2])/2;   //можно было бы пересчитывать всё заново по полиномам, но смысла мало
                    U_new[x+1] = (U_centr[x+1] + U_centr[x+2])/2;
                    
                    for (int i = 0; i < 3; i++){
                        h[i] = h[i+3];
                        h[i+3] = h[i+6];
                        h[i+6] = 0.0;
    
                        u[i] = u[i+3];
                        u[i+3] = u[i+6];
                        u[i+6] = 0.0;
                    }
                    */                    
                    //cчитаем новые значения в отдельном цикле
                    
        }
        
        for (int x = 0; x <= (X_STEPS-4); ++x) {
            H_new[x+2] = H[x+2] - dt*(H_centr[x+2]*U_centr[x+2] - H_centr[x+1]*U_centr[x+1])/DELTA_X;
            cout << "\n" << H_new[x+2] << "\n";
            U_new[x+2] = (H[x+2]*U[x+2] - dt*( H_centr[x+2]*U_centr[x+2]*U_centr[x+2] - H_centr[x+1]*U_centr[x+1]*U_centr[x+1] + G*( H_centr[x+2]*H_centr[x+2] - H_centr[x+1]*H_centr[x+1] )/2 ) / DELTA_X) / H_new[x+2];
        }
        
        
        H = move(H_new);
        U = move(U_new);

        t -= dt;

        t_cur += dt;
        add_string(PATH_H, H, t_cur);
        add_string(PATH_U, U, t_cur);
    }

    return 0;
}


