// am_option.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <math.h>              // mathematical library
#include <cmath>     // c library of math functions

using namespace std;

const double ACCURACY = 1.0e-6;

#ifndef PI 
#define PI 3.141592653589793238462643
#endif

double n(const double& z) {  // normal distribution function    
    return (1.0 / sqrt(2.0 * PI)) * exp(-0.5 * z * z);
};

double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow 
    if (z < -6.0) { return 0.0; };

    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;

    double a = fabs(z);
    double t = 1.0 / (1.0 + a * p);
    double b = c2 * exp((-z) * (z / 2.0));
    double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
    n = 1.0 - b * n;
    if (z < 0.0) n = 1.0 - n;
    return n;
};

double option_price_european_put_payout(const double& S, // spot price
    const double& K, // Strike (exercise) price,
    const double& r,  // interest rate
    const double& q,  // yield on underlying
    const double& sigma,
    const double& time) {
    double sigma_sqr = pow(sigma, 2);
    double time_sqrt = sqrt(time);
    double d1 = (log(S / K) + (r - q + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
    double d2 = d1 - (sigma * time_sqrt);
    double put_price = K * exp(-r * time) * N(-d2) - S * exp(-q * time) * N(-d1);
    return put_price;
};


double option_price_american_put_approximated_baw(const double& S,
    const double& X,
    const double& r,
    const double& b,
    const double& sigma,
    const double& time) {
    const double sigma_sqr = sigma * sigma;
    const double time_sqrt = sqrt(time);
    const double M = 2.0 * r / sigma_sqr;
    const double NN = 2.0 * b / sigma_sqr;
    const double K = 1.0 - exp(-r * time);
    double q1 = 0.5 * (-(NN - 1) - sqrt(pow((NN - 1), 2.0) + (4.0 * M / K)));

    // now find initial S value 
    double q1_inf = 0.5 * (-(NN - 1) - sqrt(pow((NN - 1), 2.0) + 4.0 * M));
    double S_star_star_inf = X / (1.0 - 1.0 / q1_inf);
    double h1 = (b * time - 2 * sigma * time_sqrt) * (X / (X - S_star_star_inf));
    double S_seed = S_star_star_inf + (X - S_star_star_inf) * exp(h1);

    double Si = S_seed;  // now do Newton iterations to solve for S**
    int no_iterations = 0;
    double g = 1;
    double gprime = 1;
    while ((fabs(g) > ACCURACY)
        && (fabs(gprime) > ACCURACY) // to avoid non-convergence
        && (no_iterations++ < 500)
        && Si > 0.0) {
        double p = option_price_european_put_payout(Si, X, r, b, sigma, time);
        double d1 = (log(Si / X) + (b + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
        g = X - Si - p + (1 - exp((b - r) * time) * N(-d1)) * Si / q1;
        gprime = (1.0 / q1 - 1.0) * (1.0 - exp((b - r) * time) * N(-d1))
            + (1.0 / q1) * exp((b - r) * time) * (1.0 / (sigma * time_sqrt)) * n(-d1);
        Si = Si - (g / gprime);
    };
    double S_star_star = Si;
    if (g > ACCURACY) {
        S_star_star = S_seed;
    };  // if not found g**
    double P = 0;
    double p = option_price_european_put_payout(S, X, r, b, sigma, time);
    if (S > S_star_star) {
        double d1 = (log(S_star_star / X)
            + (b + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
        double A1 = -(S_star_star / q1) * (1 - exp((b - r) * time) * N(-d1));
        P = p + A1 * pow((S / S_star_star), q1);
    }
    else {
        P = X - S;
    };
    return max(P, p);  // should not be lower than Black Scholes value
};

inline double phi(double S, double T, double gamma, double H, double X, double r, double b, double sigma) {
    double sigma_sqr = pow(sigma, 2);
    double kappa = 2.0 * b / sigma_sqr + 2.0 * gamma - 1.0;
    double lambda = (-r + gamma * b + 0.5 * gamma * (gamma - 1.0) * sigma_sqr) * T;
    double d1 = -(log(S / H) + (b + (gamma - 0.5) * sigma_sqr) * T) / (sigma * sqrt(T));
    double d2 = -(log((X * X) / (S * H)) + (b + (gamma - 0.5) * sigma_sqr) * T) / (sigma * sqrt(T));
    double phi = exp(lambda) * pow(S, gamma) * (N(d1) - pow((X / S), kappa) * N(d2));
    return phi;
};

double option_price_european_call_payout(const double& S, // spot price
    const double& X, // Strike (exercise) price,
    const double& r,  // interest rate
    const double& q,  // yield on underlying
    const double& sigma, // volatility
    const double& time) { // time to maturity
    double sigma_sqr = pow(sigma, 2);
    double time_sqrt = sqrt(time);
    double d1 = (log(S / X) + (r - q + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
    double d2 = d1 - (sigma * time_sqrt);
    double call_price = S * exp(-q * time) * N(d1) - X * exp(-r * time) * N(d2);
    return call_price;
};

double option_price_american_call_approximated_bjerksund_stensland(const double& S,
    const double& K,
    const double& r,
    const double& b,
    const double& sigma,
    const double& T) {

    double sigma_sqr = pow(sigma, 2);
    double B0 = max(K, (r / (r - b) * K));
    double beta = (0.5 - b / sigma_sqr) + sqrt(pow((b / sigma_sqr - 0.5), 2) + 2.0 * r / sigma_sqr);
    double Binf = beta / (beta - 1.0) * K;
    double hT = -(b * T + 2.0 * sigma * sqrt(T)) * ((K * K) / (Binf - B0));
    double XT = B0 + (Binf - B0) * (1.0 - exp(hT));
    double alpha = (XT - K) * pow(XT, -beta);
    double C = alpha * pow(S, beta);
    C -= alpha * phi(S, T, beta, XT, XT, r, b, sigma);
    C += phi(S, T, 1, XT, XT, r, b, sigma);
    C -= phi(S, T, 1, K, XT, r, b, sigma);
    C -= K * phi(S, T, 0, XT, XT, r, b, sigma);
    C += K * phi(S, T, 0, K, XT, r, b, sigma);
    double c = option_price_european_call_payout(S, K, r, b, sigma, T);
    return max(c, C);


};


double option_price_american_put_approximated_bjerksund_stensland(const double& S,
    const double& X,
    const double& r,
    const double& q,
    const double& sigma,
    const double& T) {

    return option_price_american_call_approximated_bjerksund_stensland(X, S, r - (r - q), r - q, sigma, T);
};


void test_baw_approximation_put(const double& S,
    const double& X,
    const double& r,
    const double& q,
    const double& sigma,
    const double& T) {

    //cout << " Put price using Barone-Adesi Whaley approximation = "
        //<< option_price_american_put_approximated_baw(S, X, r, q, sigma, T) << endl;
};

void approximations_examples() {
    //cout << "------------------------------------" << endl;
    //cout << "Approximations put chapter " << endl;
    //cout << "------------------------------------" << endl;
    double S = 100;   double X = 100;     double sigma = 0.20;
    double r = 0.08;  double b = -0.04;   double time = 0.25;
    double bs = option_price_american_put_approximated_bjerksund_stensland(S, X, r, b, sigma, time);

    //cout << " Put price using bjerksund_stensland approximation = "
    //    << bs << endl;

    test_baw_approximation_put(S, X, r, b, sigma, time);
};


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
