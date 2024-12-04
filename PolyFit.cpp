// ********************************************************************
// * Code PolyFit                                                     *
// * Written by Ianik Plante                                          *
// *                                                                  *
// * KBR                                                              *
// * 2400 NASA Parkway, Houston, TX 77058                             *
// * Ianik.Plante-1@nasa.gov                                          *
// *                                                                  *
// * This code is used to fit a series of n points with a polynomial  *
// * of degree k, and calculation of error bars on the coefficients.  *
// * If error is provided on the y values, it is possible to use a    *
// * weighted fit as an option. Another option provided is to fix the *
// * intercept value, i.e. the first parameter.                       *
// *                                                                  *
// * This code has been written partially using data from publicly    *
// * available sources.                                               *
// *                                                                  *  
// * The code works to the best of the author's knowledge, but some   *   
// * bugs may be present. This code is provided as-is, with no        *
// * warranty of any kind. By using this code you agree that the      * 
// * author, the company KBR or NASA are not responsible for possible *
// * problems related to the usage of this code.                      * 
// *                                                                  *   
// * The program has been reviewed and approved by export control for * 
// * public release. However some export restriction exists. Please   *    
// * respect applicable laws.                                         *
// *                                                                  *   
// ********************************************************************


#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <iomanip>
#include<opencv2\opencv.hpp>
#include"polyfit.h"

using namespace std;
using namespace cv;

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30



/*
 * zlib License
 *
 * Regularized Incomplete Beta Function
 *
 * Copyright (c) 2016, 2017 Lewis Van Winkle
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

#define STOP 1.0e-8
#define TINY 1.0e-30

 // Adapted from https://github.com/codeplea/incbeta
double incbeta(double a, double b, double x) {
    if (x < 0.0 || x > 1.0) return 1.0;

    if (a <= 0.) {
        std::cout << "Warning: a should be >0";
        return 0.;
    }

    if (b <= 0.) {
        std::cout << "Warning: b should be >0";
        return 0.;
    }


    /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
    if (x > (a + 1.0) / (a + b + 2.0)) {
        return (1.0 - incbeta(b, a, 1.0 - x)); /*Use the fact that beta is symmetrical.*/
    }

    /*Find the first part before the continued fraction.*/
    const double lbeta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);
    const double front = exp(log(x) * a + log(1.0 - x) * b - lbeta_ab) / a;

    /*Use Lentz's algorithm to evaluate the continued fraction.*/
    double f = 1.0, c = 1.0, d = 0.0;

    int i, m;
    for (i = 0; i <= 200; ++i) {
        m = i / 2;

        double numerator;
        if (i == 0) {
            numerator = 1.0; /*First numerator is 1.0.*/
        }
        else if (i % 2 == 0) {
            numerator = (m * (b - m) * x) / ((a + 2.0 * m - 1.0) * (a + 2.0 * m)); /*Even term.*/
        }
        else {
            numerator = -((a + m) * (a + b + m) * x) / ((a + 2.0 * m) * (a + 2.0 * m + 1)); /*Odd term.*/
        }

        /*Do an iteration of Lentz's algorithm.*/
        d = 1.0 + numerator * d;
        if (fabs(d) < TINY) d = TINY;
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if (fabs(c) < TINY) c = TINY;

        const double cd = c * d;
        f *= cd;

        /*Check for stop.*/
        if (fabs(1.0 - cd) < STOP) {
            return front * (f - 1.0);
        }
    }

    return 1.0; /*Needed more loops, did not converge.*/
}

double invincbeta(double y, double alpha, double beta) {

    if (y <= 0.) return 0.;
    else if (y >= 1.) return 1.;
    if (alpha <= 0.) {
        std::cout << "Warning: alpha should be >0";
        return 0.;
    }

    if (beta <= 0.) {
        std::cout << "Warning: beta should be >0";
        return 0.;
    }


    double x = 0.5;
    double a = 0;
    double b = 1;
    double precision = 1.e-8;
    double binit = y;
    double bcur = incbeta(alpha, beta, x);

    while (fabs(bcur - binit) > precision) {

        if ((bcur - binit) < 0) {
            a = x;
        }
        else {
            b = x;
        }
        x = (a + b) * 0.5;
        bcur = incbeta(alpha, beta, x);

        //std::cout << x << "\t" << bcur << "\n";


    }

    return x;


}



// Calculate the t value for a Student distribution
// Adapted from http://www.cplusplus.com/forum/beginner/216098/
// **************************************************************
double CalculateTValueStudent(const double nu, const double alpha) {

    double precision = 1.e-5;

    if (alpha <= 0. || alpha >= 1.) return 0.;

    double x = invincbeta(2. * min(alpha, 1. - alpha), 0.5 * nu, 0.5);
    x = sqrt(nu * (1. - x) / x);
    return (alpha >= 0.5 ? x : -x);


}

// Cumulative distribution for Student-t
// **************************************************************
double cdfStudent(const double nu, const double t)
{
    double x = nu / (t * t + nu);

    return 1. - incbeta(0.5 * nu, 0.5, x);
}

// Cumulative distribution for Fisher F
// **************************************************************
double cdfFisher(const double df1, const double df2, const double x) {
    double y = df1 * x / (df1 * x + df2);
    return incbeta(0.5 * df1, 0.5 * df2, y);
}

// Initialize a 2D array
// **************************************************************
double** Make2DArray(const size_t rows, const size_t cols) {

    double** array;

    array = new double* [rows];
    for (size_t i = 0; i < rows; i++) {
        array[i] = new double[cols];
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            array[i][j] = 0.;
        }
    }

    return array;

}

// Transpose a 2D array
// **************************************************************
double** MatTrans(double** array, const size_t rows, const size_t cols) {

    double** arrayT = Make2DArray(cols, rows);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            arrayT[j][i] = array[i][j];
        }
    }

    return arrayT;

}

// Perform the multiplication of matrix A[m1,m2] by B[m2,m3]
// **************************************************************
double** MatMul(const size_t m1, const size_t m2, const size_t m3, double** A, double** B) {

    double** array = Make2DArray(m1, m3);

    for (size_t i = 0; i < m1; i++) {
        for (size_t j = 0; j < m3; j++) {
            array[i][j] = 0.;
            for (size_t m = 0; m < m2; m++) {
                array[i][j] += A[i][m] * B[m][j];
            }
        }
    }
    return array;

}

// Perform the multiplication of matrix A[m1,m2] by vector v[m2,1]
// **************************************************************
void MatVectMul(const size_t m1, const size_t m2, double** A, double* v, double* Av) {


    for (size_t i = 0; i < m1; i++) {
        Av[i] = 0.;
        for (size_t j = 0; j < m2; j++) {
            Av[i] += A[i][j] * v[j];
        }
    }


}


// Calculates the determinant of a matrix 
// **************************************************************
double determinant(double** a, const size_t k) {

    double s = 1;
    double det = 0.;
    double** b = Make2DArray(k, k);
    size_t m;
    size_t n;

    if (k == 1) return (a[0][0]);

    for (size_t c = 0; c < k; c++) {

        m = 0;
        n = 0;

        for (size_t i = 0; i < k; i++) {

            for (size_t j = 0; j < k; j++) {

                b[i][j] = 0;

                if (i != 0 && j != c) {

                    b[m][n] = a[i][j];
                    if (n < (k - 2)) {
                        n++;
                    }
                    else
                    {
                        n = 0;
                        m++;
                    }
                }
            }
        }

        det = det + s * (a[0][c] * determinant(b, k - 1));
        s = -1 * s;

    }

    return (det);

}


// Perform the 
// **************************************************************
void transpose(double** num, double** fac, double** inverse, const size_t r) {

    double** b = Make2DArray(r, r);
    double deter;

    for (size_t i = 0; i < r; i++) {
        for (size_t j = 0; j < r; j++) {
            b[i][j] = fac[j][i];
        }
    }

    deter = determinant(num, r);

    for (size_t i = 0; i < r; i++) {
        for (size_t j = 0; j < r; j++) {
            inverse[i][j] = b[i][j] / deter;
        }
    }

}

// Calculates the cofactors 
// **************************************************************
void cofactor(double** num, double** inverse, const size_t f)
{

    double** b = Make2DArray(f, f);
    double** fac = Make2DArray(f, f);

    size_t m;
    size_t n;

    for (size_t q = 0; q < f; q++) {

        for (size_t p = 0; p < f; p++) {

            m = 0;
            n = 0;

            for (size_t i = 0; i < f; i++) {

                for (size_t j = 0; j < f; j++) {

                    if (i != q && j != p) {

                        b[m][n] = num[i][j];

                        if (n < (f - 2)) {
                            n++;
                        }
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
        }
    }

    transpose(num, fac, inverse, f);

}






// Perform the fit of data n data points (x,y) with a polynomial of order k
// **************************************************************
void PolyFit(const double* x, double* y, const size_t n, const size_t k, const bool fixedinter,
    const double fixedinterval, double* beta, double** Weights, double** XTWXInv) {

    // Definition of variables
    // **************************************************************
    double** X = Make2DArray(n, k + 1);           // [n,k+1]
    double** XT;                               // [k+1,n]
    double** XTW;                              // [k+1,n]
    double** XTWX;                             // [k+1,k+1]

    double* XTWY = new double[k + 1];
    double* Y = new double[n];

    size_t begin = 0;
    if (fixedinter) begin = 1;

    // Initialize X
    // **************************************************************
    for (size_t i = 0; i < n; i++) {
        for (size_t j = begin; j < (k + 1); j++) {  // begin
            X[i][j] = pow(x[i], j);
        }
    }

    // Matrix calculations
    // **************************************************************
    XT = MatTrans(X, n, k + 1);                 // Calculate XT
    XTW = MatMul(k + 1, n, n, XT, Weights);         // Calculate XT*W
    XTWX = MatMul(k + 1, n, k + 1, XTW, X);           // Calculate (XTW)*X

    if (fixedinter) XTWX[0][0] = 1.;

    cofactor(XTWX, XTWXInv, k + 1);             // Calculate (XTWX)^-1

    for (size_t m = 0; m < n; m++) {
        if (fixedinter) {
            Y[m] = y[m] - fixedinterval;
        }
        else {
            Y[m] = y[m];
        }
    }
    MatVectMul(k + 1, n, XTW, Y, XTWY);             // Calculate (XTW)*Y
    MatVectMul(k + 1, k + 1, XTWXInv, XTWY, beta);    // Calculate beta = (XTWXInv)*XTWY

    if (fixedinter) beta[0] = fixedinterval;

}


// Calculate the weights matrix
// **************************************************************
void CalculateWeights(const double* erry, double** Weights, const size_t n,
    const int type) {


    for (size_t i = 0; i < n; i++) {

        switch (type) {
        case 0:
            Weights[i][i] = 1.;
            break;
        case 1:
            Weights[i][i] = erry[i];
            break;
        case 2:
            if (erry[i] > 0.) {
                Weights[i][i] = 1. / (erry[i] * erry[i]);
            }
            else {
                Weights[i][i] = 0.;
            }
            break;
        }

    }

}



// Display the polynomial
// **************************************************************
void DisplayPolynomial(const size_t k) {

    cout << "y = ";
    for (size_t i = 0; i < (k + 1); i++) {
        cout << "A" << i;
        if (i > 0) cout << "X";
        if (i > 1) cout << "^" << i;
        if (i < k) cout << " + ";
    }
    cout << endl << endl;

}



// The main program
// **************************************************************
output Poly(vector<Point> pixels, size_t n, int width, int height) {

    output out;
    int book[8000] = { 0 };
    // Input values
    // **************************************************************
    size_t k = 2;
    // Polynomial order
    bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
    int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
    double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
    double alphaval = 0.05;                          // Critical apha value

    double* x = new double[n];
    double* y = new double[n];
    double* erry = new double[n];

    for (int i = 0; i < n; i++) {
        x[i] = pixels[i].x;
        y[i] = pixels[i].y;
        erry[i] = 0;
    }


    //    double x[] = { 490, 758, 305, 915, 746, 677, 812, 256, 973, 860,
    //448, 745, 42, 810, 71, 107, 665, 729, 721, 508,
    //106, 287, 966, 588, 756, 359, 752, 63, 994, 740,
    //967, 693, 999, 10, 67, 538, 911, 196, 533, 120,
    //23, 132, 327, 202, 395, 541, 972, 267, 699, 872,
    //686, 496, 943, 941, 88, 15, 412, 780, 249, 313,
    //460, 864, 133, 657, 609, 454, 726, 258, 2, 456,
    //241, 649, 350, 785, 906, 951, 337, 504, 790, 268,
    //360, 279, 789, 832, 205, 52, 597, 180, 168, 130,
    //489, 673, 859, 883, 86, 363, 406, 617, 271, 757 };
    //    double y[] = { 535, 339, 488, 226, 330, 979, 943, 506, 314, 681,
    //457, 284, 921, 332, 961, 865, 414, 28, 433, 784,
    //869, 115, 166, 274, 972, 860, 34, 447, 174, 246,
    //904, 65, 856, 76, 703, 803, 376, 601, 187, 232,
    //814, 409, 37, 177, 513, 286, 516, 966, 462, 491,
    //655, 411, 336, 500, 486, 621, 293, 568, 859, 123,
    //877, 745, 882, 439, 474, 193, 620, 397, 850, 393,
    //929, 361, 483, 637, 392, 191, 82, 317, 273, 518,
    //755, 477, 261, 740, 648, 693, 536, 626, 47, 576,
    //327, 509, 606, 750, 776, 870, 58, 741, 25, 145 };
    //    double erry[100] = {0};       // Data points (err on y) (if applicable)

        // Definition of other variables
        // **************************************************************

    double* coefbetax = new double[k + 1];                            // Coefficients of the polynomial
    double* coefbetay = new double[k + 1];
    double* serbeta = new double[k + 1];                             // Standard error on coefficients
    double tstudentval = 0.;                         // Student t value
    double SE = 0.;                                  // Standard error

    double** XTWXInv;                                // Matrix XTWX Inverse [k+1,k+1]
    double** Weights;                                // Matrix Weights [n,n]


    // Initialize values
    // **************************************************************

    XTWXInv = Make2DArray(k + 1, k + 1);
    Weights = Make2DArray(n, n);
    DisplayPolynomial(k);
    // Build the weight matrix
    // **************************************************************
    CalculateWeights(erry, Weights, n, wtype);

    vector<int>x_xout, x_yout, y_xout, y_yout;

    // Calculate the coefficients of the fit
    // **************************************************************
    PolyFit(x, y, n, k, fixedinter, fixedinterval, coefbetax, Weights, XTWXInv);
    //DisplayCoefs(k, nstar, tstudentval, coefbeta, serbeta);
    int j = 0;
    double sdx = 0;
    Mat resultx = Mat::zeros(Size(width, height), CV_8UC1);
    for (int i = 0; i < n; i++) {

        if (book[int(x[i])] == 0) {
            //cout << x[i] << endl;
            book[int(x[i])] = 1;
            resultx.at<uchar>(coefbetax[0] + x[i] * coefbetax[1] + pow(x[i], 2) * coefbetax[2], x[i]) = 255;
            sdx += pow(coefbetax[0] + x[i] * coefbetax[1] + pow(x[i], 2) * coefbetax[2] - y[i], 2);
            x_xout.push_back(int(x[i]));
            x_yout.push_back(int(coefbetax[0] + x[i] * coefbetax[1] + pow(x[i], 2) * coefbetax[2]));
            //out.yout.push_back();
            j++;
        }
        else {

        }

        //cout << x[i] << " ";
    }
    cout <<"sdx="<< sdx << endl;


    PolyFit(y, x, n, k, fixedinter, fixedinterval, coefbetay, Weights, XTWXInv);


    int book2[8000] = { 0 };

    Mat resulty = Mat::zeros(Size(width, height), CV_8UC1);
    double sdy = 0;
    for (int i = 0; i < n; i++) {

        if (book2[int(y[i])] == 0) {
            //cout << x[i] << endl;
            book2[int(y[i])] = 1;
            resulty.at<uchar>(y[i], coefbetay[0] + y[i] * coefbetay[1] + pow(y[i], 2) * coefbetay[2]) = 255;
            sdy += pow(coefbetay[0] + y[i] * coefbetay[1] + pow(y[i], 2) * coefbetay[2] - x[i], 2);
            y_xout.push_back(int(y[i]));
            y_yout.push_back(int(coefbetay[0] + y[i] * coefbetay[1] + pow(y[i], 2) * coefbetay[2]));
            j++;
        }


    }
    cout << "sdy=" << sdy << endl;

    if (sdx < sdy) {
        out.ind = 0;  //x is the independent var
        out.result = resultx;
        for (size_t i = 0; i < (k + 1); i++) {
            out.coef[i] = coefbetax[i];
        }
        out.xout = x_xout;
        out.yout = x_yout;
    }
    else {
        out.ind = 1;  //y is the independent var
        out.result = resulty;
        for (size_t i = 0; i < (k + 1); i++) {
            out.coef[i] = coefbetay[i];
        }
        out.xout = y_xout;
        out.yout = y_yout;
    }

    out.ord = k;


    delete[]coefbetax;
    coefbetax = NULL;
    delete[]coefbetay;
    coefbetay = NULL;
    delete[]serbeta;
    serbeta = NULL;

    return out;
}