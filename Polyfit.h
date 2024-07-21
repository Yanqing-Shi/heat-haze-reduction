#pragma once
#ifndef SECOND_FILE_H
#define SECOND_FILE_H
#include<opencv2\opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
using namespace std;
struct output {
	cv::Mat result;
	int ord;
	double coef[5] = { 0 };
	vector<int> xout;
	vector<int> yout;
	int ind;
};

double incbeta(double a, double b, double x);
double invincbeta(double y, double alpha, double beta);
double** Make2DArray(const size_t rows, const size_t cols);
double** MatTrans(double** array, const size_t rows, const size_t cols);
double** MatMul(const size_t m1, const size_t m2, const size_t m3, double** A, double** B);
void MatVectMul(const size_t m1, const size_t m2, double** A, double* v, double* Av);
double determinant(double** a, const size_t k);
void transpose(double** num, double** fac, double** inverse, const size_t r);
void cofactor(double** num, double** inverse, const size_t f);
void PolyFit(const double* x, double* y, const size_t n, const size_t k, const bool fixedinter,
	const double fixedinterval, double* beta, double** Weights, double** XTWXInv);
void CalculateWeights(const double* erry, double** Weights, const size_t n,
	const int type);
output Poly(std::vector<cv::Point> pixels, size_t n,int width,int height);

#endif // SECOND_FILE_H
