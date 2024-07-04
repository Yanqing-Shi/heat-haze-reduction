#pragma once
#include<iostream>

#include <opencv2/opencv.hpp>
#include "WindowZoom.h"
using namespace std;
using namespace cv;



bool comparePoints(const Point& p1, const Point& p2);

Vec3b average(int x, int y, Mat img, int flag);
Mat manipulation(double slope, Mat img, Mat canny, Mat post_canny, int x, int fx);

MouseParams* prep(void* param);


