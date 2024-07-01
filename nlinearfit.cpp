#include <opencv2\opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include <vector>
using namespace cv;
using namespace std;

Mat polyfit(vector<Point>& in_point, int n);

int allow()
{
	//数据输入
	Point in[6] = { Point(0,0),Point(1,0.43445),Point(2,0.99924),Point(4,2.43292),Point(6,4.77895)
		,Point(0.5,0.21723) };

	vector<Point> in_point(begin(in), end(in));

	//n:多项式阶次
	int n = 3;
	Mat mat_k = polyfit(in_point, n);


	//计算结果可视化
	Mat out(150, 500, CV_8UC3, Scalar::all(0));

	//画出拟合曲线
	for (int i = in[0].x; i < in[size(in) - 1].x; ++i)
	{
		Point2d ipt;
		ipt.x = i;
		ipt.y = 0;
		for (int j = 0; j < n + 1; ++j)
		{
			ipt.y += mat_k.at<double>(j, 0) * pow(i, j);
		}
		circle(out, ipt, 1, Scalar(255, 255, 255), FILLED);
	}

	//画出原始散点
	for (int i = 0; i < size(in); ++i)
	{
		Point ipt = in[i];
		circle(out, ipt, 3, Scalar(0, 0, 255), FILLED);
	}

	imshow("9次拟合", out);
	waitKey(0);

	return 0;
}

Mat polyfit(vector<Point>& in_point, int n)
{
	int size = in_point.size();
	//所求未知数个数
	int x_num = n + 1;
	//构造矩阵U和Y
	Mat mat_u(size, x_num, CV_64F);
	Mat mat_y(size, 1, CV_64F);

	for (int i = 0; i < mat_u.rows; ++i)
		for (int j = 0; j < mat_u.cols; ++j)
		{
			mat_u.at<double>(i, j) = pow(in_point[i].x, j);
		}

	for (int i = 0; i < mat_y.rows; ++i)
	{
		mat_y.at<double>(i, 0) = in_point[i].y;
	}

	//矩阵运算，获得系数矩阵K
	Mat mat_k(x_num, 1, CV_64F);
	mat_k = (mat_u.t() * mat_u).inv() * mat_u.t() * mat_y;
	cout << mat_k << endl;
	return mat_k;
}