//#include<opencv2\opencv.hpp>
//#include <opencv2/core.hpp>
//#include <opencv2/imgcodecs.hpp>
//#include <opencv2/highgui.hpp>
//#include<opencv2/imgproc/imgproc.hpp>
//#include<stdio.h>
//#include<algorithm>
//#include <iostream>
//#include "Polyfit.h"
//
//using namespace cv;
//using namespace std;
//
//bool drawing = false;
////int size = 5;
//vector<Point> points;
//Point temp;
//
//
//
//struct MouseParams
//{
//	Mat* img;
//	Mat mask;
//	int width;
//	int height;
//};
//
//
//
//
//void edit(void* param) {
//	MouseParams* image = reinterpret_cast<MouseParams*>(param);
//
//	Mat* img = image->img;
//	Mat mask = image->mask;
//	int height = image->height;
//	int width = image->width;
//	Mat img_canny;
//	Canny(*img, img_canny, 50, 100);
//	Mat canny_result;
//	bitwise_and(mask, img_canny, canny_result);
//
//	//imwrite("starry_night.jpg", canny_result);
//	//cout << height << " " << width << endl;
//	vector<Point> whitePixels;
//	int length = 0;
//	// Traverse the image
//	for (int i = 0; i < height; ++i) {
//		// Get pointer to the i-th row
//		uchar* row = canny_result.ptr<uchar>(i);
//
//		for (int j = 0; j < width; ++j) {
//			// Check if the pixel is white
//			if (row[j] == 255) {
//				// Store the coordinates
//				length++;
//				whitePixels.push_back(cv::Point(j, i));
//			}
//		}
//	}
//	/*for (auto elem : whitePixels) {
//		cout << elem << endl;
//	}*/
//	//Poly(whitePixels,length);
//
//	Vec3b intensity = (*img).at<Vec3b>(Point(image->width-1, image->height-1));
//	int blue = intensity.val[0];
//	int green = intensity.val[1];
//	int red = intensity.val[2];
//	//cout << blue << " " << green << " " << red << endl;
//}
//
//void onmouse(int event, int x, int y, int flag, void* param) {
//	
//	
//	MouseParams* image = reinterpret_cast<MouseParams*>(param);
//	Mat* img = image->img;
//	Mat mask = image->mask;
//	//cout << x << " " << y << endl;
//	
//	if (event==EVENT_LBUTTONDOWN&&drawing==false) {
//		drawing = true;
//	}
//	else if(event==EVENT_MOUSEMOVE && flag==EVENT_FLAG_LBUTTON){
//		
//		Point coor;
//		coor=Point(x,y);
//		if (coor != points.back()) points.push_back(coor);
//		circle(mask, Point(x, y), 5, 255, FILLED);
//		//waitKey(50);
//	}
//	else if (event == EVENT_LBUTTONUP) {
//		drawing = false;
//		
//		image->mask = mask;
//		//edit(image);
//		
//		/*for (auto elem : points) {
//			cout << elem << endl;
//		}*/
//		//cout <<endl << points.back() << endl;
//	}
//	else if (event == EVENT_MOUSEWHEEL) {
//		int xcoor = max(int(x - (*img).cols * 0.45), 0);
//		xcoor = min(xcoor, int(0.1 * (*img).cols));
//		int ycoor = max(int(y - (*img).rows * 0.45), 0);
//		xcoor = min(ycoor, int(0.1 * (*img).rows));
//		//cout << xcoor << " " << ycoor;
//		Rect span(xcoor, ycoor, (*img).cols*0.9, (*img).rows * 0.9);
//		Mat newimg = (*img)(span);
//		*img = newimg;
//		imshow("image", *img);
//	}
//	//cout << flag << endl;
//}
//
//
//
//int kk()
//{
//	MouseParams image;
//	temp=Point(-1, -1);
//	points.push_back(temp);
//	Mat img = imread("DSC_1744.jpg", IMREAD_COLOR);
//	uint8_t* myData = img.data;
//	Mat mask= Mat::zeros(Size(img.cols, img.rows), CV_8UC1);
//	Mat out;
//	image.img = &img;
//	image.width = img.cols; //1280
//	image.height = img.rows; //854
//	image.mask = mask;
//	//cout << image.height << image.width;
//	/*for (int i = 0; i < height; i++)
//	{
//		for (int j = 0; j < width; j++)
//		{
//			uint8_t val = myData[i * _stride + j];
//			printf("%d ", val);
//		}
//		printf("\n");
//	}*/
//	//imwrite("starry_night.png", img);
//
//	imshow("image", img);
//	setMouseCallback("image", onmouse, reinterpret_cast<void*>(&image));
//	waitKey(0);
//	/*
//	Vec3b intensity = img.at<Vec3b>(Point(1023,684));
//	int blue = intensity.val[0];
//	int green = intensity.val[1];
//	int red = intensity.val[2];
//	cout << blue << " " << green << " " << red << endl;
//	*/
//	//Canny(img,out,50,100);
//	//imwrite("starry_night.jpg", out);
//	//for (int j = 1; j <= size; j++) {ll
//	//	printf("%d %d\n", cl[size][0], cl[size][1]);
//	//}
//	
//	return 0;
//}