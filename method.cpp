#include"method.h"
#include"Polyfit.h"

using namespace std;
using namespace cv;


bool comparePoints(const Point& p1, const Point& p2) {
	return p1.x < p2.x;
}

Vec3b average(int x, int y, Mat img, int flag) {
	Vec3b output;
	if (flag == 1) {
		//cout << x << " " << y << endl;
		Vec3b valueLeft = img.at<Vec3b>(y, x - 2);
		Vec3b valueMidLeft = img.at<Vec3b>(y, x - 1);
		Vec3b valueCenter = img.at<Vec3b>(y, x);
		Vec3b valueMidRight = img.at<Vec3b>(y, x + 1);
		Vec3b valueRight = img.at<Vec3b>(y, x + 2);
		//Mat kernel = Mat::ones(3, 3, CV_64F) / 9;
		double avgB = (valueLeft[0] + valueMidLeft[0] + valueCenter[0] + valueMidRight[0] + valueRight[0]) / 5.0;
		double avgG = (valueLeft[1] + valueMidLeft[1] + valueCenter[1] + valueMidRight[1] + valueRight[1]) / 5.0;
		double avgR = (valueLeft[2] + valueMidLeft[2] + valueCenter[2] + valueMidRight[2] + valueRight[2]) / 5.0;
		output = Vec3b(avgB, avgG, avgR);
	}
	else if (flag == 2) {
		Vec3b valueLeft = img.at<Vec3b>(y - 2, x);
		Vec3b valueMidLeft = img.at<Vec3b>(y - 1, x);
		Vec3b valueCenter = img.at<Vec3b>(y, x);
		Vec3b valueMidRight = img.at<Vec3b>(y + 1, x);
		Vec3b valueRight = img.at<Vec3b>(y + 2, x);
		//Mat kernel = Mat::ones(3, 3, CV_64F) / 9;
		double avgB = (valueLeft[0] + valueMidLeft[0] + valueCenter[0] + valueMidRight[0] + valueRight[0]) / 5.0;
		double avgG = (valueLeft[1] + valueMidLeft[1] + valueCenter[1] + valueMidRight[1] + valueRight[1]) / 5.0;
		double avgR = (valueLeft[2] + valueMidLeft[2] + valueCenter[2] + valueMidRight[2] + valueRight[2]) / 5.0;
		output = Vec3b(avgB, avgG, avgR);
	}


	return output;
}



Mat manipulation(double slope, Mat img, Mat canny, Mat post_canny, int x, int fx) {
	//cout << slope << endl;
	if (abs(slope) >= 4.5) {
		//cout << x << " "<<fx << endl;
		for (int y = fx - 10; y < fx + 10; y++) {
			if (canny.at<uchar>(y, x) == 255) {
				int offset = fx - y;
				//cout << offset << endl;
				if (offset > 0) {
					for (int p = 0; p < offset + 5; p++) {
						//Vec3b intensity = (img).at<Vec3b>(y-p,x);

						img.at<Vec3b>(fx - p, x) = average(x, y - p, img, 1);
						//cout << img.at<Vec3b>(fx - p, x) << endl;
						//cout << red << " " << green << " " << blue << endl;

					}

				}
				else if (offset < 0) {
					for (int p = 0; p < -offset + 5; p++) {
						img.at<Vec3b>(fx + p, x) = average(x, y + p, img, 1);

						//cout << red<<" "<<green<<" "<<blue << endl;
					}
				}
				/*Rect roi(x - 1, fx - 1,3, 3);
				Mat roiImg = img(roi);
				Mat blurredRoi;
				GaussianBlur(roiImg, blurredRoi,Size(3, 3), 1);
				blurredRoi.copyTo(img(roi));*/
			}
		}
	}
	else if (1 / abs(slope) >= 4.5) {
		for (int xp = x - 10; xp < x + 10; xp++) {
			if (canny.at<uchar>(fx, xp) == 255) {
				int offset = x - xp;
				//cout << offset << endl;
				if (offset > 0) {
					for (int p = 0; p < offset + 10; p++) {
						//Vec3b intensity = (img).at<Vec3b>(y-p,x);

						img.at<Vec3b>(fx - p, x) = average(xp - p, fx, img, 2);
						//cout << img.at<Vec3b>(fx - p, x) << endl;
						//cout << red << " " << green << " " << blue << endl;

					}

				}
				else if (offset < 0) {
					for (int p = 0; p < offset + 10; p++) {
						img.at<Vec3b>(fx + p, x) = average(xp + p, fx, img, 2);

						//cout << red<<" "<<green<<" "<<blue << endl;
					}
				}
				/*Rect roi(x - 1, fx - 1, 5, 5);
				Mat roiImg = img(roi);
				Mat blurredRoi;
				GaussianBlur(roiImg, blurredRoi, Size(3, 3), 1);
				blurredRoi.copyTo(img(roi));*/

			}
		}
	}
	else if (abs(slope) > 1) {
		for (int i = -5; i < 5; i++) {
			if (canny.at<uchar>(fx + slope * i, x + i) == 255) {
				int offset = i;
				//cout << offset << endl;
				if (offset > 0) {
					for (int p = 0; p < offset + 1; p++) {
						//Vec3b intensity = (img).at<Vec3b>(y-p,x);

						img.at<Vec3b>(fx - p, x) = average(i - p, fx, img, 1);
						//cout << img.at<Vec3b>(fx - p, x) << endl;
						//cout << red << " " << green << " " << blue << endl;

					}
					//img.at<Vec3b>(x, fx) = Vec3b(0, 0, 0);
					//cout << x << " " << fx << endl;
				}
				else if (offset < 0) {
					for (int p = 0; p < offset + 1; p++) {
						img.at<Vec3b>(fx + p, x) = average(i + p, fx, img, 1);

						//cout << red<<" "<<green<<" "<<blue << endl;
					}
				}

			}
		}
	}
	else {

	}
	return img;
}






MouseParams* prep(void* param) {
	MouseParams* image = reinterpret_cast<MouseParams*>(param);

	Mat img = image->currentimg;
	imwrite("currentimg.jpg", img);
	Mat mask = image->mask;
	int height = image->height;
	int width = image->width;
	
	Mat img_canny, greyMat;
	cvtColor(img, greyMat, COLOR_BGR2GRAY);

	Mat grad_x, grad_y;
	Mat abs_grad_x, abs_grad_y;
	Sobel(greyMat, grad_x, CV_16S, 1, 0, 3, 1, 0, BORDER_DEFAULT);
	Sobel(greyMat, grad_y, CV_16S, 0, 1, 3, 1, 0, BORDER_DEFAULT);
	convertScaleAbs(grad_x, abs_grad_x);
	convertScaleAbs(grad_y, abs_grad_y);
	Mat grad;
	addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad);

	Mat img_G;
	GaussianBlur(greyMat, img_G, cv::Size(5, 5), 0);

	Canny(img_G, img_canny, 8, 15);

	Mat canny_result;
	bitwise_and(mask, img_canny, canny_result);

	//imshow("img_canny", canny_result);
	imwrite("canny.jpg", img_canny);
	imwrite("canny_result.jpg", canny_result);
	imwrite("mask.jpg", mask);


	Mat binary;
	threshold(canny_result, binary, 127, 255, THRESH_BINARY);

	// Find connected components
	Mat labels, stats, centroids;
	int num_labels = connectedComponentsWithStats(binary, labels, stats, centroids, 8, CV_32S);

	int largest_component = 1;
	int largest_area = 0;
	for (int i = 1; i < num_labels; ++i) {
		int area = stats.at<int>(i, CC_STAT_AREA);
		if (area > largest_area) {
			largest_area = area;
			largest_component = i;
		}
	}

	// Create a mask for the largest connected component
	Mat canny_clean = Mat::zeros(binary.size(), CV_8U);
	canny_clean.setTo(255, labels == largest_component);
	
	vector<Point> whitePixels; 
	imwrite("aa.jpg", canny_clean);
	int length = 0;
	// Traverse the image
	Mat book;

	
	for (int i = image->startloc.y - 5; i < image->endloc.y + 5; i++) {
		// Get pointer to the i-th row
		uchar* row = canny_result.ptr<uchar>(i);

		for (int j = image->startloc.x - 5; j < image->endloc.x + 5; j++) {
			// Check if the pixel is white
			if (row[j] == 255) {

				// Store the coordinates
				length++;
				whitePixels.push_back(Point(j, i));
				//cout << j << " ";
			}
		}
	}
	

	sort(whitePixels.begin(), whitePixels.end(), comparePoints);





	int l = 1;
	

	output output = Poly(whitePixels, length, width, height);



	Mat result_canny, result_lap;
	vector<int> xcor = output.xout;
	vector<int> ycor = output.yout;
	result_canny = output.result;
	imwrite("poly_result.jpg", result_canny);
	bitwise_or(canny_result, result_canny, result_lap);


	
	if (output.ind == 1) {
		transpose(img, img);
	}
	
	for (int i = 1; i < length - 1; i++) {
		int x = whitePixels[i].x;

		int y = whitePixels[i].y;

		if (y == whitePixels[i + 1].y) {
			l++;
		}
		else if (l > 2) {
			//cout << l << endl;
			for (int j = -2; j <= 2; j++) {
				for (int k = l - 2; k >= 1; k--) {
					Vec3b valueLeft = img.at<Vec3b>(y + j, x - k - 1);
					Vec3b valueCenter = img.at<Vec3b>(y + j, x - k);
					Vec3b valueRight = img.at<Vec3b>(y + j, x - k + 1);

					double avgB = (valueLeft[0] + valueCenter[0] + valueRight[0]) / 3.0;
					double avgG = (valueLeft[1] + valueCenter[1] + valueRight[1]) / 3.0;
					double avgR = (valueLeft[2] + valueCenter[2] + valueRight[2]) / 3.0;
					img.at<Vec3b>(y + j, x - k)[0] = avgB;
					img.at<Vec3b>(y + j, x - k)[1] = avgG;
					img.at<Vec3b>(y + j, x - k)[2] = avgR;

				}
			}
			l = 1;
		}
	}




	int a = 0;
	for (auto elem : xcor) {
		//cout << elem << endl;
		double fxp = 0;
		double nor;
		for (int i = 0; i <= output.ord; i++) {

			fxp += (i + 1) * pow(elem, i) * output.coef[i + 1];
		}
		//cout << "x=" << fx << endl;
		
		
		//if (canny_result.at<uchar>(Point(elem, fx)) == result_canny.at<uchar>(Point(elem, fx))) continue;

		//cout <<"y=" << fxp << endl;
		nor = -1 / fxp;

		//double offset = fx - nor * elem;

		img = manipulation(nor, img, canny_result, result_canny, elem, ycor[a]);
		a++;
	}
	int order = 0;
	for (auto elem : xcor) {
		for (int i = -1; i <= 1; i++) {
			//cout << elem << " " << ycoor[order] << endl;
			img.at<Vec3b>(ycor[order]+i, elem) = average(elem, ycor[order] + i, img, 1);
		}
		order++;
	}

	order = 0;

	for (auto elem : xcor) {
		img.at<Vec3b>(ycor[order], elem) = average(elem, ycor[order], img, 2);
		order++;
	}
	if (output.ind == 1) {
		transpose(img, img);
	}

	imwrite("rrsult.jpg", img);
	cout << "success" << endl;
	image->currentimg = img;

	grad.release();
	grad_x.release();
	grad_y.release();
	abs_grad_x.release();
	abs_grad_y.release();
	img.release();
	canny_clean.release();
	canny_result.release();
	result_canny.release(); 
	result_lap.release();
	binary.release();
	labels.release();
	stats.release();
	centroids.release();
	return image;
}