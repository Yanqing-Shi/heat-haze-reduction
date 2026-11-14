#include"method.h"
#include"Polyfit.h"

using namespace std;
using namespace cv;

namespace {
constexpr bool kSaveDebugImages = false;

int ClampInt(int value, int minValue, int maxValue)
{
	return max(minValue, min(value, maxValue));
}

bool InBounds(const Mat& img, int x, int y)
{
	return x >= 0 && y >= 0 && x < img.cols && y < img.rows;
}

void SaveDebugImage(const string& path, const Mat& img)
{
	if (kSaveDebugImages) {
		imwrite(path, img);
	}
}
}


bool comparePoints(const Point& p1, const Point& p2) {
	return p1.x < p2.x;
}

Vec3b average(int x, int y, const Mat& img, int flag) {
	Vec3d sum(0, 0, 0);
	if (flag == 1) {
		for (int dx = -2; dx <= 2; ++dx) {
			const Vec3b value = img.at<Vec3b>(ClampInt(y, 0, img.rows - 1), ClampInt(x + dx, 0, img.cols - 1));
			sum += Vec3d(value[0], value[1], value[2]);
		}
	}
	else if (flag == 2) {
		for (int dy = -2; dy <= 2; ++dy) {
			const Vec3b value = img.at<Vec3b>(ClampInt(y + dy, 0, img.rows - 1), ClampInt(x, 0, img.cols - 1));
			sum += Vec3d(value[0], value[1], value[2]);
		}
	}


	return Vec3b(
		saturate_cast<uchar>(sum[0] / 5.0),
		saturate_cast<uchar>(sum[1] / 5.0),
		saturate_cast<uchar>(sum[2] / 5.0));
}



Mat manipulation(double slope, Mat img, Mat canny, Mat post_canny, int x, int fx) {
	//cout << slope << endl;
	if (abs(slope) >= 4.5) {
		//cout << x << " "<<fx << endl;
		for (int y = fx - 10; y < fx + 10; y++) {
			if (!InBounds(canny, x, y)) {
				continue;
			}
			if (canny.at<uchar>(y, x) == 255) {
				int offset = fx - y;
				//cout << offset << endl;
				if (offset > 0) {
					for (int p = 0; p < offset + 5; p++) {
						//Vec3b intensity = (img).at<Vec3b>(y-p,x);

						if (InBounds(img, x, fx - p)) {
							img.at<Vec3b>(fx - p, x) = average(x, y - p, img, 1);
						}
						//cout << img.at<Vec3b>(fx - p, x) << endl;
						//cout << red << " " << green << " " << blue << endl;

					}

				}
				else if (offset < 0) {
					for (int p = 0; p < -offset + 5; p++) {
						if (InBounds(img, x, fx + p)) {
							img.at<Vec3b>(fx + p, x) = average(x, y + p, img, 1);
						}

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
			if (!InBounds(canny, xp, fx)) {
				continue;
			}
			if (canny.at<uchar>(fx, xp) == 255) {
				int offset = x - xp;
				//cout << offset << endl;
				if (offset > 0) {
					for (int p = 0; p < offset + 10; p++) {
						//Vec3b intensity = (img).at<Vec3b>(y-p,x);

						if (InBounds(img, x, fx - p)) {
							img.at<Vec3b>(fx - p, x) = average(xp - p, fx, img, 2);
						}
						//cout << img.at<Vec3b>(fx - p, x) << endl;
						//cout << red << " " << green << " " << blue << endl;

					}

				}
				else if (offset < 0) {
					for (int p = 0; p < offset + 10; p++) {
						if (InBounds(img, x, fx + p)) {
							img.at<Vec3b>(fx + p, x) = average(xp + p, fx, img, 2);
						}

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
			const int sampleY = static_cast<int>(fx + slope * i);
			const int sampleX = x + i;
			if (!InBounds(canny, sampleX, sampleY)) {
				continue;
			}
			if (canny.at<uchar>(sampleY, sampleX) == 255) {
				int offset = i;
				//cout << offset << endl;
				if (offset > 0) {
					for (int p = 0; p < offset + 1; p++) {
						//Vec3b intensity = (img).at<Vec3b>(y-p,x);

						if (InBounds(img, x, fx - p)) {
							img.at<Vec3b>(fx - p, x) = average(i - p, fx, img, 1);
						}
						//cout << img.at<Vec3b>(fx - p, x) << endl;
						//cout << red << " " << green << " " << blue << endl;

					}
					//img.at<Vec3b>(x, fx) = Vec3b(0, 0, 0);
					//cout << x << " " << fx << endl;
				}
				else if (offset < 0) {
					for (int p = 0; p < offset + 1; p++) {
						if (InBounds(img, x, fx + p)) {
							img.at<Vec3b>(fx + p, x) = average(i + p, fx, img, 1);
						}

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
	image->changed = false;

	Mat img = image->currentimg;
	SaveDebugImage("currentimg.jpg", img);
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
	SaveDebugImage("canny.jpg", img_canny);
	SaveDebugImage("canny_result.jpg", canny_result);
	SaveDebugImage("mask.jpg", mask);


	Mat binary;
	threshold(canny_result, binary, 127, 255, THRESH_BINARY);

	
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
	SaveDebugImage("aa.jpg", canny_clean);
	int length = 0;
	// Traverse the image
	Mat book;

	
	const int minY = ClampInt(image->startloc.y - 5, 0, canny_result.rows - 1);
	const int maxY = ClampInt(image->endloc.y + 5, minY, canny_result.rows - 1);
	const int minX = ClampInt(image->startloc.x - 5, 0, canny_result.cols - 1);
	const int maxX = ClampInt(image->endloc.x + 5, minX, canny_result.cols - 1);

	for (int i = minY; i <= maxY; i++) {
		// Get pointer to the i-th row
		uchar* row = canny_result.ptr<uchar>(i);

		for (int j = minX; j <= maxX; j++) {
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

	if (whitePixels.empty()) {
		cout << "No edge pixels found in the selected area." << endl;
		image->currentimg = img;
		return image;
	}





	int l = 1;
	

	output output = Poly(whitePixels, length, width, height);



	Mat result_canny, result_lap;
	vector<int> xcor = output.xout;
	vector<int> ycor = output.yout;
	result_canny = output.result;
	SaveDebugImage("poly_result.jpg", result_canny);
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
					if (!InBounds(img, x - k - 1, y + j) || !InBounds(img, x - k, y + j) || !InBounds(img, x - k + 1, y + j)) {
						continue;
					}
					Vec3b valueLeft = img.at<Vec3b>(y + j, x - k - 1);
					Vec3b valueCenter = img.at<Vec3b>(y + j, x - k);
					Vec3b valueRight = img.at<Vec3b>(y + j, x - k + 1);

					double avgB = (valueLeft[0] + valueCenter[0] + valueRight[0]) / 3.0;
					double avgG = (valueLeft[1] + valueCenter[1] + valueRight[1]) / 3.0;
					double avgR = (valueLeft[2] + valueCenter[2] + valueRight[2]) / 3.0;
					img.at<Vec3b>(y + j, x - k)[0] = saturate_cast<uchar>(avgB);
					img.at<Vec3b>(y + j, x - k)[1] = saturate_cast<uchar>(avgG);
					img.at<Vec3b>(y + j, x - k)[2] = saturate_cast<uchar>(avgR);

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
			if (InBounds(img, elem, ycor[order] + i)) {
				img.at<Vec3b>(ycor[order]+i, elem) = average(elem, ycor[order] + i, img, 1);
			}
		}
		order++;
	}

	order = 0;

	for (auto elem : xcor) {
		if (InBounds(img, elem, ycor[order])) {
			img.at<Vec3b>(ycor[order], elem) = average(elem, ycor[order], img, 2);
		}
		order++;
	}
	if (output.ind == 1) {
		transpose(img, img);
	}

	imwrite("result.jpg", img);
	cout << "success" << endl;
	image->currentimg = img;
	image->changed = true;

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
