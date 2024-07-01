#include<algorithm>
#include "OpenCVWindowExt.h"
#include "Polyfit.h"
using namespace std;
using namespace cv;
bool drawing = false;


bool comparePoints(const Point& p1, const Point& p2) {
	return p1.x < p2.x; 
}

Vec3b average(int x,int y,Mat img,int flag) {
	Vec3b output;
	if (flag == 1) {
		Vec3b valueLeft = img.at<Vec3b>(y, x - 2);
		Vec3b valueMidLeft = img.at<Vec3b>(y, x - 1);
		Vec3b valueCenter = img.at<Vec3b>(y, x);
		Vec3b valueMidRight = img.at<Vec3b>(y, x + 1);
		Vec3b valueRight = img.at<Vec3b>(y, x + 2);
		//Mat kernel = Mat::ones(3, 3, CV_64F) / 9;
		double avgB = (valueLeft[0] + valueMidLeft[0] + valueCenter[0] + valueMidRight[0] + valueRight[0]) / 5.0;
		double avgG = (valueLeft[1] + valueMidLeft[1] + valueCenter[1] + valueMidRight[1] + valueRight[1]) / 5.0;
		double avgR = (valueLeft[2] + valueMidLeft[2] + valueCenter[2] + valueMidRight[2] + valueRight[2]) / 5.0;
		output=Vec3b(avgB, avgG, avgR);
	}
	else if (flag == 2) {
		Vec3b valueLeft = img.at<Vec3b>(y-2, x);
		Vec3b valueMidLeft = img.at<Vec3b>(y-1, x);
		Vec3b valueCenter = img.at<Vec3b>(y, x);
		Vec3b valueMidRight = img.at<Vec3b>(y+1, x);
		Vec3b valueRight = img.at<Vec3b>(y+2, x);
		//Mat kernel = Mat::ones(3, 3, CV_64F) / 9;
		double avgB = (valueLeft[0] + valueMidLeft[0] + valueCenter[0] + valueMidRight[0] + valueRight[0]) / 5.0;
		double avgG = (valueLeft[1] + valueMidLeft[1] + valueCenter[1] + valueMidRight[1] + valueRight[1]) / 5.0;
		double avgR = (valueLeft[2] + valueMidLeft[2] + valueCenter[2] + valueMidRight[2] + valueRight[2]) / 5.0;
		output=Vec3b(avgB, avgG, avgR);  
	}

	
	return output;
}



Mat manipulation(double slope,Mat img, Mat canny, Mat post_canny,int x,int fx) {
	//cout << slope << endl;
	if (abs(slope) >= 4.5) {
	//cout << x << " "<<fx << endl;
		for (int y = fx - 10; y < fx+10; y++) {
			if (canny.at<uchar>(y,x) == 255) {
				int offset = fx - y;
				//cout << offset << endl;
				if (offset > 0) {
					for (int p = 0; p < offset+10; p++) { 
						//Vec3b intensity = (img).at<Vec3b>(y-p,x);
						
						img.at<Vec3b>(fx - p, x) = average(x, y - p,img,1);
						//cout << img.at<Vec3b>(fx - p, x) << endl;
						//cout << red << " " << green << " " << blue << endl;
						
					}

				}
				else if(offset<0){
					for (int p = 0; p < -offset+10; p++) {
						img.at<Vec3b>(fx + p, x) = average(x, y + p, img,1);

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
	}else if (1/abs(slope) >=4.5) {
		for (int xp = x - 10; xp < x + 10; xp++) {
			if (canny.at<uchar>(fx, xp) == 255) {
				int offset = x - xp;
				//cout << offset << endl;
				if (offset > 0) {
					for (int p = 0; p < offset + 10; p++) {
						//Vec3b intensity = (img).at<Vec3b>(y-p,x);

						img.at<Vec3b>(fx - p, x) = average(xp-p, fx, img,2);
						//cout << img.at<Vec3b>(fx - p, x) << endl;
						//cout << red << " " << green << " " << blue << endl;

					}

				}
				else if (offset < 0) {
					for (int p = 0; p < offset + 10; p++) {
						img.at<Vec3b>(fx + p, x) = average(xp+p, fx, img,2);

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
			if (canny.at<uchar>(fx+slope*i, x+i) == 255) {
				int offset = i;
				//cout << offset << endl;
				if (offset > 0) {
					for (int p = 0; p < offset + 1; p++) {
						//Vec3b intensity = (img).at<Vec3b>(y-p,x);

						img.at<Vec3b>(fx - p, x) = average(i - p, fx, img, 3);
						//cout << img.at<Vec3b>(fx - p, x) << endl;
						//cout << red << " " << green << " " << blue << endl;

					}
					//img.at<Vec3b>(x, fx) = Vec3b(0, 0, 0);
					//cout << x << " " << fx << endl;
				}
				else if (offset < 0) {
					for (int p = 0; p < offset + 1; p++) {
						img.at<Vec3b>(fx + p, x) = average(i + p, fx, img, 3);

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
	Mat img_canny,greyMat;
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
	
	//cout << height << " " << width << endl;

	vector<Point> whitePixels;
	int length = 0;
	// Traverse the image
	Mat book;
	for (int i = 0; i < height; ++i) {
		// Get pointer to the i-th row
		uchar* row = canny_result.ptr<uchar>(i);

		for (int j = 0; j < width; ++j) {
			// Check if the pixel is white
			if (row[j] == 255) {

				// Store the coordinates
				length++;
				whitePixels.push_back(Point(j, i));
				//cout << j << " ";
			}
		}
	}
	//cout << endl;
	/*for (auto elem : whitePixels) {
		cout << elem << endl;
	}*/
	sort(whitePixels.begin(), whitePixels.end(), comparePoints);
	
	//imshow("a", result_lap);
	//if (xcoor.empty()) cout << "empty";
	//cout<<xcoor.size();


	
	int l = 1;
	for (int i = 1; i < length - 1; i++) {
		int x = whitePixels[i].x;
		
		int y = whitePixels[i].y;
		
		if (y == whitePixels[i + 1].y ) {
			l++;
		}
		else if(l>2){
			//cout << l << endl;
			for (int j = -2; j <= 0; j++) {
				for (int k = l-2; k >=1; k--) {
					Vec3b valueLeft = img.at<Vec3b>(y+j, x-k-1);
					Vec3b valueCenter = img.at<Vec3b>(y+j, x-k);
					Vec3b valueRight = img.at<Vec3b>(y+j, x-k+1);

					double avgB = (valueLeft[0] + valueCenter[0] +  valueRight[0]) / 3.0;
					double avgG = (valueLeft[1]  + valueCenter[1] +  valueRight[1]) / 3.0;
					double avgR = (valueLeft[2] + valueCenter[2] + valueRight[2]) / 3.0;
					img.at<Vec3b>(y + j, x - k)[0] = avgB;
					img.at<Vec3b>(y + j, x - k)[1] = avgG;
					img.at<Vec3b>(y + j, x - k)[2] = avgR;
					
				}
			}
			l = 1;
		}
	}

	output output = Poly(whitePixels,length,width,height);



	Mat result_canny, result_lap;
	vector<int> xcoor = output.xout;
	result_canny = output.result;
	imwrite("poly_result.jpg", result_canny);
	bitwise_or(canny_result, result_canny, result_lap);

	vector<int> ycoor;
	
	for (auto elem : xcoor) {
		//cout << elem << endl;
		double fx=0;
		double fxp = 0;
		double nor;
		for (int i = 0; i <= output.ord; i++) {
			
			fx += int(pow(elem, i) * output.coef[i]);
			fxp+= (i+1)*pow(elem, i) * output.coef[i + 1];
		}
		ycoor.push_back(fx);
		//if (canny_result.at<uchar>(Point(elem, fx)) == result_canny.at<uchar>(Point(elem, fx))) continue;

		//cout << fxp << endl;
		nor = -1 / fxp;
		
		//double offset = fx - nor * elem;
		
		img=manipulation(nor, img, canny_result, result_canny,elem,fx); 

	}
	int order = 0;
	for (auto elem : xcoor) {
		for (int i = -2; i <= 2; i++) {
			cout << elem << " " << ycoor[order] << endl;
			img.at<Vec3b>(ycoor[order], elem) = average(elem, ycoor[order] + i, img, 1);
		}
		order++;
	}

	order = 0;

	for (auto elem : xcoor) {
		img.at<Vec3b>(ycoor[order], elem) = average(elem, ycoor[order], img, 2); 
	}

	imwrite("rrsult.jpg", img);
	cout << "success" << endl;
	image->currentimg = img;
	return image;
}













void MouseCall (int event, int x, int y, int flag, void* param)
{
	COpenCVWindowExt* pParent= (COpenCVWindowExt*)(param);

	MouseParams* mouse = pParent->Mouse;
	Mat mask = (*mouse).mask;
	
	
	if (event == EVENT_MOUSEWHEEL)
	{
		//pParent->SetImage();
		//cout << pParent->m_iScaleTimes << endl;
		if (getMouseWheelDelta (flag) > 0 && pParent->m_iScaleTimes != pParent->m_iMaxScaleTimes)
			pParent->m_iScaleTimes++;
		else if (getMouseWheelDelta (flag) < 0 && pParent->m_iScaleTimes != pParent->m_iMinScaleTimes)
			pParent->m_iScaleTimes--;

		if (pParent->m_iScaleTimes == 0)
			pParent->m_dCompensationX = pParent->m_dCompensationY = 0;

		//Pixel value = mouse offset (v.s. window's left top position)
		//But using x, y is not correct in Wheel Event. So, use pre recorded value
		x = pParent->ptLButtonDown.x;
		y = pParent->ptLButtonDown.y;
		double dPixelX = (pParent->m_iHorzScrollBarPos + x + pParent->m_dCompensationX) / pParent->m_dNewScale;
		double dPixelY = (pParent->m_iVertScrollBarPos + y + pParent->m_dCompensationY) / pParent->m_dNewScale;

		pParent->m_dNewScale = pParent->m_dInitialScale * pow (pParent->m_dScaleRatio, pParent->m_iScaleTimes);

        int iW = pParent->m_iOrgW;
        int iH = pParent->m_iOrgH;
        pParent->m_iHorzScrollBarRange_Max = int(pParent->m_dNewScale * iW - pParent->m_dInitialScale * iW);
        pParent->m_iVertScrollBarRange_Max = int(pParent->m_dNewScale * iH - pParent->m_dInitialScale * iH);
        int iBarPosX = int(dPixelX * pParent->m_dNewScale - x + 0.5);
        int iBarPosY = int(dPixelY * pParent->m_dNewScale - y + 0.5);
        pParent->SetHorzBarPos(iBarPosX);
        pParent->SetVertBarPos(iBarPosY);
        pParent->m_dCompensationX = -iBarPosX + (dPixelX * pParent->m_dNewScale - x);
        pParent->m_dCompensationY = -iBarPosY + (dPixelY * pParent->m_dNewScale - y);

		pParent->RefreshImage ();
	}
	else if (event == EVENT_RBUTTONDOWN)
	{
		pParent->ptRButtonDown.x = x;
		pParent->ptRButtonDown.y = y;
	}
	else if (flag == EVENT_FLAG_RBUTTON)
	{
		int iRButtonOffsetX = x - pParent->ptRButtonDown.x;
		int iRButtonOffsetY = y - pParent->ptRButtonDown.y;


		int iBarPosX = pParent->m_iHorzScrollBarPos_copy - iRButtonOffsetX;
		pParent->SetHorzBarPos (iBarPosX);
		

		int iBarPosY = pParent->m_iVertScrollBarPos_copy - iRButtonOffsetY;
		pParent->SetVertBarPos (iBarPosY);

		
		
		pParent->RefreshImage ();
	}
	else if (event == EVENT_MOUSEMOVE)
	{
		pParent->ptLButtonDown.x = x;
		pParent->ptLButtonDown.y = y;
		pParent->m_iHorzScrollBarPos_copy = pParent->m_iHorzScrollBarPos;
		pParent->m_iVertScrollBarPos_copy = pParent->m_iVertScrollBarPos;
		if (flag == EVENT_FLAG_LBUTTON) {
			//cout << x << " " << y << endl;
			circle(mask, Point(x / (pParent->m_dNewScale * 1280) * 6016, y / (pParent->m_dNewScale * 854) * 4016),5, 255, FILLED);
		}

	}else if (event == EVENT_LBUTTONDOWN && drawing == false) {
		drawing = true;
		mouse->height = mouse->currentimg.rows;
		mouse->width = mouse->currentimg.cols;
		Mat mask = Mat::zeros(Size(mouse->width, mouse->height), CV_8UC1);
		mouse->mask = mask; 
		//cout << x << " " << y << endl;
	}
	
	else if (event == EVENT_LBUTTONUP) {
		drawing = false;
		//imwrite("starry_night.jpg", mask);
		pParent->Mouse=prep(mouse);
		Rect org(Point(pParent->m_iHorzScrollBarPos / (pParent->m_dNewScale * 1280) * 6016, 
			pParent->m_iVertScrollBarPos / (pParent->m_dNewScale * 854) * 4016), Size(6016 / pParent->m_dNewScale, 4016 / pParent->m_dNewScale));

		pParent->Mouse->currentimg.copyTo(pParent->m_matSrc(org));


		
		pParent->RefreshImage();

	
		/*for (auto elem : points) {
			cout << elem << endl;
		}*/
		//cout <<endl << points.back() << endl;
	}
}
COpenCVWindowExt::COpenCVWindowExt (String strWindowName, int iFlag)
{
	namedWindow (strWindowName, iFlag);
	Mouse = new MouseParams;
	m_strWindowName = strWindowName;
	//MouseParams image;
	
	
	setMouseCallback(strWindowName, MouseCall, this);
	

	

	//Initial values

	m_iScaleTimes = 0;
	m_iMaxScaleTimes = 12;
	m_iMinScaleTimes = 0;
	m_vecMatResize.resize (m_iMaxScaleTimes + 1);
	m_dCompensationX = 0;
	m_dCompensationY = 0;
	m_dInitialScale = 2;
	m_dNewScale = 1;
	m_dScaleRatio = 1.25;

	m_iOrgW = 0;
	m_iOrgH = 0;

	m_iHorzScrollBarPos = 0;
	m_iVertScrollBarPos = 0;

	m_iHorzScrollBarRange_Min = 0;
	m_iHorzScrollBarRange_Max = 1;
	m_iVertScrollBarRange_Min = 0;
	m_iVertScrollBarRange_Max = 1;
	

}

COpenCVWindowExt::~COpenCVWindowExt ()
{
	setMouseCallback (m_strWindowName, NULL);
}

bool COpenCVWindowExt::ImRead (String strFileName)
{
	m_matSrc = imread (strFileName);
	if (m_matSrc.empty ())
		return false;
	m_iOrgW = m_matSrc.cols;
	m_iOrgH = m_matSrc.rows;
	//Size sizeInitial (m_matSrc.cols * m_dInitialScale, m_matSrc.rows * m_dInitialScale);
	Size sizeInitial (1280, 854);
	//cout << m_matSrc.cols * m_dInitialScale << " " << m_matSrc.rows * m_dInitialScale << endl;
	resize (m_matSrc, m_vecMatResize[0], sizeInitial);
	//imwrite("aa.jpg", img_canny);
	if (!m_vecMatResize[0].empty()) {
		imshow(m_strWindowName, m_vecMatResize[0]);
		Mouse->img = m_vecMatResize[0];
		Mouse->width = m_vecMatResize[0].cols; //1280
		Mouse->height = m_vecMatResize[0].rows; //854
		//cout << m_vecMatResize[0].cols << " " << m_vecMatResize[0].rows;
		Mouse->currentimg = m_vecMatResize[0];
	}
		

	return !m_vecMatResize[0].empty ();
}

void COpenCVWindowExt::SetInitailScale (double dScale)
{
	if (dScale <= 0)
		return;

	m_dInitialScale = dScale;
	m_dNewScale = dScale;
}

void COpenCVWindowExt::RefreshImage ()
{
	
	if (m_matSrc.empty ())
		return;
	
	Size size (int (m_dNewScale * 1280), int (m_dNewScale * 854));
		//cout << m_dNewScale << endl;
		//cout << int(m_dNewScale * m_matSrc.cols) << " " << int(m_dNewScale * m_matSrc.rows) << endl;
	resize (m_matSrc, m_vecMatResize[m_iScaleTimes], size);
		//imshow(m_strWindowName, m_vecMatResize[m_iScaleTimes]);
	
	
	int iW = m_vecMatResize[0].cols, iH = m_vecMatResize[0].rows;//-1: for bar size
	

	Rect rectShow (Point (m_iHorzScrollBarPos, m_iVertScrollBarPos), Size (iW, iH));
	//cout << rectShow.tl() << " " << rectShow.br();
	imshow (m_strWindowName, m_vecMatResize[m_iScaleTimes] (rectShow));
	Rect org(Point(m_iHorzScrollBarPos / (m_dNewScale * 1280) * 6016, m_iVertScrollBarPos / (m_dNewScale * 854) * 4016),Size(6016/ m_dNewScale,4016/ m_dNewScale));
	Mouse->currentimg = m_matSrc(org);//m_vecMatResize[m_iScaleTimes](rectShow);
	
	imwrite("aa.jpg", Mouse->currentimg);
	
}

//void COpenCVWindowExt::SetImage() {
//	int iW = m_vecMatResize[0].cols, iH = m_vecMatResize[0].rows;//-1: for bar size
//	Rect rectShow(Point(m_iHorzScrollBarPos, m_iVertScrollBarPos), Size(iW, iH));
//	Mouse->currentimg.copyTo(m_vecMatResize[m_iScaleTimes](rectShow));
//}



void COpenCVWindowExt::SetHorzBarPos (int iPos)
{
	if (iPos > m_iHorzScrollBarRange_Max)
		m_iHorzScrollBarPos = m_iHorzScrollBarRange_Max;
	else if (iPos < 0)
		m_iHorzScrollBarPos = m_iHorzScrollBarRange_Min;
	else
		m_iHorzScrollBarPos = iPos;

}

void COpenCVWindowExt::SetVertBarPos (int iPos)
{
	if (iPos > m_iVertScrollBarRange_Max)
		m_iVertScrollBarPos = m_iVertScrollBarRange_Max;
	else if (iPos < 0)
		m_iVertScrollBarPos = m_iVertScrollBarRange_Min;
	else
		m_iVertScrollBarPos = iPos;
}
