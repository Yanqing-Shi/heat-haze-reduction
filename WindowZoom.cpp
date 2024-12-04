#include<algorithm>
#include "WindowZoom.h"
#include "Polyfit.h"
#include"method.h"
using namespace std;
using namespace cv;
bool drawing = false;



void MouseCall (int event, int x, int y, int flag, void* param)
{
	COpenCVWindowExt* pParent= (COpenCVWindowExt*)(param);
	MouseParams* mouse = pParent->Mouse;



	Mat mask = (*mouse).mask;
	//cout << x << " " << y << endl;
	
	if (event == EVENT_MOUSEWHEEL)
	{
		//pParent->SetImage();
		//cout << pParent->m_iScaleTimes << endl;
		if (getMouseWheelDelta (flag) > 0 && pParent->m_iScaleTimes != pParent->m_iMaxScaleTimes)
			pParent->m_iScaleTimes++;
		else if (getMouseWheelDelta (flag) < 0 && pParent->m_iScaleTimes > pParent->m_iMinScaleTimes)
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
	else if (event == EVENT_RBUTTONDOWN && pParent->m_iScaleTimes!=0)
	{
		pParent->ptRButtonDown.x = x;
		pParent->ptRButtonDown.y = y;
	}
	else if (flag == EVENT_FLAG_RBUTTON && pParent->m_iScaleTimes != 0)
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
			int xcurrent = x / (pParent->m_dNewScale * 1280) * 6016;
			int ycurrent = y / (pParent->m_dNewScale * 854) * 4016;
			circle(mask, Point(xcurrent, ycurrent),4, 255, FILLED);
			if (xcurrent < mouse->startloc.x) mouse->startloc.x = xcurrent;

			if (ycurrent < mouse->startloc.y) mouse->startloc.y = ycurrent;
			if (xcurrent > mouse->endloc.x) mouse->endloc.x = xcurrent;
			if (ycurrent > mouse->endloc.y) mouse->endloc.y = ycurrent;
			
		}
		//pParent->OverlapImage(x, y);
	}else if (event == EVENT_LBUTTONDOWN && drawing == false) {
		drawing = true;
		mouse->height = mouse->currentimg.rows;
		mouse->width = mouse->currentimg.cols;
		Mat mask = Mat::zeros(Size(mouse->width, mouse->height), CV_8UC1);
		mouse->mask = mask; 

		mouse->startloc = Point(x / (pParent->m_dNewScale * 1280) * 6016, y / (pParent->m_dNewScale * 854) * 4016);
		mouse->endloc = Point(x / (pParent->m_dNewScale * 1280) * 6016, y / (pParent->m_dNewScale * 854) * 4016);
		//cout << x << " " << y << endl;
		mask.release();
	}
	
	else if (event == EVENT_LBUTTONUP) {
		drawing = false;
		
		//imwrite("starry_night.jpg", mask);
		pParent->Mouse=prep(mouse);
		Rect org(Point(pParent->m_iHorzScrollBarPos / (pParent->m_dNewScale * 1280) * 6016, 
			pParent->m_iVertScrollBarPos / (pParent->m_dNewScale * 854) * 4016), Size(6016 / pParent->m_dNewScale, 4016 / pParent->m_dNewScale));

		pParent->Mouse->currentimg.copyTo(pParent->m_matSrc(org));

		pParent->history.push(pParent->m_matSrc.clone());
		
		pParent->RefreshImage();

	
		/*for (auto elem : points) {
			cout << elem << endl;
		}*/
		//cout <<endl << points.back() << endl;
	}
	else if (event == EVENT_MBUTTONDOWN) {
		if (pParent->history.size() > 1) {
			pParent->history.pop();                  // Remove the current state
			pParent->m_matSrc = pParent->history.top().clone();  // Revert to previous state
			pParent->RefreshImage();

		}

		else {
			std::cout << "No more undos available!" << std::endl;
		}
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
	history.push(m_matSrc.clone());

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

void COpenCVWindowExt::OverlapImage(int x,int y) 
{

	Mat m (Mouse->height, Mouse->width, CV_8UC4, Scalar(0, 0, 0, 100));
	
	//circle(m, Point(x, y), 4, Scalar(255, 255, 255, 100), -1);
	Mat overlapped;
	Mat mat;
	int iW = m_vecMatResize[0].cols, iH = m_vecMatResize[0].rows;
	Rect rectShow(Point(m_iHorzScrollBarPos, m_iVertScrollBarPos), Size(iW, iH));
	cvtColor(m_vecMatResize[m_iScaleTimes](rectShow), mat, COLOR_BGR2BGRA);
	overlapped=mat+m;
	imshow(m_strWindowName, overlapped);
	mat.release();
	m.release();
	overlapped.release();
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
	
	
	int iW = m_vecMatResize[0].cols, iH = m_vecMatResize[0].rows;
	

	Rect rectShow (Point (m_iHorzScrollBarPos, m_iVertScrollBarPos), Size (iW, iH));
	//cout << rectShow.tl() << " " << rectShow.br();
	imshow (m_strWindowName, m_vecMatResize[m_iScaleTimes] (rectShow));
	Rect org(Point(m_iHorzScrollBarPos / (m_dNewScale * 1280) * 6016, m_iVertScrollBarPos / (m_dNewScale * 854) * 4016),Size(6016/ m_dNewScale,4016/ m_dNewScale));
	Mouse->currentimg = m_matSrc(org);//m_vecMatResize[m_iScaleTimes](rectShow);
	

}



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

