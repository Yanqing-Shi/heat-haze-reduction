#include<algorithm>
#include <cmath>
#include "WindowZoom.h"
//#include "Polyfit.h"
#include"method.h"
using namespace std;
using namespace cv;
bool drawing = false;

namespace {
constexpr int kMaxViewWidth = 1280;
constexpr int kMaxViewHeight = 854;
constexpr size_t kMaxUndoSteps = 10;

int ClampInt(int value, int minValue, int maxValue)
{
	return max(minValue, min(value, maxValue));
}
}



static void HandleMouse (int event, int x, int y, int flag, void* param)
{
	COpenCVWindowExt* pParent= (COpenCVWindowExt*)(param);
	MouseParams* mouse = pParent->Mouse;


	Mat mask = (*mouse).mask;
	//cout << x << " " << y << endl;
	
	if (event == EVENT_MOUSEWHEEL)
	{
		//pParent->SetImage();
		//cout << pParent->m_iScaleTimes << endl;
		if (getMouseWheelDelta(flag) > 0 && pParent->m_iScaleTimes != pParent->m_iMaxScaleTimes)
			pParent->m_iScaleTimes++;
		else if (getMouseWheelDelta(flag) < 0 && pParent->m_iScaleTimes != pParent->m_iMinScaleTimes)
			pParent->m_iScaleTimes--;

		if (pParent->m_iScaleTimes == 0)
			pParent->m_dCompensationX = pParent->m_dCompensationY = 0;

		//Pixel value = mouse offset (v.s. window's left top position)
		//But using x, y is not correct in Wheel Event. So, use pre recorded value
		x = pParent->ptLButtonDown.x;
		y = pParent->ptLButtonDown.y;
		double dPixelX = (pParent->m_iHorzScrollBarPos + x + pParent->m_dCompensationX) / pParent->m_dNewScale;
		double dPixelY = (pParent->m_iVertScrollBarPos + y + pParent->m_dCompensationY) / pParent->m_dNewScale;

		pParent->m_dNewScale = pParent->m_dInitialScale * pow(pParent->m_dScaleRatio, pParent->m_iScaleTimes);

		int iW = pParent->m_viewSize.width;
		int iH = pParent->m_viewSize.height;
		pParent->m_iHorzScrollBarRange_Max = int(pParent->m_dNewScale * iW - pParent->m_dInitialScale * iW);
		pParent->m_iVertScrollBarRange_Max = int(pParent->m_dNewScale * iH - pParent->m_dInitialScale * iH);
		int iBarPosX = int(dPixelX * pParent->m_dNewScale - x + 0.5);
		int iBarPosY = int(dPixelY * pParent->m_dNewScale - y + 0.5);
		pParent->SetHorzBarPos(iBarPosX);
		pParent->SetVertBarPos(iBarPosY);
		pParent->m_dCompensationX = -iBarPosX + (dPixelX * pParent->m_dNewScale - x);
		pParent->m_dCompensationY = -iBarPosY + (dPixelY * pParent->m_dNewScale - y);

		pParent->RefreshImage();
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
		int iBarPosY = pParent->m_iVertScrollBarPos_copy - iRButtonOffsetY;
		
		pParent->SetHorzBarPos(iBarPosX);
		pParent->SetVertBarPos(iBarPosY);

		
		
		pParent->RefreshImage ();
	}
	else if (event == EVENT_MOUSEMOVE)
	{
		pParent->ptLButtonDown.x = x;
		pParent->ptLButtonDown.y = y;
		pParent->m_iHorzScrollBarPos_copy = pParent->m_iHorzScrollBarPos;
		pParent->m_iVertScrollBarPos_copy = pParent->m_iVertScrollBarPos;
		
		if (flag == EVENT_FLAG_LBUTTON && drawing && !mask.empty()) {
			Point imagePoint = pParent->WindowToImagePoint(x, y);
			circle(mask, imagePoint,4, 255, FILLED);
			if (imagePoint.x < mouse->startloc.x) mouse->startloc.x = imagePoint.x;

			if (imagePoint.y < mouse->startloc.y) mouse->startloc.y = imagePoint.y;
			if (imagePoint.x > mouse->endloc.x) mouse->endloc.x = imagePoint.x;
			if (imagePoint.y > mouse->endloc.y) mouse->endloc.y = imagePoint.y;
			
		}
	}else if (event == EVENT_LBUTTONDOWN && drawing == false) {
		drawing = true;
		mouse->height = mouse->currentimg.rows;
		mouse->width = mouse->currentimg.cols;
		Mat mask = Mat::zeros(Size(pParent->m_iOrgW, pParent->m_iOrgH), CV_8UC1);
		mouse->mask = mask; 

		Point imagePoint = pParent->WindowToImagePoint(x, y);
		mouse->startloc = imagePoint;
		mouse->endloc = imagePoint;
		//cout << x << " " << y << endl;
		mask.release();
	}
	
	else if (event == EVENT_LBUTTONUP) {
		const bool wasDrawing = drawing;
		drawing = false;

		// A release without a press in this window (or the second release of a
		// double-click) has no fresh mask to process.
		if (!wasDrawing || mouse->mask.empty())
			return;

		pParent->Mouse=prep(mouse);
		if (!pParent->Mouse->changed)
			return;

		pParent->m_matSrc = pParent->Mouse->currentimg;

		pParent->history.push_back(pParent->m_matSrc.clone());
		if (pParent->history.size() > kMaxUndoSteps + 1)
			pParent->history.pop_front();           // Drop the oldest state

		pParent->RefreshImage();
	}
	else if (event == EVENT_MBUTTONDOWN) {
		if (pParent->history.size() > 1) {
			pParent->history.pop_back();             // Remove the current state
			pParent->m_matSrc = pParent->history.back().clone();  // Revert to previous state
			pParent->RefreshImage();

		}

		else {
			cout << "No more undos available!" << endl;
		}
	}
}
// OpenCV calls this from its window message loop, where an escaping exception
// would terminate the program, so every failure is caught and reported here.
void MouseCall (int event, int x, int y, int flag, void* param)
{
	try {
		HandleMouse(event, x, y, flag, param);
	}
	catch (const std::exception& e) {
		cerr << "Operation failed: " << e.what() << endl;
		drawing = false;

		// prep() edits the image in place, so go back to the last complete state.
		COpenCVWindowExt* pParent = (COpenCVWindowExt*)(param);
		try {
			if (!pParent->history.empty()) {
				pParent->m_matSrc = pParent->history.back().clone();
				pParent->RefreshImage();
			}
		}
		catch (const std::exception& e2) {
			cerr << "Could not restore the image: " << e2.what() << endl;
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
	m_dCompensationX = 0;
	m_dCompensationY = 0;
	m_dInitialScale = 2;
	m_dNewScale = 1;
	m_dScaleRatio = 1.25;

	m_iOrgW = 0;
	m_iOrgH = 0;
	m_viewSize = Size(kMaxViewWidth, kMaxViewHeight);

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
	delete Mouse;
	Mouse = nullptr;
}

bool COpenCVWindowExt::ImRead (String strFileName)
{
	m_matSrc = imread (strFileName);
	if (m_matSrc.empty ())
		return false;
	history.clear();
	history.push_back(m_matSrc.clone());

	m_iOrgW = m_matSrc.cols;
	m_iOrgH = m_matSrc.rows;
	const double fitScale = min(kMaxViewWidth / static_cast<double>(m_iOrgW),
		kMaxViewHeight / static_cast<double>(m_iOrgH));
	m_viewSize = Size(max(1, static_cast<int>(round(m_iOrgW * fitScale))),
		max(1, static_cast<int>(round(m_iOrgH * fitScale))));
	//cout << m_matSrc.cols * m_dInitialScale << " " << m_matSrc.rows * m_dInitialScale << endl;
	resize (m_matSrc, m_matView, m_viewSize);
	if (!m_matView.empty()) {
		imshow(m_strWindowName, m_matView);
		Mouse->img = m_matView;
		Mouse->width = m_matView.cols; //1280
		Mouse->height = m_matView.rows; //854
		Mouse->currentimg = m_matSrc;
	}

	return !m_matView.empty ();
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

	// Render only the visible part of the source, never the whole zoomed image.
	// The margin keeps the interpolation neighbours inside the cropped region.
	Rect rectSrc = CurrentSourceRect();
	if (rectSrc.area() <= 0)
		return;
	rectSrc = Rect(rectSrc.x - 2, rectSrc.y - 2, rectSrc.width + 4, rectSrc.height + 4)
		& Rect(0, 0, m_matSrc.cols, m_matSrc.rows);

	// View pixel (x, y) shows the source point below; this is the same mapping
	// WindowToImagePoint() uses, so the brush stays aligned with the display.
	const double dStepX = m_iOrgW / max(1.0, m_dNewScale * m_viewSize.width);
	const double dStepY = m_iOrgH / max(1.0, m_dNewScale * m_viewSize.height);
	const Matx23d viewToSrc(
		dStepX, 0, dStepX * (m_iHorzScrollBarPos + 0.5) - 0.5 - rectSrc.x,
		0, dStepY, dStepY * (m_iVertScrollBarPos + 0.5) - 0.5 - rectSrc.y);

	warpAffine (m_matSrc (rectSrc), m_matView, viewToSrc, m_viewSize,
		INTER_LINEAR | WARP_INVERSE_MAP, BORDER_REPLICATE);
	imshow (m_strWindowName, m_matView);
	Mouse->currentimg = m_matSrc;
}

Point COpenCVWindowExt::WindowToImagePoint(int x, int y) const
{
	if (m_iOrgW <= 0 || m_iOrgH <= 0 || m_viewSize.width <= 0 || m_viewSize.height <= 0) {
		return Point(0, 0);
	}

	const double scaledWidth = max(1.0, m_dNewScale * m_viewSize.width);
	const double scaledHeight = max(1.0, m_dNewScale * m_viewSize.height);
	const int imageX = static_cast<int>((m_iHorzScrollBarPos + x) * m_iOrgW / scaledWidth);
	const int imageY = static_cast<int>((m_iVertScrollBarPos + y) * m_iOrgH / scaledHeight);

	return Point(ClampInt(imageX, 0, m_iOrgW - 1), ClampInt(imageY, 0, m_iOrgH - 1));
}

Rect COpenCVWindowExt::CurrentSourceRect() const
{
	if (m_iOrgW <= 0 || m_iOrgH <= 0 || m_viewSize.width <= 0 || m_viewSize.height <= 0) {
		return Rect();
	}

	const double scaledWidth = max(1.0, m_dNewScale * m_viewSize.width);
	const double scaledHeight = max(1.0, m_dNewScale * m_viewSize.height);

	const int x0 = ClampInt(static_cast<int>(m_iHorzScrollBarPos * m_iOrgW / scaledWidth), 0, m_iOrgW - 1);
	const int y0 = ClampInt(static_cast<int>(m_iVertScrollBarPos * m_iOrgH / scaledHeight), 0, m_iOrgH - 1);
	const int x1 = ClampInt(static_cast<int>((m_iHorzScrollBarPos + m_viewSize.width) * m_iOrgW / scaledWidth), x0 + 1, m_iOrgW);
	const int y1 = ClampInt(static_cast<int>((m_iVertScrollBarPos + m_viewSize.height) * m_iOrgH / scaledHeight), y0 + 1, m_iOrgH);

	return Rect(Point(x0, y0), Point(x1, y1));
}



void COpenCVWindowExt::SetHorzBarPos (int iPos)
{
	if (iPos > m_iHorzScrollBarRange_Max)
		m_iHorzScrollBarPos = m_iHorzScrollBarRange_Max;
	else if (iPos < 0)
		m_iHorzScrollBarPos = m_iHorzScrollBarRange_Min;
	else
		m_iHorzScrollBarPos = iPos;
	//cout << m_iHorzScrollBarPos << " ";

}

void COpenCVWindowExt::SetVertBarPos (int iPos)
{
	if (iPos > m_iVertScrollBarRange_Max)
		m_iVertScrollBarPos = m_iVertScrollBarRange_Max;
	else if (iPos < 0)
		m_iVertScrollBarPos = m_iVertScrollBarRange_Min;
	else
		m_iVertScrollBarPos = iPos;

	//cout << m_iVertScrollBarPos << endl;
}
