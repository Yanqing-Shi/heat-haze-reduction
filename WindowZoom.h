#include <opencv2/opencv.hpp>
#include<iostream>
#include<stack>
using namespace std;
using namespace cv;
#pragma once


typedef struct MouseParams
{
	Mat img;
	Mat mask;
	int width;
	int height;
	Mat currentimg;
	Point startloc, endloc;
};

class COpenCVWindowExt
{

public:
	COpenCVWindowExt (String strWindowName, int iFlag = 1);
	~COpenCVWindowExt ();
	stack<Mat> history;
	bool ImRead (String strFileName);
	void SetInitailScale (double dScale);
	MouseParams* Mouse;

	void RefreshImage ();
	void OverlapImage(int x,int y);
	Mat m_matSrc;

	vector<Mat> m_vecMatResize;
	String m_strWindowName;

	double m_dInitialScale;
	double m_dNewScale;
	double m_dScaleRatio;
	double m_dCompensationX;
	double m_dCompensationY;

	int m_iScaleTimes;
	int m_iMaxScaleTimes;
	int m_iMinScaleTimes;
	
	int m_iOrgW;
	int m_iOrgH;
	

	Point ptLButtonDown;
	Point ptRButtonDown;

	void SetHorzBarPos (int iPos);
	void SetVertBarPos (int iPos);
	int m_iHorzScrollBarPos;
	int m_iVertScrollBarPos;
	int m_iHorzScrollBarPos_copy;
	int m_iVertScrollBarPos_copy;


	int m_iHorzScrollBarRange_Min;
	int m_iHorzScrollBarRange_Max;
	int m_iVertScrollBarRange_Min;
	int m_iVertScrollBarRange_Max;


	
};



