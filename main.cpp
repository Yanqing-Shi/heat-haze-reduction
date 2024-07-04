

#include <iostream>
#include "WindowZoom.h"
int main(int argc, char* argv[])
{
    //std::cout << "Hello World!\n";
	COpenCVWindowExt window ("Src");
	window.SetInitailScale (1);
	
	window.ImRead ("Dsc_1744.jpg");
	waitKey (0);
	destroyAllWindows();
}
