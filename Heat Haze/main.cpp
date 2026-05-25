#include "WindowZoom.h"

#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	COpenCVWindowExt window("Src");
	window.SetInitailScale(1);

	const std::string imagePath = argc > 1 ? argv[1] : "DSC_1744.jpg";
	if (!window.ImRead(imagePath)) {
		std::cerr << "Could not open image: " << imagePath << std::endl;
		std::cerr << "Usage: Heat Haze.exe [image path]" << std::endl;
		return 1;
	}

	std::cout << "Left-drag: select haze line | mouse wheel: zoom | right-drag: pan | middle-click: undo" << std::endl;
	waitKey(0);
	destroyAllWindows();

	return 0;
}
