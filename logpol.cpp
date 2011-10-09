#include <iostream>
#include <math.h>
#include <cv.h>
#include <highgui.h>

//Use math constants defined in math.h
#define _USE_MATH_DEFINES

//Global variables
int mouse_x = -1,mouse_y = -1, sampling_angles = 64;
cv::Mat src;
int write_count=0;
std::string file_open;
//This function computes the Log-Polar Transformm and its inverse
void computeLogPol(){

  cv::imshow("Inverse-Retinal", src);
}

void help(){
    std::cout   << "This program demonstrates the logpolartransform and"
		<< "the inverse into euclidean sapace or the retinal image.\n"
            	<< "Usage:\n"
		<< "./assignment1 <num_sampling_angles> <image_name>"
		<< "Default is 64 Cameraman.tif.\n" << std::endl;
}

void mouseHandler(int event, int y, int x, int flags, void *param){
  if( event == CV_EVENT_LBUTTONDOWN ){
    mouse_x = x;
    mouse_y = y;
    std::cout << "Point " << x << " , " << y << " selected" << std::endl;
    computeLogPol();
  }
}

int main(int argc, char** argv)
{
  sampling_angles      = argc >= 2 ? atoi(argv[1]): 64;
  const char *filename = argc >= 3 ? argv[2]      : "Cameraman.tif";
  file_open = filename;

  src = cv::imread(filename, 0);
  if(src.empty()){
    help();
    std::cout << "can not open " << filename << std::endl;
    return -1;
  }

  //Load image for display
  IplImage *img0;
  img0 = cvLoadImage(filename, CV_LOAD_IMAGE_GRAYSCALE);

  //Create main window and set mouse call back if mouse event and display
  cvNamedWindow("Source", CV_WINDOW_AUTOSIZE);
  cvSetMouseCallback( "Source", mouseHandler, NULL );
  cvShowImage("Source", img0);

  //Create the default display with logpol at the center of the image
  computeLogPol();
  
  while(1) {
    //Wait for keyboard input
    char key = cvWaitKey(0);
    //Key pressed, quit the program
    if ( key ) break;
  }

  return 0;
}

