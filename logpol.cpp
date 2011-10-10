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
  //Get Image Dimensions
  int rows,cols;
  rows = src.rows;
  cols = src.cols;

  //Some variable declarations
  const int selPtx =( mouse_x <= 0 || mouse_y <= 0 ) ? (cols/2) : mouse_x;
  const int selPty =( mouse_x <= 0 || mouse_y <= 0 ) ? (rows/2) : mouse_y;
  cv::Mat logp, ilogp;

  //First call to function by default computes the logpoltransform wrt the
  //center of the image
  if( mouse_x <= 0 || mouse_y <= 0 ){
    std::cout<<"Creating default logpol transform\n";
  }
  //Subsequent calls are due to mouse press events
  else{
    std::cout<<"Creating logpol transform after mouse event\n";
    cvDestroyWindow("Log-Polar");
    cvDestroyWindow("Inverse-Retinal");
  }

  //Compute normalization factor so that you sample as much as possible
  //within the image bounds
  double norm_fact = cols>rows?(cols/2):(rows/2);
  norm_fact = norm_fact/log(norm_fact);

  //The size of the log-polar image is always im_size/2
  //You can chage it relative to the point selected but eh, I'm lazy
  int logp_h = rows>cols ? (rows/2):(cols/2);
  logp = cvCreateMat(logp_h,sampling_angles,src.type());

  //Interpolate through all the pixels in the polar image
  for( int i=0; i<sampling_angles; ++i ){ // i == theta
    //reduce the angle to (0,pi/2]
    double red_ang;
    red_ang = 2*M_PI*((double)i)/((double)sampling_angles);
    int count = 0;
    while( M_PI/2 < red_ang ){
      red_ang -= (M_PI/2);
      ++count;
    }
    //Get the sine and cosine of the reduced angle
    double red_ang_sin, red_ang_cos;
    if( count==1 || count==3 ){
      red_ang_sin = cos( red_ang );
      red_ang_cos = sin( red_ang );
    }
    else{
      red_ang_sin = sin( red_ang );
      red_ang_cos = cos( red_ang );
    }

    for( int j=1; j<=logp_h; ++j){ //j == rho
      //Get the x,y position
      double dist = exp(((double)j)/norm_fact); //rho on the x-y plane
      double pos_x = red_ang_cos*dist;
      double pos_y = red_ang_sin*dist;

      if(count==0){
	pos_x = selPtx + pos_x;
	pos_y = selPty + pos_y;
      }
      else if(count==1){ 
	pos_x = selPtx - pos_x;
	pos_y = selPty + pos_y;
      }
      else if(count==2){ 
	pos_x = selPtx - pos_x;
	pos_y = selPty - pos_y;
      }
      else if(count==3){ 
	pos_x = selPtx + pos_x;
	pos_y = selPty - pos_y;
      }

      //Bounds check and set the output pixel to 0 if out of bounds
      if( pos_x >= cols || pos_x < 0 || pos_y >= rows || pos_y < 0 ){
	logp.at<unsigned char>(i,(j-1)) = 0.0;
      }
      //Use bilinear interpolation and set rest of the pixels
      else{
	//Following the convention on the wiki page
	double x2 = ceil( pos_x ), x1 = floor( pos_x );
	double y2 = ceil( pos_y ), y1 = floor( pos_y );
	double r1 = (x2-pos_x)/(x2-x1)*src.at<unsigned char>(x1,y1) 
		    + (pos_x-x1)/(x2-x1)*src.at<unsigned char>(x2,y1);
	double r2 = (x2-pos_x)/(x2-x1)*src.at<unsigned char>(x1,y2) 
		    + (pos_x-x1)/(x2-x1)*src.at<unsigned char>(x2,y2);
 	double p = (y2-pos_y)/(y2-y1)*r1 + (pos_y-y1)/(y2-y1)*r2;
	if( p>254 )
	  logp.at<unsigned char>(i,(j-1)) = 255;
	else if( p<0 )
	  logp.at<unsigned char>(i,(j-1)) = 0;
	else
	  logp.at<unsigned char>(i,(j-1)) = (unsigned char)p;
      }
    } 
  }
  cv::imshow("Log-Polar", logp);

  //Create the inverse image
  ilogp = cvCreateMat(cols,rows,src.type());

  for( int i=0; i<cols; ++i ){
    double x_diff = i-selPtx;
    for( int j=0; j<rows; ++j ){
      //Find the co-ordinates of (i,j) in the log-polar co-ordinates
      double y_diff = j-selPty;
      double arc_tan;
      if( x_diff != 0 )
	arc_tan = atan(fabs((y_diff/x_diff)));
      double theta = x_diff==0 ? 0 : arc_tan;
      if( y_diff >= 0 && x_diff < 0 )
	theta = M_PI - theta;
      else if( y_diff < 0 && x_diff < 0 )
	theta = M_PI + theta;
      else if( y_diff < 0 && x_diff > 0 )
	theta = 2*M_PI - theta;
      theta = theta*sampling_angles/(2*M_PI);
      double distance = (sqrt(x_diff*x_diff+y_diff*y_diff));
      //Use bilinear interpolation and set the pixels
      //Following the convention on the wiki page
      double x2 = ceil( theta ), x1 = floor( theta );
      double y2 = ceil( distance ), y1 = floor( distance );
      if( x1<0 || x1>=sampling_angles || y1<0 || y2>=logp_h ){
	ilogp.at<unsigned char>(i,j) = 0;
      }
      else{
	double r1 = (x2-theta)/(x2-x1)*logp.at<unsigned char>(x1,y1) 
		    + (theta-x1)/(x2-x1)*logp.at<unsigned char>(x2,y1);
	double r2 = (x2-theta)/(x2-x1)*logp.at<unsigned char>(x1,y2) 
		    + (theta-x1)/(x2-x1)*logp.at<unsigned char>(x2,y2);
 	double p = (y2-distance)/(y2-y1)*r1 + (distance-y1)/(y2-y1)*r2;
	ilogp.at<unsigned char>(i,j) = p;
      }
    }
  }
  cv::imshow("Inverse-Retinal", ilogp);
  if(1){
    std::stringstream ss;
    ss<<write_count;
    std::string out_str1 = "log_pol_" + ss.str() + "_" + file_open;
    std::string out_str2 = "inv_ret_" + ss.str() + "_" + file_open;
    cv::imwrite(out_str1.c_str(),logp);
    cv::imwrite(out_str2.c_str(),ilogp);
    ++write_count;
  }
  return;
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

