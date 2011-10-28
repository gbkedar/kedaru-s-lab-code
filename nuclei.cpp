#include <iostream>
#include <math.h>
#include <cv.h>
#include <highgui.h>
#include <omp.h>
#include <utility>
#include <limits.h>

//Use math constants defined in math.h
#define _USE_MATH_DEFINES

//A indexing convenience
#define BIP(a,b) BinImPad.at<unsigned char>(a,b)

//Global variables
std::string file_open;
bool write_files = true;
int rows,cols;

typedef struct { int index; double bin_count; } hist_bin;
bool binny_sort( hist_bin i, hist_bin j ) { return ((i.bin_count)<(j.bin_count)); }
//Computes a threshold based on the histogram of the image
int GetThresh(cv::Mat src){
  //Assume read image is of type unsigned char
  //Compute histogram
  std::vector<float> bins(256,0);
  int threshold;
  rows = src.rows;
  cols = src.cols;
  for( int y=0; y<src.rows; ++y )
    for( int x=0; x<src.cols; ++x )
      ++bins[src.at<unsigned char>(x,y)];

  //Smooth the histogam
  for( int i=0; i<256; ++i ){
    int low  = (i-2)<0   ? 0   : (i-2);
    int high = (i+2)>255 ? 255 : (i+2);
    double sum = 0;
    for( int j=low; j<=high; ++j ) sum += bins[j];
    sum /= (high-low+1);
    bins[i] = sum;
  }

  std::vector<hist_bin> binny;
  for( int i=0; i<256; ++i ) {
    hist_bin temp;
    temp.index = i; temp.bin_count = bins[i];
    binny.push_back( temp );
  }
  sort(binny.begin(), binny.end(), binny_sort);

  //Since they peaks are well separated we find the modes
  int mode1, mode2;
  mode1 = binny.back().index;

  bool binny2f = false;
  while( !binny2f ){
    binny.pop_back();
    mode2 = binny.back().index; binny2f = true;
    //Check if it is the highest in a 5-neighborhood
    int low  = (mode2-2)<0   ? 0   : (mode2-2);
    int high = (mode2+2)>255 ? 255 : (mode2+2);
    for( int j=low; j<=high; ++j ) if( bins[j]>binny.back().bin_count ) binny2f = false;
  }
  return ((mode1+mode2)/2);
}

//This function binarizes the nuclei in an image and colors connected components
void computebin(cv::Mat src, cv::Mat BinIm ){
  //Get Image Dimensions
  int thresh = GetThresh(src);
  int nthreads, tid;

  #pragma omp parallel private(tid,nthreads)
  {
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    int start, end;
    start = tid*(rows/nthreads);
    end   = (tid+1)*(rows/nthreads);
    if( tid == (nthreads-1) ) end = rows;
    for( int j=start; j<rows; ++j )
      for( int i=0; i<cols; ++i )
	if( src.at<unsigned char>(i,j) > thresh  ){
	  BinIm.at<unsigned char>(i,j)  = 0;
	}
	else{
	  BinIm.at<unsigned char>(i,j)	= 255;
	}
  }

  if(write_files){
    std::string out_str1 = "binary_" + file_open;
    cv::imwrite(out_str1.c_str(),BinIm);
  }
  std::cout<<"The threshold is: "<<thresh<<std::endl;
}


void help(){
    std::cout   << "This program demonstrates a naieve nucleus segmentation"
            	<< "Usage:\n"
		<< "./assignment2 <image_name>" << std::endl;
}

int main(int argc, char** argv)
{
  const char *filename = argc > 1 ? argv[2] : "cells.png";
  file_open = filename;

  cv::Mat src;
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
  cvShowImage("Source", img0);

  rows = src.rows;
  cols = src.cols;
  cv::Mat BinIm;
  BinIm = cvCreateMat(cols,rows,src.type());

  //Create the default display with logpol at the center of the image
  computebin(src,BinIm);
  while(1) {
    //Wait for keyboard input
    char key = cvWaitKey(0);
    //Key pressed, quit the program
    if ( key ) break;
  }

  return 0;
}

