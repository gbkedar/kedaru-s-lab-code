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
//[(b*(cols+2)+a)]

//Global variables
std::string file_open;
bool write_files = true;
int rows,cols;

//Computes the area and the centroid for a given label
bool GetStats(cv::Mat LabelIm, int i, int j, cv::Mat BinImPad, std::vector<int> &one_stat){
  //Seed with the first pixel
  std::pair<int,int> seed(i,j);
  std::vector< std::pair<int,int> > indices, queueue;
  indices.push_back(seed); queueue.push_back(seed);
  BIP(i,j) = 0;
  double x_men = 0, y_men = 0; //centroid

  //Find indices of all pixels with the same label
  while( !queueue.empty() ){
    std::pair<int, int> index = queueue.back();
    queueue.pop_back();
    //Check 8-neighborhood
    for( int k=-1; k<2; ++k )
      for( int l=-1; l<2; ++l )
	if( (k==0) && (l==0) ) continue;
	  else{
	    if( BIP( (index.first+k),(index.second+l) ) ){
	      std::pair<int,int> sed( (index.first+k),(index.second+l) );
	      queueue.push_back( sed );
	      indices.push_back( sed );
	      LabelIm.at<unsigned char>( sed.first, sed.second ) = BIP( sed.first,sed.second );
	      BIP( sed.first, sed.second ) = 0;
	      x_men += sed.first;
	      y_men += sed.second;
	    }
          }
  }
  x_men = round(x_men/((double)indices.size()));
  y_men = round(y_men/((double)indices.size()));
  one_stat[1] = indices.size();
  one_stat[2] = (int)x_men;
  one_stat[3] = (int)y_men;

  if( (indices.size())>=15 ){ return true; }
  else{
    std::vector< std::pair<int,int> >::iterator it;
    for( it=indices.begin() ; it < indices.end(); ++it )
      LabelIm.at<unsigned char>(it->first, it->second) = 0;
    return false;
  }
}

//Color one blob for one thread
void ColorBlob( cv::Mat BinImPad, std::vector< std::pair<int,int> > *indices,
	const int &start, const int &end, unsigned char &current_index, unsigned char &replace_ind ){
  //Check the 8-connected neighborhood and set them if they are within bounds
  std::pair<int,int> CI = indices->back();
  indices->pop_back();
  if( (CI.second-1)>=start ){//Bounds checking for current thread
    for( int i=-1; i<2; ++i ) if( BIP((CI.first+i),(CI.second-1)) ==replace_ind ){
      BIP( (CI.first+i),(CI.second-1) ) = current_index;
      std::pair<int,int> PushInd( (CI.first+i),(CI.second-1) );
      indices->push_back( PushInd );
    }
  }
  if( (CI.second+1)<end ){//Bounds checking for next thread
    for( int i=-1; i<2; ++i ) if( BIP((CI.first+i),(CI.second+1)) ==replace_ind ){
      BIP( (CI.first+i),(CI.second+1) ) = current_index;
      std::pair<int,int> PushInd( (CI.first+i),(CI.second+1) );
      indices->push_back( PushInd );
    }
  }
  if( BIP((CI.first-1),(CI.second))==replace_ind ){
    BIP( (CI.first-1),(CI.second) ) = current_index;
    std::pair<int,int> PushInd( (CI.first-1),(CI.second) );
    indices->push_back( PushInd );

  }
  if( BIP((CI.first+1),(CI.second))==replace_ind ){
    BIP( (CI.first+1),(CI.second) ) = current_index;
    std::pair<int,int> PushInd( (CI.first+1),(CI.second) );
    indices->push_back( PushInd );
  }
  if(!indices->empty()){
    ColorBlob( BinImPad, indices, start, end, current_index, replace_ind );
  }
}

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
void computebin(cv::Mat src, cv::Mat BinImPad, cv::Mat BinIm ){
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
	  BIP((i+1),(j+1))	 	= 0;
	}
	else{
	  BinIm.at<unsigned char>(i,j)	= 255;
	  BIP((i+1),(j+1)) 		= UCHAR_MAX;
	}
  }

  //Set the borders of the padded image
  for( int i=0; i<(cols+2); ++i ){
    BIP(i,0)		= 0;
    BIP(i,(rows+1))	= 0;
  }
  for( int i=1; i<(rows+2); ++i ){
    BIP(0,i)		= 0;
    BIP((cols+1),i)	= 0;
  }

  if(write_files){
    std::string out_str1 = "binary_" + file_open;
    cv::imwrite(out_str1.c_str(),BinIm);
  }
  std::cout<<"The threshold is: "<<thresh<<std::endl;
}

void computeLabels(cv::Mat src, cv::Mat BinImPad){
  cv::Mat LabelIm  = cvCreateMat(cols,rows,src.type());
  for( int j=0; j<rows; ++j )
    for( int i=0; i<cols; ++i ) LabelIm.at<unsigned char>(i,j)=0;
  int nthreads, tid;
  //omp_set_num_threads(4);
  bool redo_coloring = true;
  #pragma omp parallel private(tid,nthreads)
  {
    //Color blob recursively as soon as a pixel in FG is found
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    int end1;
    unsigned char start_index, label_index;
    label_index = UCHAR_MAX * tid / nthreads;
    ++label_index; start_index = label_index;
    const int start = tid*(rows/nthreads) + 1;
    end1 = (tid+1)*(rows/nthreads) + 1;
    if( tid==(nthreads-1) ) ++end1; ++end1;
    const int end = end1;
    for( int j=start; j<end; ++j )
      for( int i=1; i<=cols; ++i ){
	if( BIP(i,j)==UCHAR_MAX ){
	  std::vector< std::pair<int,int> > indices;
	  std::pair<int,int> one_index; one_index.first = i; one_index.second = j;
	  indices.push_back(one_index);
          BIP(i,j) = label_index;
	  unsigned char temp = UCHAR_MAX; 
	  ColorBlob( BinImPad, &indices, start, end, label_index, temp );
	  ++label_index;
	}
      }
    #pragma omp barrier
    //Re-color blob if the previous thread has assigned it
    //to a different color
    if( tid==0 && write_files ){
      std::string out_str1 = "label_intermediate.png";
      cv::imwrite(out_str1.c_str(),BinImPad);
    }
    while( redo_coloring ){
      #pragma omp barrier
      if( tid == 0 ) redo_coloring = false;
      for( int i=1; i<=cols; ++i ){
	if( (BIP((i-1),(start-1)) && BIP(i,start)) &&
	    (BIP((i-1),(start-1))!=BIP(i,start)) ){
	  redo_coloring = true;
	  unsigned char lund=BIP(i,start);
	  BIP(i,start)=BIP((i-1),(start-1));
	  std::pair<int,int> one_index; one_index.first = i; one_index.second = start;
	  std::vector< std::pair<int,int> > indices;
	  indices.push_back(one_index);
	  unsigned char temp = BIP(i,start);
	  ColorBlob( BinImPad, &indices, start, end, temp, lund );
	}
	else if( (BIP(i,(start-1)) && BIP(i,start)) &&
	    	(BIP(i,(start-1))!=BIP(i,start)) ){
	  redo_coloring = true;
	  unsigned char lund=BIP(i,start);
	  BIP(i,start)=BIP(i,(start-1));
	  std::pair<int,int> one_index; one_index.first = i; one_index.second = start;
	  std::vector< std::pair<int,int> > indices;
	  indices.push_back(one_index);
	  unsigned char temp = BIP(i,start);
	  ColorBlob( BinImPad, &indices, start, end, temp, lund );
	}
	else if( BIP((i+1),(start-1)) && BIP(i,start) &&
	   	 (BIP((i+1),(start-1))!=BIP(i,start)) ){
	  redo_coloring = true;
	  unsigned char lund=BIP(i,start);
	  BIP(i,start)=BIP((i+1),(start-1));
	  std::pair<int,int> one_index; one_index.first = i; one_index.second = start;
	  std::vector< std::pair<int,int> > indices;
	  indices.push_back(one_index);
	  unsigned char temp = BIP(i,start);
	  ColorBlob( BinImPad, &indices, start, end, temp, lund );
	}
      }
      #pragma omp barrier
      if( tid==0 && write_files ){
	std::string out_str1 = "label_intermediate1.png";
	cv::imwrite(out_str1.c_str(),BinImPad);
      }
    }
  }

  //Output the centroids and areas
  std::vector< std::vector<int> > statistics;
  for( int j=1; j<=rows; ++j )
    for( int i=1; i<=cols; ++i )
      if( ((int)BIP(i,j))>0){
	std::vector<int> one_stat(4); 
	one_stat[0] = BIP(i,j);
	bool bigEnough = GetStats(LabelIm,i,j,BinImPad,one_stat);
	if(bigEnough) statistics.push_back( one_stat );
	//else std::cout<<one_stat[0]<<"\t"<<one_stat[1]<<"\t"<<one_stat[2]<<"\t"<<one_stat[3]<<"\n";
      }

  std::string out_str1 = "label_final.png";
  cv::imwrite(out_str1.c_str(),LabelIm);
  out_str1 = "binary_" + file_open;
  IplImage *img = cvLoadImage(out_str1.c_str(), CV_LOAD_IMAGE_COLOR);
  CvFont font;
  cvInitFont(&font, CV_FONT_HERSHEY_SIMPLEX, 0.3, 0.3, 0, 1);
  std::vector< std::vector<int> >::iterator it; int i=1;
  for( it=statistics.begin() ; it < statistics.end(); ++it ){
    std::stringstream ss;
    ss<<i<<","<<it->at(1);++i;
    std::string sss=ss.str();
    cvPutText(img, sss.c_str(), cvPoint((it->at(3)),(it->at(2))), &font, cvScalar(255, 255, 100, 0));
    cvCircle(img, cvPoint((it->at(3)),(it->at(2))), 3, cvScalar(0,255,0), 1);
  }
  cvShowImage("Binary Image", img);

  out_str1 = "final.png";
  cvSaveImage(out_str1.c_str(),img);

  return;
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
  cv::Mat BinImPad;
  BinImPad = cvCreateMat((rows+2),(cols+2),src.type());
  cv::Mat BinIm;
  BinIm = cvCreateMat(cols,rows,src.type());

  //Create the default display with logpol at the center of the image
  computebin(src, BinImPad, BinIm);
  computeLabels(src, BinImPad);
  //cv::imshow("Binary Image", BinIm);
  while(1) {
    //Wait for keyboard input
    char key = cvWaitKey(0);
    //Key pressed, quit the program
    if ( key ) break;
  }

  return 0;
}

