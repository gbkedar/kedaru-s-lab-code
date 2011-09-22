#ifndef _FDDIJKSTRA_CXX
#define _FDDIJKSTRA_CXX

#include<iostream>
#include<math.h>
#include<fstream>
#include<vector>
#include<algorithm>
#include<limits.h>
#include<queue>

#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionConstIterator.h"

#define IM(i,j) image.at(i*size2+j)
#define DIM(i,j) dist_im.at(i*size2+j)
#define N(i,j) neighbor[i*num_cells+j]
#define D(i,j) Distance[i*num_cells+j]

class Edge{
	public:
	double weight;
	int x; //End vertex
	int y; //End vertex
};

bool edgecomparator( Edge a, Edge b ){ return( a.weight < b.weight ); }

double change_weights_to_geodesic( int cell_class, const char *channel_filename, float scaling, int *neighbor,
                                        double *Distance, int num_cells, int *x, int *y, int*cl ){
	typedef unsigned short PixelType;
	const   unsigned int Dimension = 2;
	typedef itk::Image< PixelType, Dimension > ImageType;

	typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleUsUsType;
	typedef itk::ImageFileReader< ImageType >  ReaderType;

	typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( channel_filename );

	try{
		reader->Update();
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return 0;
	}

	RescaleUsUsType::Pointer rescale = RescaleUsUsType::New();
	rescale->SetOutputMinimum( 0 );
	rescale->SetOutputMaximum( scaling );
	rescale->SetInput( reader->GetOutput() );   
	try{
		rescale->Update();
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return 0;
	}

	ImageType::Pointer scaledimg = rescale->GetOutput();
	ConstIteratorType pix_buf( scaledimg,scaledimg->GetRequestedRegion() );

	int size1 = scaledimg->GetLargestPossibleRegion().GetSize()[0];                                                             
	int size2 = scaledimg->GetLargestPossibleRegion().GetSize()[1];

	std::vector<int> image;

	for( int i=0; i<size1; ++i ){
		for( int j=0; j<size2; ++j ){
			ImageType::IndexType curpixind;
			curpixind[0]=i;
			curpixind[1]=j;
			pix_buf.SetIndex(curpixind);
			int curpix = pix_buf.Get();
			int temp;
			temp = curpix;
			image.push_back( temp );
		}
	}

	std::cout<<"Starting Geodesic weight computation\n";
	double mean=0;
	int countt=0;
	for( int i=0; i<num_cells; ++i ){

		//std::cout<<i<<"\t";

		int min_x = x[i], max_x = x[i];
		int min_y = y[i], max_y = y[i];
		int test_count = 0;
		for( int j=i+1; j<num_cells; ++j ){
			if( N(i,j) ){
				if( min_y > y[j] )
					min_y = y[j];
				if( max_y < y[j] )
					max_y = y[j];
				if( min_x > x[j] )
					min_x = x[j];
				if( max_x < x[j] )
					max_x = x[j];
				++test_count;
			}
		}
		if( min_y == max_y || max_x==min_x ){}
		else
		{
		//std::cout<<i<<" "<<test_count<<std::endl;
		min_y = ((min_y-50)>0) ? (min_y-50) : 0;
		min_x = ((min_x-50)>0) ? (min_x-50) : 0;
		max_y = ((max_y+50)<(size2-1)) ? (max_y+50) : (size2-1);
		max_x = ((max_x+50)<(size1-1)) ? (max_x+50) : (size1-1);

		int exx=x[i], why=y[i];
		if( exx==(size1-1) ) --exx;
		if( why==(size2-1) ) --why;
		if( exx==0 ) ++exx;
		if( why==0 ) ++why;

		std::queue<Edge> edges;
		for( int j=-1; j<2; ++j ){
			for( int k=-1; k<2; ++k ){
				if( j==0 && k==0 ){}
				else{
					Edge temp;
					int exxpj = exx+j, whypk = why+k;
					temp.x = exxpj;
					temp.y = whypk;
					double wtt;
					if( IM(exx,why) > IM(exxpj,whypk) )
						wtt = IM(exx,why) - IM(exxpj,whypk);
					else
						wtt = IM(exxpj,whypk) - IM(exx,why);
					temp.weight = wtt;
					edges.push( temp );
				}
			}
		}

		std::vector<double> dist_im((size1*size2),INT_MAX);

		DIM(exx,why) = 0;

		while( !edges.empty() ){
			Edge temp = edges.front();
			edges.pop();

			if( DIM(temp.x,temp.y) >= temp.weight ){
				DIM(temp.x,temp.y) = temp.weight;
				if( temp.x>min_x && temp.y>min_y && temp.x<max_x  && temp.y<max_y ){
					for( int j=-1; j<2; ++j ){
						for( int k=-1; k<2; ++k ){
							if( j==0 && k==0 ){}
							else{
								double new_weight;
								int xpj = temp.x+j, ypk = temp.y+k;
								if( IM(temp.x,temp.y) > IM(xpj,ypk) )
									new_weight = DIM(temp.x,temp.y) + 
										(IM(temp.x,temp.y)-IM(xpj,ypk))+0.01;
								else
									new_weight = DIM(temp.x,temp.y) + 
										(IM(xpj,ypk)-IM(temp.x,temp.y))+0.01;
								if( new_weight < DIM(xpj,ypk) ){
									Edge temp1;
									temp1.x = xpj;
									temp1.y = ypk;
									temp1.weight = new_weight;
									edges.push(temp1);
									DIM(xpj,ypk) = new_weight;
								}
							}
						}
					}
				}
			}
		}

		for( int j=i+1; j<num_cells; ++j ){
			if( N(i,j) ){
				std::cout<< DIM(x[j],y[j]) << " "<<D(i,j)<<"\t";
				mean += DIM(x[j],y[j]); ++countt;
				D(i,j) = DIM(x[j],y[j]);
				if( DIM(x[j],y[j])==INT_MAX )
					std::cout<<"Check Neighs: "<<i<<"\t"<<j<<std::endl;
			}
		}
		std::cout<< "\n";
		dist_im.clear();
		}
	}

	mean = mean/countt;
	double dummy=0;
	for( int i=0; i<num_cells; ++i )
		for( int j=i+1; j<num_cells; ++j )
			if( N(i,j) )
				dummy += (mean-D(i,j))*(mean-D(i,j))/mean;
	double stddev = sqrt(dummy);

	std::cout<<std::endl<<"Geodesic distance done\n";

	image.clear();

	return (mean+stddev);
}

#endif
