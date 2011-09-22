#include<iostream>
#include<sstream>
#define _USE_MATH_DEFINES
#include<math.h>
#include<fstream>
#include<vector>
#include<algorithm>
#include<limits.h>
#include<iomanip>

#include "FDDijkstra.cxx"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkLineIterator.h"

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#define N(i,j) neighbor[i*num_cells+j]
#define D(i,j) Distance[i*num_cells+j]
#define OD(i,j) out_dist[i*new_rows+j]
#define OLL(i,j) out_list[i*new_rows+j]
#define EDC(i,j) edge_costs[(num_vars-(num_cells-i+(num_cells-i)*(num_cells-i+1)/2)+j-(i+1))]
#define E(i,j) EP[(num_vars-(num_cells-i+(num_cells-i)*(num_cells-i+1)/2)+j-(i+1))]
#define pi M_PI
#define TESTING_MIP true


class distcontainer{
public:
int cell1;
int cell2;
double dist;
//double a1;
//double a2;
};

bool comparator( distcontainer a, distcontainer b){ return( a.dist < b.dist ); }

void cast_rays(double *radius, int *x, int *y, int num_rays, int num_cells, int *neighbor, double *Distance){
//	int counter=0;
//	double md=0;
	for(int i=0; i<num_cells; ++i){
		for(int j=0; j<num_rays; ++j){
			std::vector<distcontainer> potential_neighs;
			int k;
			for(k=0; k<num_cells; ++k){
				//Ref wiki on raycasting
				//The index i is the point source 's' and j is the center
				//of the sphere 'c'
				double v[2];
				v[0] = (x[k]-x[i]);
				v[1] = (y[k]-y[i]);
				if( fabs(v[0])<0.1 ) v[0]=0.00001*v[0];
				if( fabs(v[1])<0.1 ) v[1]=0.00001*v[1];
				//Computing the direction of the ray 'd'
				double d_angle;
				if( j <= (num_rays/2) )
					d_angle = pi*j/(num_rays/2);
				else
					d_angle = pi + pi*(j-(num_rays/2))/(num_rays/2);
				double v_angle = atan(fabs((double)v[1])/fabs((double)v[0]));
				//Compute the direction of v
				if( v[0] < 0 && v[1] <= 0 )
					v_angle = pi - v_angle;
				else if ( v[0] <= 0 && v[1] > 0 )
					v_angle = pi + v_angle;
				else if ( v[0] > 0  && v[1] > 0 )
					v_angle = 2*pi - v_angle;
//				if( i == 0 && j==0 )
//					std::cout<<v[0]<<"\t"<<v[1]<<"\t"<<v_angle<<std::endl;
				double ang_diff;
				if( v_angle>d_angle )
					ang_diff = v_angle-d_angle;
				else
					ang_diff = d_angle-v_angle;
				if( ang_diff > pi )
					ang_diff = 2*pi-ang_diff;
//				if( i==0 && k==3 ) std::cout<< ang_diff << std::endl;
				if( ang_diff > (pi/2) )
					continue;
				double VdotD = sqrt(v[0]*v[0]+v[1]*v[1])*cos(ang_diff);
				VdotD = VdotD*VdotD;
//				if( VdotD < 0 ) std::cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<VdotD<<std::endl;
				double VminR = v[0]*v[0]+v[1]*v[1]-(double)(radius[k]*radius[k]);
				if( VdotD > VminR ){
					distcontainer temp;
					temp.cell1 = i;
					temp.cell2 = k;
					double distance = sqrt(v[0]*v[0]+v[1]*v[1]);
					temp.dist  = distance;
//					temp.a1 = v_angle;
//					temp.a2 = d_angle;
//					if( distance < 0 )
//						std::cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<distance<<std::endl;
					potential_neighs.push_back(temp);
//					if( i == 1 )
//						std::cout<<"br hr\n";
					//potential_neighs
//					++N(i,k)=1;
//					++N(k,i)=1;
//					if( i==1 && (k==0||k==3) ){
//						++counter;
//						std::cout<<"Ray: "<<d_angle<<"\tv_angle: "<<v_angle<<
//						"\tangle_diff: "<<ang_diff<<"\tVDD: "<<VdotD<<"\tVMR: "<<VminR<<"\t"<<counter<<std::endl;
//					}
//					if( k == 65 )
//						std::cout<<"br hr\n";
				}
			}
			//Add code to check first intersection with ray
			if(!potential_neighs.empty()){
				sort( potential_neighs.begin(), potential_neighs.end(), comparator );
				++N( potential_neighs.at(0).cell1, potential_neighs.at(0).cell2 );// = 1;
				++N( potential_neighs.at(0).cell2, potential_neighs.at(0).cell1 );// = 1;
				D( potential_neighs.at(0).cell2, potential_neighs.at(0).cell1 )= potential_neighs.at(0).dist;
				D( potential_neighs.at(0).cell1, potential_neighs.at(0).cell2 )= potential_neighs.at(0).dist;
//				if(potential_neighs.size()>1)
//				if( potential_neighs.at(0).cell1==1 &&  (potential_neighs.at(0).cell2==0||potential_neighs.at(0).cell2==3) ) {
//					std::cout<<"Ray: "<<potential_neighs.at(0).a2<<"\tv_angle: "<<potential_neighs.at(0).a1<<"\tangle_diff: "
//					<<"\t\tVDD: "<<"\t\tVMR: "<<"\tN:"<<N( potential_neighs.at(0).cell2, potential_neighs.at(0).cell1 )
//					<<"\tD:"<< D( potential_neighs.at(0).cell2, potential_neighs.at(0).cell1 )<<std::endl;
//				}

//				if( md < potential_neighs.at(0).dist )
//					md = potential_neighs.at(0).dist;
			}
//			if( j == 33 )
//				std::cout<<"break here\n";
			potential_neighs.clear();
		}
	}
//	std::cout<<"md="<<md<<"\n";
//	std::cout<<counter<<std::endl;
}


/*
void change_weights_to_geodesic( int cell_class, const char *channel_filename, float scaling, int *neighbor,
					double *Distance, int num_cells, int *x, int *y, int*cl ){
	typedef unsigned short PixelType;
	const   unsigned int Dimension = 2;
	typedef itk::Image< PixelType, Dimension > ImageType;

	typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleUsUsType;
	typedef itk::ImageFileReader< ImageType >  ReaderType;
	typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractImageType;
	typedef itk::FastMarchingImageFilter< ImageType, ImageType > FastMarchFilterType;

	typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( channel_filename );

	try{
		reader->Update();
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
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
		return;
	}

	ImageType::Pointer scaledimg = rescale->GetOutput();

	ConstIteratorType pix_buf( scaledimg,scaledimg->GetRequestedRegion() );
	pix_buf.GoToBegin();
	std::cout<<"First index in the image: "<<pix_buf.GetIndex()[0]<<"\t"<< pix_buf.GetIndex()[1]<<std::endl;
	++pix_buf;
	std::cout<<"Second index in the image: "<<pix_buf.GetIndex()[0]<<"\t"<< pix_buf.GetIndex()[1]<<std::endl;

	int size1 = scaledimg->GetLargestPossibleRegion().GetSize()[0];
	int size2 = scaledimg->GetLargestPossibleRegion().GetSize()[1];


	for( int i=0; i<num_cells; ++i ){
		;
	}


	for( int i=0; i<num_cells; ++i ){
		int min_x = x[i], max_x = x[i];
		int min_y = y[i], max_y = y[i];
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
			}
		}
		min_y = ((min_y-40)>0) ? (min_y-40) : 0;
		min_x = ((min_x-40)>0) ? (min_x-40) : 0;
		max_y = ((max_y+40)<size2) ? (max_y+40) : size2;
		max_x = ((max_x+40)<size1) ? (max_x+40) : size1;

		int sizee1 = max_x-min_x+1;
		int sizee2 = max_y-min_y+1;
		const int V=size1*size2;

		typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;
 
		typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
			boost::no_property, EdgeWeightProperty > Graph;
        
		typedef boost::graph_traits < Graph >::vertex_descriptor vertex_descriptor;
		typedef boost::graph_traits < Graph >::edge_descriptor edge_descriptor;
		typedef std::pair<int, int> Edge;

		// Create a graph
		Graph g;

		std::vector<Graph::vertex_descriptor> vertices;
		for( int k=0; k<V; ++k ){
			Graph::vertex_descriptor v0 = boost::add_vertex(g);
			vertices.push_back(v0);
		}

		// Add weighted edges
		for( int h=sizee1; h<=max_x; ++h ){
			//std::cout<<i<<"\t";
			for( int j=sizee2; j<=max_y; ++j ){
				ImageType::IndexType curpixind;
				curpixind[0]=h;
				curpixind[1]=j;
				pix_buf.SetIndex(curpixind);
				int curpix = pix_buf.Get();
				for( int k=-1; k<1; ++k ){
					for( int l=-1; l<1; ++l ){
						if( k==0 && l==0 ){}
						else{
							curpixind[0]=h+k;
							curpixind[1]=j+l;
							pix_buf.SetIndex(curpixind);
							int curpix1 = pix_buf.Get();
							int diff;
							if( curpix1 > curpix )
								diff = curpix1-curpix;
							else
								diff = curpix-curpix1;
							EdgeWeightProperty weight0(diff);
							boost::add_edge(vertices.at(h*size2+j), vertices.at((h+k)*size2+(j+l)), weight0, g);
						}
					}
				}
			}
		}

		// Compute shortest paths from v0 to all vertices, and store the output in parents and distances
		std::vector<int> distances; // To store distances
		distances.resize(V);
		std::vector<vertex_descriptor> parents(V); // To store parents
		std::cout<<"\tstarting dijkstra "<<i;
		boost::dijkstra_shortest_paths(g, vertices.at(x[i]*size2+y[i]), boost::predecessor_map(&parents[0]).distance_map(&distances[0]));
		std::cout<<" DD";

		for( int j=i+1; j<num_cells; ++j ){
			if( D(i,j) > 0.00001 ){
				double distance;
				distance=distances[vertices.at(x[j]*size2+y[j])];
				D(i,j)=distance;
				D(j,i)=distance;
			}
		}
		parents.clear();
		distances.clear();
		std::cout<<" JJ"<<std::endl;
	}
*/
/*
	//Start: Johnsons all pairs shortest path
	typedef boost::adjacency_list < boost::listS, boost::vecS, boost::directedS, boost::no_property,
		boost::property< boost::edge_weight_t, int, boost::property< boost::edge_weight2_t, int > > > Graph;

	typedef std::pair<int, int> Edge;
	
	std::vector<Edge> edge_array;
	std::vector<int> weights;
	for( int i=1; i<size1; ++i )
		for( int j=1; j<size2; ++j ){
			ImageType::IndexType curpixind;
			curpixind[0]=i;
			curpixind[1]=j;
			pix_buf.SetIndex(curpixind);
			int curpix = pix_buf.Get();
			for( int k=-1; k<1; ++k )
				for( int l=-1; l<1; ++l ){
					if( k==0 && l==0 ){}
					else{
						curpixind[0]=i+k;
						curpixind[1]=j+l;
						pix_buf.SetIndex(curpixind);
						int curpix1 = pix_buf.Get();
						int diff;
						if( curpix1 > curpix )
							diff = curpix1-curpix;
						else
							diff = curpix-curpix1;
						Edge edge1=Edge(i*size1+j,(i+k)*size1+j+l);
						Edge edge2=Edge((i+k)*size1+j+l,i*size1+j);
						edge_array.push_back(edge1);
						edge_array.push_back(edge2);
						weights.push_back( diff );
						weights.push_back( diff );
					}
				}
		}
	const std::size_t E = edge_array.size();
	const int V=size1*size2;

	Graph g(V);
	for (std::size_t j = 0; j < E; ++j)
		boost::add_edge(edge_array.at(j).first, edge_array.at(j).second, g);

	boost::property_map < Graph, boost::edge_weight_t >::type w = get(boost::edge_weight, g);

	int ind=0;
	boost::graph_traits < Graph >::edge_iterator e, e_end;
	for (boost::tie(e, e_end) = edges(g); e != e_end; ++e,++ind)
		w[*e]=weights.at(ind);

	std::vector<int>d(V,(std::numeric_limits < int >::max)());
	int **D;
	D = (int**) malloc(V*sizeof(int*));
	for( int i=0; i<V; ++i )
		D[i] = (int*) malloc(V*sizeof(int));
	std::cout<<"HERE\n";
	boost::johnson_all_pairs_shortest_paths(g, D, boost::distance_map(&d[0]));
	std::cout<<"HERE1\n";
	//End:Johnsons all pairs shortest path
*/
/*
	//Start: Dijkstra
	typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;
 
	typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
		boost::no_property, EdgeWeightProperty > Graph;
        
	typedef boost::graph_traits < Graph >::vertex_descriptor vertex_descriptor;
	typedef boost::graph_traits < Graph >::edge_descriptor edge_descriptor;
	typedef std::pair<int, int> Edge;

	// Create a graph
	Graph g;

	std::vector<Graph::vertex_descriptor> vertices;
	for( int i=0; i<(size1*size2); ++i ){
		Graph::vertex_descriptor v0 = boost::add_vertex(g);
		vertices.push_back(v0);
	}

	// Add weighted edges
	for( int i=1; i<size1; ++i ){
		std::cout<<i<<"\t";
		for( int j=1; j<size2; ++j ){
			ImageType::IndexType curpixind;
			curpixind[0]=i;
			curpixind[1]=j;
			pix_buf.SetIndex(curpixind);
			int curpix = pix_buf.Get();
			for( int k=-1; k<1; ++k ){
				for( int l=-1; l<1; ++l ){
					if( k==0 && l==0 ){}
					else{
						curpixind[0]=i+k;
						curpixind[1]=j+l;
						pix_buf.SetIndex(curpixind);
						int curpix1 = pix_buf.Get();
						int diff;
						if( curpix1 > curpix )
							diff = curpix1-curpix;
						else
							diff = curpix-curpix1;
						EdgeWeightProperty weight0(diff);
						boost::add_edge(vertices.at(i*size2+j), vertices.at((i+k)*size2+(j+l)), weight0, g);
					}
				}
			}
		}
	}

	const int V=size1*size2;

	std::vector<int>d(V,(std::numeric_limits < int >::max)());
	int **D;
	D = (int**) malloc(V*sizeof(int*));
	for( int i=0; i<V; ++i )
		D[i] = (int*) malloc(V*sizeof(int));
	std::cout<<"HERE\n";

	// Compute shortest paths from v0 to all vertices, and store the output in parents and distances
	std::cout<<"Starting computation of geodesic distances\n";
	std::cout<<num_cells<<" Cells. Processing cell number:\n ";
	for( int i=0; i<num_cells; ++i ){
		std::vector<int> distances; // To store distances
		distances.resize(V);
		std::vector<vertex_descriptor> parents(V); // To store parents
		std::cout<<"\t"<<i;
		boost::dijkstra_shortest_paths(g, vertices.at(x[i]*size2+y[i]), boost::predecessor_map(&parents[0]).distance_map(&distances[0]));
		for( int j=i+1; j<num_cells; ++j ){
			if( D(i,j) > 0.00001 ){
				double distance;
				distance=distances[vertices.at(x[j]*size2+y[j])];
				D(i,j)=distance;
				D(j,i)=distance;
			}
		}
		parents.clear();
		distances.clear();
	}
	std::cout<<std::endl;
	//End: Dijkstra
*/

//}

void compute_ip_for_ducts( int *neighbor, double *Distance, int num_cells, double dum_val ){
	IloEnv env;
	try{
		int num_vars = (num_cells*(num_cells+1))/2+num_cells;
		IloNumVar::Type varType   = ILOINT;

		//Declare limits for Edges present -- as binary variables
		IloNumArray edge_min(env,num_vars);
		for( int i=0; i<num_vars; ++i )
			edge_min[i] = 0;
		IloNumArray edge_max(env,num_vars);
		for( int i=0; i<num_vars; ++i )
			edge_max[i] = 1;
		//Add edge costs
		IloNumArray edge_costs(env, num_vars);

		for( int i=0; i<(num_cells-1); ++i ){
			//Add edges
			for( int j=i+1; j<num_cells; ++j )
				EDC(i,j) = D(i,j);
			//Add dummy edges
			EDC(i,num_cells) = dum_val;
			EDC(i,(num_cells+1)) = dum_val;
		}

		//Setup Model and binary variables for Edgepresent
		IloNumVarArray	EP(env);
		IloModel	mod(env);
		EP.clear();
		IloNumVarArray tmp(env, edge_min, edge_max, varType);
		EP.add(tmp);
		tmp.end();

		//Set up minimzation problem
		mod.add(IloMinimize(env, IloScalProd(EP,edge_costs)));

		//Add constraints
		for( int i=0; i<num_cells; ++i ){
			IloExpr expr(env);
			for( int j=0; j<i; ++j )
				expr += E(j,i);
			for( int j=(i+1); j<num_cells; ++j )
				expr += E(i,j);
			expr += E(i,num_cells);
			expr += E(i,(num_cells+1));
			mod.add( expr == 2 );
			expr.end();
		}
		for( int i=0; i<num_cells; ++i ){
			IloExpr expr(env);
			expr += E(i,num_cells);
			mod.add( expr == E(i,(num_cells+1)));
			expr.end();
		}

		//Solve
		IloCplex cplex(mod);
		//cplex.exportModel("Ducts.lp");
		cplex.solve();
		cplex.out() << "solution status = " << cplex.getStatus() << endl;
		cplex.out() << "Objective = " << cplex.getObjValue() << endl;

		for( int i=0; i<num_cells; ++i ){
			for( int j=i+1; j<num_cells; ++j ){
				if( cplex.getValue(E(i,j)) ){
					N(i,j)=1;
					N(j,i)=1;
				}
				else{
					N(i,j)=0;
					N(j,i)=0;
				}
			}
		}
	}
	catch (IloException& ex) {
		std::cerr << "Error: " << ex << endl;
	}   
	catch (...) {
		std::cerr << "Error: " << endl;
	}
}


int main(int argc, char *argv[]){

	//Simple test to check if CPLEX is linked properly
	IloEnv env;
	try{
		IloNumVar::Type varType   = ILOINT;

		IloNumArray	xy_min(env,2); xy_min[0] = 0; xy_min[1] = 0;
		std::stringstream num_feed(stringstream::in | stringstream::out);
		IloNumArray	xy_max(env,2); xy_max[0] = 10; xy_max[1]= 10;
		IloNumArray	weight(env,6);
		weight[0] = 0; weight[1] = 1; weight[2] = 1; weight[3] = -1; weight[4] = -1; weight[5] = -1;
		IloNumArray	lims(env,3); lims[0] = 4; lims[1] = 0; lims[2] = -1.1;
		IloNumArray	xy_cost(env,2); xy_cost[0] = 1; xy_cost[1] = 1;

		IloModel	mod(env);
		IloNumVarArray	xy(env);

		xy.clear();
		IloNumVarArray tmp(env, xy_min, xy_max, varType);
		xy.add(tmp);
		tmp.end();

		int m = 3;
		int n = 2;

		mod.add(IloMinimize(env, IloScalProd(xy,xy_cost)));

		for( int i=0; i<m; ++i ){
			IloExpr expr(env);
			for( int j=0; j<n; ++j )
				expr += xy[j]*weight[((i*n)+j)];
			mod.add( expr <= lims[i] );
			expr.end();
		}

		IloCplex cplex(mod);
		//cplex.exportModel("xy_test_ip.lp");
		cplex.solve();
		cplex.out() << "solution status = " << cplex.getStatus() << endl;
	/*	cplex.out() << "x+y = " << cplex.getObjValue() << endl;
		cplex.out() << "x" << " = " << cplex.getValue(xy[0]) << endl;
		cplex.out() << "y" << " = " << cplex.getValue(xy[1]) << endl;*/
		if( cplex.getValue(xy[1])>0.9 )
			std::cout<<"CPLEX MIP Present and working\n";
		else{
			std::cerr<<"ERROR IN CPLEX MIP\n";
			throw(-1);
		}
	}
	catch (IloException& ex) {
		std::cerr << "Error: " << ex << endl;
	}   
	catch (...) {
		std::cerr << "Error" << endl;
	}

	if( argc < 3 ){
		std::cout<<"Usage "<<argv[0]<<" <input_text_file> <output_text_file> <OPTIONAL_Channel_class> "<<
		"<OPTIONAL_image_for_geodesic_distance> <OPTIONAL_scaling>\n";
		return 0;
	}

	std::ifstream inpf(argv[1]);
	int rows=0;
	std::string dummy;

	std::getline( inpf, dummy ); //Headers
	std::cout<<dummy<<std::endl;

	while( inpf.good() ){
		std::getline( inpf,dummy );
		++rows;
	}

	std::cout<<"Last line: "<<dummy<<"\nNumber of rows:"<<rows<<std::endl;
//	--rows;
//	--rows;
	inpf.close();

	std::ifstream inpfil(argv[1]);
	std::vector<int> values;
	int ex;
//	std::getline( inpfil, dummy ); //Headers
//	std::cout<<"HEADERS:\n"<<dummy<<std::endl;
	while( inpfil>>ex ){
		values.push_back(ex);
	}
	if( floor(values.size()/(double)rows) != values.size()/(double)rows ){
		std::cout<<"Check input file format\n";
		std::cout<<floor(values.size()/(double)rows)<<std::endl<<values.size()<<"\t"<<(double)rows<<std::endl;
	}
	inpfil.close();

	int num_features = (int)floor(values.size()/(double)rows);

	std::cout<<"Num Features:"<<num_features<<std::endl;

	int *x = new int[rows];
	int *y = new int[rows];
	double *radius = new double[rows];
	int *neighbor = new int[(rows*rows)];
	int *cl = new int[rows];
	double *Distance = new double[(rows*rows)];

	for( int i=0; i<(rows*rows); ++i){
		Distance[i]=0;
		neighbor[i]=0;
	}

	for( int i=0; i<rows; ++i ){
		y[i]=(int)values.at(i*num_features);
		x[i]=(int)values.at(i*num_features+1);
		radius[i]=sqrt((values.at(i*num_features+2))/(pi));
//	(int) values.at(i*num_features+2);(int)sqrt((values.at(i*num_features+2))/(2*pi));
		cl[i]=(int)values.at(i*num_features+3);
	}
//	std::cout<<"Last row: "<<x[rows-1]<<"\t"<<y[rows-1]<<"\t"<<radius[rows-1]<<"\t"<<cl[rows-1]<<"\n";
	int num_cells = rows;
	for( int i=0; i<rows; ++i )
		for( int j=0; j<rows; ++j )
			N(i,j)=0;

/*
	std::cout<<"JA\n";
	num_cells = 10;
	int num_vars = (num_cells*(num_cells+1))/2+num_cells;
	std::cout<<"NV:"<<num_vars<<"\n";
		for( int i=0; i<num_cells; ++i ){ //(num_cells-1)
			for( int j=i+1; j<(num_cells+2); ++j ){
				int nnn_vars = num_vars - (num_cells-i+(num_cells-i)*(num_cells-i+1)/2) + j-(i+1) ;
				std::cout<<nnn_vars<<"\t";
			}
		std::cout<<std::endl;
		}
	std::cout<<std::endl;
	return 0;
*/

	cast_rays( radius, x, y, 360, rows, neighbor, Distance );

//	x[0] = 10; y[0]=10; x[1] = 15; y[1]=5; x[2] = 5; y[2]=15; x[3] = 15; y[3]=15; x[4] = 5; y[4]=5;
//	radius[0]=2.5; radius[1]=2.5; radius[2]=2.5; radius[3]=2.5; radius[4]=2.5;
//	cast_rays( radius, x, y, 360, 5, neighbor, distance ); //should be ~60 intersections

	int dirty_class=3;

	int checksumss=0,cc=0;
	for( int i=0; i<rows; ++i )
		if( cl[i]==dirty_class )
			++cc;
		else
			for( int j=i+1; j<rows; ++j )
				if( N(i,j)>5 && cl[i]==cl[j] && cl[j]!=dirty_class ){
					++checksumss; 
				}
				else{
					N(i,j)=0;
					N(j,i)=0;
					D(i,j)=1e8;
					D(j,i)=1e8;
				}

	int *out_list;
	double *out_dist;
	int new_rows=rows;
	out_list = new int[new_rows*new_rows];
	out_dist = new double[new_rows*new_rows];
	double dum_val=0;

	if( argc > 3 ){
		int num_channels_provided = (int)argc-3;
		if( (double)((int)(double(num_channels_provided)/3)) != (double(num_channels_provided)/3) ){
			std::cout<<"usage:"<<argv[0]<<" <input_text_file> <output_text_file> <optinal_channel_class> "<<
			"<optional_image_for_geodesic_distance> <optional_scaling>\n";
			return 0;
		}
		num_channels_provided/=3;
		for(int i=0; i<1; ++i){
			for( int j=0; j<rows; ++j ){
				for( int k=j+1; k<rows; ++k ){
					if( cl[j]==atoi(argv[(3+i*3)]) && cl[k]==atoi(argv[(3+i*3)]) ){
						OLL(j,k)=N(j,k);
						OLL(k,j)=N(j,k);
						OD(j,k)=D(j,k);
						OD(k,j)=D(j,k);
					}
					else{
						OLL(j,k)=0;
						OLL(k,j)=0;
						OD(j,k)=1e8;
						OD(k,j)=1e8;
					}
				}
				//if( cl[j]==atoi(argv[(3+i*3)]) ){ std::cout<<j<<"\t"; }
			}
			//std::cout<<"\n";
		dum_val = change_weights_to_geodesic( atoi(argv[(3+i*3)]), argv[(4+i*3)], atof(argv[(5+i*3)]), out_list, 
											out_dist, num_cells, x, y, cl );//neighbor, Distance,

		int countt=0;
		double summD=0;
		for( int j=0; j<rows; ++j )
			for( int k=j+1; k<rows; ++k )
				if( OLL(j,k) ){
					summD += OD(j,k);
					++countt;
				}

/*		double mean = summD/countt;
		for( int j=0; j<rows; ++j )
			for( int k=j+1; k<rows; ++k )
				if( cl[j]==atoi(argv[(3+i*3)]) && cl[k]==atoi(argv[(3+i*3)]) )
					dum_val +=( (OD(j,k)-mean)*(OD(j,k)-mean)/countt );

		dum_val = sqrt(dum_val);
		dum_val = mean + 2*dum_val;
*/

		std::cout<<"dum_val "<<dum_val<<std::endl;

		compute_ip_for_ducts( out_list, out_dist, num_cells, dum_val);
		}
	}

	typedef unsigned short PixelType;
	const   unsigned int Dimension = 2;
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleUsUsType;
	typedef itk::ImageFileReader< ImageType >  ReaderType;
	typedef itk::ImageFileWriter< ImageType >  WriterType;

	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();
	reader->SetFileName( argv[4] );
	writer->SetFileName( "out.tif" );
	ImageType::Pointer image = ImageType::New();
	image = reader->GetOutput();
	reader->Update();
	ImageType::RegionType region = image->GetLargestPossibleRegion();
	ImageType::IndexType corner1,corner2;
	for( int i=0; i<rows; ++i ){
		for( int j=i+1; j<rows; ++j ){
			if( OLL(i,j) ){
				corner1[0] = x[i];
				corner1[1] = y[i];
				corner2[0] = x[j];
				corner2[1] = y[j];
				itk::LineIterator<ImageType> it(image, corner1, corner2);
				it.GoToBegin();
				while (!it.IsAtEnd()){
					it.Set(255);
					++it;
				}
			}
		}
	}
	writer->SetInput( image );
	writer->Update();

	new_rows=rows-cc;
	int ind1=0,ind2=0,ind=0;
	double distA=0;
	for( int i=0; i<rows; ++i ){
		ind2=0;
		if( cl[i]!=dirty_class ){
			for( int j=0; j<rows; ++j ){
				if(cl[j]!=dirty_class){
					OD(ind1,ind2)=D(i,j);
					OLL(ind1,ind2)=N(i,j);
					if( j>i && N(i,j) ){
						distA += D(i,j);
						++ind;
//						if( ind1==0 && 100 < D(i,j) )
//							std::cout<<i<<(pi*(double)radius[i]*(double)radius[i])
//							<<"\t"<<j<<"\t"<<(pi*(double)radius[j]*(double)radius[j])<<std::endl;
					}
					++ind2;
				}
			}
			++ind1;
		}
		if( cl[i]!=dirty_class && ind2!=new_rows )
			std::cout<<"I am bad. STOP!\n";
			return 0;
	}
	distA /= ind;

	double var=0;
	for( int i=0; i<new_rows; ++i )
		for( int j=i+1; j<new_rows; ++j )
			var +=  (distA-OD(i,j))*(distA-OD(i,j))/ind;
	var /= ind;
	double bsval = exp(var*-1.0);

	std::cout<<"Num Neigh:"<<checksumss<<"\nNum non-class cells: "<<cc<<"\nClass cells:"
	<<new_rows<<" "<<ind1<<" mean="<<distA<<" var="<<var<<" bsval="<<bsval<<std::endl;

	std::ofstream output;
	output.open(argv[2]);
//	output.precision(2);
	output << "param cells := "<<new_rows<<";"<<std::endl;
	output << "param wearethisclose:\t";
	for( int i=0; i<=(new_rows+1); ++i )
		output<<i+1<<"\t";
	output<<":=\n";
	for( int i=0; i<=(new_rows+1); ++i ){
		output<<i+1<<"\t";
		if( i< new_rows ){
			for( int j=0; j<new_rows; ++j)
				if( OLL(i,j) )
					output<<OD(i,j)<<"\t";
				else
					output<<"0\t";
			output<<bsval<<"\t"<<bsval<<std::endl;
		}
		else{
			for( int j=0;j<=(new_rows+1);++j )
				output<<bsval<<"\t";
			output<<bsval;
		}
	}
	output<<" ;";


	output.close();

	delete[] x;
	delete[] y;
	delete[] radius;
	delete[] neighbor; 
	delete[] cl;
	delete[] Distance;
	delete[] out_list;
	delete[] out_dist;

	return 0;
}

