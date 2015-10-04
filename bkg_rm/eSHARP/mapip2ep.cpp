/* MAP IP TO EP
 *
 *	Maps internal points (IP) to external points (EP)
 * 
 *
 *
 * topfer@ualberta.ca 
 *
 * created Feb. 2014
 */

#include <cassert>

#include <vector>
using std::vector;

#include <cstddef>
using std::size_t;

#include <cstdlib>
using std::exit;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ios;

#include <stdio.h>

#include <fstream>
using std::ifstream;
using std::ostream;
using std::ofstream;

#include <iomanip>
using std::setw;
using std::left;

#include <string>
using std::string;

#include <cmath>
using std::sqrt;
using std::fabs;

//--------------------
// Global variables
//
// defaults:
string DEFAULT_SAVEFLDR(".") ;

size_t gridSizeImg[3] ;
size_t numVoxelsImg ;
size_t radiiSphere[3] ;



class subs3
{
public:
	// member data
	int i, j, k ;

	// constructor
	subs3() : i(0), j(0), k(0) {}
	subs3(int i_, int j_, int k_) : i(i_), j(j_), k(k_) {}
	
	friend subs3 operator+(const subs3 &u, const subs3 &v)
	{
		return subs3(u.i + v.i, u.j +v.j, u.k + v.k);
	}

	friend subs3 operator-(const subs3 &u, const subs3 &v)
	{
		return subs3( u.i - v.i, u.j - v.j, u.k - v.k);
	}

	friend ostream& operator<<(ostream& os, const subs3 &v)
	{
		os.precision(4);
		os << setw(20) << v.i << setw(20) << v.j << setw(20) << v.k ;
		return os;
	}
	
};

class displacementsIP2EP
{
public:
	// member data
	// NB: type INT is not necessarily ideal here but used for simplicity.
	vector<int> indexEP, row, col, dx, dy, dz ;

	// constructor
	displacementsIP2EP( ) 
	{ 	
		vector<int> indexEP;
		vector<int> row; 
		vector<int> col; 
		vector<int> dx;
		vector<int> dy;
		vector<int> dz;
	}

	void writevector2file( const vector<int> &v, string filename )
	{
		FILE * fout ;
		cerr << filename.c_str( ) << endl ;
		fout = fopen(filename.c_str(), "wb" )  ;
		
		if (!fout)
		{
			cerr << "Error: could not write to " << filename.c_str( ) << endl ;
			exit(1) ;
		}
		
		fwrite( &v[0], sizeof(int), v.size(), fout ) ;
		fclose( fout ) ;
	}

	void write2file( string saveFldr ) 
	{
		cerr << endl << "Saving..." << endl ;

		writevector2file( dx, saveFldr + string("displacementsX.bin") ) ;
		writevector2file( dy, saveFldr + string("displacementsY.bin") ) ;
		writevector2file( dz, saveFldr + string("displacementsZ.bin") ) ;
		writevector2file( row, saveFldr + string("rows.bin") ) ;
		writevector2file( col, saveFldr + string("columns.bin") ) ;
		writevector2file( indexEP, saveFldr + string("indexEP.bin") ) ;
	}	


};


//--------------------
// Function prototypes
void usage(void) ;

inline size_t subs2ind(size_t gridSize[3], subs3 subs) ;

inline subs3 ind2subs(size_t gridSize[3], size_t ind) ;

inline bool iswithingrid(size_t gridSize[3], subs3 subs) ;

void loadarray(vector<bool> &array, string filenameLoad) ;

vector<subs3> circumscribeneighbourhood( size_t radius[3] ) ; 

displacementsIP2EP assignEP2IP( vector<bool> &maskIP, vector<bool> &maskEP, vector<subs3> &sphericalNeighbourhood ) ;






















int main(int argc, char* argv[])
{

	if ( ( argc < 9 ) || (argc > 10 ) ) 
		usage() ;
	else
	{
		string saveFldr ;
		if(argc == 9)
			saveFldr = DEFAULT_SAVEFLDR ;
		else if (argc == 10)
			saveFldr = string( argv[9] ) ;

	gridSizeImg[0] = atoi( argv[1] ) ;
	gridSizeImg[1] = atoi( argv[2] ) ;
	gridSizeImg[2] = atoi( argv[3] ) ;

	numVoxelsImg = gridSizeImg[0] * gridSizeImg[1] * gridSizeImg[2] ;

	radiiSphere[0] = atoi( argv[4] ) ;
	radiiSphere[1] = atoi( argv[5] ) ;
	radiiSphere[2] = atoi( argv[6] ) ;

	// load internal points (mask/Boolean array)
	vector<bool> maskIP, maskEP ;
	
	loadarray( maskIP, string ( argv[7] )  ) ;
	loadarray( maskEP, string ( argv[8] )  ) ;

	//initialize spherical neighbourhood
	vector<subs3> sphericalNeighbourhood = circumscribeneighbourhood( radiiSphere ) ;

	//assign EP to IP 
	displacementsIP2EP D = assignEP2IP( maskIP, maskEP, sphericalNeighbourhood ) ;
	
	//write to file
	D.write2file( saveFldr ) ;

	}
	return 0;
}























void usage(void)
{
	cerr << endl << "Usage: ./mapip2ep nX nY nZ radX radY radZ filePathIP filePathEP fileSaveFldr" << endl 
                     << endl << "nX/Y/Z" << endl 
                     << endl << "\t" << "grid dimensions (i.e. num freq. encodes/phase encodes/slices)" << endl
                     << endl << "radX/Y/Z" << endl 
                     << endl << "\t" << "radii of sphere (i.e. ellipsoid - radii are in units of voxels to accommodate anisotropic resolutions)" << endl
                     << endl << "filePathIP/EP" << endl
                     << endl << "\t" << "location of 'internal/external point' (IP/EP) .bin files" << endl
                     << endl << "[optional: ]" << endl
                     << endl << "fileSaveFldr" << endl
                     << endl << "\t" << "where to output .bin files" << endl  
                     << endl ;

	exit(1);
}

inline size_t subs2ind(size_t gridSize[3], subs3 subs)
{//subscripts (i,j,k) to linear indices

	return subs.k*(gridSize[0]*gridSize[1]) + subs.j*gridSize[0] + subs.i ;
}

inline subs3 ind2subs(size_t gridSize[3], size_t ind)
{//linear indices to subscripts (i,j,k)

	subs3 subs ;	
		subs.i = ind % gridSize[0] ;
		subs.j = floor( double( ind % (gridSize[0] * gridSize[1]) )/ double( gridSize[0]) ) ;
		subs.k = floor( double(ind) / double(gridSize[0] * gridSize[1]) ) ;    

	return subs ;
}

inline bool iswithingrid(size_t gridSize[3], subs3 subs)
{//check if subscripts (i,j,k) are within [0,gridSize)

	if( ( subs.i < 0 ) || ( subs.i >= gridSize[0] ) ||
		( subs.j < 0 ) || ( subs.j >= gridSize[1] ) ||
		( subs.k < 0 ) || ( subs.k >= gridSize[2] ) ) 
		return false ; 
	else 
		return true ; 

}

void loadarray(vector<bool> &array, string filenameLoad)
{
	cerr  << endl << "Loading file: " << filenameLoad << endl ;

	ifstream fin( filenameLoad.c_str( ), std::ifstream::binary) ;

	if (!fin.is_open())
	{
		cerr << "Error: could not open " << filenameLoad << endl;
		exit(1);
	}

	fin.seekg( 0, fin.end ) ;
	size_t length = fin.tellg( ) ;
	fin.seekg( 0, fin.beg ) ;

	char * buffer = new char [length] ;

	cerr << "Size (bytes): " << setw(10) << length << endl ;
	assert( length == numVoxelsImg ) ;

	// read data as a block:
	fin.read( buffer, length );

	if (fin)
		cout << "Load successful." << endl ;
	else
		cout << "Error: only " << fin.gcount() << " could be read" << endl;

	fin.close();

	for (size_t n = 0; n < numVoxelsImg ; ++n)
	{
		array.push_back( bool( buffer[n] ) );
	}

}

vector<subs3> circumscribeneighbourhood( size_t radius[3] )
{// circumscribe (ellipsoidal) neighbourhood & return vector subscripts within it

	size_t gridSize[3] ;
		gridSize[0] = radius[0]*2 + 1 ;
		gridSize[1] = radius[1]*2 + 1 ;
		gridSize[2] = radius[2]*2 + 1 ;

	vector<subs3> neighbourhood ;

	subs3 midPointSubs ;
		midPointSubs.i = ( gridSize[0] - 1 )/2 ;
		midPointSubs.j = ( gridSize[1] - 1 )/2 ;
		midPointSubs.k = ( gridSize[2] - 1 )/2 ;

	double radius2[3] ;
		radius2[0] = double( radius[0]*radius[0] ) ;
		radius2[1] = double( radius[1]*radius[1] ) ;
		radius2[2] = double( radius[2]*radius[2] ) ;

	for (size_t i = 0; i < gridSize[0]; ++i)
		for (size_t j = 0; j < gridSize[1]; ++j)
			for (size_t k = 0; k < gridSize[2]; ++k)
			{
				if( 1.0 >= (double( (i - midPointSubs.i )*( i - midPointSubs.i ) )/radius2[0]
                                          + double( (j - midPointSubs.j )*( j - midPointSubs.j ) )/radius2[1]
                                          + double( (k - midPointSubs.k )*( k - midPointSubs.k ) )/radius2[2] ) )
				{
					neighbourhood.push_back( subs3(i, j, k ) - midPointSubs ) ;
				}
			}

	return neighbourhood ;

}

displacementsIP2EP assignEP2IP( vector<bool> &maskIP, vector<bool> &maskEP, vector<subs3> &sphericalNeighbourhood )
{

	displacementsIP2EP D ;

	size_t rowEP = 0 ; 

	for( size_t voxel = 0; voxel < numVoxelsImg; ++ voxel )	
		if( maskEP[ voxel ] == true) 
		{
			subs3 subsEP = ind2subs( gridSizeImg, voxel ) ;

			for( size_t pointInNeighbourhood = 0; pointInNeighbourhood < sphericalNeighbourhood.size( ); ++pointInNeighbourhood )
			{
				subs3 subsTmp = sphericalNeighbourhood[ pointInNeighbourhood ] + subsEP ;

				if( iswithingrid( gridSizeImg, subsTmp ) )
				{	
					if( maskIP[ subs2ind( gridSizeImg, subsTmp ) ] )
					{//IP is within neighbourhood of EP
					// subsTmp therefore refers to IP subscripts

						size_t colIP = subs2ind( gridSizeImg, subsTmp ) ;

						subs3 deltaIP2EP = subsEP - subsTmp ;

						D.indexEP.push_back( voxel ) ;	
						D.row.push_back( rowEP ) ;
						D.col.push_back( colIP ) ;
						D.dx.push_back( deltaIP2EP.i ) ;	
						D.dy.push_back( deltaIP2EP.j ) ;	
						D.dz.push_back( deltaIP2EP.k ) ;	
						
					}
				}
			}
		
			++rowEP ; 
		}

	cerr << endl << "num EP/num voxels: " << double(rowEP)/double(numVoxelsImg) << endl ;

	return D ;

}
