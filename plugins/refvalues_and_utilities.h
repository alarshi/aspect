#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include <string>
#include <map>

#include "nr3.h"
//#include "hash_mod.h"
//#include "ran_mod.h"

using namespace std;


namespace UniversalConst
{
	const double PI = 3.141592653589793;
	//const long double PI = 3.14159265358979323846264338328L;
	
	// Conversion from degree to radians
	const double DEG2RAD = PI/180.0;
	
	// Conversion from radians to degrees
	const double RAD2DEG = 180.0/PI;
	
	// Gravitational constant
	// G=(6.67384+-0.00080)x10^11 m3\kg s2
	// [Mohr et al., 2012], Williams et al. 2014
	const double G_gravconst = 6.67384e-11;
}



namespace EarthConst
{
	// Mean radius:
	// Various (2000). David R. Lide, ed. Handbook of Chemistry and Physics (81st ed.). CRC. ISBN 0-8493-0481-4. -- Wikipedia [6]
	const double mean_radius = 6.3710e6;  // meters [6]
	
	
	// Equatorial radius:
	//"Selected Astronomical Constants, 2011". The Astronomical Almanac. Archived from the original on 26 August 2013. Retrieved 25 February 2011. -- Wikipedia [7]
	// World Geodetic System (WGS-84). Available online from National Geospatial-Intelligence Agency. -- Wikipedia [8]
	const double equatorial_radius = 6.3781e6;  // meters [7][8]
	
	
	// Polar radius:
	// Cazenave, Anny (1995). "Geoid, Topography and Distribution of Landforms" (PDF). In Ahrens, Thomas J. Global Earth Physics: A Handbook of Physical Constants. Washington, DC: American Geophysical Union. ISBN 0-87590-851-9. Archived from the original (PDF) on 16 October 2006. Retrieved 3 August 2008. -- Wikipedia [9]
	const double polar_radius = 6.3568e6;  // meters  [9]
	
	
	// Flattening:
	// International Earth Rotation and Reference Systems Service (IERS) Working Group (2004). "General Definitions and Numerical Standards" (PDF). In McCarthy, Dennis D.; Petit, Gerard. IERS Conventions (2003) (PDF). IERS Technical Note No. 32. Frankfurt am Main: Verlag des Bundesamts fur Kartographie und Geodasie. p. 12. ISBN 3-89888-884-3. Retrieved 29 April 2016. -- Wikipedia [10]
	const double flattening = 0.0033528;  // [10]
	//earth_flat = 1.0/298.257222101 # (ETRS89)
	
	
	// Average Atmospheric Pressure on the surface of the Earth [Pa]
	// !!! CITATION NEED !!!
	const double surface_pressure = 101325.0;
	
	
	// PREM Outer Core Radius
	const double PREM_CoreRad = 3480000.000000;

}


namespace MyUtilities
{


	/**
	 * Some templated functions from:
	 * Numerical Recipes. The Art of Scientific Computing. 3rd ed. Press et. al.
	 */
	template<class T>
	inline const T &MAX(const T &a, const T &b)
			{return b > a ? (b) : (a);}
	
	inline float MAX(const double &a, const float &b)
			{return b > a ? (b) : float(a);}
	
	inline float MAX(const float &a, const double &b)
			{return b > a ? float(b) : (a);}
	
	
	
	/**
	* Function to calculate the normalized difference
	*/
	double normdiff (double newval, double refval)
	{
		if (refval == 0.0) return (newval - refval);
		else return ((newval - refval) / refval);
	}
	
	
	
	/**
	 * Split a string by a delimiter into a vector of string parts.
	 */
	std::vector<std::string> split(const std::string& s, char delimiter)
	{
		// Vector of string to save tokens 
		std::vector<std::string> tokens;
		
		std::string token;
		
		// Turn input string into a stringstream
		std::istringstream tokenStream(s);
		
		// Tokenizing w.r.t. delimiter
		while (std::getline(tokenStream, token, delimiter))
		{
			tokens.push_back(token);
		}
		return tokens;
	}
	
	
	
	/**
	 * Simple function that takes in a vector as input and returns
	 * a sum over elements starting from element a to element b,
	 * including the endpoints elements a and b.
	 * 
	 * This function assumes that element a < element b,
	 * and that every element inbetween is summed and counted,
	 * (consecutive).
	 * 
	 * There are two versions of this function to deal with
	 * a simple 1-d vector, and 2-d vector.
	 * 
	 * Specified in the input arguments are:
	 * 		input vector
	 * 		element a
	 * 		element b
	 * 		vector row (for 2-d vectors only)
	 */
	double sumvec_consecutive (std::vector<double> vector_in, int element_a, int element_b)
	{
		double sum = 0;
		
		for (unsigned int i=element_a, count=0; i<(element_b+1); ++i, ++count)
		{
			sum = sum + vector_in[i];
		}
		return sum;
	}
	int sumvec_consecutive (std::vector<int> vector_in, int element_a, int element_b)
	{
		double sum = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			sum = sum + vector_in[i];
		}
		return sum;
	}
	double sumvec_consecutive (std::vector<std::vector<double> > vector_in, int element_a, int element_b, int sumaxis, int fixed_axis_element)
	{
		double sum = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			if (sumaxis==0) sum = sum + vector_in[i][fixed_axis_element];
			else if (sumaxis==1) sum = sum + vector_in[fixed_axis_element][i];
		}
		return sum;
	}
	int sumvec_consecutive (std::vector<std::vector<int> > vector_in, int element_a, int element_b, int sumaxis, int fixed_axis_element)
	{
		double sum = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			if (sumaxis==0) sum = sum + vector_in[i][fixed_axis_element];
			else if (sumaxis==1) sum = sum + vector_in[fixed_axis_element][i];
		}
		return sum;
	}
	
	
	
	/**
	 * Simple function that takes in a vector as input,
	 * sums over elements starting from element a to element b,
	 * including the endpoints elements a and b, and then
	 * returns the arithmitic average.
	 * 
	 * This function assumes that element a < element b,
	 * and that every element inbetween is summed and counted,
	 * (consecutive).
	 * 
	 * There are two versions of this function to deal with
	 * a simple 1-d vector, and 2-d vector.
	 * 
	 * Specified in the input arguments are:
	 * 		input vector
	 * 		element a
	 * 		element b
	 * 		vector row (for 2-d vectors only)
	 */
	double avgvec_consecutive (std::vector<double> vector_in, int element_a, int element_b)
	{
		double sum = 0;
		int count = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			sum = sum + vector_in[i];
			++count;
		}
		return (sum/count);
	}
	int avgvec_consecutive (std::vector<int> vector_in, int element_a, int element_b)
	{
		double sum = 0;
		int count = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			sum = sum + vector_in[i];
			++count;
		}
		return (sum/count);
	}
	double avgvec_consecutive (std::vector<std::vector<double> > vector_in, int element_a, int element_b, int sumaxis, int fixed_axis_element)
	{
		double sum = 0;
		int count = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			if (sumaxis==0) sum = sum + vector_in[i][fixed_axis_element];
			else if (sumaxis==1) sum = sum + vector_in[fixed_axis_element][i];
			++count;
		}
		return (sum/count);
	}
	double avgvec_consecutive (std::vector<std::vector<int> > vector_in, int element_a, int element_b, int sumaxis, int fixed_axis_element)
	{
		double sum = 0;
		int count = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			if (sumaxis==0) sum = sum + vector_in[i][fixed_axis_element];
			else if (sumaxis==1) sum = sum + vector_in[fixed_axis_element][i];
			++count;
		}
		return (sum/count);
	}
	
	
	
	/**
	 * Function to create and return a linearly
	 * increasing, evenly spaced vector.
	 * 
	 * Inputs:
	 * 			a = "start point" of type T
	 * 			b = "end point" of type T
	 * 			N = "number of increments" of type int
	 * 
	 * Output:	std::vector of type T
	 *
	 * Modified from https://gist.github.com/lorenzoriano/5414671
	 *	20181005
	 */
	template <typename T>
	std::vector<T> linspace(T a, T b, size_t N)
	{
		T h = (b - a) / static_cast<T>(N-1);
		std::vector<T> xs(N);
		typename std::vector<T>::iterator x;
		T val;
		for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) *x = val;
		
		
		return xs;
	}
	
	
	
	/**
	 * Function to create and return a linearly
	 * increasing, evenly spaced vector on log scale.
	 * 
	 * Inputs:
	 * 			a = "start point" of type T
	 * 			b = "end point" of type T
	 * 			N = "number of increments" of type int
	 * 			base = "the log base value" of type T
	 * 					the default value is 10.0
	 * 
	 * Output:	std::vector of type T
	 *
	 * Modified from https://gist.github.com/lorenzoriano/5414671
	 *	20181005
	 */
	template <typename T>
	std::vector<T> logspace(T a, T b, size_t N, T base=10)
	{
		T h = (b - a) / static_cast<T>(N-1);
		std::vector<T> xs(N);
		typename std::vector<T>::iterator x;
		T val;
		for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) *x = pow(base,val);
		
		
		return xs;
	}
	
	
	
	/**
	 * Function to calculate the mass of a spherical
	 * shell with a given density, where density is
	 * constant within the shell.
	 */
	template <class T>
	struct CalcCumulativeMass
	{
		double density;
		
		CalcCumulativeMass (double den) : density(den) {}
		double operator()(double radius)
		{
			
			return (4.0*UniversalConst::PI*density*radius*radius); //Needs to be integrated over radius
			
		}
	};
	
	
	
	/**
	 * !!!	UNFINISHED	!!!
	 * Funciton to calculate gravity at radii
	 * with given density.
	 * 
	 * !!!	MAY NOT NEED THIS FUNCTION AT ALL !!!
	 * !!!	THIS FUNCTION DOES NOTHING AT THE MOMENT	!!!
	 */
	double compute_gravity(std::vector <double> density, std::vector <double> radius, int irad, const double Gconst)
	{
		// Version based on Mark Panning's
		// Matlab code for Mars.
		
		const double PI = 3.141592653589793;
		return 4.0*PI;
		//return (4.0*PI*density*radius*radius);
		
		
		//gravity_grav[1:] = Gconst*mass_at_radii[1:]/(radii_grav[1:]**2.0)
		//gravity_grav[0] = 0.0
	}
	
	
	
	/**
	 * Funciton to calculate pressure at radii
	 * with given density and gravity.
	 * 
	 */
	template <class T>
	struct ComputePressure
	{
		std::vector<double> density, gravity, depth_increasing;
		
		ComputePressure (std::vector<double> den, std::vector<double> grav, std::vector<double> dep) : density(den), gravity(grav), depth_increasing(dep) {}
		double operator()(unsigned int irad, unsigned int idep1, unsigned int idep2)
		{
			return (density[irad]*gravity[irad] * (depth_increasing[idep2] - depth_increasing[idep1]) );
		}
	};
	
	
	
	/**
	 * Compute the average density of a series of spherical
	 * shells from the surface down to a current depth.
	 */
	template <class T>
	struct ComputeAverageDensity
	{
		std::vector<double> density, radius;
		
		ComputeAverageDensity (std::vector<double> den, std::vector<double> rad) : density(den), radius(rad) {}
		double operator()(int r1, int r2)
		{
			double avgdensity = 0.0; //= ( radius[r1] - radius[r1-1] ) * density[r1];
			
			for (unsigned int i=r1; i<r2; ++i)
			{
				avgdensity += ( radius[i+1] - radius[i] ) * density[i];
			}
			/*
			for (unsigned int i=r1+1; i==r2; --i)
			{
				avgdensity += ( radius[i] - radius[i+1] ) * density[i];
				//avgdensity = density[i];
				//cout << i << "\n";
			}
			*/
			
			return ( avgdensity / (radius[r2] - radius[r1]) );
		}
	};
	
	
	
	/**
	 * Function to compute the average density of several material components.
	 * 
	 * The general scheme is to determine the total mass and total volume
	 * of the combined material components. This is done by determining
	 * the fraction of volume and mass for each component, and summing
	 * them up. The average density is then simply:
	 * 
	 * 	AVGDENSITY = TOTMASS / TOTVOL
	 * 
	 * This particular function is valid for non-interactive solids.
	 * In the case of interactive solids and liquids, it is possible
	 * that the combined volume is not the same as the sum of
	 * individual volumes. Thus, more information would be needed to
	 * calculate the average density.
	 */
	double avg_density_multmatcomp (std::vector<double> densities, std::vector<double> volumes, std::vector<double> fractions)
	{
		double total_mass, total_volume;
		double avg_density;
		
		total_mass = total_volume = 0.0;
		
		for (unsigned int i=0; i<densities.size(); ++i)
		{
			total_mass = total_mass + fractions[i]*densities[i]*volumes[i];
			total_volume = total_volume + fractions[i]*volumes[i];
		}
		avg_density = total_mass / total_volume;
		
		return avg_density;
	}
	
	
	
	/**
	 * Function to calculate and populate a vector of depth increments
	 * from a vector of radii. The depth vector is populated such that
	 * the first element is the surface, and succesive elements move
	 * deeper. In other words, depth increases with increasing radii
	 * elements.
	 * 
	 * It is assumed the input radii are increasing from some depth
	 * onward toward the surface with increasing elements, and that
	 * the final radius element is the surface.
	 */
	std::vector<double> make_increasing_depth(std::vector<double> radius_in)
	{
		std::vector<double> depth_increasing(radius_in.size());
		for (unsigned int irad=0, idepth=radius_in.size()-1; irad<radius_in.size(); ++irad, --idepth)
		{
			depth_increasing[idepth] = EarthConst::mean_radius - radius_in[irad];
		}
		return depth_increasing;
	}
	
	
	
	/**
         * The following set of routines solve a set of
	 * linear equations (A %cdot x = b)
	 * using Lower Upper (LU) decomposition.
         */
	//  --------------------------------------------------
	
	struct LUdcmp
	{
	  // Object for solving linear equations A*x=b using LU decomposition,  and related functions.
	  
	  int n;
	  MatDoub lu;                                              // Stores the decomposition
	  VecInt indx;                                             // Stores the permutation
	  double d;                                                // Used by det
	  
	  LUdcmp(MatDoub_I &a);                                    // Constructor. Argument is the matrix A
	  
	  void solve(VecDoub_I &b,  VecDoub_O &x);                 // Solve for a single right-hand side
	  void solve(MatDoub_I &b,  MatDoub_O &x);                 // Solve for multiple right-hand sides
	  
	  void inverse(MatDoub_O &ainv);                           // Calculate matrix inverse A^(-1)
	  
	  double det();                                            // Return determinant of A
	  
	  void mprove(VecDoub_I &b,  VecDoub_IO &x);               // Discussed in SS2.5
	  MatDoub_I &aref;                                         // Used only by mprove
	  
	};
	
	
	LUdcmp::LUdcmp(MatDoub_I &a) : n(a.nrows()), lu(a),  aref(a),  indx(n)
	{
	  //  Given a matrix a[0...n-1][0...n-1], this routine replaces it by the LU decomposition of a
	  //  rowwise permutation of itself. a is input. On output, it is arranged as in equation (2.3.14)
	  //  above; indx[0...n-1] is an output vector that records the row permutation effected by the
	  //  partial pivoting; d is output as +-1 depending on whether the number of row interchanges
	  //  was even or odd,  respectively. This routine is used in combination with "solve" to solve linear
	  //  equations or invert a matrix.
	  
	  const double TINY = 1.0e-40;                             // A small number
	  int i, imax, j, k;
	  double big, temp;
	  std::vector<double> vv(n);                               // vv stores the implicit scaling of each row
	  
	  d = 1.0;                                                 // No row interchanges yet
	  
	  for (i = 0; i<n; i++)                                    // Loop over rows to get the implicit scaling information
	  {
	    big = 0.0;
	    
	    for (j=0; j<n; j++) if ( (temp = abs(lu[i][j])) > big ) big = temp;
	    
	    if (big == 0.0) throw("Singular matrix in LUdcmp");
	    
	    // No nonzero largest element
	    vv[i] = 1.0/big;                                       // Save the scaling
	  }
	  
	  for (k = 0; k<n; k++)                                    // This is the outermost kij loop
	  {
	    big = 0.0;                                             // Initialize for the search for largest pivot element
	    
	    imax = k;
	    
	    for (i = k; i<n; i++)
	    {
	      temp = vv[i] * abs(lu[i][k]);
	      if (temp>big)
	      {                                                    // Is the figure of merit for the pivot better than the best so far?
		big = temp;
		imax = i;
	      }
	    }
	    
	    if (k != imax)
	    {							// Do we need to interchange rows?
	      for (j=0; j<n; j++)
	      {							// Yes, do so...
		temp = lu[imax][j];
		lu[imax][j] = lu[k][j];
		lu[k][j] = temp;
	      }
	      
	      d = -d;						// ...and change the parity of d.
	      vv[imax] = vv[k];					// Also interchange the scale factor
	    }
	    
	    indx[k] = imax;
	    
	    if (lu[k][k] == 0.0) lu[k][k]=TINY;
	    // If the pivot element is zero, the matrix is singular (at least to the precision of
	    // the algorithm). For some applications on singular matrices, it is desirable to
	    // substitute TINY for zero.
	    
	    for (i=k+1; i<n; i++)
	    {
	      temp = lu[i][k] /= lu[k][k];			// Divide by the pivot element
	      
	      for (j=k+1; j<n; j++) lu[i][j] -= temp*lu[k][j];	// Innermost loop: reduce remaining submatrix
	    }
	  }
	  
	}
	
	
	void LUdcmp::solve(VecDoub_I &b, VecDoub_O &x)
	{
	  // Solves the set of n linear equations A*x=b using the stored LU decomposition of A.
	  // b[0..n-1] is input as the right-hand side vector b, while x returns the solution vector x;
	  // b and x may reference the same vector, in which case the solution overwrites the input.
	  // This routine takes into account the possibility that b will begin with many zero elements,
	  // so it is efficient for use in matrix inversion.
	  
	  int i, ii=0, ip, j;
	  double sum;
	  
	  if (b.size() != n || x.size() != n) throw("LUdcmp::solve bad sizes");
	  
	  for (i=0; i<n; i++) x[i] = b[i];
	  
	  // When ii is set to a positive value, it will become the index of the
	  // first nonvanishing elements of b. We now do the forward substitution,
	  // equation (2.3.6). The only new wrinkle is to unscramble the permutation
	  // as we go.
	  for (i=0; i<n; i++)
	  {
	    ip = indx[i];
	    sum = x[ip];
	    x[ip] = x[i];
	    
	    if (ii != 0) for (j=ii-1; j<i; j++) sum -= lu[i][j]*x[j];
	    
	    // A nonzero element was encountered, so from now on we
	    // will have to do the sums in the loop above.
	    else if (sum != 0.0) ii = i+1;
	    
	    x[i] = sum;
	  }
	  
	  for (i=n-1; i>=0; i--)
	  { 				// Now we do the backsubstitution, equation (2.3.7)
	    sum = x[i];
	    for (j=i+1; j<n; j++) sum -= lu[i][j]*x[j];
	    x[i] = sum/lu[i][i];	// Store a component of the solution vector X
	  }				// All done!
	  
	}
	
	
	void LUdcmp::solve(MatDoub_I &b, MatDoub_O &x)
	{
	  // Solves m sets of n linear equations A*X=B using the stored LU decompostion of A.
	  // The matrix b[0..n-1][0..m-1] inputs the right-hand sides, while x[0..n-1][0..m-1]
	  // returns the solution A^(-1)*B. b and x may reference the same matrix, in which
	  // case the solution overwrites the input.
	  
	  int i, j, m=b.ncols();
	  
	  if ( b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols() ) throw("LUdcmp::solve bad sizes");
	  
	  VecDoub xx(n); //std::vector<double> xx(n);
	  
	  for (j=0; j<m; j++)
	  {				// Copy and solve each column in turn
	    for (i=0; i<n; i++) xx[i] = b[i][j];
	    solve(xx,xx);
	    for (i=0; i<n; i++) x[i][j] = xx[i];
	  }
	}
	
	
	void LUdcmp::inverse(MatDoub_O &ainv)
	{
	  // Using the stored LU decomposition, return in ainv the matrix inverse A^(-1)
	  
	  int i, j;
	  ainv.resize(n,n);
	  
	  for (i=0; i<n; i++)
	  {
	    for (j=0; j<n; j++) ainv[i][j] = 0.0;
	    ainv[i][i] = 1.0;
	  }
	  
	  solve(ainv, ainv);
	}
	
	
	double LUdcmp::det()
	{
	  // Using the stored LU decomposition, return the determinant of the matrix A
	  
	  double dd = d;
	  for (int i=0; i<n; i++) dd *= lu[i][i];
	  return dd;
	}
	
	
	void LUdcmp::mprove(VecDoub_I &b, VecDoub_IO &x)
	{
	  // Improves a solution vector x[0..n-1] of the linear set of equations A*x=b.
	  // The vectors b[0..n-1] and x[0..n-1] are input. On output, x[0..n-1] is
	  // modified, to an improved set of values.
	  
	  int i, j;
	  VecDoub r(n); //std::vector<double> r(n);
	  
	  // Calculate the right-hand side, accumulating the residual in higher precision
	  for (i=0; i<n; i++)
	  {
	    Ldoub sdp = -b[i];
	    
	    for (j=0; j<n; j++) sdp += (Ldoub)aref[i][j] * (Ldoub)x[j];
	    r[i] = sdp;
	  }
	  
	  // Solve for the error term, and subtract it from the old solution
	  solve(r,r);
	  for (i=0; i<n; i++) x[i] -= r[i];
	}
	//  --------------------------------------------------
}



namespace PTUtilities
{
	
	// Declare a class with the necessary functions to complete the work
	//template
	// !!! THESE FUNCTIONS CAN LIKELY BE MERGED A BIT !!!
	class PTLookups
	{
	public:
		std::vector<std::vector<double> > load_T_profile(std::string Tfilename);
		std::vector<std::vector<double> > load_entire_orig_T_file(std::string Tfilename);
		std::vector<std::vector<double> > load_P_profile(std::string Pfilename);
		std::vector<std::vector<double> > load_entire_P_file(std::string Pfilename);
	};
	
	
	
	std::vector< std::vector<double> >
	PTLookups::load_T_profile(std::string Tfilename)
	{
		
		std::ifstream Tfile;
		//std::string Tfilename = "../geotherms/geotherm_Data_Fig1b_in_GlisovicEtAl_2015GRL_mod.dat";
		std::vector<std::vector <double> > data;
		
		Tfile.open(Tfilename);
		
		
		// Skip header
		std::string temp;
		for (unsigned int ihead=0; ihead<2; ++ihead)
		{
			std::getline(Tfile,temp); //throw these lines away
			//std::cout << temp << "\n";
		}
		
		
		
		// Read in lines until the right coordinate is reached
		// then read in the information and exit
		
		double dummy;
		double radius_in, G_in, M_in, N_in;
		//double z_orig_in, z_mod_in, radius_in, G_in, M_in, N_in;
		
		
		// Collect coordinates and data
		while (std::getline(Tfile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data Point (coordinate and data: Vp,Vs,Rho)
			// z(m)_orig z(m)_mod radius(m) G(K) M(K) N(K)
			//line >> z_orig_in >> z_mod_in >> radius_in >> G_in >> M_in >> N_in;
			line >> dummy >> dummy >> radius_in >> G_in; // >> M_in >> N_in;
			
			values.emplace_back (radius_in);
			values.emplace_back (G_in);
			
			data.emplace_back (values);
			
		}
		
		
		
		Tfile.close();
		//std::cout << "END OF TLOOKUP...temperature is "<< temperature;
		
		return data;
	}
	
	
	
	/**
	 * This function reads in the original temperture profile
	 * file from Fig1b_in_GlisovicEtAl_2015GRL. The profile is
	 * shifted in rdepth by 3.0km deeper than the surface. So
	 * a depth of 0.0km in the file is actually 3.0km deep.
	 * The profile also starts at the surface and proceeds to
	 * the CMB.
	 * This function shifts the profile 3.0km deeper and
	 * reverses the profile so that radius=0.0km is at the top.
	 * The radius is also converted from km to meters, standard
	 * SI units.
	 * This function returns all the columns found in the file.
	 */
	std::vector< std::vector<double> >
	PTLookups::load_entire_orig_T_file(std::string Tfilename)
	{
		
		std::ifstream Tfile;
		//std::string Tfilename = "../geotherms/geotherm_Data_Fig1b_in_GlisovicEtAl_2015GRL.dat";
		
		Tfile.open(Tfilename);
		
		
		// Skip header
		std::string temp;
		for (unsigned int ihead=0; ihead<1; ++ihead)
		{
			std::getline(Tfile,temp); //throw these lines away
			//std::cout << temp << "\n";
		}
		
		
		
		// Read in lines until the right coordinate is reached
		// then read in the information and exit
		
		std::vector<std::vector <double> > data_in, data_out;
		double depth_in, G_in, M_in, N_in;
		
		
		// Collect coordinates and data
		while (std::getline(Tfile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			//depth(km) G(K) M(K) N(K)
			line >> depth_in >> G_in >> M_in >> N_in;
			
			values.emplace_back (depth_in);
			values.emplace_back (G_in);
			values.emplace_back (M_in);
			values.emplace_back (N_in);
			
			data_in.emplace_back (values);
			
		}
		
		
		// Reverse and shift the profile, and also
		// convert the depth to units of meters.
		for (unsigned int k=data_in.size(); k>0; --k)
		{
			unsigned int i = k-1; // The vectors are minus 1 of the count
			
			std::vector<double> values;
			
			// Original depth, then shifted depth, then equivalent Earth radius
			values.emplace_back (data_in[i][0]*1000.0); // convert to meters
			values.emplace_back ((data_in[i][0]+3.0)*1000.0); // convert to meters and shift
			values.emplace_back ( EarthConst::mean_radius - ( (data_in[i][0]+3.0)*1000.0) ); // covert to meters, shift, and convert to radius
			
			// Temperature values
			values.emplace_back (data_in[i][1]);  // G(K)
			values.emplace_back (data_in[i][2]);  // M(K)
			values.emplace_back (data_in[i][3]);  // N(K)
			
			data_out.emplace_back (values);
		}
		
		
		
		Tfile.close();
		
		return data_out;
	}
	
	
	
	std::vector< std::vector<double> >
	PTLookups::load_P_profile(std::string Pfilename)
	{
		
		std::ifstream Pfile;
		//std::string Pfilename = "../PRESSURE_PROFILES/PressureProfile_STW105.dat";
		Pfile.open(Pfilename);
		std::vector<std::vector <double> > data;
		
		// Skip the header lines
		std::string temp;
		for (unsigned int ihead=0; ihead<1; ++ihead)
		{
			std::getline (Pfile, temp); //throw these lines away
		}
		
		
		
		// Read in lines until the right coordinate is reached
		// then read in the information and exit
		
		double dummy;
		double radius_in, pressure_in;
		//double radius_in, depth_increasing_in, density_in, mass_at_radii_in, gravity_in, pressure_in;
		
		
		while (std::getline(Pfile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data Point (coordinate and data: Vp,Vs,Rho)
			// radius(m) depth(m) density(kg/m3) cumulativemass(kg) gravity(kg*m/s2) pressure(Pa)
			//line >> radius_in >> depth_increasing_in >> density_in >> mass_at_radii_in >> gravity_in >> pressure_in;
			line >> radius_in >> dummy >> dummy >> dummy >> dummy >> pressure_in;
			
			values.emplace_back (radius_in);
			values.emplace_back (pressure_in);
			
			data.emplace_back (values);
			
		}
		
		
		
		Pfile.close();
		
		return data;
	}
	
	
	
	std::vector< std::vector<double> >
	PTLookups::load_entire_P_file(std::string Pfilename)
	{
		
		std::ifstream Pfile;
		//std::string Pfilename = "../PRESSURE_PROFILES/PressureProfile_STW105.dat";
		Pfile.open(Pfilename);
		std::vector<std::vector <double> > data;
		
		// Skip the header lines
		std::string temp;
		for (unsigned int ihead=0; ihead<1; ++ihead)
		{
			std::getline (Pfile, temp); //throw these lines away
		}
		
		
		
		// Read in lines until the right coordinate is reached
		// then read in the information and exit
		// radius(m) depth(m) density(kg/m3) cumulativemass(kg) gravity(kg*m/s2) pressure(Pa)
		double radius, depth, density, cumulativemass, gravity, pressure;
		
		
		while (std::getline(Pfile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data Points
			// radius(m) depth(m) density(kg/m3) cumulativemass(kg) gravity(kg*m/s2) pressure(Pa)
			line >> radius >> depth >> density >> cumulativemass >> gravity >> pressure;
			
			values.emplace_back (radius);
			values.emplace_back (depth);
			values.emplace_back (density);
			values.emplace_back (cumulativemass);
			values.emplace_back (gravity);
			values.emplace_back (pressure);
			
			data.emplace_back (values);
			
		}
		
		
		
		Pfile.close();
		
		return data;
	}
}



namespace SearchUtilities
{
  // A place to put file search utilities
  
  
  /**
   * Functions to define Points and Boxes
   */
  //  --------------------------------------------------
  template <int dim>
  struct Point
  {
    // Simple structure to represent a point in dim dimensions
    
    double x[dim];		// the coordinates
    
    // Copy constructor
    Point(const Point &p)
    {
      for (int i=0; i<dim; i++) x[i] = p.x[i];
    }
    
    // Assignment operator
    Point& operator= (const Point &p)
    {
      for (int i=0; i<dim; i++) x[i] = p.x[i];
      return *this;
    }
    
    bool operator== (const Point &p) const
    {
      for (int i=0; i<dim; i++) if (x[i] != p.x[i]) return false;
      return true;
    }
    
    Point(double x0=0.0, double x1=0.0, double x2=0.0)
    {
      // Constructor by coordinate values. Arguments beyond the required
      // number are not used and can be omitted.
      
      x[0] = x0;
      if (dim > 1) x[1] = x1;
      if (dim > 2) x[2] = x2;
      if (dim > 3) throw("Point not implemented for dim > 3");
    }
  };
  
  
  template <int dim>
  double dist(const Point<dim> &p, const Point<dim> &q)
  {
    // Returns the distance between two points in dim dimensions
    
    double dd = 0.0;
    
    for (int j=0; j<dim; j++) dd += (q.x[j]-p.x[j])*(q.x[j]-p.x[j]);
    
    return sqrt(dd);
  }
  
  
  template <int dim>
  struct Box
  {
    // Structure to represent a Cartesian box in dim dimensions
    
    // Diagonally opposite corners (min of all coodinates and max of all coordinates) are stored as two points
    Point<dim> lo, hi;
    
    Box() {}
    Box(const Point<dim> &mylo, const Point<dim> &myhi) : lo(mylo), hi(myhi) {}
  };
  
  
  template <int dim>
  double dist(const Box<dim> &b, const Point<dim> &p)
  {
    // If point p lies outside box b, the distance to the nearest point on b is returned.
    // If p is inside b or on its surface, zero is returned.
    
    double dd = 0;
    
    for (int i=0; i<dim; i++)
    {
      if (p.x[i]<b.lo.x[i]) dd += (p.x[i]-b.lo.x[i])*(p.x[i]-b.lo.x[i]);
      if (p.x[i]>b.hi.x[i]) dd += (p.x[i]-b.lo.x[i])*(p.x[i]-b.hi.x[i]);
    }
    return sqrt(dd);
  }
  
  
  template <int dim>
  struct Boxnode : Box<dim>
  {
    // Node in a binary tree of boxes containing points. See text for details.
    
    int mom, dau1, dau2, ptlo, pthi;
    
    Boxnode() {}
    Boxnode(Point<dim> mylo, Point<dim> myhi, int mymom, int myd1, int myd2, int myptlo, int mypthi) : Box<dim>(mylo, myhi), mom(mymom), dau1(myd1), dau2(myd2), ptlo(myptlo), pthi(mypthi) {}
  };
  //  --------------------------------------------------
  
  
  
  /**
   * Setup the KDTree search method:
   * 	search a set of n-dimensional coordinates
   * 	with m-dimensional data to locate specific
   * 	points and nearest neighbors
   */
  //  --------------------------------------------------
  
  template <int dim>
  struct KDtree
  {
    // Structure for implementing a KD tree
    
    static const double BIG;			// Size of the root box, value set below
    int nboxes, npts;				// Number of boxes, number of points
    std::vector< Point<dim> > &ptss;		// Reference to the vector of points in the KD tree
    Boxnode<dim> *boxes;			// the array of Boxnodes that form the tree
    std::vector<int> ptindx, rptindx;		// Index of points (see text), and reverse index
    double *coords;				// Point coordinates rearranged contiguously
    
    KDtree(std::vector< Point<dim> > &pts);	// Constructor
    ~KDtree() {delete [] boxes;}
    
    // Next, utility functions for use after the tree is constructed. See below.
    double disti(int jpt, int kpt);
    int locate(Point<dim> pt);
    int locate(int jpt);
    
    // Next, applications that use the KD tree. See text
    int nearest(Point<dim> pt);
    void nnearest(int jpt, int *nn, double *dn, int n);
    static void sift_down(double *heap, int *ndx, int nn);	// Used by nnearest
    int locatenear(Point<dim> pt, double r, int *list, int nmax);
  };
  
  
  template <int dim>
  const double KDtree<dim>::BIG(1.0e99);
  
  int selecti(const int k, int *indx, int n, double *arr)
  {
    // Permutes indx[0..n-1] to make arr[indx[0..k-1]] <= arr[indx[k+1..n-1]].
    // The array arr is not modified. See comments in the routine select.
    
    int i, ia, ir, j, l, mid;
    double a;
    
    l = 0;
    ir = n-1;
    
    for (;;)
    {
      if (ir <= l+1)
      {
	if (ir == l+1 && arr[indx[ir]] < arr[indx[l]]) SWAP(indx[l], indx[ir]);
	return indx[k];
      }
      else
      {
	mid = (l+ir) >> 1;
	SWAP(indx[mid], indx[l+1]);
	
	if (arr[indx[l]] > arr[indx[ir]]) SWAP(indx[l], indx[ir]);
	if (arr[indx[l+1]] > arr[indx[ir]]) SWAP(indx[l+1], indx[ir]);
	if (arr[indx[l]] > arr[indx[l+1]]) SWAP(indx[l], indx[l+1]);
	
	i = l+1;
	j = ir;
	ia = indx[l+1];
	a = arr[ia];
	
	for (;;)
	{
	  do i++; while (arr[indx[i]] < a);
	  do j--; while (arr[indx[j]] > a);
	  
	  if (j<i) break;
	  SWAP(indx[i], indx[j]);
	}
	
	indx[l+1] = indx[j];
	indx[j] = ia;
	if (j >= k) ir = j-1;
	if (j <= k) l = i;
      } 
    }
  }
  
  
  template <int dim>
  KDtree<dim>::KDtree(std::vector< Point<dim> > &pts) : ptss(pts), npts(pts.size()), ptindx(npts), rptindx(npts)
  {
    // Construct a KDtree from a vector of points
    
    int ntmp, m, k, kk, j, nowtask, jbox, np, tmom, tdim, ptlo, pthi;
    int *hp;
    double *cp;
    int taskmom[50], taskdim[50];		// Enough stack for 2^50 points!
    
    // Initialize the index of points
    for (k=0; k<npts; k++) ptindx[k] = k;
    
    // Calculate the number of boxes and allocate memory for them
    m = 1;
    for (ntmp = npts; ntmp; ntmp >>= 1) m <<= 1;
    
    nboxes = 2*npts - (m >> 1);
    if (m < nboxes) nboxes = m;
    nboxes--;
    
    boxes = new Boxnode<dim>[nboxes];
    
    // Copy the point coordinates into a contiguous array
    coords = new double[dim*npts];
    for (j=0, kk=0; j<dim; j++, kk += npts)
    {
      for (k=0; k<npts; k++) coords[kk+k] = pts[k].x[j];
    }
    
    // Initialize the root box and put it on the task list for subdivision
    Point<dim> lo(-BIG, -BIG, -BIG), hi(BIG, BIG, BIG);	// Syntax OK for 2-D too
    
    boxes[0] = Boxnode<dim>(lo, hi, 0, 0, 0, 0, npts-1);
    
    jbox = 0;
    taskmom[1] = 0;		// Which box
    taskdim[1] = 0;		// Which dimension
    nowtask = 1;
    
    // Main loop over pending tasks
    while (nowtask)
    {
      tmom = taskmom[nowtask];
      tdim = taskdim[nowtask--];
      ptlo = boxes[tmom].ptlo;
      pthi = boxes[tmom].pthi;
      hp = &ptindx[ptlo];		// Points to left end of subdivision
      cp = &coords[tdim*npts];		// Points to coordinate list for current dim
      np = pthi - ptlo + 1;		// Number of points in the subdivision
      kk = (np-1)/2;			// Index of last point on left (boundary point)
      
      // Here is where all the work is done
      (void) selecti(kk, hp, np, cp);
      
      // Now create the daughters and push them onto the task list if they need further subdividing
      hi = boxes[tmom].hi;
      lo = boxes[tmom].lo;
      hi.x[tdim] = lo.x[tdim] = coords[tdim*npts + hp[kk]];
      
      boxes[++jbox] = Boxnode<dim>(boxes[tmom].lo, hi, tmom, 0, 0, ptlo, ptlo+kk);
      boxes[++jbox] = Boxnode<dim>(lo, boxes[tmom].hi, tmom, 0, 0, ptlo+kk+1, pthi);
      boxes[tmom].dau1 = jbox-1;
      boxes[tmom].dau2 = jbox;
      
      if (kk > 1)
      {
	taskmom[++nowtask] = jbox-1;
	taskdim[nowtask] = (tdim+1) % dim;
      }
      if (np - kk > 3)
      {
	taskmom[++nowtask] = jbox;
	taskdim[nowtask] = (tdim+1) % dim;
      }
    }
    
    for (j=0; j<npts; j++) rptindx[ptindx[j]] = j;	// Create reverse index
    
    delete [] coords;					// Don't need them anymore
    
  }
  
  
  template <int dim>
  double KDtree<dim>::disti(int jpt, int kpt)
  {
    // Returns the distance between two points in the KDtree given their indices
    // in the array of points, but returns a large value if the points are identical
    
    if (jpt == kpt) return BIG;
    else return dist(ptss[jpt], ptss[kpt]);
  }
  
  
  template <int dim>
  int KDtree<dim>::locate( Point<dim> pt )
  {
    // Given an arbitrary point pt, return the index of which KDtree box it is in
    
    int nb, d1, jdim;
    
    nb = jdim = 0;		// Start with the root box
    
    // As far as possible down the tree
    while (boxes[nb].dau1)
    {
      d1 = boxes[nb].dau1;
      
      if (pt.x[jdim] <= boxes[d1].hi.x[jdim]) nb = d1;
      else nb = boxes[nb].dau2;
      
      // Increment the dimension cyclically
      jdim = ++jdim % dim;
    }
    
    return nb;
  }
  
  
  template <int dim>
  int KDtree<dim>::locate( int jpt )
  {
    // Given the index of a point in the KDtree, return the index of which box it is in
    
    int nb, d1, jh;
    
    // The reverse index tells where the point lies in the index of points
    jh = rptindx[jpt];
    
    nb = 0;
    
    while (boxes[nb].dau1)
    {
      d1 = boxes[nb].dau1;
      
      if (jh <= boxes[d1].pthi) nb = d1;
      else nb = boxes[nb].dau2;
    }
    
    return nb;
  }
  
  
  template <int dim>
  int KDtree<dim>::nearest( Point<dim> pt)
  {
    // Given an arbitrary location pt, return the index of the nearest point in the KDtree
    
    int i, k, nrst, ntask;
    int task[50];		// Stack for boxes waiting to be opened
    double dnrst = BIG, d;
    
    // First stage, we find the nearest KDtree point in same box as pt
    k = locate(pt);		// Which box is pt in?
    
    // Find nearest
    for (i=boxes[k].ptlo; i<=boxes[k].pthi; i++)
    {
      d = dist(ptss[ptindx[i]], pt);
      
      if (d < dnrst)
      {
	nrst = ptindx[i];
	dnrst = d;
      }
    }
    
    // Second stage, we traverse the tree opening only possibly better boxes
    task[1] = 0;
    ntask = 1;
    while (ntask)
    {
      k = task[ntask--];
      
      // Distance to closest point in box
      if (dist(boxes[k], pt) < dnrst)
      {
	// If not an end node, put on task list
	if (boxes[k].dau1)
	{
	  task[++ntask] = boxes[k].dau1;
	  task[++ntask] = boxes[k].dau2;
	}
	else
	{
	  // Check the 1 or 2 points in the box
	  
	  for (i=boxes[k].ptlo; i<=boxes[k].pthi; i++)
	  {
	    d = dist(ptss[ptindx[i]], pt);
	    
	    if (d < dnrst)
	    {
	      nrst = ptindx[i];
	      dnrst = d;
	    }
	  }
	}
      }
    }
    
    return nrst;
  }
  
  
  template <int dim>
  void KDtree<dim>::nnearest(int jpt, int *nn, double *dn, int n)
  {
    // Given the index jpt of a point in a KDtree, return a list nn[0..n-1] of indices of the
    // n points in the tree nearest to point jpt, and a list dn[0..n-1] of their distances.
    
    int i, k, ntask, kp;
    int task[50];		// Stack for boxes to be opened
    double d;
    
    if (n > npts-1) throw("too many neighbors requested");
    
    for (i=0; i<n; i++) dn[i] = BIG;
    
    // Find smallest mother box with enough points to initialize the heap
    kp = boxes[locate(jpt)].mom;
    
    while (boxes[kp].pthi - boxes[kp].ptlo < n) kp = boxes[kp].mom;
    
    // Examine its points and save the n closest
    for (i=boxes[kp].ptlo; i<=boxes[kp].pthi; i++)
    {
      if (jpt == ptindx[i]) continue;
      d = disti(ptindx[i], jpt);
      
      if (d < dn[0])
      {
	dn[0] = d;
	nn[0] = ptindx[i];
	if (n>1) sift_down(dn, nn, n);		// Maintain the heap structure
      }
    }
    
    // Now we traverse the tree opening only possibly better boxes
    task[1] = 0;
    ntask = 1;
    while (ntask)
    {
      k = task[ntask--];
      
      if (k == kp) continue;		// Don't redo the box used to initialize
      
      if (dist(boxes[k], ptss[jpt]) < dn[0])
      {
	// If not an end node, put it on task list
	if (boxes[k].dau1)
	{
	  task[++ntask] = boxes[k].dau1;
	  task[++ntask] = boxes[k].dau2;
	}
	// Check the 1 or 2 points in the box
	else
	{
	  for (i=boxes[k].ptlo; i<=boxes[k].pthi; i++)
	  {
	    d = disti(ptindx[i], jpt);
	    if (d < dn[0])
	    {
	      dn[0] = d;
	      nn[0] = ptindx[i];
	      if (n>1) sift_down(dn, nn, n);	// Maintain the heap
	    }
	  }
	}
      }
    }
    
    return;
  }
  
  
  template <int dim>
  void KDtree<dim>::sift_down(double *heap, int *ndx, int nn)
  {
    // Fix heap[0..nn-1], whose first element (only) may be wrongly filed. Make a
    // corresponding permutation in ndx[0..nn-1]. The algorithm is identical
    // to that used by sift_down in hpsort.
    
    int n = nn -1;
    int j, jold, ia;
    double a;
    
    a = heap[0];
    ia = ndx[0];
    jold = 0;
    
    j = 1;
    while (j <= n)
    {
      if (j < n && heap[j] < heap[j+1]) j++;
      if (a >= heap[j]) break;
      
      heap[jold] = heap[j];
      ndx[jold] = ndx[j];
      
      jold = j;
      j = 2*j + 1;
    }
    
    heap[jold] = a;
    ndx[jold] = ia;
  }
  
  
  template <int dim>
  int KDtree<dim>::locatenear( Point<dim> pt, double r, int *list, int nmax)
  {
    // Given a point pt and radius r, returns a value nret such that list[0..nret-1] is a list of
    // all KDtree points within a radius r of pt, up to a user-specified maximum of nmax points.
    
    int k, i, nb, nbold, nret, ntask, jdim, d1, d2;
    int task[50];
    
    nb = jdim = nret = 0;
    
    if (r < 0.0) throw("radius must be nonnegative");
    
    // Find the smallest box that contains the "ball" of radius r
    while (boxes[nb].dau1)
    {
      nbold = nb;
      d1 = boxes[nb].dau1;
      d2 = boxes[nb].dau2;
      
      // Only need to check the dimension that divides the daughters
      if (pt.x[jdim] + r <= boxes[d1].hi.x[jdim]) nb = d1;
      else if (pt.x[jdim] - r >= boxes[d2].lo.x[jdim]) nb = d2;
      
      jdim = ++jdim % dim;
      
      // Neither daughter encloses the ball
      if (nb == nbold) break;
    }
    
    // Now traverse the tree below the starting box only as needed
    task[1] = nb;
    ntask = 1;
    
    while (ntask)
    {
      k = task[ntask--];
      
      if (dist(boxes[k], pt) > r) continue;		// Box and ball are disjoint
      
      // Expand box further when possible
      if (boxes[k].dau1)
      {
	task[++ntask] = boxes[k].dau1;
	task[++ntask] = boxes[k].dau2;
      }
      // Otherwise process points in the box
      else
      {
	for (i=boxes[k].ptlo; i<=boxes[k].pthi; i++)
	{
	  if (dist(ptss[ptindx[i]], pt) <= r && nret < nmax) list[nret++] = ptindx[i];
	  if (nret == nmax) return nmax;		// Not enough space!
	}
      }
    }
    
    return nret;
  }
  //  --------------------------------------------------
}


namespace InterpolationRoutines
{
	// A place to put general interpolation stuff
	
	struct Base_interp
	{
		Int n, mm, jsav, cor, dj;
		const Doub *xx, *yy;
		Base_interp(VecDoub_I &x, const Doub *y, Int m)
			: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
			dj = MyUtilities::MAX(1,(int)pow((Doub)n,0.25));
		}

		Doub interp(Doub x) {
			Int jlo = cor ? hunt(x) : locate(x);
			return rawinterp(jlo,x);
		}

		Int locate(const Doub x);
		Int hunt(const Doub x);
		
		Doub virtual rawinterp(Int jlo, Doub x) = 0;

	};
	Int Base_interp::locate(const Doub x)
	{
		Int ju,jm,jl;
		if (n < 2 || mm < 2 || mm > n) throw("locate size error");
		Bool ascnd=(xx[n-1] >= xx[0]);
		jl=0;
		ju=n-1;
		while (ju-jl > 1) {
			jm = (ju+jl) >> 1;
			if (x >= xx[jm] == ascnd)
				jl=jm;
			else
				ju=jm;
		}
		cor = abs(jl-jsav) > dj ? 0 : 1;
		jsav = jl;
		return MyUtilities::MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
	}
	Int Base_interp::hunt(const Doub x)
	{
		Int jl=jsav, jm, ju, inc=1;
		if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
		Bool ascnd=(xx[n-1] >= xx[0]);
		if (jl < 0 || jl > n-1) {
			jl=0;
			ju=n-1;
		} else {
			if (x >= xx[jl] == ascnd) {
				for (;;) {
					ju = jl + inc;
					if (ju >= n-1) { ju = n-1; break;}
					else if (x < xx[ju] == ascnd) break;
					else {
						jl = ju;
						inc += inc;
					}
				}
			} else {
				ju = jl;
				for (;;) {
					jl = jl - inc;
					if (jl <= 0) { jl = 0; break;}
					else if (x >= xx[jl] == ascnd) break;
					else {
						ju = jl;
						inc += inc;
					}
				}
			}
		}
		while (ju-jl > 1) {
			jm = (ju+jl) >> 1;
			if (x >= xx[jm] == ascnd)
				jl=jm;
			else
				ju=jm;
		}
		cor = abs(jl-jsav) > dj ? 0 : 1;
		jsav = jl;
		return MyUtilities::MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
	}
	
	
	
	struct Linear_interp : Base_interp
	{
		Linear_interp(VecDoub_I &xv, VecDoub_I &yv)
			: Base_interp(xv,&yv[0],2)  {}
		Doub rawinterp(Int j, Doub x) {
			if (xx[j]==xx[j+1]) return yy[j];
			else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
		}
	};
	
	
	
	struct Poly_interp : Base_interp
	{
		Doub dy;
		Poly_interp(VecDoub_I &xv, VecDoub_I &yv, Int m)
			: Base_interp(xv,&yv[0],m), dy(0.) {}
		Doub rawinterp(Int jl, Doub x);
	};

	Doub Poly_interp::rawinterp(Int jl, Doub x)
	{
		Int i,m,ns=0;
		Doub y,den,dif,dift,ho,hp,w;
		const Doub *xa = &xx[jl], *ya = &yy[jl];
		VecDoub c(mm),d(mm);
		dif=abs(x-xa[0]);
		for (i=0;i<mm;i++) {
			if ((dift=abs(x-xa[i])) < dif) {
				ns=i;
				dif=dift;
			}
			c[i]=ya[i];
			d[i]=ya[i];
		}
		y=ya[ns--];
		for (m=1;m<mm;m++) {
			for (i=0;i<mm-m;i++) {
				ho=xa[i]-x;
				hp=xa[i+m]-x;
				w=c[i+1]-d[i];
				if ((den=ho-hp) == 0.0) throw("Poly_interp error");
				den=w/den;
				d[i]=hp*den;
				c[i]=ho*den;
			}
			y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
		}
		return y;
	}

	
	
	struct Spline_interp : Base_interp
	{
		VecDoub y2;
		
		Spline_interp(VecDoub_I &xv, VecDoub_I &yv, Doub yp1=1.e99, Doub ypn=1.e99)
		: Base_interp(xv,&yv[0],2), y2(xv.size())
		{sety2(&xv[0],&yv[0],yp1,ypn);}

		Spline_interp(VecDoub_I &xv, const Doub *yv, Doub yp1=1.e99, Doub ypn=1.e99)
		: Base_interp(xv,yv,2), y2(xv.size())
		{sety2(&xv[0],yv,yp1,ypn);}

		void sety2(const Doub *xv, const Doub *yv, Doub yp1, Doub ypn);
		Doub rawinterp(Int jl, Doub xv);
	};
	void Spline_interp::sety2(const Doub *xv, const Doub *yv, Doub yp1, Doub ypn)
	{
		Int i,k;
		Doub p,qn,sig,un;
		VecDoub u(n-1);
		if (yp1 > 0.99e99)
			y2[0]=u[0]=0.0;
		else {
			y2[0] = -0.5;
			u[0]=(3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
		}
		for (i=1;i<n-1;i++) {
			sig=(xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
			p=sig*y2[i-1]+2.0;
			y2[i]=(sig-1.0)/p;
			u[i]=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
			u[i]=(6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p;
		}
		if (ypn > 0.99e99)
			qn=un=0.0;
		else {
			qn=0.5;
			un=(3.0/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]));
		}
		y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
		for (k=n-2;k>=0;k--)
			y2[k]=y2[k]*y2[k+1]+u[k];
	}
	Doub Spline_interp::rawinterp(Int jl, Doub x)
	{
		Int klo=jl,khi=jl+1;
		Doub y,h,b,a;
		h=xx[khi]-xx[klo];
		if (h == 0.0) throw("Bad input to routine splint");
		a=(xx[khi]-x)/h;
		b=(x-xx[klo])/h;
		y=a*yy[klo]+b*yy[khi]+((a*a*a-a)*y2[klo]
			+(b*b*b-b)*y2[khi])*(h*h)/6.0;
		return y;
	}
	
	
	
	struct Spline2D_interp {
		Int m,n;
		const MatDoub &y;
		const VecDoub &x1;
		VecDoub yv;
		NRvector<Spline_interp*> srp;

		Spline2D_interp(VecDoub_I &x1v, VecDoub_I &x2v, MatDoub_I &ym)
			: m(x1v.size()), n(x2v.size()), y(ym), yv(m), x1(x1v), srp(m) {
			for (Int i=0;i<m;i++) srp[i] = new Spline_interp(x2v,&y[i][0]);
		}

		~Spline2D_interp(){
			for (Int i=0;i<m;i++) delete srp[i];
		}

		Doub interp(Doub x1p, Doub x2p) {
			for (Int i=0;i<m;i++) yv[i] = (*srp[i]).interp(x2p);
			Spline_interp scol(x1,yv);
			return scol.interp(x1p);
		}
	};
}


namespace RadialBasisFunctions
{
	// Radial Basis Functions -------------------------------
	
	// Virtual Base Class
	struct RBF_fn
	{
		// Abstract base class template for any particular radial basis function See specific examples below.
		virtual double rbf(double r) = 0;
	};
	
	// Interpolation routine
	struct RBF_interp
	{
		
		// Object for radial basis function interpolation using n points in dim dimensions.
		// Call constructor once, than interp as many times as desired.
		
		int dim, n;
		const MatDoub &pts; // !!! For now this uses MatDoub from nr3.h
		const VecDoub &vals; //const std::vector<double> &vals;
		VecDoub w; //std::vector<double> w;
		RBF_fn &fn;
		Bool norm;
		
		RBF_interp(MatDoub_I &ptss, VecDoub_I &valss, RBF_fn &func, Bool nrbf=false) : dim(ptss.ncols()), n(ptss.nrows()), pts(ptss), vals(valss), w(n), fn(func), norm(nrbf)
		{
			// Constructor. The n X dim matrix ptss inputs the data points, the vector valss the function values.
			// func contains the chosen radial basis function, derived from the class RBF_fn.
			// The default value of nrbf gives RBF interpolation; set it to 1 for NRBF.
			int i,  j;
			double sum;
			MatDoub rbf(n, n);
			VecDoub rhs(n); //std::vector<double> rhs(n);
			
			for (i=0; i<n; i++)                                      // Fill the matrix phi(|ri - rj|) and the rhs vector
			{
			sum = 0.0;
			
			for (j = 0; j<n; j++)
			{
				sum += (rbf[i][j] = fn.rbf(rad(&pts[i][0],  &pts[j][0]) ) );
			}
			
			if (norm) rhs[i] = sum*vals[i];
			else rhs[i] = vals[i];
			}
			
			MyUtilities::LUdcmp lu(rbf);                            // Solve the set of linear equations
			lu.solve(rhs,  w);
		
		}
		
		double interp(VecDoub_I &pt)
		{
			//  Return the interpolated function value at a dim-dimensional point pt.
			
			double fval,  sum = 0.0,  sumw = 0;
			
			if (pt.size() !=  dim) throw("RBF_interp bad pt size");
			
			for ( int i = 0; i<n; i++)
			{
			fval = fn.rbf(rad(&pt[0],  &pts[i][0]) );
			sumw += w[i] *fval;
			sum += fval;
			}
			
			return norm  ?  sumw/sum : sumw;
		}
		
		
		double rad(const double *p1,  const double *p2)
		{
			// Euclidean distance
			
			double sum = 0.0;
			
			for (int i = 0; i<dim; i++) sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
			return sqrt(sum);
		}
		
	};
	
	
	struct RBF_multiquadratic :  RBF_fn
	{
		// Instantiate this and send to RBF_interp to get multiquadratic interpolation
		
		double r02;
		
		RBF_multiquadratic(double scale = 1.0) : r02(scale*scale) {}
		
		// Constructor argument is the scale factor. See text.
		double rbf(double r) { return sqrt( (r*r) + r02 );}
	};
	
	
	struct RBF_thinplate :  RBF_fn
	{
		// Same as above,  but for thin-plate spline.
		
		double r0;
		
		RBF_thinplate(double scale = 1.0) : r0(scale) {}
		
		double rbf(double r) { return r <= 0.0  ?  0.0 : ( (r*r) * log(r/r0) ); }
	};
	
	
	struct RBF_gauss :  RBF_fn
	{
		//  Same as above,  but for Gaussian
		
		double r0;
		
		RBF_gauss(double scale = 1.0) : r0(scale) {}
		
		double rbf(double r) { return exp( -0.5*(r/r0)*(r/r0) ); }
	};
	
	
	struct RBF_inversemultiquadratic :  RBF_fn
	{
		//  Same as above,  but for inverse multiquadratic
		
		double r02;
		
		RBF_inversemultiquadratic(double scale = 1.0) : r02(scale*scale) {}
		
		double rbf(double r) { return 1.0/sqrt( (r*r) + r02 ); }
	};
	
	
	
	/**
	 * Function to retrieve the pre-determined r0 values needed
	 * for the RBF routines.
	 * 
	 * The returned values should have come from a bootstrapping
	 * process to determine the optimum r0 value (highest
	 * accuracy interpolation).
	 * 
	 * When using the sister file to the original PerpleX
	 * data file, these r0 values should be in the same row
	 * order as the original Perplex data file.
	 * !!! ADD COLUMN CHECK !!!
	 */
  	//template <int querydim, int datadim>
	std::vector<std::vector <double> >
	load_r0_data (std::string filein) { //(const MPI_Comm &comm) {
		
		
		std::vector<std::vector <double> > data;
				
		std::string temp;
		// Read data from disk and distribute among processes
		std::ifstream in;//(Utilities::read_and_distribute_file_content(filename,comm));
		in.open(filein);
		
		
			
		// Skip header line
		std::getline(in,temp); // throw these lines away
		
		
		while (getline(in,temp)) //(in.good() && !in.eof())   // ????
		{
			std::istringstream line(temp);
			
			// read the two components of the coordinates
			double P_in, T_in;
			line >> P_in >> T_in;
			
			std::vector<double> values;
			for (unsigned int i=0; i<12; ++i) // !!! NEED TO CHANGE THIS TO FIND END OF LINE INSTEAD OF HARD-CODED NUMBER OF COLUMNS !!!
			{
				double value;
				line >> value;

				values.emplace_back(value);
			}
			data.emplace_back (values);
		}
		
		in.close();
		
		return data;
	}
	
	
	
	struct Shep_interp {
		Int dim, n;
		const MatDoub &pts;
		const VecDoub &vals;
		Doub pneg;

		Shep_interp(MatDoub_I &ptss, VecDoub_I &valss, Doub p=2.)
		: dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss),
		vals(valss), pneg(-p) {}

		Doub interp(VecDoub_I &pt) {
			Doub r, w, sum=0., sumw=0.;
			if (pt.size() != dim) throw("RBF_interp bad pt size");
			for (Int i=0;i<n;i++) {
				if ((r=rad(&pt[0],&pts[i][0])) == 0.) return vals[i];
				sum += (w = pow(r,pneg));
				sumw += w*vals[i];
			}
			return sumw/sum;
		}

		Doub rad(const Doub *p1, const Doub *p2) {
			Doub sum = 0.;
			for (Int i=0;i<dim;i++) sum += SQR(p1[i]-p2[i]);
			return sqrt(sum);
		}
	};
  
}



namespace IntegrationRoutines
{
	/**
	 * Polynominal interpolation routines from:
	 * Numerical Recipes. The Art of Scientific Computing. 3rd ed. Press et. al.
	 * 
	 * Each algorithm is marked with the section number it is found in.
	 * Comments are completely preserved.
	 */
	
	// Sec 3.1
	struct Base_interp
	{
		// Abstract base class used by all interpolation routines in this
		// chapter. Only the routine interp is called directly by the user.
		
		int n, mm, jsav, cor, dj;
		const double *xx, *yy;
		
		// Constructor: Set up for interpolating on a table of
		// x's and y's of length m. Normally called by a derived class,
		// not by the user.
		Base_interp (const std::vector <double> &x, const double *y, int m) : n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y)
		{
			dj = MyUtilities::MAX(1, (int)pow((double)n, 0.25 ) );
		}
		
		// Given a value x, return an interpolated value, using data pointed
		// to by xx and yy.
		double interp(double x)
		{
			int jlo = cor ? hunt(x) : locate(x);
			return rawinterp(jlo, x);
		}
		
		int locate(const double x);  // see definitions below
		int hunt(const double x);
		
		// Derived classes provide this as the actual interpolation method.
		double virtual rawinterp(int jlo, double x) = 0;
		
	};
	
	// Sec 3.2
	struct Poly_interp : Base_interp
	{
		// Polynomial interpolation object. Construct with x and y vectors,
		// and the number M of points to be used locally (polynomial order
		// plus one), then call interp for interpolated values.
		
		double dy;
		
		Poly_interp (const std::vector <double> &xv, const std::vector <double> &yv, int m) : Base_interp (xv, &yv[0], m), dy(0.) {}
		double rawinterp (int jl, double x);
	};
	
	// Sec 3.2
	double Poly_interp::rawinterp(int jl, double x)
	{
		// Given a value x, and using pointers to data xx and yy, this routine
		// returns an interpolated value y, and stores an error estimate dy.
		// The returned value is obtained by mm-point polynomial interpolation
		// on the subrange xx[jl..jl+mm-1].
		
		int i, m, ns=0;
		double y, den, dif, dift, ho, hp, w;
		const double *xa = &xx[jl], *ya = &yy[jl];
		
		std::vector <double> c(mm), d(mm);
		
		dif = abs(x - xa[0]);
		
		// Here we find the index ns of the closest table entry,
		// and initialize the tableau of c's and d's.
		for (i=0; i<mm; i++)
		{
			if ( (dift = abs(x-xa[i]) ) < dif )
			{
				ns = i;
				dif = dift;
			}
			
			c[i] = ya[i];
			d[i] = ya[i];
		}
		
		// This is the initial approximation to y.
		y = ya[ns--];
		
		// For each column of the tableau,
		// we loop over the current c's and d's and update them.
		for (m=1; m<mm; m++)
		{
			for (i=0; i<mm-m; i++)
			{
				ho = xa[i] - x;
				hp = xa[i+m] - x;
				w = c[i+1] - d[i];
				
				// This error can occur only if two input xa's are
				// (to within roundoff) identical.
				if ( (den = ho - hp) == 0.0 ) throw("Poly_interp error");
				
				den = w / den;
				
				// Here the c's and d's are updated.
				d[i] = hp * den;
				c[i] = ho * den;
			}
			
			// After each column in the tableau is completed, we decide
			// which correction, c or d, we want to add to our accumulating
			// value of y, i.e., which path to take through the tableau
			// --- forking up or down. We do this in such a way as to take
			// the most "straight line" route through the tableau to its
			// apex, updating ns accordingly to keep track of where we are.
			// This route keeps the partial approximations centered
			// (insofar as possible) on the target x. The last dy added
			// is thus the error indication.
			y += ( dy = (2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]) );
		}
		
		return y;
	}
	
	
	
	
	/**
	 * Integration routines from:
	 * Numerical Recipes. The Art of Scientific Computing. 3rd ed. Press et. al.
	 * 
	 * Each algorithm is marked with the section number it is found in.
	 * Comments are completely preserved.
	 */
	
	// Sec 4.2
	struct Quadrature{
		// Abstract base class for elementary quadrature algorithms
		int n;  // Current level of refinement
		virtual double next() = 0;
		// Returns the value of the integral at the nth stage of refinement.
		// The function next() must be defined in the derived class.
	};
	
	
	
	// Sec 4.2
	template<class T>
	struct Trapzd : Quadrature
	{
		// Routine implementing the extended trapezoidal rule.
		double a, b, s;  // Limits of integration and current value of integral
		T &func;
		Trapzd() {};
		Trapzd(T &funcc, const double aa, const double bb) : func(funcc), a(aa), b(bb) {n=0;}
		// The constructor takes as inputs func, the function or functor
		// to be integrated between limits a and b, also input.
		double next()
		{
			// Returns the nth stage of refinement of the extended trapezoidal
			// rule. On the first call (n=1), the routine returns the crudest
			// estimate of INTERGRAL{f(x)dx} from a to h. Subsequent calls
			// set n=2,3,... and improve the accuracy by adding
			// 2^(n-2) additional interior points.
			
			double x, tnm, sum, del;
			int it, j;
			n++;
			
			if (n == 1)
			{
				return ( s = 0.5 * (b-a) * (func(a) + func(b)) );
			}
			else
			{
				for (it=1,j=1; j<n-1; j++) it <<= 1;
				tnm = it;
				del = (b-a) / tnm;  // This is the spacing of the points to be added
				x = a + 0.5 * del;
				for (sum=0.0,j=0; j<it; j++,x+=del) sum += func(x);
				s = 0.5 * (s + (b-a) * sum / tnm);  // This replaces s by its refined value
				return s;
			}
		}
	};
	
	
	// Sec 4.2
	template<class T>
	struct Trapzd_int : Quadrature
	{
		// Routine implementing the extended trapezoidal rule.
		double s;  // Limits of integration and current value of integral
		int a, b;
		T &func;
		Trapzd_int() {};
		Trapzd_int(T &funcc, const int aa, const int bb) : func(funcc), a(aa), b(bb) {n=0;}
		// The constructor takes as inputs func, the function or functor
		// to be integrated between limits a and b, also input.
		double next()
		{
			// Returns the nth stage of refinement of the extended trapezoidal
			// rule. On the first call (n=1), the routine returns the crudest
			// estimate of INTERGRAL{f(x)dx} from a to h. Subsequent calls
			// set n=2,3,... and improve the accuracy by adding
			// 2^(n-2) additional interior points.
			
			int x;
			double tnm, sum, del;
			int it, j;
			n++;
			cout << n << "\n";
			if (n == 1)
			{
				return ( s = 0.5 * (b-a) * (func(a) + func(b)) );
			}
			else
			{
				for (it=1,j=1; j<n-1; j++) it <<= 1;
				tnm = it;
				del = (b-a) / tnm;  // This is the spacing of the points to be added
				x = int(a + 0.5 * del);
				for (sum=0.0,j=0; j<it; j++,x+=int(del)) sum += func(x);
				s = 0.5 * (s + (b-a) * sum / tnm);  // This replaces s by its refined value
				cout << s << "\n";
				return s;
			}
		}
	};
	
	
	
	// Sec 4.3
	template <class T>
	double qromb (T &func, double a, double b, const double eps=1.0e-10)
	{
		// Returns the integral of the function or functor func from a to b.
		// Integration is performed by Romberg's method of order 2K, where
		// , e.g., K=2 is Simpson's rule.
		
		
		// Here EPS is the fractional accuracy desired, as determined by the
		// extrapolation error estimate; JMAX limits the total number of steps;
		// K is the number of points used in the extrapolation.
		const int JMAX=20, JMAXP=JMAX+1, K=5;
		
		
		// These store the successive trapezoidal approximations and their
		// relative stepsizes.
		std::vector <double> s(JMAX), h(JMAXP);
		
		Poly_interp polint(h, s, K);
		
		h[0] = 1.0;
		
		Trapzd<T> t(func, a, b);
		
		for (int j=1; j<=JMAX; j++)
		{
			s[j-1] = t.next();
			
			if (j >= K)
			{
				double ss = polint.rawinterp (j-K, 0.0);
				
				if ( abs(polint.dy) <= eps*abs(ss) ) return ss;
			}
			
			h[j] = 0.25 * h[j-1];
			// This is a key step: The factor is 0.25 even though the stepsize
			// is decreased by only 0.5. This makes the extrapolation a
			// polynomial in h^2 as allowed by NumericalRecipes_eq (4.2.1),
			// not just a polynomial in h.
		}
		
		throw("Too many steps in routine qromb");
	}
	
	
	// Sec 4.3 -- Overloaded Class -- Modified to accept integers a and b
	template <class T>
	double qromb (T &func, int a, int b, const double eps=1.0e-10)
	{
		// Returns the integral of the function or functor func from a to b.
		// Integration is performed by Romberg's method of order 2K, where
		// , e.g., K=2 is Simpson's rule.
		
		
		// Here EPS is the fractional accuracy desired, as determined by the
		// extrapolation error estimate; JMAX limits the total number of steps;
		// K is the number of points used in the extrapolation.
		const int JMAX=20, JMAXP=JMAX+1, K=5;
		
		
		// These store the successive trapezoidal approximations and their
		// relative stepsizes.
		std::vector <double> s(JMAX), h(JMAXP);
		
		Poly_interp polint(h, s, K);
		
		h[0] = 1.0;
		
		Trapzd_int<T> t(func, a, b);
		
		for (int j=1; j<=JMAX; j++)
		{
			s[j-1] = t.next();
			
			if (j >= K)
			{
				double ss = polint.rawinterp (j-K, 0.0);
				
				if ( abs(polint.dy) <= eps*abs(ss) ) return ss;
			}
			
			h[j] = 0.25 * h[j-1];
			// This is a key step: The factor is 0.25 even though the stepsize
			// is decreased by only 0.5. This makes the extrapolation a
			// polynomial in h^2 as allowed by NumericalRecipes_eq (4.2.1),
			// not just a polynomial in h.
		}
		
		throw("Too many steps in routine qromb");
	}
}



namespace PerpleXUtilities
{
	template <int querydim, int datadim>
	class PerpleXLookup
	{
		/** Things and steps this function should do:
		 * ------------------------------------------
		 * -- Load the necessary n-dim coordinate and m-dim data from a table
		 * -- Read in a target point that needs to be assessed
		 * -- Read in the interpolation type to find the corresponding
		 * -- Setup the KDTree from dealii
		 *    -- Load only coordinate points into a KDTree
		 *    -- Use the KDTree to return the n-number of closest points to a target point using:
		 *    -- KDTree< dim >::get_closest_points (const Point< dim > &target, const unsigned int n_points) const
		 *    data that goes with the target point
		 * -- Return data values for use
		 */
		
		
		public:
		std::vector<double> get_data (std::string filein, double T, double P); //(const MPI_Comm &comm); // Load the data table into memory
		//void closest_points; // Setup the KDTree and return the specified number of closest points
		//void data_interpolation; // Interpolate the data values based on data at the closest points
		std::string get_header (std::string filein);
		std::vector< SearchUtilities::Point<querydim> > load_coord (std::string filein);
		std::vector<std::vector <double> > load_data (std::string filein);
	};
	
	
	
	template <int querydim, int datadim>
	std::vector<double>
	PerpleXLookup<querydim,datadim>::get_data (std::string filein, double T, double P) { //(const MPI_Comm &comm) {
		
		
		std::vector<double> coords(querydim); // !!!	CREATE A POINT CLASS TO DEAL WITH COORDINATES	!!!
		std::vector<double> data;
		double dummy, rho, vp, vs;
		double mindist;
		int counter;
		
		std::vector<std::string> keys = {"T(K)","P(bar)","rho,kg/m3","alpha,1/K","cp,J/K/kg","vp,km/s","vs,km/s","h,J/kg"};
		//P(bar) T(K) rho,kg/m3 alpha,1/K beta,1/bar Ks,bar Gs,bar v0,km/s vp,km/s vs,km/s s,J/K/kg h,J/kg cp,J/K/kg V,J/bar/mol
		
		
		std::string temp;
		// Read data from disk and distribute among processes
		std::ifstream in;//(Utilities::read_and_distribute_file_content(filename,comm));
		in.open(filein);
		
		
			
		// Skip header
		for (unsigned int i=0; i<12; ++i)
		{
			std::getline(in,temp); // throw these lines away
		}
		
		double P_found, T_found;  // !!!	TEMPORARY STATEMENT FOR CHECKING CODE	!!!
		
		counter = 0;
		// Read in the full table
		while (getline(in,temp)) //(in.good() && !in.eof())   // ????
		{
			std::istringstream line(temp);
			
			// read the two components of the coordinates
			double P_in, T_in;
			line >> P_in >> T_in;
			
			// !!!	TEMPORARY CONVERSION -- BAR -> Pa	!!!
			P_in = P_in * 1.0e5;
			
			
			// Find the distance
			double Pdiff = P_in - P;
			double Tdiff = T_in - T;
			if (counter == 0)
			{
				mindist = std::sqrt( Pdiff*Pdiff + Tdiff*Tdiff );
				line >> rho >> dummy >> dummy >> dummy >> dummy >> dummy >> vp >> vs;
				//P_found = P_in;
				//T_found = T_in;
			}
			
			double currentdist = std::sqrt( Pdiff*Pdiff + Tdiff*Tdiff );
			
			if (currentdist < mindist)
			{
				mindist = currentdist;
				line >> rho >> dummy >> dummy >> dummy >> dummy >> dummy >> vp >> vs;
				//P_found = P_in;
				//T_found = T_in;
			}
			
			++counter;
			
		}
		
		// Populate the data vector
		data.emplace_back (rho);
		data.emplace_back (vp*1.0e3);
		data.emplace_back (vs*1.0e3);
		
		//std::cout << " T " << T << "  T_found " << T_found << "  P " << P << "  P_found " << P_found;
		//std::cout << "  RHO " << rho << "  VP " << vp << "  VS " << vs << "\n";
		
		in.close();
		
		return data;
	}
	
	
	template <int querydim, int datadim>
	std::string
	PerpleXLookup<querydim,datadim>::get_header (std::string filein) { //(const MPI_Comm &comm) {
		
				
		std::string temp;
		// Read data from disk and distribute among processes
		std::ifstream in;//(Utilities::read_and_distribute_file_content(filename,comm));
		in.open(filein);
		
		
		// Skip header
		for (unsigned int i=0; i<12; ++i)
		{
			std::getline(in,temp); // throw these lines away
		}
        
        
		// Read the header that describes the data columns
		std::getline(in,temp);
		
		in.close();
		
		
		return temp;
	}
	
	
	
	template <int querydim, int datadim>
	std::vector< SearchUtilities::Point<querydim> >
	PerpleXLookup<querydim,datadim>::load_coord (std::string filein) { //(const MPI_Comm &comm) {
		
		
		std::vector< SearchUtilities::Point<querydim> > coords;
		//std::vector<double> coords(querydim); // !!!	CREATE A POINT CLASS TO DEAL WITH COORDINATES	!!!
				
		std::string temp;
		// Read data from disk and distribute among processes
		std::ifstream in;//(Utilities::read_and_distribute_file_content(filename,comm));
		in.open(filein);
		
		
			
		// Skip header
		for (unsigned int i=0; i<13; ++i)
		{
			std::getline(in,temp); // throw these lines away
		}
		
		
		while (getline(in,temp)) //(in.good() && !in.eof())   // ????
		{
			std::istringstream line(temp);
			SearchUtilities::Point<querydim> coordinates;
			// read the two components of the coordinates
			double P_in, T_in;
			line >> P_in >> T_in;
			
			//std::get<0>(coordinates) = P_in;
			coordinates.x[0] = P_in;
			coordinates.x[1] = T_in;
			//std::cout << "P " << P_in << " T " << T_in << " Coord " << coordinates.x[0] << "," << coordinates.x[1] << "\n"; 
			////std::vector<double> coordinates(querydim);
			//line >> coordinates;
			coords.emplace_back(coordinates);	
			//coords.emplace_back(coordinates.x[0],coordinates.x[1]);
		}
		
		in.close();
		
		//std::cout << "A point: " << coords[2].x[0] << "," << coords[2].x[1] << "\n";
		
		return coords;
	}
	
	
	
	template <int querydim, int datadim>
	std::vector<std::vector <double> >
	PerpleXLookup<querydim,datadim>::load_data (std::string filein) { //(const MPI_Comm &comm) {
		
		
		std::vector<std::vector <double> > data;
		
		std::string temp;
		// Read data from disk and distribute among processes
		std::ifstream in;//(Utilities::read_and_distribute_file_content(filename,comm));
		in.open(filein);
		
		
			
		// Skip header
		for (unsigned int i=0; i<13; ++i)
		{
			std::getline(in,temp); // throw these lines away
		}
		
		
		while (getline(in,temp)) //(in.good() && !in.eof())   // ????
		{
			std::istringstream line(temp);
			
			// read the two components of the coordinates
			double P_in, T_in;
			line >> P_in >> T_in;
			
			std::vector<double> values;
			for (unsigned int i=0; i<datadim; ++i) // !!! NEED TO CHANGE THIS TO FIND END OF LINE INSTEAD OF HARD-CODED NUMBER OF COLUMNS !!!
			{
				double value;
				line >> value;

				values.emplace_back(value);
			}
			data.emplace_back (values);
		}
		
		in.close();
		
		return data;
	}
	
}



namespace SeismicUtilities
{
	/**
	 * The two following functions compute the Voigt averages of the
	 * elastic moduli, Bulk Modulus and Shear Modulus.
	 * 
	 * Voigt averages are the "upper" bound estimate of the true average.
	 * It can be interpreted as the ratio of average strain within a
	 * composite material. In other words, the Voigt average is found by
	 * assuming that the strain is uniform throughout the composite
	 * material (iso-strain), which yields an upper bound estimate.
	 */
	double BulkMod_voigt (std::vector<double> bulkmods, std::vector<double> fractions)
	{
		double avg_bulkmod = 0.0;
		
		for (unsigned int i=0; i<bulkmods.size(); ++i)
		{
			avg_bulkmod = avg_bulkmod + fractions[i]*bulkmods[i];
		}
		
		return avg_bulkmod;
	}
	double ShearMod_voigt (std::vector<double> shearmods, std::vector<double> fractions)
	{
		double avg_shearmod = 0.0;
		
		for (unsigned int i=0; i<shearmods.size(); ++i)
		{
			avg_shearmod = avg_shearmod + fractions[i]*shearmods[i];
		}
		
		return avg_shearmod;
	}
	
	
	
	/**
	 * The two following functions compute the Reuss averages of the
	 * elastic moduli, Bulk Modulus and Shear Modulus.
	 * 
	 * Reuss averages are the "lower" bound estimate of the true average.
	 * It can be interpreted as the ratio of average stress within a
	 * composite material. In other words, the Reuss average is found by
	 * assuming that the stress is uniform throughout the composite
	 * material (iso-stress), which yields a lower bound estimate.
	 */
	double BulkMod_reuss (std::vector<double> bulkmods, std::vector<double> fractions)
	{
		double avg_bulkmod = 0.0;
		
		for (unsigned int i=0; i<bulkmods.size(); ++i)
		{
			avg_bulkmod = avg_bulkmod + ( fractions[i] / bulkmods[i] );
		}
		
		avg_bulkmod = 1.0 / avg_bulkmod;
		
		return avg_bulkmod;
	}
	double ShearMod_reuss (std::vector<double> shearmods, std::vector<double> fractions)
	{
		double avg_shearmod = 0.0;
		
		for (unsigned int i=0; i<shearmods.size(); ++i)
		{
			avg_shearmod = avg_shearmod + ( fractions[i] / shearmods[i] );
		}
		
		avg_shearmod = 1.0 / avg_shearmod;
		
		return avg_shearmod;
	}
	
	
	
	/**
	 * The two following functions compute the Voigt-Reuss-Hill averages
	 * of the elastic moduli, Bulk Modulus and Shear Modulus, from the
	 * upper (Voigt) and lower (Reuss) bound estimates.
	 * 
	 * These funcions call the Voigt and Reuss averaging functions,
	 * then compute the average of the results to provide an estimate
	 * of the true average.
	 */
	double bulkmod_voigtreusshill (std::vector<double> bulkmods, std::vector<double> fractions)
	{
		double bulkmod_voigtavg = BulkMod_voigt(bulkmods, fractions);
		double bulkmod_reussavg = BulkMod_reuss(bulkmods, fractions);
		
		return ( (bulkmod_voigtavg + bulkmod_reussavg) / 2.0);
	}
	double shearmod_voigtreusshill (std::vector<double> shearmods, std::vector<double> fractions)
	{
		double shearmod_voigtavg = ShearMod_voigt(shearmods, fractions);
		double shearmod_reussavg = ShearMod_reuss(shearmods, fractions);
		
		return ( (shearmod_voigtavg + shearmod_reussavg) / 2.0);
	}
	
	
	
	/**
	 * The two following functions compute average seismic velocities
	 * given average density, bulk and shear moduli.
	 */
	double vp_avg (double BulkMod_avg, double ShearMod_avg, double Density_avg)
	{
		return sqrt((BulkMod_avg+(4.0/3.0)*ShearMod_avg)/Density_avg);
	}
	double vs_avg (double ShearMod_avg, double Density_avg)
	{
		return sqrt(ShearMod_avg/Density_avg);
	}
	
	
	
	/**
	 * The following functions calculate bulk sound speed, vphi,
	 * given the isotropic Vp, Vs, and density.
	 */
	double
	calc_vphi_from_vpvs (double vp, double vs)
	{
		return ( sqrt( (vp*vp) - ((4.0/3.0)*(vs*vs)) ) );
	}
	double
	calc_vphi_from_bulk_den (double bulkmod, double density)
	{
		return ( sqrt( bulkmod/density ) );
	}
	
	
	
	/**
	 * The following two functions calculate the
	 * isotropic bulk and shear moduli given the
	 * isotropic Vp, Vs, and density.
	 */
	double
	calc_bulkmod_from_vpvsden (double vp, double vs, double density)
	{
		return ( ( (vp*vp) - ( (4.0/3.0)*(vs*vs) ) ) * density );
	}
	double
	calc_shearmod_from_vsden (double vs, double density)
	{
		return ( ( (vs*vs) * density ) );
	}
	
	
	
	
	/**
	 * The following functions convert the anisotropic Pv, Ph, Sv, and Sh values
	 * to the Voigt average isotropic Vp, Vs.
	 * 
	 * Conversion equations are from Appendix A of:
	 *      Panning, Mark, and Barbara Romanowicz.
	 *      "A three-dimensional radially anisotropic model of shear velocity in the whole mantle."
	 *      Geophysical Journal International 167.1 (2006): 361-379.
	 *
	 * Eqn (A8): Vs*Vs = (2*Vsv*Vsv + Vsh*Vsh) / 3
	 * 
	 * Eqn (A9): Vp*Vp = (Vpv*Vpv + 4*Vph*Vph) / 5
	 * 
	 * 
	 * These are approximations of:
	 * Eqn (A5): rho*Vp*Vp = (1/15) * (3C + (8 + 4*eta)*A + 8*(1 − eta)*L)
	 * 
	 * and
	 * 
	 * Eqn (A6): rho*Vs*Vs = (1/15) * (C + (1 + 2*eta)*A + (6 + 4*eta)*L + 5*N)
	 * where eta = F / (A - 2*L), and A, C, L, N, F are the Love parameters (Love 1927).
	 * 
	 * For eta ~ 1, eqn's A5 and A6 simplify to A8 and A9, respectively.
	 * See reference for more detail.
	 * 
	 * The eta value in the anisotropic PREM is equal to 1.0 except in the upper
	 * mantle, where it is >0.9. A more accurate calculation would be to use eqn's
	 * A5 and A6 for these depths, and should be revisited.
	 * 
	 * 
	 * Written 21 Sep 2018. PMBremner
	 **/
	double convert_PvPh_to_isoVp (double Vpv, double Vph)
	{
		double isoVp=0;
		
		// Calculate the Voigt average isotropic velocity
		
		// Eqn (A9):
		isoVp = sqrt( (Vpv*Vpv + 4.0*Vph*Vph) / 5.0 );
		
		return isoVp;
	}
	double convert_SvSh_to_isoVs (double Vsv, double Vsh)
	{
		double isoVs=0;
		
		// Calculate the Voigt average isotropic velocity
		
		// Eqn (A8)
		isoVs = sqrt( (2.0*Vsv*Vsv + Vsh*Vsh) / 3.0 );
		
		return isoVs;
	}
	
}



namespace PREMUtilities
{
	
	class PREMLookup
	{
	public:
		std::vector<double> load_coord (std::string filename, bool IncludeCore);
		std::vector<double> load_avg_den_vp_vs_overdepthrange (std::string filename, double mindepth, double maxdepth);
		std::vector< std::vector<double> > load_den_vp_vs (std::string filename, bool IncludeCore);
		std::vector< std::vector<double> > load_coord_and_den_vp_vs (std::string filename, bool IncludeCore);
		std::vector< std::vector<double> > load_coord_and_all_data (std::string filename, bool IncludeCore);
	};
	
	
	
	std::vector<double>
	PREMLookup::load_coord (std::string filename, bool IncludeCore=true)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<double> radii;
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			
			// Read in radius
			double radius;
			line >> radius;
			
			if (radius <= EarthConst::PREM_CoreRad && IncludeCore==true) radii.emplace_back (radius);
			else if (radius > EarthConst::PREM_CoreRad) radii.emplace_back (radius);
			
		}
		
		infile.close();
		
		return radii;
	}
	
	
	
	std::vector<double>
	PREMLookup::load_avg_den_vp_vs_overdepthrange (std::string filename, double mindepth, double maxdepth)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<double> data;
		double density_avg, Vp_avg, Vs_avg;
		double density_sum=0.0, Vp_sum=0.0, Vs_sum=0.0;
		int depthcount = 0;
		
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data Point (radius and data: Vp,Vs,Rho)
			double radius;
			double density, Vp, Vs;
			line >> radius >> density >> Vp >> Vs;
			
			double PREMdepth = EarthConst::mean_radius - radius;
			if (PREMdepth>=mindepth && PREMdepth<=maxdepth)
			{
				density_sum = density_sum + density;
				Vp_sum = Vp_sum + Vp;
				Vs_sum = Vs_sum + Vs;
				depthcount = ++depthcount;
			}
			
		}
		
		infile.close();
		
		density_avg = density_sum / depthcount;
		Vp_avg = Vp_sum / depthcount;
		Vs_avg = Vs_sum / depthcount;
		
		data.emplace_back (density_avg);
		data.emplace_back (Vp_avg);
		data.emplace_back (Vs_avg);
		
		
		return data;
	}
	
	
	
	std::vector< std::vector<double> >
	PREMLookup::load_den_vp_vs (std::string filename, bool IncludeCore=true)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<std::vector <double> > data;
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data Point (radius and data: Vp,Vs,Rho)
			double radius;
			double density, Vp, Vs;
			line >> radius >> density >> Vp >> Vs;
			
			if (radius <= EarthConst::PREM_CoreRad && IncludeCore==true)
			{
				
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				
				data.emplace_back (values);
			}
			else if (radius > EarthConst::PREM_CoreRad)
			{
				
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				
				data.emplace_back (values);
			}
		}
		
		infile.close();
		
		return data;
	}
	
	
	std::vector< std::vector<double> >
	PREMLookup::load_coord_and_den_vp_vs (std::string filename, bool IncludeCore=true)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<std::vector <double> > data;
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data Point (radius and data: Vp,Vs,Rho)
			double radius;
			double density, Vp, Vs;
			line >> radius >> density >> Vp >> Vs;
			
			if (radius <= EarthConst::PREM_CoreRad && IncludeCore==true)
			{
				values.emplace_back (radius);
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				
				data.emplace_back (values);
			}
			else if (radius > EarthConst::PREM_CoreRad)
			{
				values.emplace_back (radius);
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				
				data.emplace_back (values);
			}
		}
		
		infile.close();
		
		return data;
	}
	
	
	std::vector< std::vector<double> >
	PREMLookup::load_coord_and_all_data (std::string filename, bool IncludeCore=true)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<std::vector <double> > data;
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double radius;
			double density, Vp, Vs, Q_kappa, Q_mu, eta;
			line >> radius >> density >> Vp >> Vs >> Q_kappa >> Q_mu >> eta;
			
			if (radius <= EarthConst::PREM_CoreRad && IncludeCore==true)
			{
				values.emplace_back (radius);
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				values.emplace_back (Q_kappa);
				values.emplace_back (Q_mu);
				values.emplace_back (eta);
				
				
				data.emplace_back (values);
			}
			else if (radius > EarthConst::PREM_CoreRad)
			{
				values.emplace_back (radius);
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				values.emplace_back (Q_kappa);
				values.emplace_back (Q_mu);
				values.emplace_back (eta);
				
				
				data.emplace_back (values);
			}
		}
		
		infile.close();
		
		return data;
	}
}



namespace GypsumUtilities
{
	
	//template <class T>
	class GyPSuMLookup
	{
	public:
		std::vector<double> load_minmaxdepth (std::string filename);
		std::vector< std::vector<double> > load_coord (std::string filename);
		std::vector< std::vector<double> > load_data_perturb (std::string filename);
		std::vector< std::vector<double> > load_data_absolute (std::string filename, std::string reffile);
		std::vector< std::vector<double> > load_coord_data_absolute (std::string filename, std::string reffile);
		double perturb_to_absolute (double pertval, double refval);
		std::vector<double> load_temperature_variation (std::string filename);
		std::vector<double> load_pressure_variation (std::string filename);
		void make_combined_data (std::string densityfile, std::string vpfile, std::string vsfile, std::string outfile);
	};
	
	
	
	std::vector<double>
	GyPSuMLookup::load_minmaxdepth (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<double> depths;
		
		
		in_f.open(filename);
		
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		depths.emplace_back (mindepth);
		depths.emplace_back (maxdepth);
		
		in_f.close();
		
		
		return depths;
	}
	
	
	
	std::vector< std::vector<double> >
	GyPSuMLookup::load_coord (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<std::vector <double> > coordinates;
		
		
		in_f.open(filename);
		
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			line >> lat >> lon;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			
			coordinates.emplace_back (values);
		}
		
		in_f.close();
		
		
		return coordinates;
	}
	
	
	
	std::vector< std::vector<double> >
	GyPSuMLookup::load_data_perturb (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<std::vector <double> > data;
		
		
		in_f.open(filename);
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double density, vp, vs;
			line >> lat >> lon >> density >> vp >> vs;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (density);
			values.emplace_back (vp);
			values.emplace_back (vs);
			
			data.emplace_back (values);
		}
		
		in_f.close();
		
		
		return data;
	}
	
	
	
	std::vector< std::vector<double> >
	GyPSuMLookup::load_data_absolute (std::string filename, std::string reffile)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		PREMUtilities::PREMLookup lookup;
		GyPSuMLookup pert_abs;
		std::string temp;
		std::vector<std::vector <double> > data;
		
		
		in_f.open(filename);
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		// Retrieve the avg PREM values over the considered depth range
		// *** BE SURE TO FEED IN MINDEPTH AND MAXDEPTH IN METERS ***
		std::vector<double> ref_den_vp_vs = lookup.load_avg_den_vp_vs_overdepthrange(reffile, mindepth*1000.0, maxdepth*1000.0);
		std::cout << ref_den_vp_vs[0] << " " << ref_den_vp_vs[1] << " " << ref_den_vp_vs[2] << "\n";
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double density, vp, vs;
			line >> lat >> lon >> density >> vp >> vs;
			
			// Covert to absolute values
			// *** BE SURE TO CONVERT dlnrho, dlnVp, dlnVs from PERCENT PERTURBATION TO DECIMAL PERTURBATION ***
			density = pert_abs.perturb_to_absolute(density/100.0, ref_den_vp_vs[0]);
			vp = pert_abs.perturb_to_absolute(vp/100.0, ref_den_vp_vs[1]);
			vs = pert_abs.perturb_to_absolute(vs/100.0, ref_den_vp_vs[2]);
			
			//values.emplace_back (lat);
			//values.emplace_back (lon);
			values.emplace_back (density);
			values.emplace_back (vp);
			values.emplace_back (vs);
			
			data.emplace_back (values);
		}
		
		in_f.close();
		
		
		return data;
	}
	
	
	
	std::vector< std::vector<double> >
	GyPSuMLookup::load_coord_data_absolute (std::string filename, std::string reffile)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		PREMUtilities::PREMLookup lookup;
		GyPSuMLookup pert_abs;
		std::string temp;
		std::vector<std::vector <double> > data;
		
		
		in_f.open(filename);
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		// Retrieve the avg PREM values over the considered depth range
		// *** BE SURE TO FEED IN MINDEPTH AND MAXDEPTH IN METERS ***
		std::vector<double> ref_den_vp_vs = lookup.load_avg_den_vp_vs_overdepthrange(reffile, mindepth*1000.0, maxdepth*1000.0);
		std::cout << ref_den_vp_vs[0] << " " << ref_den_vp_vs[1] << " " << ref_den_vp_vs[2] << "\n";
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double density, vp, vs;
			line >> lat >> lon >> density >> vp >> vs;
			
			// Covert to absolute values
			// *** BE SURE TO CONVERT dlnrho, dlnVp, dlnVs from PERCENT PERTURBATION TO DECIMAL PERTURBATION ***
			density = pert_abs.perturb_to_absolute(density/100.0, ref_den_vp_vs[0]);
			vp = pert_abs.perturb_to_absolute(vp/100.0, ref_den_vp_vs[1]);
			vs = pert_abs.perturb_to_absolute(vs/100.0, ref_den_vp_vs[2]);
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (density);
			values.emplace_back (vp);
			values.emplace_back (vs);
			
			data.emplace_back (values);
		}
		
		in_f.close();
		
		
		return data;
	}
	
	
	
	double
	GyPSuMLookup::perturb_to_absolute (double pertval, double refval)
	{
		return ( (pertval*refval) + refval);
	}
	
	
	
	std::vector<double>
	GyPSuMLookup::load_temperature_variation (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<double> temperature_v;
		
		
		in_f.open(filename);
		
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			
			// Read In Input Data
			double lat, lon, Tvar;
			line >> lat >> lon >> Tvar;
			
			temperature_v.emplace_back (Tvar);
		}
		
		in_f.close();
		
		
		return temperature_v;
	}
	
	
	
	std::vector<double>
	GyPSuMLookup::load_pressure_variation (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<double> pressure_v;
		
		
		in_f.open(filename);
		
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon, Pvar;
			line >> lat >> lon;
			
			pressure_v.emplace_back (Pvar);
		}
		
		in_f.close();
		
		
		return pressure_v;
	}
	
	
	
	void
	GyPSuMLookup::make_combined_data (std::string densityfile, std::string vpfile, std::string vsfile, std::string outfile)
	{
		// Open Input File
		std::ifstream density_f, vp_f, vs_f;
		
		
		std::string temp;
		std::vector<std::vector <double> > densities, vpvals, vsvals;
		
		
		// -------------------------------------------
		//
		//  READ IN THE DENSITY FILE
		//
		// -------------------------------------------
		
		density_f.open(densityfile);
		
		// Read in the min and max depths  (in km)
		double density_mindepth, density_maxdepth;
		std::getline(density_f,temp);
		std::istringstream line(temp);
		
		line >> density_mindepth >> density_maxdepth;
		
		
		// Collect coordinates and data
		while (getline(density_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double value;
			line >> lat >> lon >> value;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (value);
			
			densities.emplace_back (values);
		}
		
		density_f.close();
		
		
		// -------------------------------------------
		//
		//  READ IN THE VP FILE
		//
		// -------------------------------------------
		
		vp_f.open(vpfile);
		
		// Read in the min and max depths  (in km)
		double vp_mindepth, vp_maxdepth;
		std::getline(vp_f,temp);
		std::istringstream line2(temp);
		
		line2 >> vp_mindepth >> vp_maxdepth;
		
		
		// Collect coordinates and data
		while (getline(vp_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double value;
			line >> lat >> lon >> value;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (value);
			
			vpvals.emplace_back (values);
		}
		
		vp_f.close();
		
		
		// -------------------------------------------
		//
		//  READ IN THE VS FILE
		//
		// -------------------------------------------
		
		vs_f.open(vsfile);
		
		// Read in the min and max depths  (in km)
		double vs_mindepth, vs_maxdepth;
		std::getline(vs_f,temp);
		std::istringstream line3(temp);
		
		line3 >> vs_mindepth >> vs_maxdepth;
		
		
		// Collect coordinates and data
		while (getline(vs_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double value;
			line >> lat >> lon >> value;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (value);
			
			vsvals.emplace_back (values);
		}
		
		vs_f.close();
		
		// -------------------------------------------
		//
		//  CHECK THE (DEPTH, LAT, LON) AND COMBINE
		//
		// -------------------------------------------
		
		// Check that all depths are the same
		if (density_mindepth != vp_mindepth)
		{
			std::cout << density_mindepth << " " << vp_mindepth << "\n";
			throw("Min depth ranges don't match, check input files for density and vp!");
		}
		if (density_maxdepth != vp_maxdepth)
		{
			std::cout << density_maxdepth << " " << vp_maxdepth << "\n";
			throw("Max depth ranges don't match, check input files for density and vp!");
		}
		
		if (vs_mindepth != vp_mindepth)
		{
			std::cout << vs_mindepth << " " << vp_mindepth << "\n";
			throw("Min depth ranges don't match, check input files for vp and vs!");
		}
		if (vs_maxdepth != vp_maxdepth)
		{
			std::cout << vs_maxdepth << " " << vp_maxdepth << "\n";
			throw("Max depth ranges don't match, check input files for vp and vs!");
		}
		
		
		// Check the number of data points
		if (densities.size() != vpvals.size()) throw("Data lengths don't match! Check input files for density and vp.");
		if (vsvals.size() != vpvals.size()) throw("Data lengths don't match! Check input files for vp and vs.");
		
		
		// Name and open the output file
		std::ofstream out_f;
		out_f.open(outfile);
		
		out_f << density_mindepth << " " << density_maxdepth << "\n";
		
		
		// Check coordinates...if good, paste the columns together
		for (unsigned int i=0; i<densities.size(); ++i)
		{
			if (densities[i][0]==vpvals[i][0] && densities[i][1]==vpvals[i][1])
			{
				if (vsvals[i][0]==vpvals[i][0] && vsvals[i][1]==vpvals[i][1])
				{
					out_f << densities[i][0] << " " << densities[i][1] << " " << densities[i][2] << " " << vpvals[i][2] << " " << vsvals[i][2] << "\n";
				}
			}
		}
		
		out_f.close();
		
		
		return;
	}
	
	
	
}


