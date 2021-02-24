//
//     Algorithms.hpp
//
//     Version 2.3  
//
//     Subroutines for Data Structures and Algorithms in C/C++
//     Copyright@ Yuan-Hsiang Chang, Ph.D
//     Department of Information and Computer Engineering
//     Chung Yuan Christian University
//
//     Update: May 2017
//
//
//     -------------------------------------------------------------------------------------------------------
//     Instructions
//     -------------------------------------------------------------------------------------------------------
//     1.  All subroutines were either implemented or collected from the given references below.
//     2.  The subroutines were originally compiled and tested using Visual C/C++.NET on Microsoft Windows.
//     3.  Templates are used for various KeyTypes (e.g., char, int, float, double, etc.).
//     4.  All subroutines may contain bugs, therefore are still subject to further revision.
//     5.  Some subtroutines may be incomplete (not optimal), and thus require further improvement.
//     6.  If you have any suggestions or would like to make contributions, please let me know 
//         (Extra bonus will be given accordingly).
//
//
//     -------------------------------------------------------------------------------------------------------
//     References:
//     -------------------------------------------------------------------------------------------------------
//     Introduction to Algorithms, 3rd Edition (Cormen)
//     Numerical Recipes in C++ (Press)
//     Data Structures and the Standard Template Library (Collins)
//     Data Structures and Algorithms in C++, 2nd Edition (Drozdek)
//     Data Structures and Algorithm Analysis in C, 2nd Ed (Weiss)
//     Introduction to the Design & Analysis of Algorithms (Levitin)
//     Fundamentals of Data Structures in C++ (Horowitz)
//     名題精選百則使用C語言 (冼鏡光)
//     資料結構使用C++ (蔡明志)
//     C 語言於演算法與資料結構之實習應用 (河西朝雄)
//
//
//     -------------------------------------------------------------------------------------------------------
//     C/C++ Data Type   Bytes                           Range
//     -------------------------------------------------------------------------------------------------------
//     char                1                          -128 ~ 127 
//     unsigned char       1                             0 ~ 255 
//     short int           2                       -32,768 ~ 32,767 
//     unsigned short      2                             0 ~ 65,535 
//     int                 4                -2,147,483,648 ~ 2,147,483,647 
//     unsigned int        4                             0 ~ 4,294,967,295 
//     long int            4                -2,147,483,648 ~ 2,147,483,647 
//     unsigned long       4                             0 ~ 4,294,967,295 
//     long long           8    -9,223,372,036,854,775,808 ~ 9,223,372,036,854,775,807 
//     unsigned long long  8                             0 ~ 18,446,744,073,709,551,615 
//     float               4                      -3.4e-38 ~ 3.4e+38 
//     double              8                     -1.7e-308 ~ 1.7e+308 
//
//
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <complex>

using namespace std;

typedef double DP;        // Numerical Recipes for double precision

// Numerical Recipes standard error handler
inline void nrerror( const string error_text )
{
	cout << "Numerical Recipes run-time error..." << endl;
	cout << "...now exiting to system..." << endl;
	exit(1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Date Class
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Given month / day / year, the month must be an integer between 1 and 12, inclusive;
//  the day must be an integer between 1 and the maximum number of days for the given
//  month and year, inclusive; the year must be an integer between 1800 and 2200, inclusive.
//
//  Reference: Data Structures and the Standard Template Library (Collins)
//
#define MONTHS_IN_YEAR    12  // 12 months a year 
#define DAYS_IN_WEEK       7  // 7 days a week
#define MIN_YEAR        1800  // Minimum year 
#define MAX_YEAR        2200  // Maximum year
#define MIN_FIRST_DAY      3  // first day of 1/1/1800 was Wednesday

class Date
{
public:
	Date();                                   // Date constructor: 1/1/1800 if not specified
	Date( int month, int day, int year );     // Date constructor: month/day/year
	void Previous();                          // Previous date
	void Next();                              // Next date
	void Display();                           // Display date 

	bool IsValid();                           // Is the date valid (or legal)?
	bool IsLeapYear( int year );              // Is the year a leap year?   
	int  DayOfWeek();                         // Day of week
	int  DaysInMonth( int month, int year );  // Days in month
	int  DaysLeftInYear();                    // Days left in the year
	int  FirstDay();                          // First day

protected:
    int dayIn;
    int monthIn;
    int yearIn;
}; 


Date::Date()
{
	monthIn = 1;  dayIn = 1;  yearIn = 1800;
}


Date::Date(int month, int day, int year)
{
	monthIn = month;  dayIn = day;  yearIn = year;
}


void Date::Previous()
{
	if ( dayIn > 1 )
		dayIn--;
	else if ( dayIn == 1 && monthIn != 1 )
	{
		monthIn--;
		dayIn = DaysInMonth( monthIn, yearIn );
	}
	else
	{
		monthIn = MONTHS_IN_YEAR;  dayIn = 31;  yearIn--;
	}
}


void Date::Next()
{
	if ( dayIn < DaysInMonth( monthIn, yearIn ) )
		dayIn++;
	else if ( dayIn == DaysInMonth( monthIn, yearIn ) && monthIn != MONTHS_IN_YEAR )
	{
		monthIn++;  dayIn = 1;
	}	
	else
	{
		monthIn = 1;  dayIn = 1;  yearIn++;
	}
}


void Date::Display()
{
	char *monthNames[] = { "January", "February", "March", "April", "May", "June", "July","August", 
						   "September", "October", "November", "December" };
	char *dayNames[] = { "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday" };
	
	if ( IsValid() )
	{
		cout << monthNames[monthIn-1] << " " << dayIn << ", " << yearIn;
		cout << " (" << dayNames[DayOfWeek()] << ")" << endl;
	}
	else
		cout << "Date is not valid" << endl;
}


bool Date::IsValid()
{
	if ( monthIn < 1 || monthIn > MONTHS_IN_YEAR || dayIn < 1 || yearIn < MIN_YEAR || yearIn > MAX_YEAR )
		return false;
	return dayIn <= DaysInMonth( monthIn, yearIn );
}


bool Date::IsLeapYear( int year )
{
	if ( year % 4 == 0 && ( year % 100 != 0 || year % 400 == 0 ) )
		return true;
	return false;
}


int Date::DayOfWeek()
{
	int thisDay;    // Take days from Jan 1 to this day, add day of week for Jan 1, % 7.
	if ( IsLeapYear ( yearIn ) )
		thisDay = ( 365 - DaysLeftInYear() + FirstDay() ) % DAYS_IN_WEEK;
	else
		thisDay = ( 364 - DaysLeftInYear() + FirstDay() ) % DAYS_IN_WEEK;
	return thisDay;
}


int Date::DaysInMonth( int month, int year )
{
	// 30 days for April, June, September, and November
	if ( month == 4 || month == 6 || month == 9 || month == 11 )
		return 30;
	// 31 days for January, March, May, July, August, October, December
	if ( month == 1 || month == 3 || month == 5 || month == 7 || month  == 8 || month == 10 || month == 12 )
		return 31;
	// One year is slightly less than 365.25 days
	if ( year % 4 != 0 || ( year % 100 == 0 && year % 400 != 0 ) )
		return  28;
	return 29;
}


int Date::DaysLeftInYear()
{
	int daysLeft = DaysInMonth( monthIn, yearIn ) - dayIn;
	for ( int i = monthIn + 1; i <= MONTHS_IN_YEAR; i++ )
		daysLeft += DaysInMonth ( i, yearIn );
	return daysLeft;
} // method daysLeftInYear


int Date::FirstDay()
{
	int temp = MIN_FIRST_DAY + ( yearIn - MIN_YEAR ) + ( yearIn - 1 - MIN_YEAR ) / 4;
	for ( int i = MIN_YEAR + 100; i < yearIn; i = i + 100 )
		if ( !IsLeapYear ( i ) )  temp--;
	return temp % DAYS_IN_WEEK;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Random Number Generators 
//
////////////////////////////////////////////////////////////////////////////////////////////////////
// 
//  Random number generator with uniform distribution [0,1].   
//  Set or reset idum to any integer value to initialize the sequence; 
//  idum must not be altered between calls for successive deviates in a sequence.
//  This routine is to replace the system-supplied rand().
//
//  Reference: Numerical Recipes in C++ (Press)
//
double ran0( int &idum )
{
	const int IA = 16807, IM = 2147483647, IQ = 127773;
	const int IR = 2836, MASK = 123459876;
	const DP AM = 1.0 / DP( IM );
	int k;
	DP ans;

	idum ^= MASK;
	k = idum / IQ;
	idum = IA * ( idum - k * IQ ) - IR * k;
	if ( idum < 0 ) idum += IM;
	ans = AM * idum;
	idum ^= MASK;
	return ans;
}


// 
//  Random number generator with uniform distribution [0,1].   
//  Set or reset idum to any integer value to initialize the sequence; 
//  idum must not be altered between calls for successive deviates in a sequence.
//  This routine is to replace the system-supplied rand().
//
//  Reference: Numerical Recipes in C++ (Press)
//
double ran1( int &idum )
{
	const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, NTAB = 32;
	const int NDIV = ( 1 + ( IM - 1 ) / NTAB );
	const DP EPS = 3.0e-16, AM = 1.0 / IM, RNMX = ( 1.0 - EPS );
	static int iy = 0;
	static int iv[NTAB];
	int j, k;
	DP temp;

	if ( idum <= 0 || !iy ) {
		if ( -idum < 1 ) idum = 1;
		else idum = -idum;
		for ( j = NTAB + 7 ; j >= 0 ; j-- ) {
			k = idum / IQ;
			idum = IA * ( idum - k * IQ ) - IR * k;
			if ( idum < 0 ) idum += IM;
			if ( j < NTAB ) iv[j] = idum;
		}
		iy = iv[0];
	}
	k = idum / IQ;
	idum = IA * ( idum - k * IQ ) - IR * k;
	if ( idum < 0 ) idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = idum;
	if ( ( temp = AM * iy ) > RNMX ) return RNMX;
	else return temp;
}


// 
//  Random number generator with uniform distribution [0,1].   
//  Set or reset idum to any integer value to initialize the sequence; 
//  idum must not be altered between calls for successive deviates in a sequence.
//  This routine is to replace the system-supplied rand().
//
//  Reference: Numerical Recipes in C++ (Press)
//
double ran2( int &idum )
{
	const int IM1 = 2147483563, IM2 = 2147483399;
	const int IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774;
	const int IR1 = 12211, IR2 = 3791, NTAB = 32, IMM1 = IM1 - 1;
	const int NDIV = 1 + IMM1 / NTAB;
	const DP EPS = 3.0e-16, RNMX = 1.0 - EPS, AM = 1.0 / DP( IM1 );
	static int idum2 = 123456789, iy = 0;
	static int iv[NTAB];
	int j,k;
	DP temp;

	if ( idum <= 0 ) 
	{
		idum = ( idum == 0 ? 1 : -idum );
		idum2 = idum;
		for ( j = NTAB + 7 ; j >= 0 ; j-- ) 
		{
			k = idum / IQ1;
			idum = IA1 * ( idum - k * IQ1 ) - k * IR1;
			if ( idum < 0 ) idum += IM1;
			if ( j < NTAB ) iv[j] = idum;
		}
		iy = iv[0];
	}

	k = idum / IQ1;
	idum = IA1 * ( idum - k * IQ1 ) - k * IR1;
	if ( idum < 0 ) idum += IM1;
	k = idum2 / IQ2;
	idum2 = IA2 * ( idum2 - k * IQ2 ) - k * IR2;
	if ( idum2 < 0 ) idum2 += IM2;
	j = iy / NDIV;
	iy = iv[j] - idum2;
	iv[j] = idum;
	if ( iy < 1 ) iy += IMM1;
	if ( ( temp = AM * iy ) > RNMX ) return RNMX;
	else return temp;
}


// 
//  Random number generator with uniform distribution [0,1].   
//  Set or reset idum to any integer value to initialize the sequence; 
//  idum must not be altered between calls for successive deviates in a sequence.
//  This routine is to replace the system-supplied rand().
//
//  Reference: Numerical Recipes in C++ (Press)
//
double ran3( int &idum )
{
	static int inext,inextp;
	static int iff = 0;
	const int MBIG = 1000000000, MSEED = 161803398, MZ = 0;
	const DP FAC = ( 1.0 / MBIG );
	static int ma[56];
	int i, ii, k, mj, mk;

	if ( idum < 0 || iff == 0 ) {
		iff = 1;
		mj = labs( MSEED - labs( idum ) );
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for ( i = 1 ; i <= 54 ; i++ ) {
			ii = ( 21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if ( mk < int( MZ ) ) mk += MBIG;
			mj = ma[ii];
		}
		for ( k = 0 ; k < 4 ; k++ )
			for ( i = 1 ; i <= 55 ; i++ ) {
				ma[i] -= ma[1+(i+30) % 55];
				if ( ma[i] < int( MZ ) ) ma[i] += MBIG;
			}
		inext = 0;
		inextp = 31;
		idum = 1;
	}
	if ( ++inext == 56 ) inext=1;
	if ( ++inextp == 56 ) inextp=1;
	mj = ma[inext] - ma[inextp];
	if ( mj < int( MZ ) ) mj += MBIG;
	ma[inext] = mj;
	return mj * FAC;
}


//
//  Exponential random number generator with exponential distribution.
//  Using ran1(idum) as the source of uniform distribution.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double expdev( int &idum )
{
	DP dum;

	do {
		dum = ran1( idum );
	} while ( dum == 0.0 );
	return -log( dum );
}


//
//  Gaussian random number generator with zero mean and unit variance.
//  Using ran1(idum) as the source of uniform distribution.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double gasdev( int &idum )
{
	static int iset = 0;
	static DP gset;
	DP fac, rsq, v1, v2;

	if ( idum < 0 ) iset=0;
	if ( iset == 0 ) {
		do {
			v1 = 2.0 * ran1( idum ) - 1.0;
			v2 = 2.0 * ran1( idum ) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while ( rsq >= 1.0 || rsq == 0.0 );
		fac = sqrt( -2.0 * log( rsq ) / rsq );
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	} else {
		iset = 0;
		return gset;
	}
}


//
//  Gamma random number generator with integer order.
//  Using ran1(idum) as the source of uniform distribution.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double gamdev( const int ia, int &idum )
{
	int j;
	DP am,e,s,v1,v2,x,y;

	if ( ia < 1 ) nrerror( "Error in routine gamdev." );
	if ( ia < 6 ) {
		x = 1.0;
		for ( j = 1 ; j <= ia ; j++ ) x *= ran1( idum );
		x = -log( x );
	} else {
		do {
			do {
				do {
					v1 = ran1( idum );
					v2 = 2.0 * ran1( idum ) - 1.0;
				} while ( v1 * v1 + v2 * v2 > 1.0 );
				y = v2 / v1;
				am = ia - 1;
				s = sqrt( 2.0 * am + 1.0 );
				x = s * y + am;
			} while ( x <= 0.0 );
			e = ( 1.0 + y * y ) * exp( am * log( x / am ) - s * y );
		} while ( ran1( idum ) > e );
	}
	return x;
}


//
//  Generation of Random Bits
//  Returns as an integer a random bit, based on the 18 low-significance bits in iseed
//  (which is modified for the next call).
//
//  Reference: Numerical Recipes in C++ (Press)
// 
int irbit1( unsigned long &iseed )
{
	unsigned long newbit;

	newbit =  ( ( iseed >> 17 ) & 1 ) ^ ( ( iseed >> 4 ) & 1) ^ ( ( iseed >> 1 ) & 1 ) ^ ( iseed & 1 );
	iseed = ( iseed << 1 ) | newbit;
	return int( newbit );
}


//
//  Generation of Random Bits
//  Returns as an integer a random bit, based on the 18 low-significance bits in iseed
//  (which is modified for the next call).
//
//  Reference: Numerical Recipes in C++ (Press)
// 
int irbit2( unsigned long &iseed )
{
	const unsigned long IB1 = 1, IB2 = 2, IB5 = 16, IB18 = 131072;
	const unsigned long MASK = IB1 + IB2 + IB5;

	if ( iseed & IB18 ) {
		iseed = ( ( iseed ^ MASK ) << 1 ) | IB1;
		return 1;
	} else {
		iseed <<= 1;
		return 0;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Combinations and Permutations
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Subsets
//  Given an array A[0..n-1], list all possible subsets in lexical order.
//  There are 2^n possible subsets.
//  Note: Input array A must be already in lexical order.
//
//  Reference: 名題精選百則使用C語言 (冼鏡光)
//
template <class KeyType>
void Subsets( KeyType *A, int n )
{
	int i, position;
	vector<int> set( n );
	
	cout << "{}";                
	position = 0;          
    set[position] = 1;                     
    while ( 1 ) 
	{                       
		cout << "\n{" << A[set[0]-1];
		for ( i = 1 ; i <= position ; i++ )
			cout << "," << A[set[i]-1];
		cout << "}";

        if ( set[position] < n ) { 
			set[position+1] = set[position] + 1;   
            position++;                    
        }
        else if ( position != 0 )  
			set[--position]++;  
        else                
            break;         
    }
	cout << endl;
}


//
//  K-Subsets
//  Given an array A[0..n-1], list all possible subsets with k elements in lexical order.
//  There are C(n, k) possible k-subsets.
//  Note: Input array A must be already in lexical order.
//
//  Reference: 名題精選百則使用C語言 (冼鏡光)
//
template <class KeyType>
void K_Subsets( KeyType *A, int n, int k )
{
	int i, j, position;
	vector<int> set( n );

	for ( i = 0 ; i < k ; i++ )  
		set[i] = i + 1;    

	cout << "{" << A[set[0]-1];
    for ( j = 1 ; j < k ; j++ )
		cout << "," << A[set[j]-1];
	cout << "}";

	position = k - 1;        
    while ( 1 ) 
	{         
		if ( set[k-1] == n )  
			position--;    
        else                
            position = k - 1; 
        set[position]++;    
        for ( i = position + 1 ; i < k; i++ ) 
			set[i] = set[i-1] + 1;

		cout << "\n{" << A[set[0]-1];
		for ( j = 1 ; j < k ; j++ )
			cout << "," << A[set[j]-1];
		cout << "}";

        if ( set[0] >= n - k + 1 ) break; 
	}
	cout << endl;
}


//
//  Permutations
//  Given an array A[0..n-1], list all possible permutations in lexical order.
//  There are n! possible permutations. 
//  Note: Input array A must be already in lexical order.
//
//  Reference: Data Structures and the Standard Template Library (Collins)
//
template <class KeyType>
void PermutationRec( KeyType *A, int i, int n )
{
	KeyType temp;
	int j, k;
	if ( i < n )
	{
		for ( j = i ; j < n ; j++ )
		{
			temp = A[j];
			for ( k = j ; k > i ; k-- )
				A[k] = A[k-1];
			A[i] = temp;

			PermutationRec( A, i + 1, n );

			for ( k = i ; k < j ; k++ )
				A[k] = A[k+1];
			A[j] = temp;
		}
	}
	else
	{
		for ( j = 0 ; j < n ; j++ )
			cout << A[j] << " ";
		cout << endl;
	}
}


template <class KeyType>
void Permutations( KeyType *A, int n )
{
	PermutationRec( A, 0, n );
}


//
//  K-Permutations
//  Given an array A[0..n-1], list all possible k-permutations.
//  This routine generate the k-Subsets first & then permutation.
//  There are C(n, k) * k! possible k-permutations.
//
//  Reference: 名題精選百則使用C語言 (冼鏡光)
//
template <class KeyType>
void K_Permutations( KeyType *A, int n, int k )
{
	int i, j, position;
	vector<int> set( n );
	KeyType *B;

	B = new KeyType[k + 1];   // For temporary K-Subsets

	for ( i = 0 ; i < k ; i++ )  
		set[i] = i + 1;    

	B[0] = A[set[0]-1];
	for ( j = 1 ; j < k ; j++ )
		B[j] = A[set[j]-1];

	Permutations( B, k );

	position = k - 1;        
    while ( 1 ) 
	{         
		if ( set[k-1] == n )  
			position--;    
        else                
            position = k - 1; 
        set[position]++;    
        for ( i = position + 1 ; i < k; i++ ) 
			set[i] = set[i-1] + 1;

		B[0] = A[set[0]-1];
		for ( j = 1 ; j < k ; j++ )
			B[j] = A[set[j]-1];

		Permutations( B, k );

        if ( set[0] >= n - k + 1 ) break; 
	}
	cout << endl;

	delete [] B;
}


//
//  Random Permutation
//  Given an array A[0..n-1], randomly permutate the array.
//  This routine is useful for rearrange the A array in random order.
// 
template <class KeyType>
void RandomPermutations( KeyType *A, int n )
{
	int i, j, seed = 3000;
	for ( i = 0 ; i < n ; i++ )
	{
		j = (int)( ran1( seed ) * (double) n );
		if ( j >= 0 && j < n )
			swap( A[i], A[j] );
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Recursive & Iterative Algorithms (Divide & Conquer)
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Factorial Number 
//  Compute the Factorial number n! 
//  (Recursive Approach)
//
double Factorial( const int n )
{
	if ( n <= 1 ) return 1;
	else	      return n * Factorial( n - 1 );
}


//
//  Factorial Number 
//  Compute the Factorial number n! 
//  (Iterative Approach)
//
double Factorial_Iterative( const int n )  
{
	double f;
	int i;

	if ( n <= 1 )
		return 1;
	else
	{
		f = 1;
		for ( i = 2 ; i <= n ; i++ )
			f *= (double) i;
		return f;
	}
}


//
//  Fibonacci Number 
//  Compute the Fibonacci number F(n) = F(n-1) + F(n-2)
//  (Recursive Approach)
//
double Fibonacci( const int n )
{
	if ( n == 0 ) return 0;
	if ( n <= 2 ) return 1;
	else          return ( Fibonacci( n - 1 ) + Fibonacci( n - 2 ) );
}


//
//  Fibonacci Number
//  Compute the Fibonacci number F(n) = F(n-1) + F(n-2)
//  (Iterative Approach)
//
double Fibonacci_Iterative( const int n )  
{
	double f0, f1, temp;
	int i;

	if ( n == 0 )  return 0;
	if ( n <= 2 )  return 1;
	else
	{
		for ( f0 = f1 = 1, i = 3 ; i <= n ; i++ )
		{
			temp = f0 + f1;
			f0 = f1;
			f1 = temp;
		}
		return f1;
	}
}


//
//  Binomial Coefficient 
//  Compute the binomial coefficient n and k, C(n, k)
//  (Recursive Approach)
//
double Cnk( const int n, const int k )
{
	if ( n == k || k == 0 )
		return 1;
	else
		return Cnk( n - 1, k ) + Cnk( n - 1, k - 1 );
}


//
//  Binomial Coefficient 
//  Compute the binomial coefficient n and k, C(n, k)
//  (Iterative Approach)
//
double Cnk_Iterative( const int n, const int k )  
{
	int i, j;
	double temp;

	vector<double> c( k + 1 );
	for ( i = 0 ; i <= k ; i++ )
		c[i] = 1;
	for ( i = 1 ; i <= n - k ; i++ )
		for ( j = 1 ; j <= k ; j++ )
			c[j] += c[j-1];
	temp = c[k];
	return temp;
}


//
//  Catalan Number 
//  Compute the nth catalan number
//  (Recursive Approach)
//
double Catalan( const int n )
{
	if ( n == 0 ) return 1;
	else
	{
		int k;
		double c = 0;
		for ( k = 0 ; k <= n - 1 ; k++ )
			c += ( Catalan( k ) * Catalan( n - 1 - k ) );
		return c;
	}
}


//
//  Catalan Number
//  Compute the nth catalan number
//  (Iterative Approach)
//
double Catalan_Iterative( const int n )
{
	if ( n == 0 || n == 1 ) return 1;
	else
	{
		vector<double> c( n + 1 );
		int i, k;
		c[0] = c[1] = 1;
		for ( i = 2 ; i <= n ; i++ )
		{
			c[i] = 0;
			for ( k = 0 ; k <= i - 1 ; k++ )
				c[i] += c[k] * c[i-1-k];
		}
		return c[n];
	}
}


//
//  Ackerman's Function 
//  Compute the Ackerman's function
//  (Recursive Approach)
//
double Ackerman( const int m, const int n )
{
	if ( m == 0 )      return n+1;
	else if ( n == 0 ) return Ackerman( m - 1, 1 );
	else               return Ackerman( m - 1, Ackerman( m, n - 1 ) );
}


//
//  Tower of Hanoi Problem 
//  Given three pegs, i.e., a: source  b: temporary  c: destination, 
//  move the disks originally stacked on the source peg to the destination peg
//  subject to the restriction of the rules of the Hanoi tower.
//  (Recursive Approach)
//
//  Input:  a = 'A', b = 'B', c = 'C', n: number of disks
//  Output: Procedures of moving the disks
//
void Hanoi_Tower( const char a, const char b, const char c, const int n )
{
	if ( n == 1 )
		cout << "Move disk from " << a << " -> " << c << "\n";
	else 
	{
		// Move the n-1 disks from A to B via C
		Hanoi_Tower( a, c, b, n - 1 );

		// Move 1 disk from A to C
		Hanoi_Tower( a, b, c, 1 );
		
		// Move the n-1 disks from B to C via A
		Hanoi_Tower( b, a, c, n - 1 );
	}
}


//
//  Tower of Hanoi Problem 
//  (Iterative Approach)
//
//  Input:  a = 'A', b = 'B', c = 'C', n: number of disks
//  Output: Procedures of moving the disks
//
//  Reference: 名題精選百則使用C語言 (冼鏡光)
//
int which_disk( unsigned long &counter )
{
	unsigned long a, b, c;
	int i;
	a = counter;
	counter = b = a+1;
	for ( c = a^b, i = 0 ; c != 0 ; c >>= 1, i++ )
		;
	return i;
}

void Hanoi_Tower_Iterative( const char a, const char b, const char c, const int n ) 
{
	unsigned long number_of_moves, counter;
	int dir[2], i, disk, next, index;
	char start, end;

	number_of_moves = ( 0x01UL << n ) - 1;
	counter = 0;                          // counter for Gray code
	
	if ( n & 0x01 )                       // setup direction
		dir[0] = 0, dir[1] = 1;
	else
		dir[0] = 1, dir[1] = 0;

	vector<int> pin( n + 1 );             // setup location
	for ( i = 1 ; i <= n ; i++ )
		pin[i] = 1;

	for ( i = 1 ; i <= number_of_moves ; i++ )       
	{
		disk = which_disk( counter );     // get disk #
		index = disk & 0x01;              // compute direction index  
		next = ( pin[disk] + dir[index] ) % 3 + 1;
		
		if ( pin[disk] == 1 ) start = a;
		if ( pin[disk] == 2 ) start = b;
		if ( pin[disk] == 3 ) start = c;
		if ( next == 1 ) end = a;
		if ( next == 2 ) end = b;
		if ( next == 3 ) end = c;
		
		cout << "Move disk from " << start << " -> " << end << endl;
		pin[disk] = next;
	}
}


//
//  Maximum-Subarray Problem
//  Given an array A[1..n] of positive (or negative) integers, find the
//  maximum-subarray.
//
//  Input:  A[1..n] & n
//  Output: low, high, & sum (Maximum subarray)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Find_Max_Crossing_Subarray( int *A, int low, int mid, int high, int &cross_low, int &cross_high, int &cross_sum )
{
	int max_left, max_right, left_sum, right_sum, sum;

	left_sum = -300000;
	sum = 0;
	for ( int i = mid ; i >= low ; i-- )
	{
		sum += A[i];
		if ( sum > left_sum )
		{
			left_sum = sum;
			max_left = i;
		}
	}
	
	right_sum = -300000;
	sum = 0;
	for ( int j = mid + 1 ; j <= high ; j++ )
	{
		sum += A[j];
		if ( sum > right_sum )
		{
			right_sum = sum;
			max_right = j;
		}
	}

	cross_low  = max_left;
	cross_high = max_right;
	cross_sum  = left_sum + right_sum;
}


void Find_Maximum_Subarray( int *A, int low, int high, int &low_ans, int &high_ans, int &sum )
{
	if ( high == low )
	{
		low_ans = low;
		high_ans = high;
		sum = A[low];
		return;
	}
	else
	{
		int mid = ( low + high ) / 2;

		int left_low, left_high, left_sum;
		Find_Maximum_Subarray( A, low, mid, left_low, left_high, left_sum );

		int right_low, right_high, right_sum;
		Find_Maximum_Subarray( A, mid+1, high, right_low, right_high, right_sum );

		int cross_low, cross_high, cross_sum;
		Find_Max_Crossing_Subarray( A, low, mid, high , cross_low, cross_high, cross_sum );

		if ( left_sum >= right_sum && left_sum >= cross_sum )
		{
			low_ans  = left_low;
			high_ans = left_high;
			sum      = left_sum;
		}
		else if ( right_sum >= left_sum && right_sum >= cross_sum )
		{
			low_ans  = right_low;
			high_ans = right_high;
			sum      = right_sum;
		}
		else
		{
			low_ans  = cross_low;
			high_ans = cross_high;
			sum      = cross_sum;
		}
	}
}


void Find_Maximum_Subarray( int *A, int n )
{
	int low, high, sum;
	Find_Maximum_Subarray( A, 1, n, low, high, sum );

	cout << "Maximum-Subarray" << endl;
	cout << "Low = " << low << " High = " << high << " Sum = " << sum;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Sorting & Order Statistics 
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Selection Sort
//  Sort an array A[0..n-1] into ascending numerical order using the Selection Sort algorithm
//  The array is placed on output by its sorted rearrangement.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
void SelectionSort( KeyType *A, int n )
{
	int i, j, k;
	for ( i = 0 ; i < n ; i++ ) 
	{
		j = i;
		for ( k = i + 1 ; k < n ; k++ )   // Find the smallest in a[i+1] to a[n-1]
			if ( A[k] < A[j] ) j = k;
		swap( A[i], A[j] );               // Swap with A[j]
	}
}


//
//  Insertion Sort
//  Sorts an array A[0..n-1] into ascending numerical order using the InsertionSort algorithm
//  The array is replaced on output by its sorted rearrangement.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
void InsertionSort( KeyType *A, int n )
{
	int i, j;
	KeyType key;
	for ( j = 1 ; j < n ; j++ )  
	{
		key = A[j];
		i = j - 1;               		    
		while ( i >= 0 && A[i] > key )     // Insert A[j] into the sorted sequence A[0..j-1]
		{
			A[i+1] = A[i];
			i--;
		}
		A[i+1] = key;
	}
}


//
//  Shell Sort
//  Sorts an array A[0..n-1] into ascending numerical order using the ShellSort algorithm
//  The array is replaced on output by its sorted rearrangement.
//
//  Reference: Wikipedia the Free Encyclopedia
//
template <class KeyType>
void ShellSort( KeyType *A, int n )
{
	int i, j, inc;
	KeyType temp;

    inc = n / 2;
    while ( inc > 0 )
    {
        for (i = inc ; i < n ; i++ )
        {
            j = i;
            temp = A[i];
            while ( ( j >= inc ) && ( A[j-inc] > temp ) )
            {
                A[j] = A[j - inc];
                j = j - inc;
            }
            A[j] = temp;
        }
        inc /= 2;
    }
}


// 
//  Merge Sort 
//  Sorts an array A[0..n-1] into ascending numerical order using the MergeSort algorithm
//  The array is replaced on output by its sorted rearrangement.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
void MergeSort( KeyType *A, int n )
{
	MergeSort( A, 0, n - 1 );
}

template <class KeyType>
void MergeSort( KeyType *A, int p, int r )
{
	int q;
	if ( p < r )
	{
		q = ( p + r ) / 2;
		MergeSort( A, p, q );
		MergeSort( A, q + 1, r );
		Merge( A, p, q, r );
	}
}

template <class KeyType>
void Merge( KeyType *A, int p, int q, int r )
{
	int i, j, k;
	int n1, n2;

	n1 = q - p + 1;
	n2 = r - q;

	vector<KeyType> L(n1), R(n2);
	for ( i = 0 ; i < n1 ; i++ ) L[i] = A[p+i];
    for ( j = 0 ; j < n2 ; j++ ) R[j] = A[q+j+1];

	i = j = 0;
	for ( k = p ; k <= r ; k++ )
	{
		if ( i < n1 && j < n2 )  
		{
			if ( L[i] <= R[j] ) 
			{
				A[k] = L[i];
				i++;
			}
			else
			{
				A[k] = R[j];
				j++;
			}
		}
		else if ( i < n1 && j >= n2 )
		{
			A[k] = L[i];
			i++;
		}
		else 
		{
			A[k] = R[j];
			j++;
		}
	}
}


// 
//  Heap Sort 
//  Sorts an array A[0..n-1] into ascending numerical order using the HeapSort algorithm
//  The array is replaced on output by its sorted rearrangement.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
void HeapSort( KeyType *A, int n )
{
	int     i, heap_size;
	KeyType temp;

	Build_Max_Heap( A, n );
	heap_size = n;
	for ( i = n - 1 ; i >= 0 ; i-- ) 
	{
		swap( A[0], A[i] ); 
		heap_size--;
		Max_Heapify( A, 0, heap_size );
	}   
}


template <class KeyType>
void Build_Max_Heap( KeyType *A, int n )
{
	int i, heap_size;
	heap_size = n;

	for ( i = n / 2 - 1 ; i >= 0 ; i-- )
		Max_Heapify( A, i, heap_size );
}


template <class KeyType>
void Max_Heapify( KeyType *A, int i, int heap_size )
{
	int l, r, largest;

	l = 2 * i + 1;
	r = 2 * i + 2;
	if ( l < heap_size && A[l] > A[i] )
		largest = l;
	else
		largest = i;

	if ( r < heap_size && A[r] > A[largest] )
		largest = r;

	if ( largest != i ) 
	{
		swap( A[i], A[largest] );
		Max_Heapify( A, largest, heap_size );
	}
}


// Min-Heapify is used for the Minimum Priority Queue 
template <class KeyType>
void Min_Heapify( KeyType *A, int i, int heap_size )
{
	int l, r, smallest;

	l = 2 * i + 1;
	r = 2 * i + 2;
	if ( l < heap_size && A[l] < A[i] )
		smallest = l;
	else
		smallest = i;
	if ( r < heap_size && A[r] < A[smallest] )
		smallest = r;

	if ( smallest != i ) 
	{
		// Exchange A[i] & A[smallest]
		swap( A[i], A[smallest] );
		Min_Heapify( A, smallest, heap_size );
	}
}


// 
//  QuickSort 
//  Sorts an array A[0..n-1] into ascending numerical order using the QuickSort algorithm
//  The array is replaced on output by its sorted rearrangement. 
//  (Recursive Approach)
//  
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
int Partition( KeyType *A, int p, int r )
{
	KeyType x;
	int i, j;
	
	x = A[r];
	i = p-1;
	for ( j = p ; j <= r - 1 ; j++ )
	{
		if ( A[j] <= x )
		{
			i++;
			swap( A[i], A[j] );
		}
	}
	swap( A[i+1], A[r] );
	return i+1;
}


template <class KeyType>
void QuickSortRec1( KeyType *A, int p, int r )
{
	if ( p < r )
	{
		int q = Partition( A, p, r );
		QuickSortRec1(A, p, q-1);
		QuickSortRec1(A, q+1, r);
	}
}


template <class KeyType>
void QuickSortRec( KeyType *A, int n )
{
	QuickSortRec1( A, 0, n-1 );
}


// 
//  QuickSort 
//  Sorts an array A[0..n-1] into ascending numerical order using the QuickSort algorithm
//  The array is replaced on output by its sorted rearrangement. The iterative approach
//  is implemented.
//
//  Note: This QuickSort routine is so far the most efficient that I know of.
//
//  Reference: Numerical Recipes in C++ (Press)
//  
template <class KeyType>
void QuickSort( KeyType *A, int n )
{
	const int M = 7, NSTACK = 50;
	int i, ir, j, k, jstack = -1, l = 0;
	KeyType a;
	int istack[NSTACK];

	ir = n - 1;
	for (;;) {
		if ( ir - l < M ) {
			for ( j = l + 1 ; j <= ir ; j++ ) {
				a=A[j];
				for ( i = j - 1 ; i >= l ; i-- ) {
					if ( A[i] <= a ) break;
					A[i+1] = A[i];
				}
				A[i+1] = a;
			}
			if ( jstack < 0 ) break;
			ir = istack[jstack--];
			l = istack[jstack--];
		} else {
			k = ( l + ir ) >> 1;
			swap( A[k], A[l+1] );
			if ( A[l] > A[ir] ) {
				swap( A[l], A[ir] );
			}
			if ( A[l+1] > A[ir] ) {
				swap( A[l+1], A[ir] );
			}
			if ( A[l] > A[l+1] ) {
				swap( A[l], A[l+1] );
			}
			i = l + 1;
			j = ir;
			a = A[l+1];
			for (;;) {
				do i++; while ( A[i] < a );
				do j--; while ( A[j] > a );
				if ( j < i ) break;
				swap( A[i], A[j] );
			}
			A[l+1] = A[j];
			A[j] = a;
			jstack += 2;
			if ( jstack >= NSTACK ) nrerror("NSTACK too small in sort.");
			if ( ir - i + 1 >= j - l ) {
				istack[jstack] = ir;
				istack[jstack-1] = i;
				ir = j - 1;
			} else {
				istack[jstack] = j - 1;
				istack[jstack-1] = l;
				l = i;
			}
		}
	}
}


//
//  QuickSort2
//  Sorts an array arr[0..n-1] into ascending order using QuickSort, while making the 
//  corresponding rearrangement of the array brr[0..-1].
//
//  Reference: Numerical Recipes in C++ (Press)
//
template <class KeyType>
void QuickSort2( KeyType *arr, KeyType *brr, int n )
{
	const int M = 7, NSTACK = 50;
	int i, ir, j, k, jstack = -1, l = 0;
	DP a,b;
	int istack[NSTACK];

	ir = n - 1;
	for (;;) {
		if ( ir - l < M ) {
			for ( j = l + 1 ; j <= ir ; j++ ) {
				a = arr[j];
				b = brr[j];
				for ( i = j - 1 ; i >= l ; i-- ) {
					if ( arr[i] <= a ) break;
					arr[i+1] = arr[i];
					brr[i+1] = brr[i];
				}
				arr[i+1] = a;
				brr[i+1] = b;
			}
			if ( jstack < 0 ) break;
			ir = istack[jstack--];
			l = istack[jstack--];
		} else {
			k = ( l + ir ) >> 1;
			swap( arr[k], arr[l+1] );
			swap( brr[k], brr[l+1] );
			if ( arr[l] > arr[ir] ) {
				swap( arr[l], arr[ir] );
				swap( brr[l], brr[ir] );
			}
			if ( arr[l+1] > arr[ir] ) {
				swap( arr[l+1], arr[ir] );
				swap( brr[l+1], brr[ir] );
			}
			if ( arr[l] > arr[l+1] ) {
				swap( arr[l], arr[l+1] );
				swap( brr[l], brr[l+1] );
			}
			i = l + 1;
			j = ir;
			a = arr[l+1];
			b = brr[l+1];
			for (;;) {
				do i++; while ( arr[i] < a );
				do j--; while ( arr[j] > a );
				if ( j < i ) break;
				swap( arr[i], arr[j] );
				swap( brr[i], brr[j] );
			}
			arr[l+1] = arr[j];
			arr[j] = a;
			brr[l+1] = brr[j];
			brr[j] = b;
			jstack += 2;
			if ( jstack >= NSTACK ) nrerror("NSTACK too small in QuickSort2.");
			if ( ir - i + 1 >= j - l ) {
				istack[jstack] = ir;
				istack[jstack-1] = i;
				ir = j - 1;
			} else {
				istack[jstack] = j - 1;
				istack[jstack-1] = l;
				l = i;
			}
		}
	}
}


//
//  CountingSort
//  Sorts an integer array A[0..n-1] into ascending numerical order and the sorted output is stored
//  in B[0..n-1] using the Counting Sort algorithm. The n input elements must be in the range of 
//  0 to k (Note: this is not checked!).
//  
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void CountingSort( int *A, int *B, int n, int k )
{
	int i, j;
	vector<int> C( k + 1 );

	for ( i = 0 ; i <= k ; i++ ) 
		C[i] = 0;
	
	for ( j = 0 ; j < n ; j++ )
		C[A[j]]++;
	
	// C[j] now contains the number of elements equal to i
	for ( i = 1 ; i <= k ; i++ )
		C[i] += C[i-1];

	// C[i] now contains the number of elements less than or equal to i
	for ( j = n - 1 ; j >= 0 ; j-- )
	{
		B[C[A[j]]-1] = A[j];
		C[A[j]]--;
	}
}


//
//  Median
//  Given an array A[0..n-1], find the median of the array without sorting.
//
//  Note that the median of an array of odd number of elements is the middle element
//  in sorted order. But for an array with even number of elements, the median is the
//  average of the two elements in the middle of the array.
//
//  Reference: 名題精選百則使用C語言 (冼鏡光)
//
template <class KeyType>
KeyType Median( KeyType *A, int n )
{
	int left = 0;
	int right = n-1;
	int mid = ( left + right ) / 2;
	int split_point;

	while(1)
	{
		Split( A, left, right, &split_point );
		if ( split_point == mid )
			break;
		else if ( split_point > mid )
			right = split_point - 1;
		else
			left = split_point + 1;
	}
	return ( n & 0x01 != 0 ) ? A[mid] : ( A[mid] + A[mid+1] ) / 2;
}


template <class KeyType>
void Split( KeyType *A, int left, int right, int *split_point )
{
	KeyType split_data = A[left];
	int i;

	for ( *split_point = left, i = left + 1 ; i <= right ; i++ )
	{
		if ( A[i] < split_data )
		{
			(*split_point)++;
			swap( A[*split_point], A[i] );
		}
	}
	swap( A[left], A[*split_point] );
}


//
//  Sequential Search
//  Search if a key exists in an array A[0..n-1] (not necessarily sorted). 
//  If yes, return true. Otherwise, return false.
//  The running time is an O(n).
//
//  Reference: C 語言於演算法與資料結構之實習應用 (河西朝雄)
// 
template <class KeyType>
bool SequentialSearch( KeyType *A, int n, KeyType key )
{
	for ( int i = 0 ; i < n ; i++ )
	{
		if ( A[i] == key )
			return true;
	}
	return false;
}


//
//  Binary Search
//  Search if a key exists in an "sorted" array A[0..n-1]. 
//  If yes, return true. Otherwise, return false.
//  The running time is an O(lgn).
//
//  Reference: C 語言於演算法與資料結構之實習應用 (河西朝雄)
//
template <class KeyType>
bool BinarySearch( KeyType *A, int n, KeyType key )
{
	int low, high, mid;
	low = 0;  high = n - 1;
	while ( low <= high )
	{
		mid = ( low + high ) / 2;
		if ( A[mid] == key )
			return true;
		else if ( A[mid] < key )
			low = mid + 1;
		else
			high = mid - 1;
	}
	return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Data Structures
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Matrix
//  Allocate Memory for a Matrix (or 2D Array)
//  nr: Number of rows
//  nc: Number of columns
//
//  Example:
//  Allocate a Matrix with 10 x 10 integer elements
//     int **v;
//     v = Matrix<int> (10, 10);
//
template <class KeyType>
KeyType **Matrix( int nr, int nc )
{
	KeyType **f;
	int i;

	f = (KeyType **) malloc( nr * sizeof( KeyType * ) );
	if ( !f ) return NULL;
	for ( i = 0 ; i < nr ; i++ )
	{
		f[i] = (KeyType *) malloc( nc * sizeof( KeyType ) );
		if ( !f[i] ) return NULL;
	}
	return f;
}


template <class KeyType>
void FreeMatrix( KeyType **f, int nr, int nc )
{
	int i;
	for ( i = 0 ; i < nr ; i++ )
		free( (KeyType *) f[i] );
	free( (KeyType *) f );
}


//
//  Stack 
//
//  Operations:
//     IsEmpty  (Check if Stack is empty)
//     IsFull   (Check if Stack is full)
//     Clear    (Clear Stack)
//     Push     (Push to Stack)
//     Pop      (Pop from Stack)
//     NumKeys  (Number of Keys)
//     Top      (Return the top key without pop)
//     Display  (Display Stack)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
class Stack {
public:
	Stack();
	Stack( int n );
	~Stack();

	bool    IsEmpty();
	bool    IsFull();
	void    Clear();
	void    Push( KeyType key );
	KeyType Pop();
	int     NumKeys();
	KeyType Top();
	void    Display();

protected:
	KeyType *S;                              // Array for the stack    
	int     size;                            // Stack Size
	int     top;                             // Point to top of the stack
};	


template <class KeyType>
Stack<KeyType>::Stack()
{
	S = new KeyType[100];                    // Allocate memory for the Stack 
	size = 100;                              // Define stack size (default = 100)
	top = -1;                                // Default top of the stack 
}


template <class KeyType>
Stack<KeyType>::Stack( int n )
{
	S = new KeyType[n];                      // Allocate memory for the Stack
	size = n;                                // Define stack size
	top = -1;                                // Default top of the stack 
}


template <class KeyType>
Stack<KeyType>::~Stack()
{
	delete [] S;                             // Release memory for the stack
}


template <class KeyType>
bool Stack<KeyType>::IsEmpty()
{
	if ( top < 0 ) return(true);             // If top is at -1, the stack is empty
	else           return(false);            // Otherwise, the stack is not empty
}


template <class KeyType>
bool Stack<KeyType>::IsFull()
{
	if ( top >= size - 1 )                   // If top is at stack size -1, the stack is full
		 return( true ); 
	else                                     // Otherwise, the stack is not full
		return( false ); 
}


template <class KeyType>
void Stack<KeyType>::Clear()
{
	top = -1;
}


template <class KeyType>
void Stack<KeyType>::Push( KeyType key )
{
	if ( IsFull() )                          // If the stack is full, unable to Push 
		cout << "Push fail (Stack is full!)\n";
	else
	{
		top++;
		S[top] = key;
	}
}


template <class KeyType>
KeyType Stack<KeyType>::Pop()
{
	if ( IsEmpty() )                         // If the stack is empty, unable to Pop  
	{
		cout << "Pop fail (Stack is empty!)\n";
		return 0;
	}
	else
	{
		top--;
	    return S[top+1];
	}
}


template <class KeyType>
int Stack<KeyType>::NumKeys()
{
	return top + 1;
}


template <class KeyType>
KeyType Stack<KeyType>::Top()
{
	if(IsEmpty())
		return 0;
	else
		return S[top];
}


template <class KeyType>
void Stack<KeyType>::Display()
{
	if ( IsEmpty() )                             
		cout << "Stack is empty!\n";         // Empty stack
	else
	{
		cout << "Stack:\n";                  // Display content of the stack 
		for ( int i = 0 ; i <= top ; i++ )
			cout << S[i] << " ";
		cout << endl;
	}
}


//
//  Queue (Circular)
//
//  Operations:
//     IsEmpty  (Check if Queue is empty)
//     IsFull   (Check if Queue is full)
//     Clear    (Clear Queue)
//     Enqueue  (EnQueue to Queue)
//     Dequeue  (DeQueue from Queue)
//     NumKeys  (Number of Keys)
//     Display  (Display Queue)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
class Queue {
public:
	Queue();
	Queue( int n );
	~Queue();

	bool    IsEmpty();
	bool    IsFull();
	void    Clear();
	void    Enqueue( KeyType key );
	KeyType Dequeue();
	int     NumKeys();
	void    Display();

protected:
	KeyType *Q;
	int     size, head, tail;
}; 


template <class KeyType>
Queue<KeyType>::Queue()
{
	Q = new KeyType[100];                 // Allocate memory for Queue
	size = 100;                           // Define Queue size (default = 100)
	head = tail = -1;                     // Default for empty Queue 
}


template <class KeyType>
Queue<KeyType>::Queue( int n )
{
	Q = new KeyType[n];                   // Allocate memory for Queue
	size = n;                             // Define Queue size
	head = tail = -1;                     // Default for empty Queue 
}


template <class KeyType>
Queue<KeyType>::~Queue()
{
	delete [] Q;                          // Release memory for the stack
}


template <class KeyType>
bool Queue<KeyType>::IsEmpty()
{
	if ( head == -1 )                     // Queue is empty (head = tail = -1) 
		return true;
	else
		return false;
}


template <class KeyType>
bool Queue<KeyType>::IsFull()
{
	if ( head == 0 && tail == size - 1 )  // If head at start, tail at end location
		return true;                      // Queue is full
	else if ( tail + 1 == head )          // If head at middle && tail+1 at the same location
		return true;                      // Queue is full
	else                                  // Otherwise
		return false;                     // Queue is not full
}


template <class KeyType>
void Queue<KeyType>::Clear()
{
	head = tail = -1;
}


template <class KeyType>
void Queue<KeyType>::Enqueue( KeyType key )
{
	if ( IsFull() )                       // If Queue is full, unable to Enqueue
		cout << "Enqueue fail (Queue is full!)\n";
	else
	{
		if ( tail == size - 1 || tail == -1 ) 
		{
			Q[0] = key;
			tail = 0;
			if ( head == -1 )             // Queue is empty
				head = 0;
		}
		else                                
		{
			tail++;
			Q[tail] = key;
		}
	} 
}


template <class KeyType>
KeyType Queue<KeyType>::Dequeue()
{
	if ( IsEmpty() )                      // If Queue is empty, unable to Dequeue
	{
		cout << "Dequeue fail (Queue is empty!)\n";
		return 0;
	}
	else
	{
		KeyType x = Q[head];
		if ( head == tail )               // Only 1 left                  
			head = tail = -1;             // After Dequeue, reset location
		else if ( head == size - 1 )      // If head is at the last location
			head = 0;                     // Use 0 for circular queue     
		else
			head++;
		return x;
	}
}


template <class KeyType>
int Queue<KeyType>::NumKeys()
{
	if ( IsEmpty() ) return 0;
	else
	{
		if ( head <= tail ) return tail - head + 1;
		else                return size - head + tail + 1;
	}
}


template <class KeyType>
void Queue<KeyType>::Display()
{
	if ( IsEmpty() )
	{
		cout << "Queue is empty!\n";
	}
	else
	{
		int i;
		cout << "Queue:\n";
		if ( head <= tail )                 // Display Queue 
		{
			for ( i = head ; i <= tail ; i++ )
				cout << Q[i] << " ";
			cout << "\n";
		}
		else                               // Display Queue (Circular)
		{ 
			for ( i = head ; i < size ; i++ ) 
				cout << Q[i] << " ";
            for ( i = 0 ; i <= tail ; i++ )
				cout << Q[i] << " ";
			cout << "\n";
		}
	}
}


// 
//  Maximum Priority Queue 
//  (Also known as the Maximum Heap)
//
//  Operations:
//     Clear          (Clear Priority Queue)
//     Insert         (Insert a key to Priority Queue)
//     Maximum        (Return the Maximum of Prioruity Queue)
//     ExtractMaximum (Extract the Maximum from Priority Queue)
//     NumKeys        (Number of Keys)
//     Display        (Display Priority Queue)  
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
class MaxPriorityQueue {
public:
	MaxPriorityQueue();
	MaxPriorityQueue( int n );
	~MaxPriorityQueue();

	void    Clear();
	void    Insert( KeyType key ); 
	KeyType Maximum();
	KeyType ExtractMaximum();
	int     NumKeys();
	void    Display();

protected:
	KeyType *Q;
	int     heap_size;
	int     size;
}; 


template <class KeyType>
MaxPriorityQueue<KeyType>::MaxPriorityQueue()
{
	Q = new KeyType[100];                 // Allocate memory for priority queue
    size = 100;                           // Define priority queue size (default = 100)
	heap_size = 0;                        // Actual data size in priority queue
}


template <class KeyType>
MaxPriorityQueue<KeyType>::MaxPriorityQueue( int n )
{
	Q = new KeyType[n];                   // Allocate memory for priority queue
    size = n;                             // Define priority queue size
	heap_size = 0;                        // Actual data size in priority queue
}


template <class KeyType>
MaxPriorityQueue<KeyType>::~MaxPriorityQueue()
{
	delete [] Q;                          // Release memory for the stack
}


template <class KeyType>
void MaxPriorityQueue<KeyType>::Clear()
{
	heap_size = 0;
}


template <class KeyType>
void MaxPriorityQueue<KeyType>::Insert( KeyType key )
{
	int i, j;
	if ( heap_size == 0 )
	{
		Q[0] = key;
		heap_size++;
	}
	else
	{
		Q[heap_size] = key;
		i = heap_size;
		j = (int)( ( i - 1 ) / 2 );
		while ( i != 0 )
		{
			if ( Q[j] < Q[i] )
			{
				swap( Q[i], Q[j] );
				i = j;
				j = (int)( ( i - 1 ) / 2 );
			}
			else
				break;
		}
		heap_size++;
	}
}


template <class KeyType>
KeyType MaxPriorityQueue<KeyType>::Maximum()
{
	return Q[0];
}


template <class KeyType>
KeyType MaxPriorityQueue<KeyType>::ExtractMaximum()
{
	KeyType max;

	if ( heap_size < 1 )
		cout << "Priority Queue Underflow...\n";
	else
	{
		max = Q[0];	
		Q[0] = Q[heap_size-1];
		heap_size--;

		// Max-Heapify
		Max_Heapify( Q, 0, heap_size );
		
		return max;
	}
	return 0;
}


template <class KeyType>
int MaxPriorityQueue<KeyType>::NumKeys()
{
	return heap_size;
}


template <class KeyType>
void MaxPriorityQueue<KeyType>::Display()
{
	if ( heap_size == 0 )
		cout << "No Keys in Maximum Priority Queue\n";
	else
	{
		cout << "Maximum Priority Queue:\n";
		for ( int i = 0 ; i < heap_size ; i++ )
		{
			cout << Q[i] << " ";
		}
		cout << endl;
	}
}


// 
//  Minimum Priority Queue 
//  (Also known as the Minimum Heap)
//
//  Operations:
//     Clear          (Clear Priority Queue)
//     Insert         (Insert a key to Priority Queue)
//     Minimum        (Return the Minimum of Prioruity Queue)
//     ExtractMaximum (Extract the Minimum from Priority Queue)
//     NumKeys        (Number of Keys)
//     Display        (Display Priority Queue)  
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
class MinPriorityQueue {
public:
	MinPriorityQueue();
	MinPriorityQueue( int n );
	~MinPriorityQueue();

	void    Clear();
	void    Insert( KeyType key ); 
	KeyType Minimum();
	KeyType ExtractMinimum();
	int     NumKeys();
	void    Display();

protected:
	KeyType *Q;
	int     heap_size;
	int     size;
}; 


template <class KeyType>
MinPriorityQueue<KeyType>::MinPriorityQueue()
{
	Q = new KeyType[100];                 // Allocate memory for Queue
    size = 100;                           // Define Queue size (default = 100)
	heap_size = 0;                        // Actual data size in priority queue
}


template <class KeyType>
MinPriorityQueue<KeyType>::MinPriorityQueue(int n)
{
	Q = new KeyType[n];                   // Allocate memory for Queue
    size = n;                             // Define Queue size
	heap_size = 0;                        // Actual data size in priority queue
}


template <class KeyType>
MinPriorityQueue<KeyType>::~MinPriorityQueue()
{
	delete [] Q;                          // Release memory for the stack
}


template <class KeyType>
void MinPriorityQueue<KeyType>::Clear()
{
	heap_size = 0;
}


template <class KeyType>
void MinPriorityQueue<KeyType>::Insert( KeyType key )
{
	int i, j;
	if ( heap_size == 0 )
	{
		Q[0] = key;
		heap_size++;
	}
	else
	{
		Q[heap_size] = key;
		i = heap_size;
		j = (int)( ( i - 1 ) / 2 );
		while( i != 0 )
		{
			if ( Q[j] > Q[i] )
			{
				swap( Q[i], Q[j] );
				i = j;
				j = (int)( ( i - 1 ) / 2 );
			}
			else
				break;
		}
		heap_size++;
	}
}


template <class KeyType>
KeyType MinPriorityQueue<KeyType>::Minimum()
{
	return Q[0];
}


template <class KeyType>
KeyType MinPriorityQueue<KeyType>::ExtractMinimum()
{
	KeyType min;

	if ( heap_size < 1 )
		cout << "Priority Queue Underflow...\n";
	else
	{
		min = Q[0];	
		Q[0] = Q[heap_size-1];
		heap_size--;

		// Min-Heapify
		Min_Heapify( Q, 0, heap_size );
		return min;
	}
	return 0;
}


template <class KeyType>
int MinPriorityQueue<KeyType>::NumKeys()
{
	return heap_size;
}


template <class KeyType>
void MinPriorityQueue<KeyType>::Display()
{
	if ( heap_size == 0 )
		cout << "No Keys in Minimum Priority Queue\n";
	else
	{
		cout << "Minimum Priority Queue:\n";
		for ( int i = 0 ; i < heap_size ; i++ )
		{
			cout << Q[i] << " ";
		}
		cout << endl;
	}
}


//
//  Singly Linked List
//
//  Operations:
//     IsEmpty    (Check if Linked List is empty)
//     Clear      (Clear Linked List)
//     Seach      (Seach for a key)
//     Insert     (Insert to the head of Linked List)
//     Delete     (Delete from the head of Linked List)
//     InsertTail (Insert to the tail of Linked List)
//     DeleteTail (Delete from the tail of Linked List)
//     InsertSort (Insert in sorted order) Note: The Linked List must be already in sorted order.
//     Invert     (Invert Linked List)
//     NumKeys    (Number of Keys)
//     Display    (Display Linked List)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
struct SLLNode {                  // Node for Singly Linked List 
	KeyType key;                  // Key
	SLLNode<KeyType> *next;       // Link to next node
};


template <class KeyType>
class SinglyLinkedList {
public:
	SinglyLinkedList();
	~SinglyLinkedList();

	bool IsEmpty();
	void Clear();
	bool Search( KeyType key );
	void Insert( KeyType key );    
	void Delete();       
	void InsertTail( KeyType key );
	void DeleteTail();
	void InsertSort( KeyType key );
	void Invert();
	int  NumKeys();
	void Display();

protected:
	SLLNode<KeyType> *head;
	SLLNode<KeyType> *current;
	int n;  

	friend class AdjacencyMatrix;
	friend class AdjacencyList;
};


template <class KeyType>
SinglyLinkedList<KeyType>::SinglyLinkedList()
{
	head = new SLLNode<KeyType>;           // Create a head node
	head->next = NULL;
	n = 0;  
}


template <class KeyType>
SinglyLinkedList<KeyType>::~SinglyLinkedList()
{
	current = head->next;
	while ( current != NULL )
	{
		delete head;
		head = current;
		current = current->next;
	} 
}


template <class KeyType>
bool SinglyLinkedList<KeyType>::IsEmpty()
{
	if ( head->next == NULL ) return true;
	else                      return false;
}


template <class KeyType>
void SinglyLinkedList<KeyType>::Clear()
{
	SLLNode<KeyType> *ptr;
	ptr = head->next;
	while ( ptr != NULL )
	{
		current = ptr->next;
		delete ptr;
		ptr = current;
	}
	head->next = NULL;
	n = 0;
}


template <class KeyType>
bool SinglyLinkedList<KeyType>::Search( KeyType key )
{
	current = head->next;
	while ( current != NULL )
	{
		if ( current->key == key ) 
			return true;
		else
			current = current->next;
	}
	return false;
}


template <class KeyType>
void SinglyLinkedList<KeyType>::Insert( KeyType key )
{
	current = new SLLNode<KeyType>;        // Create a new node for insertion
	current->key = key;                    // Store key
	current->next = head->next;            // Insert new link
	head->next = current;                  // New node as the first node  
	n++;                                   // Number of keys + 1
}


template <class KeyType>
void SinglyLinkedList<KeyType>::Delete()
{
	if ( head->next == NULL )              // Empty Linked List 
		cout << "No Keys in Singly Linked List\n";
	else
	{
		current = head->next;              // Get the first node
		head->next = current->next;        // Redefine link
		delete current;                    // Delete node  
		n--;
	}
}


template <class KeyType>
void SinglyLinkedList<KeyType>::InsertTail( KeyType key )
{
	SLLNode<KeyType> *NewNode = new SLLNode<KeyType>;  // Allocate memory for a new node
	NewNode->key = key;                    // Store key 
	NewNode->next = NULL;                  // The new node is at tail of the linked list   

	if ( head->next == NULL )              // If there is no node initially
		head->next = NewNode;
	else
	{
		current = head->next;              // Find the current tail node
		while ( current->next != NULL ) 
			current = current->next;
		current->next = NewNode;           // Insert the node
	}
	n++;                                   // Number of keys + 1
}


template <class KeyType>
void SinglyLinkedList<KeyType>::DeleteTail()
{
	SLLNode<KeyType> *prev;

	if ( head->next == NULL )              // Empty linked list
		cout << "No Keys in Singly Linked List\n";
	else
	{
		prev = head;
		current = head->next;              // Find the current tail node
		while ( current->next != NULL )
		{
			prev = current;
			current = current->next;
		}
		prev->next = NULL;                 // Redefine link
		delete current;                    // Delete the tail node
		n--;
	}
}


template <class KeyType>
void SinglyLinkedList<KeyType>::InsertSort( KeyType key )
{
	SLLNode<KeyType> *prev;
	
	SLLNode<KeyType> *NewNode = new SLLNode<KeyType>;  // Allocate memory for a new node
	NewNode->key = key;                                // Store key 
	NewNode->next = NULL;                              // The new node is at tail of the linked list   

	if ( head->next == NULL )                          // If there is no node initially
		head->next = NewNode;
	else
	{
		current = head->next;              
		if ( current->next == NULL )                   // One node only   
		{
			if ( key < current->key )
			{
				head->next = NewNode;
				NewNode->next = current;
			}
			else
			{
				current->next = NewNode;
			}
		}
		else                                          // There are at least two nodes
		{
			while ( current->next != NULL )           // Search to insert
			{
				prev = current;
				current = current->next;
				if ( key < prev->key )                // Inserted key is smaller
				{
					head->next = NewNode;
					NewNode->next = prev;
					break;
				}
				else if ( key >= prev->key && key < current->key )  // Inserted key is in between
				{
					prev->next = NewNode;
					NewNode->next = current;
					break;
				}
				else                                  // Inserted key is larger      
				{
					if ( current->next == NULL )      // Already last node, simply insert, or keep searching otherwise
					{
						current->next = NewNode;
						break;
					}
				}
			}
		}
	}
	n++;                                   // Number of keys + 1
}


template <class KeyType>
void SinglyLinkedList<KeyType>::Invert()
{
	SLLNode<KeyType> *prev, *cur, *forward;
	forward = head->next;
	cur = NULL;
	while ( forward != NULL )                    
	{
		prev = cur;
		cur = forward;
		forward = forward->next;
		cur->next = prev;
	}
	head->next = cur;
}


template <class KeyType>
int SinglyLinkedList<KeyType>::NumKeys()
{
	return n;
}


template <class KeyType>
void SinglyLinkedList<KeyType>::Display()
{
	if ( head->next == NULL )
		cout << "No Keys in Singly Linked List\n";
	else
	{
		cout << "Singly Linked List:\n";
		current = head->next;
		while ( current != NULL )
		{
			cout << current->key << " ";
			current = current->next;
		}
		cout << "\n";
	}
} 


//
//  Doubly Linked list
//
//  Operations:
//     IsEmpty    (Check if Linked List is empty)
//     Clear      (Clear Linked List)
//     Seach      (Seach for a key)
//     Insert     (Insert to the head of Linked List)
//     Delete     (Delete from the head of Linked List)
//     InsertTail (Insert to the tail of Linked List)
//     DeleteTail (Delete from the tail of Linked List)
//     NumKeys    (Number of Keys)
//     Display    (Display Linked List)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
struct DLLNode {                 // Node for Doubly Linked List  
	KeyType key;                 // Key
	DLLNode<KeyType> *left;      // Left link
	DLLNode<KeyType> *right;     // Right link
};


template <class KeyType>
class DoublyLinkedList {
public:
	DoublyLinkedList();
	~DoublyLinkedList();

	bool IsEmpty();
	void Clear();
	bool Search( KeyType key );
	void Insert( KeyType key );
	void Delete();
	void InsertTail( KeyType key );
	void DeleteTail();
	int  NumKeys();
	void Display();

protected:
	DLLNode<KeyType> *ptr;
	DLLNode<KeyType> *head;
	DLLNode<KeyType> *tail;
	DLLNode<KeyType> *prev;
	DLLNode<KeyType> *current;
	int n;
};


template <class KeyType>
DoublyLinkedList<KeyType>::DoublyLinkedList()
{
	head = new DLLNode<KeyType>;           // Create a head node
	head->left = head;                     // Both left & right links to itself
	head->right = head;
	n = 0;
}


template <class KeyType>
DoublyLinkedList<KeyType>::~DoublyLinkedList()
{
	ptr = head->right;
	while ( ptr != head )
	{
		current = ptr;
		ptr = current->right;
		delete current;
	} 
	delete head;
}


template <class KeyType>
bool DoublyLinkedList<KeyType>::IsEmpty()
{
	if ( head->right == head ) 
		return true;
	else                       
		return false;
}


template <class KeyType>
void DoublyLinkedList<KeyType>::Clear()
{
	ptr = head->right;
	while ( ptr != head )
	{
		current = ptr->right;
		delete ptr;
		ptr = current;
	}
	head->right = head;
	n = 0;
}


template <class KeyType>
bool DoublyLinkedList<KeyType>::Search( KeyType key )
{
	current = head->right;
	while ( current != head )
	{
		if ( current->key == key ) 
			return true;
		else
			current = current->right;
	}
	return false;
}


template <class KeyType>
void DoublyLinkedList<KeyType>::Insert( KeyType key )
{
	ptr = new DLLNode<KeyType>;            // Create a new node for insertion
	ptr->key = key;                        // Store key

	current = head->right;                 // 1st node in the Linked List
	current->left = ptr;                   // Left link to the new node

	ptr->left = head;                      // Define links for the new node
	ptr->right = current;
	head->right = ptr;
	
	n++;
}


template <class KeyType>
void DoublyLinkedList<KeyType>::Delete()
{
	if ( head->right == head )
		cout << "No Keys in Doubly Linked List\n";
	else
	{
		current = head->right;             // The node for deletion
		head->right = current->right;      // Redirect the link 
		current->right->left = current->left;
		delete current;
		n--;
	}
}


template <class KeyType>
void DoublyLinkedList<KeyType>::InsertTail( KeyType key )
{
	ptr = new DLLNode<KeyType>;
	ptr->key = key;
	
	prev = head;                           // Previous node is the head node   
	current = head->right;                 // 1st node in the Linked List
	while ( current != head )              // Find the tail node 
	{
		prev = current;
		current = current->right;
	}
	current->left = ptr;

	ptr->left = prev;                     // Define links for the new node
	ptr->right = current;
	prev->right = ptr;
	
	n++;
}


template <class KeyType>
void DoublyLinkedList<KeyType>::DeleteTail()
{
	if ( head->right == head )
		cout << "No Keys in Doubly Linked List\n";
	else
	{
		prev = head;
		current = head->right;
		while ( current != head )           // Find the node for deletion
		{
			prev = current;
			current = current->right;
		}
		prev->left->right = head;
		delete prev;
		n--;
	}
}


template <class KeyType>
int DoublyLinkedList<KeyType>::NumKeys()
{
	return n;
}


template <class KeyType>
void DoublyLinkedList<KeyType>::Display()
{
	if ( head->right == head )
		cout << "No Keys in Doubly Linked List\n";
	else
	{
		cout << "Doubly Linked List:\n";
		current = head->right;
		while ( current != head )
		{
			cout << current->key << " ";
			current = current->right;
		}
		cout << "\n";
	}
}


//
//  Binary Search Tree (BST)
//
//  Operations:
//     Seach     (Seach for a key)
//     Insert    (Insert to BST)
//     Delete    (Delete from BST)
//     Preorder  (Preorder traversal)
//     Inorder   (Inorder traversal)
//     Postorder (Postorder traversal)
//     NumKeys   (Number of keys)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
struct BSTNode {
	KeyType key;                 // Key
	BSTNode<KeyType> *left;      // Link to left child
	BSTNode<KeyType> *right;     // Link to right child
};


template <class KeyType>
class BinarySearchTree {
public:
	BinarySearchTree();
	~BinarySearchTree();

	bool Search( KeyType key );
	void Insert( KeyType key );
	void Delete( KeyType key );
	void Preorder();
	void Inorder();
	void Postorder();
	int  NumKeys();

protected:
	BSTNode<KeyType> *root;      
	BSTNode<KeyType> *current;   
	BSTNode<KeyType> *parent;
	BSTNode<KeyType> *ptr;
	BSTNode<KeyType> *prev;
	int n;

	// Auxiliary
	void Preorder(BSTNode<KeyType> *ptr);   
	void Inorder(BSTNode<KeyType> *ptr);
	void Postorder(BSTNode<KeyType> *ptr);
};


template <class KeyType>
BinarySearchTree<KeyType>::BinarySearchTree()
{
	root = NULL;
	n = 0;
} 


template <class KeyType>
BinarySearchTree<KeyType>::~BinarySearchTree()
{
	if ( root != NULL ) delete root;
} 


template <class KeyType>
bool BinarySearchTree<KeyType>::Search( KeyType key )
{
	current = root;
	while ( current != NULL )             
	{
		if ( key == current->key ) 
			return true;
		else if ( key < current->key )
			current = current->left;
		else
			current = current->right;
	}
	return false;
}


template <class KeyType>
void BinarySearchTree<KeyType>::Insert( KeyType key )
{
	ptr = new BSTNode<KeyType>;             // Create a new BST node
	ptr->key = key;                         // Assign key
	ptr->left = ptr->right = NULL;          // No left & right children  

	if ( root == NULL )                     // Empty BST, the new node becomes the root
		root = ptr;
	else
	{
		current = root;
		while ( current != NULL )           // Search for the location to insert 
		{
			parent = current;
			if ( key < current->key )
				current = current->left;
			else
				current = current->right;
		}

		if ( key < parent->key )               // Insert the new node 
			parent->left = ptr;
		else
			parent->right = ptr;
	}
	n++;
}


template <class KeyType>
void BinarySearchTree<KeyType>::Delete(KeyType key)
{
	parent = NULL;                          
	current = root;
	while ( current != NULL )                                        // Search for the deleted node
	{
		if ( key == current->key )                                   // The deleted node is found
			break;
		parent = current;                                            // Also determine the deleted node's parent
		if ( key < current->key )               
			current = current->left;
		else                               
			current = current->right;
	}

	if ( current == NULL )                                           // Deleted node is not found
		cout << "Deleted node not found in Binary Search Tree.\n";
	else
	{
		if ( current->left == NULL && current->right == NULL )       // No left & right subtree (leaf node)
		{
			if ( current == root )                                   // The deleted node is the root with no children
				root = NULL;
			else                                                     // Otherwise, the deleted node is a leaf node
			{
				if ( current->key < parent->key )                    // Determine the deleted node's parent 
					parent->left = NULL;
				else                                            
					parent->right = NULL;
				delete current;
			}
		}
		else if ( current->left == NULL && current->right != NULL )  // Right subtree is not empty
		{
			if ( current == root )                                   // Deleted node is the root
				root = current->right;                               // Its right child becomes the root
			else                                                     // Otherwise
			{
				if ( current->key < parent->key )                    // Link its right subtree 
					parent->left = current->right;
				else
					parent->right = current->right;
			}
			delete current;
		}
		else if ( current->left != NULL && current->right == NULL )  // Left subtree is not empty
		{
			if ( current == root )                                   // Deleted node is the root
				root = current->left;
			else		                                             // Otherwise
			{
				if ( current->key < parent->key )                    // Link its left subtree
					parent->left = current->left;
				else
					parent->right = current->left;
			}
			delete current;
		}
		else                                                         // Both left & right subtrees are not empty
		{
			prev = current;
			ptr  = current->right;                                   // Search for replacing node
			while ( ptr->left != NULL )                              // Use smallest node in right subtree (default)
			{                                                        // or largest node in left subtree (not used here)
				prev = ptr;
				ptr = ptr->left;
			}

			if ( ptr->key < prev->key )                              // Link the deleted node's parent to its children
				prev->left = ptr->right;
			else
				prev->right = ptr->right;

			current->key = ptr->key;                                 // Replace key
			delete ptr;
		}
		n--;
	}
}


template <class KeyType>
void BinarySearchTree<KeyType>::Preorder()
{
	if ( root == NULL )
		cout << "No Keys in Binary Search Tree.\n";
	else
	{
		cout << "Binary Search Tree (Preorder):\n";
		Preorder(root);
		cout << "\n";
	}
}


template <class KeyType>
void BinarySearchTree<KeyType>::Preorder(BSTNode<KeyType> *ptr)
{
	if ( ptr != NULL )                      
	{
		cout << ptr->key << " ";
		Preorder(ptr->left);
		Preorder(ptr->right);
	}
}


template <class KeyType>
void BinarySearchTree<KeyType>::Inorder()
{
	if ( root == NULL )
		cout << "No Keys in Binary Search Tree.\n";
	else
	{
		cout << "Binary Search Tree (Inorder):\n";       
		Inorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void BinarySearchTree<KeyType>::Inorder( BSTNode<KeyType> *ptr )
{
	if ( ptr != NULL )
	{
		Inorder(ptr->left);
		cout << ptr->key << " ";
		Inorder(ptr->right);
	}
}


template <class KeyType>
void BinarySearchTree<KeyType>::Postorder()
{
	if ( root == NULL )
		cout << "No Keys in Binary Search Tree.\n";
	else
	{
		cout << "Binary Search Tree (Postorder):\n";
		Postorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void BinarySearchTree<KeyType>::Postorder( BSTNode<KeyType> *ptr )
{
	if ( ptr != NULL )
	{
		Postorder( ptr->left );
		Postorder( ptr->right );
        cout << ptr->key << " ";
	}
}


template <class KeyType>
int BinarySearchTree<KeyType>::NumKeys()
{
	return n;
}


//
//  AVL Tree 
//
//  Operations:
//     Insert    (Insert to AVL Tree)
//     Delete    (Delete from AVL Tree)
//     Preorder  (Preorder traversal)
//     Inorder   (Inorder traversal)
//     Postorder (Postorder traversal)
//     NumKeys   (Number of keys)
//
//  Reference: 資料結構使用C++ (蔡明志)
//
template <class KeyType>
struct AVLNode {
	KeyType key;                 // Key
	int bf;                      // Balance factor 
	AVLNode<KeyType> *left;      // Link to left child
	AVLNode<KeyType> *right;     // Link to right child
};


template <class KeyType>
class AVLTree {
public:
	AVLTree();
	~AVLTree();

	void Insert( KeyType key );
	void Delete( KeyType key );
	void Preorder();
	void Inorder();
	void Postorder();
	int  NumKeys();

protected:
	AVLNode<KeyType> *root;
	AVLNode<KeyType> *current;
	AVLNode<KeyType> *parent;
	AVLNode<KeyType> *ptr;
	AVLNode<KeyType> *prev;

	AVLNode<KeyType> *pivot;
	AVLNode<KeyType> *pivot_prev;

	int n;

	void bf_count( AVLNode<KeyType> *trees );     // Calculate balance factors for all nodes
	int  height_count( AVLNode<KeyType> *trees ); // Calculate height for subtrees
	AVLNode<KeyType> *pivot_find( void );         // Find pivot that requires rotation
	int  type_find( void );                       // Determine type of rotation
	void type_ll( void );                         // LL type
	void type_rr( void );                         // RR type
	void type_lr( void );                         // LR type
	void type_rl( void );                         // RL type

	void Preorder( AVLNode<KeyType> *ptr );
	void Inorder( AVLNode<KeyType> *ptr );
	void Postorder( AVLNode<KeyType> *ptr );
};


template <class KeyType>
AVLTree<KeyType>::AVLTree()
{
	root = NULL;
	prev = NULL;
	pivot_prev = NULL;
	n = 0;
}


template <class KeyType>
AVLTree<KeyType>::~AVLTree()
{
	if ( root != NULL ) delete root;
}


template <class KeyType>
void AVLTree<KeyType>::Insert( KeyType key )  // Based on Binary Search Tree Insertion + Balance
{
	ptr = new AVLNode<KeyType>;               // Create a new AVL Tree node
	ptr->key = key;                           // Assign key
	ptr->left = ptr->right = NULL;            // No left & right children

	if ( root == NULL )                       // Empty AVL Tree, the new node becomes the root
		root = ptr;
	else
	{
		current = root;
		while ( current != NULL )             // Search for the location to insert 
		{
			parent = current;
			if ( key < current->key )
				current = current->left;
			else
				current = current->right;
		}

		if ( key < parent->key )              // Insert the new node 
			parent->left = ptr;
		else
			parent->right = ptr;
	}

	bf_count( root );                         // Calculate balance factor
	pivot = pivot_find();                     // Find pivot for rotation
	if ( pivot != NULL )                      // Pivot exists (require rotation)
	{
		int op = type_find();                 // Determine type of rotation
		switch ( op )
		{
		case 11: type_ll();                   // LL rotation
			break;
		case 22: type_rr();                   // RR rotation
			break;
		case 12: type_lr();                   // LR rotation
			break;
		case 21: type_rl();                   // RL rotation
			break;
		}
	}
	bf_count( root );                         // Recalculate node balance factor            
	n++;
}


template <class KeyType>
void AVLTree<KeyType>::Delete( KeyType key )  // Based on Binary Search Tree Deletion + Balance
{
	AVLNode<KeyType> *del_node;

	parent = NULL;                          
	current = root;
	while ( current != NULL )                 // Search for the deleted node
	{
		if ( key == current->key )            // The deleted node is found
			break;
		parent = current;                     // Also determine the deleted node's parent
		if ( key < current->key )               
			current = current->left;
		else                               
			current = current->right;
	}

	if ( current == NULL )                                           // Deleted node is not found
		cout << "Deleted node not found in the AVL Tree.\n";
	else
	{
		if ( current->left == NULL && current->right == NULL )       // No left & right subtree (leaf node)
		{
			if ( current == root )                                   // The deleted node is the root with no children
				root = NULL;
			else                                                     // Otherwise, the deleted node is a leaf node
			{
				if ( current->key < parent->key )                    // Determine the deleted node's parent 
					parent->left = NULL;
				else                                            
					parent->right = NULL;
			}
			del_node = current;
			ptr = parent;
		}
		else if ( current->left == NULL && current->right != NULL )  // Right subtree is not empty
		{
			if ( current == root )                                   // Deleted node is the root
				root = current->right;                               // Its right child becomes the root
			else                                                     // Otherwise
			{
				if ( current->key < parent->key )                    // Link its right subtree 
					parent->left = current->right;
				else
					parent->right = current->right;
			}
			del_node = current;
			ptr = parent;
		}
		else if ( current->left != NULL && current->right == NULL )  // Left subtree is not empty
		{
			if ( current == root )                                   // Deleted node is the root
				root = current->left;
			else		                                             // Otherwise
			{
				if ( current->key < parent->key )                    // Link its left subtree
					parent->left = current->left;
				else
					parent->right = current->left;
			}
			del_node = current;
			ptr = parent;
		}
		else                                                         // Both left & right subtrees are not empty
		{
			prev = current;
			ptr  = current->right;                                   // Search for replacing node
			while ( ptr->left != NULL )                              // Use smallest node in right subtree (default)
			{                                                        // or largest node in left subtree (not used here)
				prev = ptr;
				ptr = ptr->left;
			}

			if ( ptr->key < prev->key )                              // Link the deleted node's parent to its children
				prev->left = ptr->right;
			else
				prev->right = ptr->right;

			current->key = ptr->key;                                 // Replace key
			del_node = ptr;
			ptr = prev;
		}
		delete del_node;
		n--;
	}
	
	bf_count( root );
	pivot = pivot_find();
	if ( pivot != NULL )
	{
		int op = type_find();
		switch ( op )
		{
		case 11: type_ll();
			break;
		case 22: type_rr();
			break;
		case 12: type_lr();
			break;
		case 21: type_rl();
			break;
		}
	} 
	bf_count( root );  
}


template <class KeyType>
void AVLTree<KeyType>::bf_count( AVLNode<KeyType> *trees )
{
	if ( trees != NULL )     
	{
		bf_count( trees->left );
		bf_count( trees->right );
		trees->bf = height_count( trees->left ) - height_count( trees->right );
	}
}


template <class KeyType>
int AVLTree<KeyType>::height_count(AVLNode<KeyType> *trees)
{
	if ( trees == NULL )     
		return 0;
	else if ( trees->left == NULL && trees->right == NULL )
		return 1;
	else
		return 1 + ( height_count( trees->left ) > height_count( trees->right ) ?
		             height_count( trees->left ) : height_count( trees->right ) );
}


template <class KeyType>
AVLNode<KeyType> *AVLTree<KeyType>::pivot_find( void )
{
	if ( root == ptr ) prev = NULL;
	else               prev = parent;
	
	pivot = NULL;
	current = root;
	while ( current != ptr ) 
	{
		if ( current->bf < -1 || current->bf > 1 ) // Not balanced
		{
			pivot = current;
			if ( pivot != root )
				pivot_prev = prev;
		}
		if ( ptr->key < current->key )
		{
			prev = current;
			current = current->left;
		}
		else
		{
			prev = current;
			current = current->right;
		}
	}
	return pivot;
}


template <class KeyType>
int AVLTree<KeyType>::type_find( void )
{
	int i, op_r = 0;

	current = pivot;
	for ( i = 0 ; i < 2 ; i++ )
	{
		if ( ptr->key < current->key )
		{
			current = current->left;
			if ( op_r == 0 ) op_r += 10;
			else             op_r++;
		}
		else
		{
			current = current->right;
			if ( op_r == 0 ) op_r += 20;
			else             op_r += 2;
		}
	}
	// Return 11, 22, 12, 21 for LL, RR, LR, RL types
	return op_r;
}


template <class KeyType>
void AVLTree<KeyType>::type_ll( void )
{
	AVLNode<KeyType> *pivot_next, *temp;
	pivot_next = pivot->left;
	temp = pivot_next->right;
	pivot_next->right = pivot;
	pivot->left = temp;
	if ( pivot == root )
		root = pivot_next;
	else if ( pivot_prev->left == pivot )
		pivot_prev->left = pivot_next;
	else
		pivot_prev->right = pivot_next;
}


template <class KeyType>
void AVLTree<KeyType>::type_rr( void )
{
	AVLNode<KeyType> *pivot_next, *temp;
	pivot_next = pivot->right;
	temp = pivot_next->left;
	pivot_next->left = pivot;
	pivot->right = temp;
	if ( pivot == root )
		root = pivot_next;
	else if ( pivot_prev->left == pivot )
		pivot_prev->left = pivot_next;
	else
		pivot_prev->right = pivot_next;
}


template <class KeyType>
void AVLTree<KeyType>::type_lr( void )
{
	AVLNode<KeyType> *pivot_next, *temp;
	pivot_next = pivot->left;
	temp = pivot_next->right;
	pivot->left = temp->right;
	pivot_next->right = temp->left;
	temp->left = pivot_next;
	temp->right = pivot;
	if ( pivot == root )
		root = temp;
	else if ( pivot_prev->left == pivot )
		pivot_prev->left = temp;
	else
		pivot_prev->right = temp;
}


template <class KeyType>
void AVLTree<KeyType>::type_rl( void )
{
	AVLNode<KeyType> *pivot_next, *temp;
	pivot_next = pivot->right;
	temp = pivot_next->left;
	pivot->right = temp->left;
	pivot_next->left = temp->right;
	temp->right = pivot_next;
	temp->left = pivot;
	if ( pivot == root )
		root = temp;
	else if ( pivot_prev->left == pivot )
		pivot_prev->left = temp;
	else
		pivot_prev->right = temp;
}


template <class KeyType>
void AVLTree<KeyType>::Preorder()
{
	if ( root == NULL )
		cout << "No Keys in AVL Tree.\n";
	else
	{
		cout << "AVL Tree (Preorder):\n";
		Preorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void AVLTree<KeyType>::Preorder( AVLNode<KeyType> *ptr )
{
	if ( ptr != NULL )                      
	{
		cout << ptr->key << " ";
		Preorder(ptr->left);
		Preorder(ptr->right);
	}
}


template <class KeyType>
void AVLTree<KeyType>::Inorder()
{
	if ( root == NULL )
		cout << "No Keys in AVL Tree.\n";
	else
	{
		cout << "AVL Tree (Inorder):\n";       
		Inorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void AVLTree<KeyType>::Inorder( AVLNode<KeyType> *ptr )
{
	if ( ptr != NULL )
	{
		Inorder( ptr->left );
		cout << ptr->key << " ";
		Inorder( ptr->right );
	}
}


template <class KeyType>
void AVLTree<KeyType>::Postorder()
{
	if ( root == NULL )
		cout << "No Keys in AVL Tree.\n";
	else
	{
		cout << "AVL Tree (Postorder):\n";
		Postorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void AVLTree<KeyType>::Postorder( AVLNode<KeyType> *ptr )
{
	if ( ptr != NULL )
	{
		Postorder( ptr->left );
		Postorder( ptr->right );
        cout << ptr->key << " ";
	}
}


template <class KeyType>
int AVLTree<KeyType>::NumKeys()
{
	return n;
}


//
//  Red-Black Tree (RB Tree)
//  
//  Properties:
//     1. Every node is either red or black
//     2. The root is black
//     3. Every leaf (NIL) is black
//     4. If a node is red, then both its children are black
//     5. For each node, all paths from the node to descendant leaves contain 
//        the same number of black nodes.
//
//  Operations:
//     Insert    (Insert to RB Tree)
//     Preorder  (Preorder traversal)
//     Inorder   (Inorder traversal)
//     Postorder (Postorder traversal)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
struct RBNode {
	KeyType key;               // Key
	bool red;                  // True (red) or False (black)
	RBNode<KeyType> *parent;   // Link to parent
	RBNode<KeyType> *left;     // Link to left child
	RBNode<KeyType> *right;    // Link to right child
};


template <class KeyType>
class RBTree {
public:
	RBTree();
	~RBTree();

	void Insert( KeyType key );
	void Preorder();
	void Inorder();
	void Postorder();

protected:
	RBNode<KeyType> *nil;
	RBNode<KeyType> *root;
	RBNode<KeyType> *ptr;
	RBNode<KeyType> *x, *y;

	// Auxiliary
	void Preorder( RBNode<KeyType> *ptr );   
	void Inorder( RBNode<KeyType> *ptr );
	void Postorder( RBNode<KeyType> *ptr );

	void RB_Insert_Fixup( RBNode<KeyType> *z );
	void Left_Rotate( RBNode<KeyType> *z );
	void Right_Rotate( RBNode<KeyType> *z );

	void RBTree<KeyType>::RB_Transplant( RBNode<KeyType> *u, RBNode<KeyType> *v );
};


template <class KeyType>
RBTree<KeyType>::RBTree()
{
	nil = new RBNode<KeyType>;     // Create a nil node
	nil->key = 0;
	nil->red = false;
	nil->parent = nil->left = nil->right = NULL;
	root = NULL;
};


template <class KeyType>
RBTree<KeyType>::~RBTree()
{
	delete nil;
};


template <class KeyType>
void RBTree<KeyType>::Insert( KeyType key )
{
	ptr = new RBNode<KeyType>;             
	ptr->key = key;                  
	ptr->red = true;
	ptr->parent = ptr->left = ptr->right = nil;          

	if ( root == NULL )       // Empty RB Tree
	{
		root = ptr;
		ptr->red = false;     // Root is black
	}
	else
	{
		y = nil;
		x = root;
		while ( x != nil )    
		{
			y = x;
			if ( key < x->key )
				x = x->left;
			else
				x = x->right;
		}

 		ptr->parent = y;

		if ( y == nil )
		{
			root = ptr;
		}
		else
		{
			if ( key < y->key )
				y->left = ptr;
			else
				y->right = ptr;
		}

		ptr->left = nil;
		ptr->right = nil;
		ptr->red = true;

		RB_Insert_Fixup( ptr );
	}
};


template <class KeyType>
void RBTree<KeyType>::RB_Insert_Fixup( RBNode<KeyType> *z )
{
	while ( z->parent->red )
	{
		if ( z->parent == z->parent->parent->left )
		{
			y = z->parent->parent->right;
			if ( y->red )
			{
				z->parent->red = false;          // Case 1
				y->red = false;
				z->parent->parent->red = true;
				z = z->parent->parent;
			}
			else 
			{
				if ( z == z->parent->right ) 
				{
					z = z->parent;               // Case 2       
					Left_Rotate( z );
				}
				z->parent->red = false;          // Case 3
				z->parent->parent->red = true;
				Right_Rotate( z->parent->parent );
			}
		}
		else
		{
			y = z->parent->parent->left;
			if ( y->red )
			{
				z->parent->red = false;          // Case 1
				y->red = false;
				z->parent->parent->red = true;
				z = z->parent->parent;
			}
			else 
			{
				if ( z == z->parent->left )
				{
					z = z->parent;               // Case 2
					Right_Rotate( z );
				}
				z->parent->red = false;          // Case 3 
				z->parent->parent->red = true;
				Left_Rotate( z->parent->parent );
			}			
		}
	}
	root->red = false; 
}


template <class KeyType>
void RBTree<KeyType>::Left_Rotate( RBNode<KeyType> *z )
{
	y = z->right;
	z->right = y->left;
	if ( y->left != nil )
		y->left->parent = z;
	y->parent = z->parent;
	if ( z->parent == nil )
		root = y;
	else
	{
		if ( z == z->parent->left )
			z->parent->left = y;
		else
			z->parent->right = y;
	}
	y->left = z;
	z->parent = y;
}


template <class KeyType>
void RBTree<KeyType>::Right_Rotate( RBNode<KeyType> *z )
{
	x = z->left;
	z->left = x->right;
	if ( x->right != nil )
		x->right->parent = z;
	x->parent = z->parent;
	if ( z->parent == nil )
		root = x;
	else
	{
		if ( z == z->parent->left )
			z->parent->left = x;
		else
			z->parent->right = x;
	}
	x->right = z;
	z->parent = x;
}


template <class KeyType>
void RBTree<KeyType>::RB_Transplant( RBNode<KeyType> *u, RBNode<KeyType> *v )
{
	if ( u->parent == nil )
	{
		root = v;
	}
	else if ( u == u->parent->left )
	{
		u->parent->left = v;
	}
	else
	{
		u->parent->right = v;
	}

	v->parent = u->parent;
}


template <class KeyType>
void RBTree<KeyType>::Preorder()
{
	if ( root == NULL )
		cout << "No Keys in RB Tree.\n";
	else
	{
		cout << "RB Tree (Preorder):\n";
		Preorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void RBTree<KeyType>::Preorder( RBNode<KeyType> *ptr )
{
	if ( ptr != nil )                      
	{
		if ( ptr->red ) cout << ptr->key << "(R) ";
		else            cout << ptr->key << "(B) ";
		Preorder( ptr->left );
		Preorder( ptr->right );
	}
}


template <class KeyType>
void RBTree<KeyType>::Inorder()
{
	if ( root == NULL )
		cout << "No Keys in RB Tree.\n";
	else
	{
		cout << "RB Tree (Inorder):\n";       
		Inorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void RBTree<KeyType>::Inorder( RBNode<KeyType> *ptr )
{
	if ( ptr != nil )
	{
		Inorder( ptr->left );
		if ( ptr->red ) cout << ptr->key << "(R) ";
		else            cout << ptr->key << "(B) ";
		Inorder( ptr->right );
	}
}


template <class KeyType>
void RBTree<KeyType>::Postorder()
{
	if ( root == NULL )
		cout << "No Keys in RB Tree.\n";
	else
	{
		cout << "RB Tree (Postorder):\n";
		Postorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void RBTree<KeyType>::Postorder( RBNode<KeyType> *ptr )
{
	if ( ptr != nil )
	{
		Postorder( ptr->left );
		Postorder( ptr->right );
		if ( ptr->red ) cout << ptr->key << "(R) ";
		else            cout << ptr->key << "(B) ";
	}
}


//
//  B-Tree (2-3 Tree)
//  
//  Operations:
//     Insert    (Insert to 2-3 Tree)
//     Preorder  (Preorder traversal)
//     Postorder (Postorder traversal)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
template <class KeyType>
struct TwoNode {
	int     n;                   // Number of keys (1 or 2)
	KeyType keyL, keyR;          // Keys 
	
	TwoNode<KeyType> *parent;    // Link to parent
	TwoNode<KeyType> *left;      // Link to left child
	TwoNode<KeyType> *middle;    // Link to middle child 
	TwoNode<KeyType> *right;     // Link to right child
};


template <class KeyType>
class BTree23 {
public:
	BTree23();  
	~BTree23();
 
	bool Search(KeyType x);
	void Insert( KeyType x );      
	void Preorder();      
	void Postorder();            

private:
    TwoNode<KeyType> *root;
	TwoNode<KeyType> *access( KeyType x, TwoNode<KeyType> *node );

	void putindata( KeyType key, TwoNode<KeyType> *xr, TwoNode<KeyType> *p, int k );
	void split( KeyType x_key, TwoNode<KeyType> *xr, TwoNode<KeyType> *p, int k, KeyType *y_key, TwoNode<KeyType> **yr );
	int  find_k( KeyType key, TwoNode<KeyType> *p );
	bool find_leaf( KeyType key, TwoNode<KeyType> *p, KeyType *x_key, TwoNode<KeyType> **xr );

	void Preorder( TwoNode<KeyType> *p );      
	void Postorder( TwoNode<KeyType> *p );            
};


template <class KeyType>
BTree23<KeyType>::BTree23()
{
	root = NULL;
}


template <class KeyType>
BTree23<KeyType>::~BTree23()
{
	// No need to do anything
}


template <class KeyType>
bool BTree23<KeyType>::Search(KeyType x)
{
    if ( root == NULL )
        cout << "No data in 2-3 Tree." << endl;
    else
    {
        cur = root;                     
        while ( cur != NULL )
        {
            if ( cur->n == 1 )             
            {
                if ( x < cur->keyL )       
                    cur = cur->left;
                else if ( x > cur->keyL )  
                    cur = cur->right;
                else
                    return true;        
            }
            else                        
            {
                if ( x < cur->keyL )       
                    cur = cur->left;
                else if ( ( x > cur->keyL ) && ( x < cur->keyR ) ) 
                    cur = cur->middle;
                else if ( x > cur->keyR )  
                    cur = cur->right;
                else
                    return true;        
            }
        }
        return false;                  
    }
} 


template <class KeyType>
void BTree23<KeyType>::Insert( KeyType x )
{
	root = access( x, root );
}


template <class KeyType>
void BTree23<KeyType>::putindata( KeyType x_key, TwoNode<KeyType> *xr, TwoNode<KeyType> *p, int k ) 
{
	int i;
	i = p->n;

	if ( k == 0 )             
	{
		p->keyR = p->keyL;
		p->right = p->middle;
		p->keyL = x_key;
		p->middle = xr;
		p->n++;
	}
	else if ( k == 1 ) 
	{
		p->keyR = x_key;
	    p->right = xr;
	    p->n++;
	}
}


template <class KeyType>
void BTree23<KeyType>::split( KeyType x_key, TwoNode<KeyType> *xr, TwoNode<KeyType> *p, int k, KeyType *y_key, TwoNode<KeyType> **new_node )
{
	int mid;
	if ( k <= 1 )
		mid = 1;
	else
		mid = 2;

	*new_node = new TwoNode<KeyType>;
	if ( k <= 1 )
	{
		(*new_node)->keyL = p->keyR;
		(*new_node)->middle = p->right;
	}

	(*new_node)->n = 2 - mid; 
	                         
	p->n = mid;
	if ( k <= 1 )
		putindata( x_key, xr, p, k );
	else
		putindata( x_key, xr, *new_node, k - mid );

	*y_key=p->keyR;
	(*new_node)->left = p->right;
	p->n--;
}


template <class KeyType>
int BTree23<KeyType>::find_k( KeyType key, TwoNode<KeyType> *p ) 
{
	int k = 0;

	if ( p->n == 1 )
	{
		if ( key < p->keyL )
			k = 0;
		else
			k = 1;
	}
	else if ( p->n == 2 )
	{
		if ( key < p->keyL )
			k = 0;
	    else if ( key<p->keyR )
			k = 1;
	    else 
			k = 2;
	}

	return k;
}


template <class KeyType>
bool BTree23<KeyType>::find_leaf( KeyType key, TwoNode<KeyType> *p, KeyType *x_key, TwoNode<KeyType> **xr )
{
	int k;
	TwoNode<KeyType> *ptr;
	if ( p == NULL )
	{
		*x_key = key;
		*xr = NULL;
		return true;
	}
	else
	{
		k = find_k( key, p );
 
		if ( k == 0 )
			ptr=p->left;
		else if ( k == 1 )
			ptr=p->middle;
		else
			ptr=p->right;

		if ( find_leaf( key, ptr, x_key, xr ) )
		{
			if ( p->n < 2 )
			{
				putindata( *x_key, *xr, p, k );
				return false;
			}
			else{
				split( *x_key, *xr, p, k, x_key, xr);
				return true;                    
			}
		}
		else
			return false;
	}
}


template <class KeyType>
TwoNode<KeyType> *BTree23<KeyType>::access( KeyType x,TwoNode<KeyType> *node )
{
	KeyType x_x;
	bool up=false;
	TwoNode<KeyType> *xr, *p;

	up = find_leaf( x, node, &x_x, &xr );
	if ( up ) 
	{
		p =new TwoNode<KeyType>;
		p->parent=NULL;
		p->left=NULL;
		p->middle=NULL;
		p->right=NULL;
		p->n=1;
		p->keyL=x_x;
		p->left=root;
		p->middle=xr;
		return p;
	}
	return node;
}


template <class KeyType>
void BTree23<KeyType>::Preorder()
{    
	if ( root == NULL )
		cout << "No Keys in 2-3 Tree" << "\n";
	else
	{
		cout << "2-3 Tree (Preorder):" << "\n";
		Preorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void BTree23<KeyType>::Preorder( TwoNode<KeyType> *p )
{
	if ( p != NULL )
	{
		if ( p->n > 0 )
		{
			if ( p->n == 1 )
			{
				cout << "(" << p->keyL << ")";
				Preorder( p->left );
				Preorder( p->middle );		
			}
			else if ( p->n == 2 ) 
			{
				cout << "(" << p->keyL << ",";
				cout << p->keyR << ")";
				Preorder( p->left );
				Preorder( p->middle );
				Preorder( p->right );
			}
		}
	}
}


template <class KeyType>
void BTree23<KeyType>::Postorder()
{    
	if ( root == NULL )
		cout << "No data in 2-3Tree" << "\n";
	else
	{
		cout << "2-3 Tree (Postorder):" << "\n";
		Postorder( root );
		cout << "\n";
	}
}


template <class KeyType>
void BTree23<KeyType>::Postorder( TwoNode<KeyType> *p )
{    
	if ( p != NULL )
	{
		if ( p->n > 0 )
		{
			if ( p->n == 1)
			{
				Postorder( p->left );
				Postorder( p->middle );
				cout << "(" << p->keyL << ")";
			}
			else if ( p->n == 2 )
			{
				Postorder( p->left );
				Postorder( p->middle );
				Postorder( p->right );
				cout << "(" << p->keyL << ",";
				cout << p->keyR << ")";
			}
		}
	}
}


//
//  Hash Table
//
//  Method = 1  Linear Probing
//         = 2  Quadractic Probing
//         = 3  Double Hashing
//
//  Operations:
//     Insert    (Insert to Hash Table)
//     Display   (Display Hash Table)
//
//  Note: Currently Supoort Integer Keys Only
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
class HashTable {
public:
	int method; 

	HashTable( int n );
	~HashTable();

	void Insert( int key );
	void Display();

protected:
	int *table;
	int table_size;
};


HashTable::HashTable( int n )
{
	table = new int[n];
	table_size = n;

	for ( int i = 0 ; i < n ; i++ ) 
		table[i] = 0;
}


HashTable::~HashTable()
{
	delete [] table;
}


void HashTable::Insert( int key )
{
	int i, j, flag;
	int h1, h2;

	j = (int) fmod( (double) key, (double) table_size );			
     
	if ( table[j] == 0 ) 
	{
		// If empty, simply insert
		table[j] = key;
	}
	else
	{
		// If not empty, resolve collision
		flag = 1;
		i = 0;
		while ( flag )
		{
			h1 = (int) fmod( (double) key, (double) table_size );

			if ( method == 1 )       // Linear Probing
			{
				j = (int) fmod( (double) ( h1 + i ), (double) table_size );
			}
			else if ( method == 2 )  // Quadratic Probing
			{
				cout << "i = " << i << endl;
				// Use c1 = 0, c2 = 1
				j = (int) fmod( (double) ( h1 + i + 3 * i * i ), (double) table_size );
			}
			else                     // Double Hashing
			{
				h2 = 1 + (int) fmod( (double) key, (double) ( table_size - 1 ) );
				j = (int) fmod( (double) ( h1 + i * h2 ), (double) table_size );
			}

			if ( table[j] == 0 )
			{
				table[j] = key;
				flag = 0;
			}
			else
			{
				i++;
			}
			if ( i == table_size )
			{
				cout << "Hash Table Overflow!\n";
				flag = 0;
			}
		}	
	}
}


void HashTable::Display()
{
	cout << "Hash Table:\n";
	for ( int i = 0 ; i < table_size ; i++ )
		cout << i << "  " << table[i] << endl;
}


//
//  Disjoint Set
//
//  By default, disjoint sets of {1}, {2}, ...{n} are initially formed, 
//  where n is the number of sets and each set contains 1 key only. 
//  
//  Operations:
//     Find    (Find root or set representative of a key) 
//     Union   (Union two sets)
//     Display (Display Disjoint Set)
//
//  Reference: Introduction to Algorithms (2nd Edition) Cormen
//
class DisjointSet {
public:
	DisjointSet();          
	DisjointSet( int n );     
	~DisjointSet();

	int Find( int key );
	void Union( int a, int b );    
	void Display();

protected:
	int *set;
	int set_size;
};


DisjointSet::DisjointSet()
{
	set = new int[101];
	set_size = 100;
	for ( int i = 1 ; i <= set_size ; i++ )
		set[i] = 0;
};


DisjointSet::DisjointSet( int n )
{
	set = new int[n+1];
	set_size = n;
	for ( int i = 1 ; i <= set_size ; i++ )
		set[i] = 0;
};


DisjointSet::~DisjointSet()
{
	delete [] set;
};


int DisjointSet::Find( int key )
{
	if ( set[key] == 0 ) return key;
	else                 return Find( set[key] );
}


void DisjointSet::Union( int a, int b )
{
	if ( Find( a ) == Find( b ) )
	{
		if ( a < b ) set[b] = a;
		else         set[a] = b;
	}
	else
	{
		if ( Find( a ) < Find( b ) ) set[b] = a;
		else 		                 set[a] = b;
	}
}


void DisjointSet::Display()
{
	int i, j;

	cout << "Disjoint Set:\n";   
	for ( i = 1 ; i <= set_size ; i++ )
	{
		if ( set[i] == 0 )
		{
			cout << "{" << i;
			for ( j = i + 1 ; j <= set_size ; j++ )
			{
				if ( Find(j) == i )
					cout << "," << j;
			}
			cout << "}";
		}
	}
	cout << endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Dynamic Programming
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Rod-Cutting Problem
//  Given a rod of length n inches and a table of prices pi for i = 1, 2,..., n, determine 
//  the maximum revenue rn obtainable by cutting up the rod and selling pieces.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Rod_Cutting( int *p, int n )
{
	int *r, *s;

	r = new int[ n + 1 ];    // Maximum revenue of rod length
	s = new int[ n + 1 ];    // Optimal size of first piece 


	r[0] = s[0] = 0;
	for ( int j = 1 ; j <= n ; j++ )
	{
		int q = -30000;
		for ( int i = 1 ; i <= j ; i++ )
		{
			if ( q < p[i] + r[j-i] )
			{
				q = p[i] + r[j-i];
				s[j] = i;
			}
		}
		r[j] = q;
	}

	cout << "Rod-Cutting Problem" << endl;

	for ( int i = 0 ; i <= n ; i++ )
		cout << setw( 4 ) << i;
	cout << endl;
	cout << "-----------------------------------------------------------------------------" << endl;
	for ( int i = 0 ; i <= n ; i++ )
		cout << setw( 4 ) << r[i];
	cout << endl;
	for ( int i = 0 ; i <= n ; i++ )
		cout << setw( 4 ) << s[i];
	cout << endl;

	cout << "Cuts = ";
	while ( n > 0 )
	{
		cout << s[n] << " ";
		n = n - s[n];
	}
	cout << endl;
	
	delete [] r;
	delete [] s;
}


//
//  Assembly-Line Scheduling
//  Given a manufacturing problem, find the fastest way through a factory. 
//  There are two assembly lines, each with n stations. 
//  The jth station on line i is Si,j.
//  
//  a[i][j]: Assembly time at station Si,j 
//  t[i][j]: Transfer time after station Si,j
//     e[i]: Entry time of assembly line i
//     x[i]: Exit time of assembly line i
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Assembly_Line_Scheduling( int **a, int **t, int *e, int *x, int n )
{
	int i, j, fs, ls;
	
	vector<int> f1(n+1), f2(n+1), l1(n+1), l2(n+1);
	
	f1[1] = e[1] + a[1][1];
	f2[1] = e[2] + a[2][1];

	for ( j = 2 ; j <= n ; j++ )
	{
		if ( f1[j-1] + a[1][j] <= f2[j-1] + t[2][j-1] + a[1][j] )
		{
			f1[j] = f1[j-1] + a[1][j];
			l1[j] = 1;
		}
		else
		{
			f1[j] = f2[j-1] + t[2][j-1] + a[1][j];
			l1[j] = 2;
		}
		if ( f2[j-1] + a[2][j] <= f1[j-1] + t[1][j-1] + a[2][j] )
		{
			f2[j] = f2[j-1] + a[2][j];
			l2[j] = 2;
		}
		else
		{
			f2[j] = f1[j-1] + t[1][j-1] + a[2][j];
			l2[j] = 1;
		}
	}
	if ( f1[n] + x[1] <= f2[n] + x[2] )
	{
		fs = f1[n] + x[1];
		ls = 1;
	}
	else
	{
		fs = f2[n] + x[2];
		ls = 2;
	}
	
	cout << "Cost Table:\n";
	for ( i = 1 ; i <= n ; i++ )
		cout << f1[i] << " ";
	cout << endl;
	for ( i = 1 ; i <= n ; i++ )
		cout << f2[i] << " ";
	cout << endl;

	cout << "Line Table:\n";
	for ( i = 2 ; i <= n ; i++ )
		cout << l1[i] << " ";
	cout << endl;
	for ( i = 2 ; i <= n ; i++ )
		cout << l2[i] << " ";
	cout << endl;

	cout << "Path:\n";
	i = ls;
	cout << "line " << i << " station " << n << endl;
	for ( j = n ; j >= 2 ; j-- )
	{
		if ( i == 1 ) i = l1[j];
		if ( i == 2 ) i = l2[j];
		cout << "line " << i << " station " << j-1 << endl;
	}
}


//
//  Matrix-Chain Multiplication
//  Given a sequence of A1, A2, ...An of n matrices, find the optimal parenthesization. 
//  The matrix Ai has dimensions p[i-1] x p[i] for i = 1, 2, ..n.
//
//  Input:  The dimensions p[i], i = 0..n & n: number of matrices
//  Output: The optimal parenthesization.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Print_Optimal_Parens( int **s, int i, int j )
{
	if ( i == j )
		cout << "A" << i;
	else
	{
		cout << "(";
		Print_Optimal_Parens( s, i, s[i][j] );
		Print_Optimal_Parens( s, s[i][j] + 1, j );
		cout << ")";
	}
}


void Matrix_Chain_Multiplication( int *p, int n )
{
	unsigned int **m, q;
	int **s;
	int i, j, k, l;

	m = Matrix<unsigned int> ( n + 1, n + 1 );
	s = Matrix<int> ( n + 1, n + 1 );

	for ( i = 1 ; i <= n ; i++ )
	{
		for ( j = 1 ; j <= n ; j++ )
		{
			m[i][i] = 0;
			s[i][j] = 0;
		}
	}
	for ( l = 2 ; l <= n ; l++ )
	{
		for ( i = 1 ; i <= n - l + 1 ; i++ )
		{
			j = i + l - 1;
			m[i][j] = UINT_MAX;
			for ( k = i ; k <= j - 1 ; k++ )
			{
				q = m[i][k] + m[k+1][j] + p[i-1] * p[k] * p[j];
				if ( q < m[i][j] )
				{
					m[i][j] = q;
					s[i][j] = k;
				}
			}
		}
	}

	// Results
	cout << "The m table (number of multiplications)\n";
	for ( i = 1 ; i <= n ; i++ )
	{
		for ( j = i ; j <= n ; j++ )
			cout << m[i][j] << " ";
		cout << endl;
	}

	cout << "The s table (split)\n";
	for ( i = 1 ; i < n ; i++ )
	{
		for ( j = i + 1 ; j <= n ; j++ )
			cout << s[i][j] << " ";
		cout << endl;
	}
	
	cout << "Optimal parenthesization\n";
	Print_Optimal_Parens( s, 1, n );
	cout << endl;

	FreeMatrix<unsigned int> ( m, n + 1, n + 1 );
	FreeMatrix<int> ( s, n + 1, n + 1 );
}


//
//  Longest Common Subsequence (LCS)
//  Given two strings X & Y, solve the problem of finding the LCS.
//  
//  Input:  X, Y 
//  Output: result (LCS)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
// 
void LCS( char *X, char *Y, char *result )
{
	enum direction { UPPERLEFT, UP, LEFT };
	int  i, j, m, n;
	int  **b, **c;

	m = (int) strlen( X );
	n = (int) strlen( Y );

	c = Matrix<int> ( m + 1, n + 1 );
	b = Matrix<int> ( m + 1, n + 1 );
	
	for ( i = 1 ; i <= m ; i++ )
		c[i][0] = 0;

	for ( j = 0 ; j <= n ; j++ )
		c[0][j] = 0;
	
	for ( i = 1 ; i <= m ; i++ )  // The c & b tables
	{
		for ( j = 1 ; j <= n ; j++ )
		{
			if ( X[i-1] == Y[j-1] )
			{
				c[i][j] = c[i-1][j-1] + 1;
				b[i][j] = UPPERLEFT;   
			}
			else
			{
				if ( c[i-1][j] >= c[i][j-1] )
				{
					c[i][j] = c[i-1][j];
					b[i][j] = UP;
				}
				else
				{
					c[i][j] = c[i][j-1];
					b[i][j] = LEFT;
				}
			}
		}
	}

	i = m; 
	j = n;
	result[c[i][j]] = '\0';
	while ( i !=0 && j != 0 )  // Find the LCS
	{
		if ( b[i][j] == UPPERLEFT )
		{
			result[c[i][j]-1] = X[i-1];
			i--;  
			j--;
		}
		else if ( b[i][j] == UP )
			i--;
		else
			j--;
	}

	FreeMatrix<int> ( c, m + 1, n + 1 );
	FreeMatrix<int> ( b, m + 1, n + 1 );
}


//
//  Optimal BST
//  Find the optimal binary search tree (BST) given the probabilities of n keys (nodes)
//  such that the number of nodes visited is minimized.
//  
//  p[i]: probability of key i
//  q[i]: probability of dummy key i
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Optimal_BST( double *p, double *q, int n )
{
	int    i, j, l, r;
	double **e, **w, t;
	int    **root;
	
	e = Matrix<double> ( n + 2, n + 2 );
	w = Matrix<double> ( n + 2, n + 2 );
	root = Matrix<int> ( n + 2, n + 2 );

	for ( i = 1 ; i <= n + 1 ; i++ )
	{
		e[i][i-1] = q[i-1];
		w[i][i-1] = q[i-1];
	}

	for ( l = 1 ; l <= n ; l++ )
	{
		for ( i = 1 ; i <= n - l + 1 ; i++ )
		{
			j = i + l - 1;
			e[i][j] = 1.0e300;  // big number
			w[i][j] = w[i][j-1] + p[j] + q[j];
			
			for ( r = i ; r <= j ; r++ )
			{
				t = e[i][r-1] + e[r+1][j] + w[i][j];
				if ( t < e[i][j] )
				{
					e[i][j] = t;
					root[i][j] = r;
				}
			}
		}
	}

	cout << "The e table\n";
	for ( i = 1 ; i <= n + 1 ; i++ )
	{
		for ( j = i - 1 ; j <= n ; j++ )
			cout << e[i][j] << " ";
		cout << endl;
	}

	cout << "The w table\n";
	for(i=1;i<=n+1;i++)
	{
		for(j=i-1;j<=n;j++)
			cout << w[i][j] << " ";
		cout << endl;
	}

	cout << "The root table\n";
	for ( i = 1 ; i <= n ; i++ )
	{
		for ( j = i ; j <= n ; j++ )
			cout << root[i][j] << " ";
		cout << endl;
	}

	FreeMatrix<double> ( e, n + 2, n + 2 );
	FreeMatrix<double> ( w, n + 2, n + 2 );
	FreeMatrix<int> ( root, n + 2, n + 2 );
}


//
//  Edit Distance
//  Given two strings X & Y, find the edit distance, i.e.,
//  the number of editing required to change from X to Y.
//  Possible editing: insert/delete/substitute
//
//  Input:  X, Y 
//  Output: edit distance
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
// 
#define INSERT_COST      1
#define DELETE_COST      1
#define SUBSTITUTE_COST  1

#define EDIT_NOCHANGE    0
#define EDIT_INSERT      1
#define EDIT_DELETE      2
#define EDIT_SUBSTITUTE  3

int Edit_Distance(char *X, char *Y)
{
	int i, j, m, n;
	int **d, **e, cost, edit_distance;

	m = (int) strlen( X );
	n = (int) strlen( Y );

	d = Matrix<int> ( m + 1, n + 1 );  // Edit Distance 
	e = Matrix<int> ( m + 1, n + 1 );  // Editing

	for ( i = 0 ; i <= m ; i++ )
		d[i][0] = i;

	for ( j = 1 ; j <= n ; j++ )
		d[0][j] = j;
	
	for ( i = 1 ; i <= m ; i++ )
	{
		for ( j = 1 ; j <= n ; j++ )
		{
			if ( X[i-1] == Y[j-1] )
				cost = 0;
			else
				cost = SUBSTITUTE_COST;

			d[i][j] = d[i][j-1] + INSERT_COST;
			e[i][j] = EDIT_INSERT;

			if ( d[i][j] >= d[i-1][j] + DELETE_COST ) 
			{
				d[i][j] = d[i-1][j] + 1;
				e[i][j] = EDIT_DELETE;
			}
				
			if ( d[i][j] >= d[i-1][j-1] + cost ) 
			{
				d[i][j] = d[i-1][j-1] + cost;
				if ( cost == 0 )
					e[i][j] = EDIT_NOCHANGE;
				else
					e[i][j] = EDIT_SUBSTITUTE;
			}
		}
	}
	edit_distance = d[m][n];

	cout << "The d table\n";
	for ( i = 0 ; i <= m ; i++ )
	{
		for ( j = 0 ; j <= n ; j++ )
			cout << d[i][j] << " ";
		cout << endl;
	}

	FreeMatrix<int> ( d, m + 1, n + 1 );
	FreeMatrix<int> ( e, m + 1, n + 1 );

	return edit_distance;
}


//
//  Coin Changing 
//  Given a set of base value base[0],..base[k-1] (coin change values) and and input number n
//  find the minimum number of base values needed to make change for n.
//
//  Reference: 名題精選百則使用C語言 (冼鏡光)
//
int Coin_Changing( int *base, int k, int n )
{
	int i, j, minimum;
	vector<int> money( n + 1 );

	money[0] = 0;
	money[1] = 1;
	for ( i = 2 ; i <= n ; i++ )
	{
		minimum = n;
		for ( j = 0 ; j < k ; j++ )
		{
			if ( i >= base[j] )
				minimum = min( money[ i - base[j] ] + 1, minimum );
		}
		money[i] = minimum;
	}
	return minimum;
}


//
//  0-1 Knapsack
//  A theif robbing a store find n items, where the i-th item worths value[i] dollars and 
//  weight[i] pounds, i = 1..n, he can only carry at most W pounds (knapsack). Find the 
//  maximum load. This routine implements the Dynamic Programming (DP) algorithm for solving 
//  the problem.
//
//  Reference: 名題精選百則使用C語言 (冼鏡光)
//
int DP_0_1_Knapsack( int *value, int *weight, int n, int W )
{
	int w, i, maximum;
	int **c;

	c = Matrix<int> ( n + 1, W + 1 );
	for ( w = 0 ; w <= W ; w++ )
		c[0][w] = 0;

	for ( i = 1 ; i <= n ; i++ )
	{
		c[i][0] = 0;
		for ( w = 1 ; w <= W ; w++ )
		{
			if ( weight[i] <= w )
			{
				if ( value[i] + c[i-1][w-weight[i]] > c[i-1][w] )
					c[i][w] = value[i] + c[i-1][w-weight[i]];
				else
					c[i][w] = c[i-1][w];
			}
			else
			{
				c[i][w] = c[i-1][w];
			}
		}
	}
	maximum = c[n][W];

	FreeMatrix( c, n + 1, W + 1 );
	return maximum;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Graph Algorithms
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Graph (Adjacency Matrix Representation)
//  By default, vertices are numbered as 1..n
//
//  Operations:
//     AdjacencyMatrix    (Initialize adjacency matrix)
//     ~AdjacencyMatrix   (Delete adjacency matrix)
//     SetEdge            (Set undirected edge)
//     SetDirectedEdge    (Set directed edge)
//     Display            (Display the adjacency matrix)
//     BFS                (Show the BFS sequence starting at source vertex) 
//     DFS                (Show the DFS sequence starting at source vertex)
//     Topological Sort   (Topological Sort)
//     Prim               (Prim's Algorithm)
//	   Bellman_Ford       (Bellman-Ford Algorithm)
//	   Dijkstra           (Dijkstra's Algorithm)
//     All_Pairs_Shortest_Paths (All-Pairs Shortest Paths)
//     Floyd_Warshall     (Floyd-Warshall Algorithm)
//     Transitive_Closure (Transitive Closure)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
#define INFINITY  1e30;
enum state { WHITE, GRAY, BLACK };  

class AdjacencyMatrix
{
public:
	int n_vertex;
	int n_edge;
	double **AdjMatrix;

	AdjacencyMatrix( int n );                                            
	~AdjacencyMatrix();

	void SetEdge( int start_vertex, int end_vertex );     
	void SetEdge( int start_vertex, int end_vertex, double weight );
	void SetDirectedEdge( int start_vertex, int end_vertex );
	void SetDirectedEdge( int start_vertex, int end_vertex, double weight );
	void Display();                                           

	void BFS( int source );                                    
	void DFS( int source );           
	void Topological_Sort( int source );
	void Prim( int source );
	void Bellman_Ford( int source );
	void Dijkstra( int source );
	void All_Pairs_Shortest_Paths();
	void Floyd_Warshall();
	void Transitive_Closure();

protected:
	void DFS_Visit( int u );
	void TS_DFS_Visit( int u );

	int    time;                    // visit time
	int    *color;                  // color = WHITE (not yet visit), color = GRAY (first visit), color = BLACK (done)
	double *d;                      // distance to the source vertex
	int    *pi;                     // vertex number of parent node 
	int    *f;                      // finish time

	SinglyLinkedList<int> TS_SLL;   // Linked List for the Topological Sort
};


AdjacencyMatrix::AdjacencyMatrix( int n )
{
	n_vertex = n;
	n_edge = 0;

	// Initialize the adjacency matrix
	AdjMatrix = Matrix<double> ( n_vertex + 1, n_vertex + 1 );
	for ( int i = 1 ; i <= n_vertex ; i++ )
		for ( int j = 1 ; j <= n_vertex ; j++ )
			AdjMatrix[i][j] = 0;

	color = new int[n_vertex+1];
	d     = new double[n_vertex+1];
	pi    = new int[n_vertex+1];
	f     = new int[n_vertex+1];  
}


AdjacencyMatrix::~AdjacencyMatrix()
{
	FreeMatrix( AdjMatrix, n_vertex + 1, n_vertex + 1 );

	delete [] color;	
	delete [] d;
	delete [] pi;   	
	delete [] f;
}


void AdjacencyMatrix::SetEdge( int start_vertex, int end_vertex )
{
	if ( start_vertex >= 1 && start_vertex <= n_vertex && end_vertex >= 1 && end_vertex <= n_vertex )
	{
		if ( AdjMatrix[start_vertex][end_vertex] == 0 && AdjMatrix[end_vertex][start_vertex] == 0 )
		{
			AdjMatrix[start_vertex][end_vertex] = 1;
			AdjMatrix[end_vertex][start_vertex] = 1;
			n_edge++;
		}
	}
}


void AdjacencyMatrix::SetEdge( int start_vertex, int end_vertex, double weight )
{
	if ( start_vertex >= 1 && start_vertex <= n_vertex && end_vertex >= 1 && end_vertex <= n_vertex )
	{
		if ( AdjMatrix[start_vertex][end_vertex] == 0 && AdjMatrix[end_vertex][start_vertex] == 0 )
		{
			AdjMatrix[start_vertex][end_vertex] = weight;
			AdjMatrix[end_vertex][start_vertex] = weight;
			n_edge++;
		}
	}
}


void AdjacencyMatrix::SetDirectedEdge( int start_vertex, int end_vertex )
{
	if ( start_vertex >= 1 && start_vertex <= n_vertex && end_vertex >= 1 && end_vertex <= n_vertex )
	{
		if ( AdjMatrix[start_vertex][end_vertex] == 0 )
		{
			AdjMatrix[start_vertex][end_vertex] = 1;
			n_edge++;
		}
	}
}


void AdjacencyMatrix::SetDirectedEdge( int start_vertex, int end_vertex, double weight )
{
	if ( start_vertex >= 1 && start_vertex <= n_vertex && end_vertex >= 1 && end_vertex <= n_vertex )
	{
		if ( AdjMatrix[start_vertex][end_vertex] == 0 )
		{
			AdjMatrix[start_vertex][end_vertex] = weight;
			n_edge++;
		}
	}
}


void AdjacencyMatrix::Display()
{
	int i, j;
	cout << "Graph (Adjacency Matrix Representation)\n";
	cout << "---------------------------------------\n";
	cout << "Number of Vertices = " << n_vertex << endl;
	cout << "Number of Edges = " << n_edge << endl;

	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			cout << AdjMatrix[i][j] << " ";
		}
		cout << endl;
	}
}


void AdjacencyMatrix::BFS( int source )
{
	int u, v;
	for ( u = 1 ; u <= n_vertex ; u++ )
	{
		if ( u != source )
		{
			color[u] = WHITE;
			d[u] = INFINITY;
			pi[u] = 0;
		}	
	}
	color[source] = GRAY;
	d[source] = 0;
	pi[source] = 0;

	cout << "BFS Sequences:\n";
	Queue<int> Q( n_vertex );
	Q.Enqueue( source );
	while ( !Q.IsEmpty() )
	{
		u = Q.Dequeue();
		cout << u << " ";
		for ( v = 1 ; v <= n_vertex ; v++ )
		{
			if ( AdjMatrix[u][v] )
			{
				if ( color[v] == WHITE )
				{
					color[v] = GRAY;
					d[v] = d[u] + AdjMatrix[u][v];
					pi[v] = u;
					Q.Enqueue( v );
				}
			}
		}
		color[u] = BLACK;
	}
	cout << endl;
}


void AdjacencyMatrix::DFS( int source )
{
	for ( int u = 1 ; u <= n_vertex ; u++ )
	{
		color[u] = WHITE;
		d[u] = INT_MAX;
		pi[u] = 0;
		f[u] = 0;
	}

	cout << "DFS Sequences:\n";
	time = 0;
	DFS_Visit( source );
	cout << endl;
}


void AdjacencyMatrix::DFS_Visit( int u )
{
	cout << u << " ";
	color[u] = GRAY;
	time++;
	d[u] = time;
	for ( int v = 1 ; v <= n_vertex ; v++ )
	{
		if ( AdjMatrix[u][v] )
		{
			if ( color[v] == WHITE )
			{
				pi[v] = u;
				DFS_Visit( v );
			}
		}
	}
	color[u] = BLACK;
	time++;
	f[u] = time;
}


void AdjacencyMatrix::Topological_Sort( int source )
{
	for ( int u = 1 ; u <= n_vertex ; u++ )
	{
		color[u] = WHITE;
		d[u] = INT_MAX;
		pi[u] = 0;
		f[u] = 0;
	}

	TS_SLL.Clear();
	time = 0;
	TS_DFS_Visit( source );
	
	cout << "Topological Sort\n";
	TS_SLL.current = TS_SLL.head->next;
	while ( TS_SLL.current != NULL )
	{
		cout << TS_SLL.current->key << " ";
		TS_SLL.current = TS_SLL.current->next;
	}
	cout << "\n";
}


void AdjacencyMatrix::TS_DFS_Visit( int u )
{
	color[u] = GRAY;
	time++;
	d[u] = time;
	for ( int v = 1 ; v <= n_vertex ; v++ )
	{
		if ( AdjMatrix[u][v] )
		{
			if ( color[v] == WHITE )
			{
				pi[v] = u;
				TS_DFS_Visit( v );
			}
		}
	}
	color[u] = BLACK;
	time++;
	f[u] = time;
	TS_SLL.Insert( u );
}


void AdjacencyMatrix::Prim( int source )
{
	int    i, j;
	int    *parent;    // Array to store parent node in MST
	double *key;       // Key values used to pick minimum cut
	bool   *set;       // Set of nodes not yet in MST
	double **MST;      // Adjancy matrix to store the MST 
	
	parent = new int[n_vertex + 1];
	key    = new double[n_vertex + 1];
	set    = new bool[n_vertex + 1];

	MST = Matrix<double> ( n_vertex + 1, n_vertex + 1 );

    // Initialize all vertices
	for ( i = 1 ; i <= n_vertex ; i++ )        
	{
		parent[i] = 0;
		key[i] = INFINITY;
		set[i] = false;
	}

	for ( i = 1 ; i <= n_vertex ; i++ )
		for ( j = 1 ; j <= n_vertex ; j++ )
			MST[i][j] = 0;
		
	key[source] = 0;         // Source vertex 
	parent[source] = 0;      // Source vertex is always the root
	set[source] = false;

	cout << "Prim's Algorithm (MST Sequence)\n";
	for ( i = 1 ; i <= n_vertex ; i++ )          // Add the vertices one-by-one into the MST
	{                                            // I didn't bother to use the minimum-priority-queue  
		double min = INFINITY;                   
		int min_idx = 0;
		for ( j = 1 ; j <= n_vertex ; j++ )      // Check vertices that are not yet in the MST
		{
			if ( set[j] == false && key[j] < min )
			{
				min = key[j];
				min_idx = j;
			}
		}

		set[min_idx] = true;
		cout << min_idx << " ";

		for ( j = 1 ; j <= n_vertex ; j++ )     // Add vertex and update key for minimum cut
		{
			if ( AdjMatrix[min_idx][j] != 0 && set[j] == false && AdjMatrix[min_idx][j] < key[j] )             
			{
				parent[j] = min_idx, key[j] = AdjMatrix[min_idx][j];
			}
		}
	}

	cout << endl;

	double min_cost = 0;
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		j = parent[i];
		if ( j != 0 )
		{
			MST[i][j] = MST[j][i] = 1;
			min_cost += AdjMatrix[i][j];
		}
	}
	
	cout << "Minimum Spanning Tree\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
			cout << MST[i][j] << " ";
		cout << endl;
	}

	cout << "Minimum Cost = " << min_cost << endl;

	delete [] parent;
	delete [] key;
	delete [] set;

	FreeMatrix( MST,  n_vertex + 1, n_vertex + 1 );
}


void AdjacencyMatrix::Bellman_Ford( int source )
{
	int i, j, k, u, v;
	for ( i = 1 ; i <= n_vertex ; i++ )
		d[i] = INFINITY;
	d[source] = 0;
	pi[source] = 0;

	vector<int> edge_start( n_edge );
	vector<int> edge_end( n_edge );

	k = 0;
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			if ( AdjMatrix[i][j] != 0 )
			{
				edge_start[k] = i;
				edge_end[k] = j;
				k++;
			}
		}
	}
	
	for ( i = 1 ; i <= n_vertex - 1 ; i++ )
	{
		for ( j = 0 ; j < n_edge ; j++ )
		{
			u = edge_start[j];
			v = edge_end[j];
			
			// Relax
			if ( d[v]  > d[u] + AdjMatrix[u][v] )
			{
				d[v] = d[u] + AdjMatrix[u][v];
				pi[v] = u;
			}
		}
	}

	cout << "Bellman-Ford Algorithm\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
		cout << "vertex " << i << " distance to source = " << d[i] << " parent = " << pi[i] << endl;
}


void AdjacencyMatrix::Dijkstra( int source )
{
	int i, j, u, v, n;
	double min;	

	// Initialize-Single-Source
	for ( i = 1 ; i <= n_vertex ; i++ )
		d[i] = INFINITY;
	d[source] = 0;
	pi[source] = 0;	

	vector<bool> set( n_vertex + 1 );
	for ( i = 1 ; i <= n_vertex ; i++ )
		set[i] = true;

	n = n_vertex;
	while ( n != 0 )
	{
		// Extract Minimum
		u = 0;
		min = INFINITY;
		for ( i = 1 ; i <= n_vertex ; i++ )
		{
			if ( set[i] && d[i] < min )
			{
				u = i;
				min = d[i];
			}
		}
		set[u] = false;
		n--;

		// Relax
		for ( v = 1 ; v <= n_vertex ; v++ )
		{
			if ( AdjMatrix[u][v] != 0 )
			{
				if ( d[v] > d[u] + AdjMatrix[u][v] )
				{
					d[v] = d[u] + AdjMatrix[u][v];
					pi[v] = u;
				}
			}
		}
	}

	cout << "Dijkstra's Algorithm\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
		cout << "vertex " << i << " distance to source = " << d[i] << " parent = " << pi[i] << endl;
}


void AdjacencyMatrix::All_Pairs_Shortest_Paths()
{
	int i, j, k, m;
	double **L, **W, l;
	
	L = Matrix<double> ( n_vertex + 1, n_vertex + 1 );
	W = Matrix<double> ( n_vertex + 1, n_vertex + 1 );
	
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			if ( i == j ) 
				L[i][j] = 0;
			else
			{
				if ( AdjMatrix[i][j] != 0 )
					L[i][j] = AdjMatrix[i][j];
				else
					L[i][j] = INFINITY;
			}
		}
	}

	for ( i = 1 ; i <= n_vertex ; i++ )
		for ( j = 1 ; j <= n_vertex ; j++ )
			W[i][j] = L[i][j];

	for ( m = 2 ; m <= n_vertex - 1 ; m++ )
	{
		for ( i = 1 ; i <= n_vertex ; i++ )
		{
			for ( j = 1 ; j <= n_vertex ; j++ )
			{
				l = INFINITY;
				for ( k = 1 ; k <= n_vertex ; k++ )
				{
					if ( l > L[i][j] )            l = L[i][j];
					if ( l > L[i][k] + W[k][j] )  l = L[i][k] + W[k][j]; 
				}
				L[i][j] = l;
			}
		}
	}
	
	cout << "All-Pairs Shortest Paths\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			cout << L[i][j] << " ";
		}
		cout << endl;
	}

	FreeMatrix( L, n_vertex + 1, n_vertex + 1 );
	FreeMatrix( W, n_vertex + 1, n_vertex + 1 );
}


void AdjacencyMatrix::Floyd_Warshall()
{
	int i, j, k;
	double **D, **D1;
	
	D  = Matrix<double> ( n_vertex + 1, n_vertex + 1 );
	D1 = Matrix<double> ( n_vertex + 1, n_vertex + 1 );

	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			if ( i == j ) 
				D[i][j] = 0;
			else
			{
				if ( AdjMatrix[i][j] != 0 )
					D[i][j] = AdjMatrix[i][j];
				else
					D[i][j] = INFINITY;
			}
		}
	}

	for ( k = 1 ; k <= n_vertex ; k++ )
	{
		for ( i = 1 ; i <= n_vertex ; i++ )
			for ( j = 1 ; j <= n_vertex ; j++ )
				D1[i][j] = min( D[i][j], D[i][k] + D[k][j] );

		for ( i = 1 ; i <= n_vertex ; i++ )
			for ( j = 1 ; j <= n_vertex ; j++ )
				D[i][j] = D1[i][j];
	}

	cout << "Floyd-Warshall Algorithm\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
			cout << D[i][j] << " ";
		cout << endl;
	}

	FreeMatrix( D,  n_vertex + 1, n_vertex + 1 );
	FreeMatrix( D1, n_vertex + 1, n_vertex + 1 );
}


void AdjacencyMatrix::Transitive_Closure()
{
	int i, j, k;
	double **T, **T1;
	
	T  = Matrix<double> ( n_vertex + 1, n_vertex + 1 );
	T1 = Matrix<double> ( n_vertex + 1, n_vertex + 1 );

	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			if ( i == j || AdjMatrix[i][j] != 0 ) 
				T[i][j] = 1;
			else
				T[i][j] = 0;
		}
	}

	for ( k = 1 ; k <= n_vertex ; k++ )
	{
		for ( i = 1 ; i <= n_vertex ; i++ )
			for ( j = 1 ; j <= n_vertex ; j++ )
				if ( T[i][j] == 1 || ( T[i][k] == 1 && T[k][j] == 1 ) )
					T1[i][j] = 1;

		for ( i = 1 ; i <= n_vertex ; i++ )
			for ( j = 1 ; j <= n_vertex ; j++ )
				T[i][j] = T1[i][j];
	}

	cout << "Transitive_Closure\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			cout << T[i][j] << " ";
		}
		cout << endl;
	}

	FreeMatrix( T,  n_vertex + 1, n_vertex + 1 );
	FreeMatrix( T1, n_vertex + 1, n_vertex + 1 );
}


//
//  Graph (Adjacency List Representation)
//  By default, vertices are numbered as 1..n
//
//  Operations:
//     AdjacencyList      (Initialize adjacency list)
//     ~AdjacencyList     (Delete adjacency list)
//     SetEdge            (Set undirected edge)
//     SetDirectedEdge    (Set directed edge)
//     Display            (Display the adjacency list)
//     BFS                (Show the BFS sequence starting at source vertex) 
//     DFS                (Show the DFS sequence starting at source vertex)
//     Topological Sort   (Topological Sort)
//	   Bellman_Ford       (Bellman-Ford Algorithm)
//	   Dijkstra           (Dijkstra's Algorithm)
//     All_Pairs_Shortest_Paths (All-Pairs Shortest Paths)
//     Floyd_Warshall     (Floyd-Warshall Algorithm)
//     Transitive_Closure (Transitive Closure)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
struct GraphNode {
	int vertex;
	double weight;
	GraphNode *next;
};


class AdjacencyList
{
public:
	int n_vertex;
	int n_edge;
	GraphNode **AdjList;

	AdjacencyList( int n );
	~AdjacencyList();

	void SetEdge( int start_vertex, int end_vertex );          
	void SetDirectedEdge( int start_vertex, int end_vertex ); 
	void SetEdge( int start_vertex, int end_vertex, double weight );          
	void SetDirectedEdge( int start_vertex, int end_vertex, double weight ); 
	void Display();

	void BFS( int source );
	void DFS( int source );
	void Topological_Sort( int source );
	void Prim( int source );
	void Bellman_Ford( int source );
	void Dijkstra( int source );
	void All_Pairs_Shortest_Paths();
	void Floyd_Warshall();
	void Transitive_Closure();

protected:
	GraphNode *ptr, *current;

	void DFS_Visit(int u);
	void TS_DFS_Visit(int u);

	int    time;                    // visit time
	int    *color;                  // color = WHITE (not yet visit), color = GRAY (first visit), color = BLACK (done)
	double *d;                      // distance to the source vertex
	int    *pi;                     // vertex number of parent node 
	int    *f;                      // finish time

	SinglyLinkedList<int> TS_SLL;   // Linked List for the Topological Sort
};


AdjacencyList::AdjacencyList( int n )
{
	n_vertex = n;
	n_edge = 0;

	// Initialize the adjacency list
	AdjList = (GraphNode **) new GraphNode *[n_vertex+1];
	for ( int i = 1 ; i <= n_vertex ; i++ )
		AdjList[i] = (GraphNode *) new GraphNode;

	for ( int i = 1 ; i <= n_vertex ; i++ )
	{
		AdjList[i]->vertex = i;
		AdjList[i]->weight = 0;
		AdjList[i]->next = NULL;
	}

	color = new int[n_vertex+1];
	d     = new double[n_vertex+1];
	pi    = new int[n_vertex+1];
	f     = new int[n_vertex+1];
}


AdjacencyList::~AdjacencyList()
{
	for ( int i = 1 ; i <= n_vertex ; i++ )
	{
		ptr = AdjList[i];
		current = AdjList[i]->next;
		delete ptr;
		while ( current != NULL )
		{
			ptr = current;
			current = current->next;
			delete ptr;
		}
	}

	delete [] color;	
	delete [] d;
	delete [] pi;   	
	delete [] f;
}


void AdjacencyList::SetEdge( int start_vertex, int end_vertex ) 
{
	SetDirectedEdge( start_vertex, end_vertex );
	SetDirectedEdge( end_vertex, start_vertex );
	n_edge--;
}


void AdjacencyList::SetDirectedEdge( int start_vertex, int end_vertex )
{
	if ( start_vertex >= 1 && start_vertex <= n_vertex && end_vertex >= 1 && end_vertex <= n_vertex )
	{
		ptr = new GraphNode;
		ptr->vertex = end_vertex;
		ptr->weight = 1;
		ptr->next = NULL;

		current = AdjList[start_vertex];
		while ( current->next != NULL )
			current = current->next;
		current->next = ptr;
		n_edge++;
	}
}


void AdjacencyList::SetEdge( int start_vertex, int end_vertex, double weight ) 
{
	SetDirectedEdge( start_vertex, end_vertex, weight );
	SetDirectedEdge( end_vertex, start_vertex, weight );
	n_edge--;
}


void AdjacencyList::SetDirectedEdge( int start_vertex, int end_vertex, double weight )
{
	if ( start_vertex >= 1 && start_vertex <= n_vertex && end_vertex >= 1 && end_vertex <= n_vertex )
	{
		ptr = new GraphNode;
		ptr->vertex = end_vertex;
		ptr->weight = weight;
		ptr->next = NULL;

		current = AdjList[start_vertex];
		while ( current->next != NULL )
			current = current->next;
		current->next = ptr;
		n_edge++;
	}
}


void AdjacencyList::Display()
{
	int i;
	cout << "Graph (Adjacency List Representation)  \n";
	cout << "---------------------------------------\n";
	cout << "Number of Vertices = " << n_vertex << endl;
	cout << "Number of Edges = " << n_edge << endl;

	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		current = AdjList[i]->next;
		cout << i << " -> ";
		while ( current != NULL )
		{
			cout << current->vertex << "(" << current->weight << ") -> ";
			current = current->next;
		}
		cout << "NULL\n";
	}		
}


void AdjacencyList::BFS( int source )
{
	int u, v;
	for ( u = 1 ; u <= n_vertex ; u++ )
	{
		if ( u != source )
		{
			color[u] = WHITE;
			d[u] = INFINITY;
			pi[u] = 0;
		}	
	}
	color[source] = GRAY;
	d[source] = 0;
	pi[source] = 0;

	cout << "BFS Sequences:\n";
	Queue<int> Q( n_vertex );
	Q.Enqueue( source );
	while ( !Q.IsEmpty() )
	{
		u = Q.Dequeue();
		cout << u << " ";
		current = AdjList[u]->next;
		while ( current != NULL )
		{
			v = current->vertex;
			if ( color[v] == WHITE )
			{
				color[v] = GRAY;
				d[v] = d[u] + current->weight;
				pi[v] = u;
				Q.Enqueue( v );
			}
			current = current->next;
		}
		color[u] = BLACK;
	}
	cout << endl;
}


void AdjacencyList::DFS( int source )
{
	for ( int u = 1 ; u <= n_vertex ; u++ )
	{
		color[u] = WHITE;
		d[u] = INT_MAX;
		pi[u] = 0;
		f[u] = 0;
	}

	cout << "DFS Sequences:\n";
	time = 0;
	DFS_Visit( source );
	cout << endl;
}


void AdjacencyList::DFS_Visit( int u )
{
	GraphNode *ptr;
	int v;

	cout << u << " ";
	color[u] = GRAY;
	time++;
	d[u] = time;

	ptr = AdjList[u]->next;
	while ( ptr != NULL )
	{
		v = ptr->vertex;
		if ( color[v] == WHITE )
		{
			pi[v] = u;
			DFS_Visit(v);
		}
		ptr = ptr->next;
	}
	color[u] = BLACK;
	time++;
	f[u] = time;
}


void AdjacencyList::Topological_Sort( int source )
{
	for ( int u = 1 ; u <= n_vertex ; u++ )
	{
		color[u] = WHITE;
		d[u] = INT_MAX;
		pi[u] = 0;
		f[u] = 0;
	}

	TS_SLL.Clear();
	time = 0;
	TS_DFS_Visit( source );
	
	cout << "Topological Sort\n";
	TS_SLL.current = TS_SLL.head->next;
	while ( TS_SLL.current != NULL )
	{
		cout << TS_SLL.current->key << " ";
		TS_SLL.current = TS_SLL.current->next;
	}
	cout << "\n";
}


void AdjacencyList::Prim( int source )
{
	int    i, j;
	int    *parent;    // Array to store parent node in MST
	double *key;       // Key values used to pick minimum cut
	bool   *set;       // Set of nodes not yet in MST
	double **MST;      // Adjancy matrix to store the MST 
	
	parent = new int[n_vertex + 1];
	key    = new double[n_vertex + 1];
	set    = new bool[n_vertex + 1];

	MST = Matrix<double> ( n_vertex + 1, n_vertex + 1 );

    // Initialize all vertices
	for ( i = 1 ; i <= n_vertex ; i++ )        
	{
		parent[i] = 0;
		key[i] = INFINITY;
		set[i] = false;
	}

	for ( i = 1 ; i <= n_vertex ; i++ )
		for ( j = 1 ; j <= n_vertex ; j++ )
			MST[i][j] = 0;
		
	key[source] = 0;         // Source vertex 
	parent[source] = 0;      // Source vertex is always the root
	set[source] = false;

	cout << "Prim's Algorithm (MST Sequence)\n";
	for ( i = 1 ; i <= n_vertex ; i++ )          // Add the vertices one-by-one into the MST
	{                                            // I didn't bother to use the minimum-priority-queue  
		double min = INFINITY;                   
		int min_idx = 0;
		for ( j = 1 ; j <= n_vertex ; j++ )      // Check vertices that are not yet in the MST
		{
			if ( set[j] == false && key[j] < min )
			{
				min = key[j];
				min_idx = j;
			}
		}

		set[min_idx] = true;
		cout << min_idx << " ";

		current = AdjList[min_idx]->next;
		while ( current != NULL )
		{
			if ( current->weight != 0 && set[current->vertex] == false && current->weight < key[current->vertex] )
			{
				parent[current->vertex] = min_idx, key[current->vertex] = current->weight;
			}
			current = current->next;
		}
	}

	cout << endl;

	double min_cost = 0;
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		j = parent[i];
		if ( j != 0 )
		{
			MST[i][j] = MST[j][i] = 1;
			
			current = AdjList[i]->next;
			while ( current != NULL )
			{
				if ( current->vertex == j )
				{
					min_cost += current->weight;
				}
				current = current->next;
			}
		}
	}
 	
	cout << "Minimum Spanning Tree\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
			cout << MST[i][j] << " ";
		cout << endl;
	}

	cout << "Minimum Cost = " << min_cost << endl;

	delete [] parent;
	delete [] key;
	delete [] set;

	FreeMatrix( MST,  n_vertex + 1, n_vertex + 1 );
}


void AdjacencyList::TS_DFS_Visit( int u )
{
	GraphNode *ptr;
	int v;

	color[u] = GRAY;
	time++;
	d[u] = time;

	ptr = AdjList[u]->next;
	while ( ptr != NULL )
	{
		v = ptr->vertex;
		if ( color[v] == WHITE )
		{
			pi[v] = u;
			TS_DFS_Visit( v );
		}
		ptr = ptr->next;
	}
	color[u] = BLACK;
	time++;
	f[u] = time;
	TS_SLL.Insert( u );
}


void AdjacencyList::Bellman_Ford( int source )
{
	int i, j, k, u, v;
	for ( i = 1 ; i <= n_vertex ; i++ )
		d[i] = INFINITY;
	d[source] = 0;
	pi[source] = 0;

	vector<int> edge_start( n_edge );
	vector<int> edge_end( n_edge );

	k = 0;
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		current = AdjList[i]->next;
		while ( current != NULL )
		{
			edge_start[k] = i;
			edge_end[k] = current->vertex;
			k++;
			current = current->next;
		}
	}		

	for ( i = 1 ; i <= n_vertex - 1 ; i++ )
	{
		for ( j = 0 ; j < n_edge ; j++ )
		{
			u = edge_start[j];
			v = edge_end[j];
			
			double weight;
			current = AdjList[u]->next;
			while ( current != NULL )
			{
				if ( current->vertex == v )
				{
					weight = current->weight;
					break;
				}
				else
					current = current->next;
			}

			// Relax
			if ( d[v]  > d[u] + weight )
			{
				d[v] = d[u] + weight;
				pi[v] = u;
			}
		}
	}

	cout << "Bellman-Ford Algorithm\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
		cout << "vertex " << i << " distance to source = " << d[i] << " parent = " << pi[i] << endl;
}


void AdjacencyList::Dijkstra(int source)
{
	int i, j, u, v, n;
	double min;	

	// Initialize-Single-Source
	for ( i = 1 ; i <= n_vertex ; i++ )
		d[i] = INFINITY;
	d[source] = 0;
	pi[source] = 0;	

	vector<bool> set( n_vertex + 1 );
	for ( i = 1 ; i <= n_vertex ; i++ )
		set[i] = true;

	n = n_vertex;
	while ( n != 0 )
	{
		// Extract Minimum
		u = 0;
		min = INFINITY;
		for ( i = 1 ; i <= n_vertex ; i++ )
		{
			if ( set[i] && d[i] < min )
			{
				u = i;
				min = d[i];
			}
		}
		set[u] = false;
		n--;

		// Relax
		current = AdjList[u]->next;
		while ( current != NULL )
		{
			v = current->vertex;				

			if ( d[v] > d[u] + current->weight )
			{
				d[v] = d[u] + current->weight;
				pi[v] = u;
			}
			current = current->next;
		}
	}

	cout << "Dijkstra's Algorithm\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
		cout << "vertex " << i << " distance to source = " << d[i] << " parent = " << pi[i] << endl;
}


void AdjacencyList::All_Pairs_Shortest_Paths()
{
	int i, j, k, m;
	double **L, **W, l;
	
	L = Matrix<double> ( n_vertex + 1, n_vertex + 1 );
	W = Matrix<double> ( n_vertex + 1, n_vertex + 1 );
	
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			if ( i == j ) L[i][j] = 0;
			else	      L[i][j] = INFINITY;
		}
	}

	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		current = AdjList[i]->next;
		while ( current != NULL )
		{
			L[i][current->vertex] = current->weight;
			current = current->next;
		}
	}

	for ( i = 1 ; i <= n_vertex ; i++ )
		for ( j = 1 ; j <= n_vertex ; j++ )
			W[i][j] = L[i][j];
	
	for ( m = 2 ; m <= n_vertex - 1 ; m++ )
	{
		for ( i = 1 ; i <= n_vertex ; i++ )
		{
			for ( j = 1 ; j <= n_vertex ; j++ )
			{
				l = INFINITY;
				for ( k = 1 ; k <= n_vertex ; k++ )
				{
					if ( l > L[i][j] )            l = L[i][j];
					if ( l > L[i][k] + W[k][j] )  l = L[i][k] + W[k][j]; 
				}
				L[i][j] = l;
			}
		}
	}
	
	cout << "All-Pairs Shortest Paths\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			cout << L[i][j] << " ";
		}
		cout << endl;
	}

	FreeMatrix( L, n_vertex + 1, n_vertex + 1 );
	FreeMatrix( W, n_vertex + 1, n_vertex + 1 );
}


void AdjacencyList::Floyd_Warshall()
{
	int i, j, k;
	double **D, **D1;
	
	D  = Matrix<double> ( n_vertex + 1, n_vertex + 1 );
	D1 = Matrix<double> ( n_vertex + 1, n_vertex + 1 );

	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			if ( i == j )  D[i][j] = 0;
			else	       D[i][j] = INFINITY;
		}
	}

	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		current = AdjList[i]->next;
		while ( current != NULL )
		{
			D[i][current->vertex] = current->weight;
			current = current->next;
		}
	}		

	for ( k = 1 ; k <= n_vertex ; k++ )
	{
		for ( i = 1 ; i <= n_vertex ; i++ )
			for ( j = 1 ; j <= n_vertex ; j++ )
				D1[i][j] = min( D[i][j], D[i][k]+D[k][j] );

		for ( i = 1 ; i <= n_vertex ; i++ )
			for ( j = 1 ; j <= n_vertex ; j++ )
				D[i][j] = D1[i][j];
	}

	cout << "Floyd-Warshall Algorithm\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
			cout << D[i][j] << " ";
		cout << endl;
	}

	FreeMatrix( D,  n_vertex + 1, n_vertex + 1 );
	FreeMatrix( D1, n_vertex + 1, n_vertex + 1 );
}


void AdjacencyList::Transitive_Closure()
{
	int i, j, k;
	double **T, **T1;
	
	T  = Matrix<double> ( n_vertex + 1, n_vertex + 1 );
	T1 = Matrix<double> ( n_vertex + 1, n_vertex + 1 );

	for ( i = 1 ; i <= n_vertex ; i++ )
		T[i][i] = 1;

	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		current = AdjList[i]->next;
		while ( current != NULL )
		{
			T[i][current->vertex] = 1;
			current = current->next;
		}
	}	

	for ( k = 1 ; k <= n_vertex ; k++ )
	{
		for ( i = 1 ; i <= n_vertex ; i++ )
			for ( j = 1 ; j <= n_vertex ; j++ )
				if ( T[i][j] == 1 || ( T[i][k] == 1 && T[k][j] == 1 ) )
					T1[i][j] = 1;

		for ( i = 1 ; i <= n_vertex ; i++ )
			for ( j = 1 ; j <= n_vertex ; j++ )
				T[i][j] = T1[i][j];
	}

	cout << "Transitive_Closure\n";
	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			cout << T[i][j] << " ";
		}
		cout << endl;
	}

	FreeMatrix( T,  n_vertex + 1, n_vertex + 1 );
	FreeMatrix( T1, n_vertex + 1, n_vertex + 1 );
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Flow Network
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  This subroutine implements the algorithm for solving the flow network problem.
//  The algorithm is based on the Ford-Fulkerson method, in which the BFS is used 
//  for the augmenting path (i.e., Edmonds-Karp algorithm).
//
//  By default, vertices are numbered as 1..n (source = 1, sink = n)
//
//  Operations:
//     FlowNetwork        (Initialize the flow network)
//     ~FlowNetwork       (Delete the flow network)
//     SetEdge            (Set directed edge)
//     Display            (Display the flow network)
//     BFS                (BFS to determine if the augmenting path exists) 
//     MaxFlow			  (Find the maximum flow)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
class FlowNetwork
{
public:
	int n_vertex;          // Number of vertices
	int n_edge;            // Number of edges
	int source, sink;      // Source (default = 1) & sink (default = n_vertex)
	
	int **flow, **cap;     // Flow & Capacity
	int *color, *pi;       // For the BFS

	FlowNetwork( int n );                                            
	~FlowNetwork();

	void SetEdge( int start_vertex, int end_vertex, int capacity );
	void Display();

	bool BFS();
	int  MaxFlow();
};


FlowNetwork::FlowNetwork( int n )
{
	n_vertex = n;
	n_edge = 0;
	source = 1;
	sink = n;

	// Initialize the capacity to 0
	cap = Matrix<int> ( n_vertex + 1, n_vertex + 1 );  
	for ( int i = 1 ; i <= n_vertex ; i++ )
		for ( int j = 1 ; j <= n_vertex ; j++ )
			cap[i][j] = 0;

	// Initialize the flow to 0
	flow = Matrix<int> ( n_vertex + 1, n_vertex + 1 );
	for ( int i = 1 ; i <= n_vertex ; i++ )
		for ( int j = 1 ; j <= n_vertex ; j++ )
			flow[i][j] = 0;

	// For the BFS
	color = new int[n_vertex+1];
	pi    = new int[n_vertex+1];
}


FlowNetwork::~FlowNetwork()
{
	FreeMatrix( cap, n_vertex + 1, n_vertex + 1 );
	FreeMatrix( flow, n_vertex + 1, n_vertex + 1 );
}


void FlowNetwork::SetEdge( int start_vertex, int end_vertex, int capacity )
{
	if ( start_vertex >= 1 && start_vertex <= n_vertex && end_vertex >= 1 && end_vertex <= n_vertex )
	{
		if ( cap[start_vertex][end_vertex] == 0 )
		{
			cap[start_vertex][end_vertex] = capacity;
			n_edge++;
		}
	}
}


void FlowNetwork::Display()
{
	int i, j;
	cout << "Flow Network\n";
	cout << "---------------------------------------\n";
	cout << "Number of Vertices = " << n_vertex << endl;
	cout << "Number of Edges = " << n_edge << endl;

	for ( i = 1 ; i <= n_vertex ; i++ )
	{
		for ( j = 1 ; j <= n_vertex ; j++ )
		{
			cout << cap[i][j] << " ";
		}
		cout << endl;
	}
}


bool FlowNetwork::BFS() 
{
	int u, v;
	for ( u = 1 ; u <= n_vertex ; u++ )
	{
		if ( u != source )
		{
			color[u] = WHITE;
		}	
	}
	color[source] = GRAY;
	pi[source] = 0;

	Queue<int> Q( n_vertex );
	Q.Enqueue( source );
	while ( !Q.IsEmpty() )
	{
		u = Q.Dequeue();
		for ( v = 1 ; v <= n_vertex ; v++ )
		{
			if ( color[v] == WHITE && cap[u][v] - flow[u][v] > 0 )  
			{
				color[v] = GRAY;
				Q.Enqueue( v );
				pi[v] = u;
			}
		}
		color[u] = BLACK;
	}

	if ( color[sink] == BLACK )  // If the sink is reachable, then the augmenting path exists. 
		return true;
	else
		return false;
}


int FlowNetwork::MaxFlow()
{
	int u, max_flow = 0;

	while ( BFS() )   // Determine if the augmenting path exists
	{
		int increment = 10000000;  
		for ( u = sink ; pi[u] >= source ; u = pi[u] ) 
		{
			if ( increment > cap[pi[u]][u] - flow[pi[u]][u] )  
			{
				increment = cap[pi[u]][u] - flow[pi[u]][u];
			}
		}
    
		for ( u = sink ; pi[u] >= source ; u = pi[u] )  // Update the flows given the path
		{
			flow[pi[u]][u] += increment;
			flow[u][pi[u]] -= increment;
		}

		max_flow += increment; 
    }
    
	return max_flow;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Matrix Operations
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Matrix Addition
//  Add two matrices A & B, each with dimension m x n
//  
//  Input:  A[0..m-1][0..n-1] & B[0..m-1][0..n-1]
//  Output: C[0..m-1][0..n-1] 
//
void Matrix_Addition( double **A, double **B, double **C, int m, int n )
{
	for ( int i = 0 ; i < m ; i++ )
		for ( int j = 0 ; j < n ; j++ )
			C[i][j] = A[i][j] + B[i][j];
}


//
//  Matrix Subtraction
//  Subtract two matrices A & B, each with dimension m x n
//  
//  Input:  A[0..m-1][0..n-1] & B[0..m-1][0..n-1]
//  Output: C[0..m-1][0..n-1] 
//
void Matrix_Subtraction( double **A, double **B, double **C, int m, int n )
{
	for ( int i = 0 ; i < m ; i++ )
		for ( int j = 0 ; j < n ; j++ )
			C[i][j] = A[i][j] - B[i][j];
}


//
//  Matrix Multiplication
//  Mutiply two matrices A & B, each with dimension n x n
//  
//  Input:  A[0..n-1][0..n-1] & B[0..n-1][0..n-1]
//  Output: C[0..n-1][0..n-1]
//
void Matrix_Multiplication( double **A, double **B, double **C, int n )
{
	for ( int i = 0 ; i < n ; i++ )
	{
		for ( int j = 0 ; j < n ; j++ )
		{
			C[i][j] = 0.0;
			for ( int k = 0 ; k < n ; k++ ) 
				C[i][j] += A[i][k] * B[k][j];
		}
	}
}


//
//  Matrix Multiplication
//  Mutiply two matrices A & B, with dimensions p x q and q x r
//  
//  Input:  A[0..p-1][0..q-1] & B[0..q-1][0..r-1]
//  Output: C[0..p-1][0..r-1]
//
void Matrix_Multiplication( double **A, double **B, double **C, int p, int q, int r )
{
	for ( int i = 0 ; i < p ; i++ )
	{
		for ( int j = 0 ; j < r ; j++ )
		{
			C[i][j] = 0.0;
			for ( int k = 0 ; k < q ; k++ ) 
				C[i][j] += A[i][k] * B[k][j];
		}
	}
}


//
//  Matrix Multiplication
//  Mutiply two matrices A & B, with dimensions n x n and n x 1
//  Generally used for the linear equations such as Ax = b.
//
//  Input:  A[0..n-1][0..n-1] & B[0..n-1]
//  Output: C[0..n-1]
//
void Matrix_Multiplication( double **A, double *B, double *C, int n )
{
	for ( int i = 0 ; i < n ; i++ )
	{
		C[i] = 0.0;
		for ( int j = 0 ; j < n ; j++ ) 
			C[i] += A[i][j] * B[j];
	}
}


//
//  LU Decomposition
//  Given an n x n matrix A[1..n][1..n], compute an LU Decomposition of A such that A = LU.
//  L: Lower-triangular matrix (unit)
//  U: Upper-triangular matrix
//
//  Input: A matrix (dimension: n x n)
//  Output: A matrix is replaced with L & U matrices 
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void LU_Decomposition( double **A, int n )
{
	int i, j, k;
	for ( k = 0 ; k < n ; k++ )
	{
		for ( i = k + 1 ; i < n ; i++ )
		{
			if ( A[k][k] == 0.0 ) cout << "LU decomposition error - divide by 0\n";
			else
				A[i][k] = A[i][k] / A[k][k];
		}

		for ( i = k + 1 ; i < n ; i++ )
		{
			for ( j = k + 1 ; j < n ; j++ )
			{
				A[i][j] -= A[i][k] * A[k][j];
			}
		}
	}
}


//
//  LUP Decomposition
//  Given an n x n matrix A[0..n-1][0..n-1], compute an LUP Decomposition of A such that PA = LU.
//  L: Lower-triangular matrix (unit)
//  U: Upper-triangular matrix
//  P: Permutation index (representing the Permutation matrix)
//
//  Input: A matrix (dimension: n x n)
//  Output: A matrix is replaced with L & U matrices
//          P is the permutation index 
//          d is output as 1 or -1 depending on whether the number of row interchanges was even or odd
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void LUP_Decomposition( double **A, int *P, double &d, int n )
{
	int i, j, k, kp;
	double max;

	for ( i = 0 ; i < n ; i++ )
		P[i] = i+1;

	d = 1.0;
	for ( k = 0 ; k < n ; k++ )
	{
		max = (double) 0.0;
		kp = k;
		for ( i = k ; i < n ; i++ )
		{
			if ( fabs(A[i][k]) > max )
			{
				max = fabs( A[i][k] );
				kp = i;
			}
		}
		if ( max == 0.0 ) cout << "Error: A is a singular matrix" << endl;
		
		if ( kp != k )
		{
			swap( P[k], P[kp] );
			for ( i = 0 ; i < n ; i++ ) 
				swap( A[k][i], A[kp][i] );
			d = -d;
		}
		
		for ( i = k + 1 ; i < n ; i++ )
		{
			if ( A[k][k] != 0.0 ) 
				A[i][k] /= (double) A[k][k];
		}

		for ( i = k + 1 ; i < n ; i++ )
		{
			for ( j = k + 1 ; j < n ; j++ )
			{
				A[i][j] -= A[i][k]*A[k][j];
			}
		}
	}
}


//
//  Solve Linear System using Forward and Back Substitution
//  Given the matrix A[0..n-1][0..n-1], solve the set of n linear equations Ax = b.
//
//  Input: A is the matrix after LUP Decomposition (not the original matrix A)
//         P is the permutation index 
//         b is the right-hand side vector
//  Output: b vector is replaced with the solution vector
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void LUP_Solve( double **A, int *P, double *b, int n )
{
	int i, j;
	double *x, *y, sum;
	
	x = new double[n];
	y = new double[n];
	
	y[0] = b[P[0]-1];
	for ( i = 1 ; i < n ; i++ )                  // Back substitution
	{
		sum = 0.0;
		for ( j = 0 ; j < i ; j++ )
			sum += A[i][j] * y[j];
		y[i] = b[P[i]-1] - sum;
	}
	x[n-1] = y[n-1] / A[n-1][n-1];

	for ( i = n - 2 ; i >= 0 ; i-- )             // Forward substitution
	{
		sum = 0.0;
		for ( j = i + 1 ; j < n ; j++ )
			sum += A[i][j] * x[j];
		x[i] = ( y[i] - sum ) / A[i][i];
	}

	for ( i = 0 ; i < n ; i++ )
		b[i] = x[i];

	delete [] x;
	delete [] y;
}


//
//  Solve Linear System using LUP Decomposition
//  Given the matrix A[0..n-1][0..n-1] and b[0..n-1], solve the set of n linear equations A x = b.
//  
//  Input:  A is the input matrix (A is replaced with the LUP decomposition)
//          b is the input vector
//  Output: b is replaced with the solution vector x
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Solve_Linear_System( double **A, double *b, int n )
{
	int *P;
	double d;
	P = new int[n];
	LUP_Decomposition( A, P, d, n );
	LUP_Solve( A, P, b, n );
	delete [] P;
}


//
//  Matrix Inverse
//  Given a matrix A[0..n-1][0..n-1], find the inverse of the matrix A.
//  Upon return, the matrix A is replaced with its inverse.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Matrix_Inverse( double **A, int n )
{
	int i, j;
	int *P;
	double **temp, *b, d;

	P = new int[n];
	b = new double[n];
	temp = Matrix<double>( n, n );

	LUP_Decomposition( A, P, d, n );         // LUP Decomposition of A (only once)
	for ( j = 0 ; j < n ; j++ )              // Solve for each column
	{
		for ( i = 0 ; i < n ; i++ ) b[i] = 0.0;
		b[j] = 1.0;
		LUP_Solve( A, P, b, n );
		for ( i = 0 ; i < n ; i++ ) temp[i][j] = b[i];
	}

	for ( i = 0 ; i < n ; i++ )              // Inverse of A
		for ( j = 0 ; j < n ; j++ )
			A[i][j] = temp[i][j];

	delete [] P;
	delete [] b;
	FreeMatrix( temp, n, n );
}


//
//  Determinant
//  Given a matrix A[0..n-1][0..n-1], find the detrminant of A.
//  This routine is intended to keep the matrix A intact.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
double Determinant( double **A, int n )
{
	int i, j;
	double **temp, d, ans;
	int *P;

	temp = Matrix<double>( n, n );
	P = new int[n];

	for ( i = 0 ; i < n ; i++ )             // Temporary for LUP Decomposition
		for ( j = 0 ; j < n ; j++ )
			temp[i][j] = A[i][j];	

	LUP_Decomposition( temp, P, d, n );      
	ans = 1.0;                              // Compute the determinant 
	for ( i = 0 ; i < n ; i++ )
		ans *= temp[i][i];
	ans *= (double) d;

	FreeMatrix( temp, n, n );
	delete [] P;

	return ans;
}


// 
//  Gauss-Jordan Elimination
//  Linear equation soluton by Gauss-Jordan elimination. The input matrix is
//  a[0..n-1][0..n-1]. b[0..n-1][0..m-1] is input containing the m right-hand
//  side vectors. On output, a is replaced by its matrix inverse, and b is 
//  replaced by the corresponding set of solution vectors.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void gaussj( double **a, double **b, int n, int m )
{
	int i,icol,irow,j,k,l,ll;
	DP big,dum,pivinv;

	vector<int> indxc( n ), indxr( n ), ipiv( n );

	for ( j = 0 ; j < n ; j++ ) 
		ipiv[j] = 0;

	for ( i = 0 ; i < n ; i++) 
	{
		big=0.0;
		for ( j = 0 ; j < n ; j++ )
			if ( ipiv[j] != 1 )
				for ( k = 0 ; k < n ; k++ ) 
				{
					if ( ipiv[k] == 0 ) 
					{
						if ( fabs(a[j][k]) >= big ) 
						{
							big = fabs( a[j][k] );
							irow = j;
							icol = k;
						}
					}
				}

		++( ipiv[icol] );
		if ( irow != icol ) 
		{
			for ( l = 0 ; l < n ; l++ ) swap( a[irow][l], a[icol][l] );
			for ( l = 0 ; l < m ; l++ ) swap( b[irow][l], b[icol][l] );
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if ( a[icol][icol] == 0.0 ) nrerror( "gaussj: Singular Matrix" );
		pivinv=1.0 / a[icol][icol];
		a[icol][icol] = 1.0;
		for ( l = 0 ; l < n ; l++ ) a[icol][l] *= pivinv;
		for ( l = 0 ; l < m ; l++ ) b[icol][l] *= pivinv;
		for ( ll = 0 ; ll < n ; ll++ )
			if ( ll != icol ) 
			{
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l = 0 ; l < n ; l++ ) a[ll][l] -= a[icol][l] * dum;
				for (l = 0 ; l < m ; l++ ) b[ll][l] -= b[icol][l] * dum;
			}
	}

	for ( l = n - 1 ; l >= 0 ; l-- ) 
	{
		if ( indxr[l] != indxc[l] )
			for ( k = 0 ; k < n ; k++ )
				swap( a[k][indxr[l]], a[k][indxc[l]] );
	}
}
 

//
//  LU Decomposition
//  Given a matrix a[0..n-1][0..n-1], this routine replaces it by LU decomposition of a
//  rowwise permutation of itself. a is input. On output, it is re-arranged with the L
//  and U matrices. index[0..n-1] is an output vector that records the row permutation
//  effected by the partial pivoting; d is output as +-1 depending on whether the number
//  of row interchanges was even or odd, respectively. This routine is used in combination
//  with lubksb to solve linera equations or invert a matrix.
//
//  Remark: This routine uses the Crout's method for pivoting. This routine is essentially
//  the same as the LUP decomposition as described in Introduction to Algorithms (Cormen).
//
//  Reference: Numerical Recipes in C++ (Press)
//
void ludcmp( double **a, int *indx, double &d, int n )
{
	const DP TINY = 1.0e-20;
	int i, imax, j, k;
	DP big, dum, sum, temp;

	vector<double> vv( n );
	d = 1.0;
	for (i = 0 ; i < n ; i++ ) 
	{
		big=0.0;
		for ( j = 0 ; j < n ; j++ )
			if ( ( temp=fabs( a[i][j] ) ) > big ) big=temp;
		if ( big == 0.0 ) nrerror( "Singular matrix in routine ludcmp" );
		vv[i] = 1.0 / big;
	}

	for ( j = 0 ; j < n ; j++ ) 
	{
		for ( i = 0 ; i < j ; i++ ) 
		{
			sum = a[i][j];
			for ( k = 0 ; k < i ; k++ ) 
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}

		big=0.0;
		for ( i = j ; i < n ; i++ ) 
		{
			sum = a[i][j];
			for ( k = 0 ; k < j ; k++ ) 
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ( ( dum = vv[i] * fabs(sum) ) >= big ) 
			{
				big = dum;
				imax = i;
			}
		}
		if ( j != imax ) 
		{
			for ( k = 0 ; k < n ; k++ ) 
			{
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if ( a[j][j] == 0.0 ) a[j][j] = TINY;
		if ( j != n - 1 ) 
		{
			dum = 1.0 / ( a[j][j] );
			for ( i = j + 1 ; i < n ; i++ ) a[i][j] *= dum;
		}
	}
}


//
//  Solves the set of n linear equations AX = b. Here a[0..n-1][0..n-1] is input, not
//  as the matrix A but rather as its LU decomposition, determined by the routine ludcmp.
//  index[0..n-1] is input as the permutation vector returned by ludcmp. b[0..n-1] is 
//  input as the right-hand side vector b, and returns with the solution vector X. a and
//  indx are not modified by this routine and can be left in place for successive calls 
//  with the different right-hand sides b. This routine takes into account the probability
//  that b will begin with many zero elements, so it is efficient for use in matrix
//  inversion.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void lubksb( double **a, int *indx, double *b, int n )
{
	int i, ii = 0, ip, j;
	DP sum;

	for ( i = 0 ; i < n ; i++ ) 
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if ( ii != 0 )
			for ( j = ii - 1 ; j < i ; j++ ) sum -= a[i][j] * b[j];
		else if ( sum != 0.0 )
			ii = i + 1;
		b[i] = sum;
	}
	for ( i = n - 1 ; i >= 0 ; i-- ) 
	{
		sum = b[i];
		for ( j = i + 1 ; j < n ; j++ ) 
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}


//
//  Singular Value Decomposition
//  Given a matrix a[0..m-1][0..n-1], this routine compute its singular value decomposition.
//  A = U W V^T. The matrix U replaces a on output. The diagonal matrix of singular values
//  W is output as a vector w[0..n-1]. The matrix V (not the transpost V^T) is output as
//  v[0..n-1][0..n-1].
//
//  Reference: Numerical Recipes in C++ (Press)
//
template <class T>
inline const T SQR( const T a ) { return a*a; }

template <class T>
inline const T MAX( const T &a, const T &b ) { return b > a ? (b) : (a); }

template <class T>
inline const T MIN( const T &a, const T &b ) { return b < a ? (b) : (a); }

template <class T>
inline const T SIGN( const T &a, const T &b ) { return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }


double pythag( const DP a, const DP b )
{
	DP absa,absb;

	absa = fabs( a );
	absb = fabs( b );
	if ( absa > absb ) return absa * sqrt( 1.0 + SQR( absb / absa ) );
	else return ( absb == 0.0 ? 0.0 : absb * sqrt( 1.0 + SQR( absa / absb ) ) );
}


void svdcmp( double **a, double *w, double **v, int m, int n )
{
	bool flag;
	int i, its, j, jj, k, l, nm;
	DP anorm, c, f, g, h, s, scale, x, y, z;

	vector<double> rv1( n );
	g = scale = anorm = 0.0;
	for ( i = 0 ; i < n ; i++ ) 
	{
		l = i + 2;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if ( i < m ) 
		{
			for ( k = i ; k < m ; k++ ) 
				scale += fabs( a[k][i] );

			if ( scale != 0.0 ) 
			{
				for ( k = i ; k < m ; k++ ) 
				{
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN( sqrt( s ), f );
				h = f * g - s;
				a[i][i] = f - g;
				for ( j = l - 1 ; j < n ; j++ ) 
				{
					for ( s = 0.0, k = i ; k < m ; k++ ) 
						s += a[k][i] * a[k][j];
					f=s / h;
					for ( k = i ; k < m ; k++ ) a[k][j] += f * a[k][i];
				}
				for ( k = i ; k < m ; k++ ) 
					a[k][i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if ( i + 1 <= m && i != n ) 
		{
			for ( k = l - 1 ; k < n ; k++ ) 
				scale += fabs( a[i][k] );
			
			if ( scale != 0.0 ) 
			{
				for ( k = l - 1 ; k < n ; k++ ) 
				{
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l-1];
				g = -SIGN( sqrt( s ), f );
				h = f * g - s;
				a[i][l-1] = f - g;
				for ( k = l - 1 ; k < n ; k++ ) rv1[k] = a[i][k] / h;
				for ( j = l - 1 ; j < m ; j++) 
				{
					for ( s = 0.0, k = l - 1 ; k < n ; k++ ) 
						s += a[j][k] * a[i][k];
					for ( k = l - 1 ; k < n ; k++ )          
						a[j][k] += s * rv1[k];
				}
				for ( k = l - 1 ; k < n ; k++ ) 
					a[i][k] *= scale;
			}
		}
		anorm=MAX( anorm, ( fabs( w[i] ) + fabs( rv1[i] ) ) );
	}

	for ( i = n - 1 ; i >= 0 ; i-- ) 
	{
		if ( i < n - 1 ) 
		{
			if ( g != 0.0 ) 
			{
				for ( j = l ; j < n ; j++ )
					v[j][i] = ( a[i][j] / a[i][l] ) / g;
				
				for ( j = l ; j < n ; j++ ) 
				{
					for ( s = 0.0, k = l ; k < n ; k++ ) s += a[i][k] * v[k][j];
					for ( k = l ; k < n ; k++ )          v[k][j] += s * v[k][i];
				}
			}

			for ( j = l ; j < n ; j++ ) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}

	for ( i = MIN( m, n ) - 1 ; i >= 0 ; i-- ) 
	{
		l = i + 1;
		g = w[i];
		for ( j = l ; j < n ; j++ ) 
			a[i][j] = 0.0;
		
		if (g != 0.0) 
		{
			g = 1.0 / g;
			for ( j = l ; j < n ; j++ ) 
			{
				for ( s = 0.0, k = l ; k < m ; k++ ) 
					s += a[k][i] * a[k][j];
				f = ( s / a[i][i] ) * g;
				for ( k = i ; k < m ; k++ ) 
					a[k][j] += f * a[k][i];
			}
			for ( j = i ; j < m ; j++ ) 
				a[j][i] *= g;
		} 
		else 
			for ( j = i ; j < m ; j++ ) 
				a[j][i]=0.0;
		++a[i][i];
	}

	for ( k = n - 1 ; k >= 0 ; k-- ) 
	{
		for ( its = 0 ; its < 30 ; its++ ) 
		{
			flag = true;
			for ( l = k ; l >= 0 ; l-- ) 
			{
				nm = l - 1;
				if ( fabs( rv1[l] ) + anorm == anorm ) 
				{
					flag = false;
					break;
				}
				if ( fabs( w[nm] ) + anorm == anorm ) 
					break;
			}

			if ( flag ) 
			{
				c = 0.0;
				s = 1.0;
				for ( i = l - 1 ; i < k + 1 ; i++ ) 
				{
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ( fabs( f ) + anorm == anorm ) 
						break;
					g = w[i];
					h = pythag( f, g );
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for ( j = 0 ; j < m ; j++ ) 
					{
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}
				}
			}

			z=w[k];
			if (l == k) 
			{
				if ( z < 0.0 ) 
				{
					w[k] = -z;
					for ( j = 0 ; j < n ; j++ ) 
						v[j][k] = -v[j][k];
				}
				break;
			}
			if ( its == 29 ) nrerror( "no convergence in 30 svdcmp iterations" );
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0 * h * y );
			g = pythag( f, 1.0 );
			f = ( ( x - z ) * ( x + z ) + h * ( ( y / ( f + SIGN( g,f ) ) ) - h ) ) / x;
			c = s = 1.0;
			for ( j = l ; j <= nm ; j++ ) 
			{
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag( f, h );
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for ( jj = 0 ; jj < n ; jj++ ) 
				{
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z=pythag( f, h );
				w[j] = z;
				if ( z ) 
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for ( jj = 0 ; jj < m ; jj++ ) 
				{
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
}


//
//  SVD Multiplication 
//  Given matrix a[0..m-1][0..n-1] (i.e., U matrix), w[0..n-1], and v[0..n-1][0..n-1] returned from svdcmp 
//  Compute A = U W V^T and store the result in the matrix a.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void svdmul( double **a, double *w, double **v, int m, int n )
{
	int i, j, k;
	double **temp;

	temp = Matrix<double> ( m, n );
	for ( i = 0 ; i < m ; i++ )
		for ( j = 0 ; j < n ; j++ )
			temp[i][j] = 0.0;
	
	for ( i = 0 ; i < m ; i++ )
	{
		for ( j = 0 ; j < n ; j++ )
		{
			for ( k = 0 ; k < n ; k++ )
			{
				if ( k == j )
					temp[i][j] += a[i][k] * w[k];
			}
		}
	}

	for ( i = 0 ; i < m ; i++ )
		for ( j = 0 ; j < n ; j++ )
			a[i][j] = 0.0;

	for ( i = 0 ; i < m ; i++ )
	{
		for ( j = 0 ; j < n ; j++ )
		{
			for ( k = 0 ; k < n ; k++ )
			{
				a[i][j] += temp[i][k] * v[j][k];
			}
		}
	}
	FreeMatrix( temp, m, n );
}


//
//  Toeplitz System
//  Solve the Toeplitz system. The Toeplitz matrix need not be symmetric.
//  y[0..n-1] and r[0..2*n-2] are input array; x[0..n-1] is the output array.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void toeplz( double *r, double *x, double *y, int n )
{
	int j, k, m, m1, m2, n1;
	DP pp, pt1, pt2, qq, qt1, qt2, sd, sgd, sgn, shn, sxn;

	n1=n-1;
	if (r[n1] == 0.0) nrerror("toeplz-1 singular principal minor");
	x[0]=y[0]/r[n1];
	if (n1 == 0) return;
	vector<double> g(n1), h(n1);
	g[0]=r[n1-1]/r[n1];
	h[0]=r[n1+1]/r[n1];
	for (m=0;m<n;m++) {
		m1=m+1;
		sxn = -y[m1];
		sd = -r[n1];
		for (j=0;j<m+1;j++) {
			sxn += r[n1+m1-j]*x[j];
			sd += r[n1+m1-j]*g[m-j];
		}
		if (sd == 0.0) nrerror("toeplz-2 singular principal minor");
		x[m1]=sxn/sd;
		for (j=0;j<m+1;j++)
			x[j] -= x[m1]*g[m-j];
		if (m1 == n1) return;
		sgn = -r[n1-m1-1];
		shn = -r[n1+m1+1];
		sgd = -r[n1];
		for (j=0;j<m+1;j++) {
			sgn += r[n1+j-m1]*g[j];
			shn += r[n1+m1-j]*h[j];
			sgd += r[n1+j-m1]*h[m-j];
		}
		if (sgd == 0.0) nrerror("toeplz-3 singular principal minor");
		g[m1]=sgn/sgd;
		h[m1]=shn/sd;
		k=m;
		m2=(m+2) >> 1;
		pp=g[m1];
		qq=h[m1];
		for (j=0;j<m2;j++) {
			pt1=g[j];
			pt2=g[k];
			qt1=h[j];
			qt2=h[k];
			g[j]=pt1-pp*qt2;
			g[k]=pt2-pp*qt1;
			h[j]=qt1-qq*pt2;
			h[k--]=qt2-qq*pt1;
		}
	}
	nrerror("toeplz - should not arrive here!");
}


//
//  Cholesky Decomposition
//  Given a positive-definite symmetric matrix a[0..n-1][0..n-1], this routine constructs its
//  Cholesky decomposition, A = L L^T. On input, only the upper triangle of a need be given;
//  it is not modified. The Cholesky factor L is returned in the lower triangle of a, except for its
//  diagonal elements which are returned in p[0..n-1].
//
//  Reference: Numerical Recipes in C++ (Press)
//
void choldc(double **a, double *p, int n)
{
	int i,j,k;
	DP sum;

	for (i=0;i<n;i++) {
		for (j=i;j<n;j++) {
			for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0)
					nrerror("choldc failed");
				p[i]=sqrt(sum);
			} else a[j][i]=sum/p[i];
		}
	}
}


//
//  Solve Linear System using Cholesky Decomposition
//  Solve the set of n linear equations A x = b, where a is a positive-definite symmetric matrix.
//  a[0..n-1][0..n-1] and p[0..n-1] are input as the output of the routine choldc. Only the
//  lower subdiagonal portion of a is accessed. b[0..n-1] is input as the right-hand side
//  vector. The solution vector is returned in x[0..n-1]. a, n, and p are not modified and can
//  be left in place for successive calls with different right-hand sides b. b is not modified.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void cholsl(double **a, double *p, double *b, double *x, int n)
{
	int i,k;
	DP sum;

	for (i=0;i<n;i++) {
		for (sum=b[i],k=i-1;k>=0;k--) sum -= a[i][k]*x[k];
		x[i]=sum/p[i];
	}
	for (i=n-1;i>=0;i--) {
		for (sum=x[i],k=i+1;k<n;k++) sum -= a[k][i]*x[k];
		x[i]=sum/p[i];
	}
}


//
//  QR Decomposition
//  Given a matrix a[0..n-1][0..n-1], compute its QR decomposition. The upper triangular matrix R
//  is returned in the upper triangle of a, except for the diagonal elements of R which are returned
//  in d[0..n-1]. The orthogonal matrix Q is represented as a product of n - 1 Householder matrices
//  Q_0..Q_n-2, where Q_j = 1 - u_j * u_j / c_j. The ith component of u_j is zero for i = 0,..,j-1
//  while the nonzero components are returned in a[i][j] for i = j,..,n-1. sing returns as true
//  if singularity is encountered during the decomposition, but the decomposition is still completed
//  in this case; otherwise it returns false.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void qrdcmp(double **a, double *c, double *d, bool &sing, int n)
{
	int i,j,k;
	DP scale,sigma,sum,tau;

	sing=false;
	for (k=0;k<n-1;k++) {
		scale=0.0;
		for (i=k;i<n;i++) scale=MAX(scale,fabs(a[i][k]));
		if (scale == 0.0) {
			sing=true;
			c[k]=d[k]=0.0;
		} else {
			for (i=k;i<n;i++) a[i][k] /= scale;
			for (sum=0.0,i=k;i<n;i++) sum += SQR(a[i][k]);
			sigma=SIGN(sqrt(sum),a[k][k]);
			a[k][k] += sigma;
			c[k]=sigma*a[k][k];
			d[k] = -scale*sigma;
			for (j=k+1;j<n;j++) {
				for (sum=0.0,i=k;i<n;i++) sum += a[i][k]*a[i][j];
				tau=sum/c[k];
				for (i=k;i<n;i++) a[i][j] -= tau*a[i][k];
			}
		}
	}
	d[n-1]=a[n-1][n-1];
	if (d[n-1] == 0.0) sing=true;
}


//
//  Solves the set of n linear equations R x = b, where R is the upper triangular matrix stored in a
//  and d. a[0..n-1][0..n-1] and d[0..n-1] are input as the output of the routien qrdcmp and are
//  not modified. b[0..n-1] is input as the right-hand side vector, and is overwritten with the
//  solution vector on output.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void rsolv(double **a, double *d, double *b, int n)
{
	int i,j;
	DP sum;

	b[n-1] /= d[n-1];
	for (i=n-2;i>=0;i--) {
		for (sum=0.0,j=i+1;j<n;j++) sum += a[i][j]*b[j];
		b[i]=(b[i]-sum)/d[i];
	}
}


//
//  Solve Linear System using QR Decomposition
//  Solves the set of n linear equations A x = b. a[0..n-1][0..n-1], c[0..n-1], and d[0..n-1] are input 
//  as the output of the routine qrdcmp and are not modified. b[0..n-1] is input as the right-hand side 
//  vector, and is overwritten with the solution vector on output.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void qrsolv(double **a, double *c, double *d, double *b, int n)
{
	int i,j;
	DP sum,tau;

	for (j=0;j<n-1;j++) {
		for (sum=0.0,i=j;i<n;i++) sum += a[i][j]*b[i];
		tau=sum/c[j];
		for (i=j;i<n;i++) b[i] -= tau*a[i][j];
	}
	rsolv(a,d,b,n);
}


//
//  Jacobi Transformations of a Symmetric Matrix
//  Compute all eigenvalues and eigenvectors of a real symmetric matrix a[0..n-1][0..n-1].
//  On output, elements of a above the diagonal are destroyed. d[0..n-1] returns the eigenvalues
//  of a. v[0..n-1][0..n-1] is a matrix whose columns contain, on output, the normalized
//  eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
//
//  Reference: Numerical Recipes in C++ (Press)
//
inline void rot(double **a, const DP s, const DP tau, const int i, const int j, const int k, const int l)
{
	DP g,h;

	g=a[i][j];
	h=a[k][l];
	a[i][j]=g-s*(h+g*tau);
	a[k][l]=h+s*(g-h*tau);
}


void jacobi(double **a, double *d, double **v, int &nrot, int n)
{
	int i,j,ip,iq;
	DP tresh,theta,tau,t,sm,s,h,g,c;

	vector<double> b(n), z(n);
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0)
			return;
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
						a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip;j++)
						rot(a,s,tau,j,ip,j,iq);
					for (j=ip+1;j<iq;j++)
						rot(a,s,tau,ip,j,j,iq);
					for (j=iq+1;j<n;j++)
						rot(a,s,tau,ip,j,iq,j);
					for (j=0;j<n;j++)
						rot(v,s,tau,j,ip,j,iq);
					++nrot;
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");
}


//
//  Given the eigenvalues d[0..n-1] and eigenvectors v[0..n-1][0..n-1] as output from
//  jacobi or tqli, this routine sorts the eigenvalues into descending order, and
//  rearranges the columns of v correspondingly. The method is straight insertion.
//
void eigsrt(double *d, double **v, int n)
{
	int i,j,k;
	DP p;

	for (i=0;i<n-1;i++) {
		p=d[k=i];
		for (j=i;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Linear Programming
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Simplex Algorithm
//  
void simp1( double **a, const int mm, int *ll, const int nll, const int iabf, int &kp, DP &bmax )
{
	int k;
	DP test;

	if ( nll <= 0 )
		bmax = 0.0;
	else {
		kp = ll[0];
		bmax = a[mm][kp];
		for ( k = 1 ; k < nll ; k++ ) {
			if ( iabf == 0 )
				test = a[mm][ll[k]] - bmax;
			else
				test = fabs( a[mm][ll[k]] ) - fabs( bmax );
			if ( test > 0.0 ) {
				bmax = a[mm][ll[k]];
				kp = ll[k];
			}
		}
	}
}


void simp2( double **a, const int m, const int n, int &ip, const int kp )
{
	const DP EPS = 1.0e-14;
	int k,i;
	DP qp,q0,q,q1;

	ip = 0;
	for ( i = 0 ; i < m ; i++ )
		if ( a[i+1][kp] < -EPS ) break;
	if ( i + 1 > m ) return;
	q1 = -a[i+1][0] / a[i+1][kp];
	ip = i + 1;
	for ( i = ip ; i < m ; i++ ) {
		if ( a[i+1][kp] < -EPS ) {
			q = -a[i+1][0] / a[i+1][kp];
			if ( q < q1 ) {
				ip = i + 1;
				q1 = q;
			} else if ( q == q1 ) {
				for ( k = 0 ; k < n ; k++ ) {
					qp = -a[ip][k+1] / a[ip][kp];
					q0 = -a[i][k+1] / a[i][kp];
					if ( q0 != qp ) break;
				}
				if ( q0 < qp ) ip = i + 1;
			}
		}
	}
}


void simp3( double **a, const int i1, const int k1, const int ip, const int kp )
{
	int ii,kk;
	DP piv;

	piv = 1.0 / a[ip][kp];
	for ( ii = 0 ; ii < i1 + 1 ; ii++ )
		if ( ii != ip ) {
			a[ii][kp] *= piv;
			for ( kk = 0 ; kk < k1 + 1 ; kk++ )
				if ( kk != kp )
					a[ii][kk] -= a[ip][kk] * a[ii][kp];
		}
	for ( kk = 0 ; kk < k1 + 1 ; kk++ )
		if ( kk != kp ) a[ip][kk] *= -piv;
	a[ip][kp]=piv;
}


//
//  Simplex Algorithm
//  Solve linear program as follows:
//
//  Maximize    z = a00 x0 +a01 x1 + ... + a0,N-1 x_N-1
//  Subject to: ai0 x0 + ai1 x1 + ... ai,N-1 x_N-1 <= bi  i = 1..m1
//              aj0 x0 + aj1 x1 + ... aj,N-1 x_N-1 >= bj  j = m1+1...m1+m2
//              ak0 x0 + ak1 x1 + ... ak,N-1 x_N-1 = bk   k = m1+m2+1...m1+m2+m3
//              x0 >= 0, x1 >= 0, ...x_N-1 >= 0
//  
//  Given the matrix a[i][k], i=0..m+1,k=0..n
//  There are m1 <= constraints, m2 >= constrants, and m3 = constraints
//  and m = m1 + m2 + m3.
//             
//  icase =  0  a finite solution is found
//        = +1  unbounded
//        = -1  no solution
//
//  Reference: Numerical Recipes in C++ (Press)
//
void simplx( double **a, const int m1, const int m2, const int m3, int &icase, int *izrov, int *iposv, int n )
{
	const DP EPS = 1.0e-14;
	int i,k,ip,is,kh,kp,nl1;
	DP q1,bmax;

	int m = m1 + m2 + m3;
	int *l1, *l3;
	l1 = new int[ n + 1 ];
	l3 = new int[ m ];
	nl1 = n;
	for ( k = 0 ; k < n ; k++ ) {
		l1[k] = k + 1;
		izrov[k] = k;
	}
	for ( i = 1 ; i <= m ; i++ ) {
		if ( a[i][0] < 0.0 ) nrerror("Bad input tableau in simplx");
		iposv[i-1] = n + i - 1;
	}
	if ( m2 + m3 != 0 ) {
		for ( i = 0 ; i < m2 ; i++ ) l3[i]=1;
		for ( k = 0 ; k < (n + 1) ; k++ ) {
			q1=0.0;
			for ( i = m1 + 1 ; i < m + 1 ; i++ ) q1 += a[i][k];
			a[m+1][k] = -q1;
		}
		for (;;) {
			simp1( a, m + 1, l1, nl1, 0, kp, bmax );
			if ( bmax <= EPS && a[m+1][0] < -EPS ) {
				icase = -1;
				return;
			} else if ( bmax <= EPS && a[m+1][0] <= EPS ) {
				for ( ip = m1 + m2 + 1 ; ip < m + 1 ; ip++ ) {
					if ( iposv[ip-1] == (ip + n - 1 ) ) {
						simp1( a, ip, l1, nl1, 1, kp, bmax );
						if ( bmax > EPS )
							goto one;
					}
				}
				for (i = m1 + 1 ; i <= m1 + m2 ; i++ ) 
					if ( l3[i - m1 - 1] == 1)
						for ( k = 0 ; k < n + 1 ; k++ )
							a[i][k] = -a[i][k];
				break;
			}
			simp2( a, m, n, ip, kp );
			if ( ip == 0 ) {
				icase = -1;
				return;
			}
	one:	simp3( a, m + 1, n, ip, kp );
			if ( iposv[ip-1] >= (n + m1 + m2 ) ) {
				for ( k = 0 ; k < nl1 ; k++ )
					if ( l1[k] == kp ) break;
				--nl1;
				for ( is = k ; is < nl1 ; is++ ) l1[is] = l1[is+1];
			} else {
				kh = iposv[ip-1] - m1 - n + 1;
				if ( kh >= 1 && l3[kh-1] ) {
					l3[kh-1]=0;
					++a[m+1][kp];
					for ( i = 0 ; i < m + 2 ; i++ )
						a[i][kp]= -a[i][kp];
				}
			}
			swap( izrov[kp-1], iposv[ip-1] );
		}
	}
	for (;;) {
		simp1( a, 0, l1, nl1, 0, kp, bmax );
		if ( bmax <= EPS ) {
			icase=0;
			return;
		}
		simp2( a, m, n, ip, kp );
		if ( ip == 0 ) {
			icase = 1;
			return;
		}
		simp3( a, m, n, ip, kp );
		swap( izrov[kp - 1],iposv[ip - 1] );
	}
	delete [] l1;
	delete [] l3;
}


//
//  Simplex Algorithm
//  Solve linear program using the simplex algorithm (Numerical Recipes)
//
//  Linear programs must be given in standard form:
//  A[0..m-1][0..n-1], b[0..m-1] and c[0..n-1]
// 
//  Example:
//  Maximize    3x1 + x2 + 2x3
//  Subject to   x1 +  x2 + 3x3 <= 30
//              2x1 + 2x2 + 5x3 <= 24
//              4x1 +  x2 + 2x3 <= 36
//              x1, x2, x3 >= 0
//  
//  Set A = 1 1 3   b = 30    c = 3 
//          2 2 5       24        1
//          4 1 2       36        2
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void simplex( double **A, double *b, double *c, int m, int n )
{
	double **a;
	int icase, *izrov, *iposv;
	int i, j, k, l, m1, m2;

	a = Matrix<double> ( m + 3, n + 1 );
	izrov = new int[n+1];
	iposv = new int[m+1];
	
	// Determine number of <= and >= constraints
	m1 = m2 = 0;
	for ( i = 0 ; i < m ; i++ )
	{
		if ( b[i] >= 0 ) m1++;
		else             m2++;
	}

	// Define c (maximize z)
	a[0][0] = 0;
	for ( j = 0 ; j < n ; j++ )
		a[0][j+1] = c[j];

	// Define b
	k = 1; l = m1 + 1;
	for ( i = 0 ; i < m ; i++ )
	{
		if ( b[i] >= 0 )
			a[k++][0] = b[i];
		else
			a[l++][0] = -b[i];
	}
		
	// Define A
	k = 1; l = m1 + 1;
	for ( i = 0 ; i < m ; i++ )
	{
		if ( b[i] >= 0 )
		{
			for ( j = 0 ; j < n ; j++ )
				a[k][j+1] = -A[i][j];
			k++;
		}
		else
		{
			for ( j = 0 ; j < n ; j++ )
				a[l][j+1] = -A[i][j];
			l++;
		}
	}

	// Simplex
	simplx( a, m1, m2, 0, icase, izrov, iposv, n );

	// Solution vector
	for ( i = 0 ; i < n ; i++ ) c[i] = 0;
	for ( i = 0 ; i < m ; i++ )
	{
		if ( iposv[i] < n )
			c[iposv[i]] = a[i+1][0];
	}

	// Report results
	if ( icase > 0 ) 
		cout << "The linear program is infeasible (no solution).\n";
	else if ( icase < 0 )
		cout << "The linear program is unbounded (infinity).\n";
	else
	{
		cout << "The objective value = " << a[0][0] << endl;
		cout << "The solution is: (";
		for ( i = 0 ; i < n ; i++ )
		{
			cout << c[i];
			if ( i < n - 1 )
				cout << ", ";
		}
		cout << ")" << endl;
	}

	FreeMatrix( a, m + 3, n + 1 );
	delete [] izrov;
	delete [] iposv;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Fourier Transforms
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Fast Fourier Transform
//  Replaces data[0,..2*nn-1] by its discrete Fourier transform, if isign is input as 1; or
//  replaces data[0,..2*nn-1] by its inverse discrete Fourier transform, if isign is input as
//  -1. data is a complex array of length nn or, equivalently, a real array of length 2*nn.
//  nn MUST be an integer power of 2 (this is not checked for!) 
//
//  Reference: Numerical Recipes in C++ (Press)
//
void four1( double *data, int nn, int isign )
{
	int n, mmax, m, j, istep, i;
	DP wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;

	n = nn << 1;
	j = 1;
	for ( i = 1 ; i < n ; i += 2 )            // This is the bit-reversal section of the routine.
	{                      
		if ( j > i ) 
		{                        
			swap( data[j-1], data[i-1] );       
			swap( data[j], data[i] );
		}
		m = nn;
		while ( m >= 2 && j > m ) 
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax = 2;
	while ( n > mmax )                               // Outer loop executed lg nn times
	{                        
		istep = mmax << 1;              
		theta = isign * ( 6.28318530717959 / mmax);  // Initialize the trigonometric recurrence.
		wtemp = sin( 0.5 * theta );
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin( theta );
		wr = 1.0;
		wi = 0.0;
		for ( m = 1 ; m < mmax ; m += 2 )            // Here are the two nested inner loops.
		{               
			for ( i = m ; i <= n ; i += istep )      // This is the Danielson-Lanczos formula.
			{           
				j = i + mmax;
				tempr = wr * data[j-1] - wi * data[j];
				tempi = wr * data[j] + wi * data[j-1];
				data[j-1] = data[i-1] - tempr;
				data[j] = data[i] - tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr = ( wtemp = wr ) * wpr - wi * wpi + wr;      // Trigonometric recurrence. 
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}

	if ( isign == -1 )                           // Inverse FFT (Revised) 
	{ 
		for ( i = 0 ; i < n ; i++ )
			data[i] /= (double) nn;
	}
}


//
//  Given two real input arrays data1[0..n-1] and data2[0..n-1], this routine calls 
//  four1 and returns two complex output arrays, fft1[0..2n-1] and fft2[0..2n-1], each
//  of complex length n (i.e., real length 2*n), which contains the dicrete Fourier
//  transforms of the respective data arrays. n MUST be an integer power of 2.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void twofft(double *data1, double *data2, int n, double *fft1, double *fft2)
{
	int nn3,nn2,jj,j;
	DP rep,rem,aip,aim;

	nn3=1+(nn2=n+n);
	for (j=0,jj=0;j<n;j++,jj+=2) {
		fft1[jj]=data1[j];
		fft1[jj+1]=data2[j];
	}
	four1(fft1,n,1);
	fft2[0]=fft1[1];
	fft1[1]=fft2[1]=0.0;
	for (j=2;j<n+1;j+=2) {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j]= -aim;
		fft2[j]=aip;
		fft2[j+1]= -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	}
}


//
//  Fourier Transform of Real-Valued Data
//  Calculates the Fourier transform of a set of n real-valued data points. Replaces this data (which
//  is stored in array data[0..n-1] by the positive frequency half of its complex Fourier transform.
//  The real-valued first and last components of the complex transform are returned as elements
//  data[0] and data[1], respectively. n must be a power of 2. This routine also calculates the
//  inverse transform of a complex data array if it is the transform of real data.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void realft(double *data, int n, const int isign)
{
	int i,i1,i2,i3,i4;
	DP c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793238/DP(n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n/2,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	for (i=1;i<(n>>2);i++) {
		i2=1+(i1=i+i);
		i4=1+(i3=n-i1);
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r= -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4]= -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[0] = (h1r=data[0])+data[1];
		data[1] = h1r-data[1];
	} else {
		data[0]=c1*((h1r=data[0])+data[1]);
		data[1]=c1*(h1r-data[1]);
		four1(data,n/2,-1);
	}
}


//
//  Sine Transform
//  Calculates the sine transform of a set of n real-valued data points stored in array y[0..n-1].
//  The number n must be a power of 2. On exit y is replaced by its transform. This program,
//  without changes, also calculates the inverse sine transform, but in this case the output array
//  should be multiplied by 2/n.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void sinft(double *y, int n)
{
	int j;
	DP sum,y1,y2,theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;

	theta=3.141592653589793238/DP(n);
	wtemp=sin(0.5*theta);
	wpr= -2.0*wtemp*wtemp;
	wpi=sin(theta);
	y[0]=0.0;
	for (j=1;j<(n>>1)+1;j++) {
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
		y1=wi*(y[j]+y[n-j]);
		y2=0.5*(y[j]-y[n-j]);
		y[j]=y1+y2;
		y[n-j]=y1-y2;
	}
	realft(y,n,1);
	y[0]*=0.5;
	sum=y[1]=0.0;
	for (j=0;j<n-1;j+=2) {
		sum += y[j];
		y[j]=y[j+1];
		y[j+1]=sum;
	}
}


//  
//  Cosine Transform
//  Calculates the "staggered" cosine transform of a set y[0..n-1] of real-valued data points.
//  The transformed data replace the original data in array y. n must be a power of 2. Set isign
//  to +1 for a transform, and to -1 for an inverse transform. 
//
//  Reference: Numerical Recipes in C++ (Press)
//
void cosft(double *y, int n, const int isign)
{
	const DP PI=3.141592653589793238;
	int i;
	DP sum,sum1,y1,y2,ytemp,theta,wi=0.0,wi1,wpi,wpr,wr=1.0,wr1,wtemp;

	theta=0.5*PI/n;
	wr1=cos(theta);
	wi1=sin(theta);
	wpr = -2.0*wi1*wi1;
	wpi=sin(2.0*theta);
	if (isign == 1) {
		for (i=0;i<n/2;i++) {
			y1=0.5*(y[i]+y[n-1-i]);
			y2=wi1*(y[i]-y[n-1-i]);
			y[i]=y1+y2;
			y[n-1-i]=y1-y2;
			wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
			wi1=wi1*wpr+wtemp*wpi+wi1;
		}
		realft(y,n,1);
		for (i=2;i<n;i+=2) {
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
			y1=y[i]*wr-y[i+1]*wi;
			y2=y[i+1]*wr+y[i]*wi;
			y[i]=y1;
			y[i+1]=y2;
		}
		sum=0.5*y[1];
		for (i=n-1;i>0;i-=2) {
			sum1=sum;
			sum += y[i];
			y[i]=sum1;
		}
	} else if (isign == -1) {
		ytemp=y[n-1];
		for (i=n-1;i>2;i-=2)
			y[i]=y[i-2]-y[i];
		y[1]=2.0*ytemp;
		for (i=2;i<n;i+=2) {
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
			y1=y[i]*wr+y[i+1]*wi;
			y2=y[i+1]*wr-y[i]*wi;
			y[i]=y1;
			y[i+1]=y2;
		}
		realft(y,n,-1);
		for (i=0;i<n/2;i++) {
			y1=y[i]+y[n-1-i];
			y2=(0.5/wi1)*(y[i]-y[n-1-i]);
			y[i]=0.5*(y1+y2);
			y[n-1-i]=0.5*(y1-y2);
			wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
			wi1=wi1*wpr+wtemp*wpi+wi1;
		}
	}
}


//
//  Convolution using FFT
//  Convolves or deconvolves a real data set data[0..n-1] (including any user-supplied
//  zero padding) with a response function respns[0..m-1], where m is an odd integer <= n.
//  The response function must be stored in wrap-around order: the first half of the array
//  respns contains the impulse response function at positive times, while the second half
//  of the array contains the impulse response function at negative times, counting down
//  from the highest element respns[m-1]. On input isign is +1 for convolution, -1 for
//  deconvolution. The answer is returned in ans[0..n-1]. n MUST be an integer power of two.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void convlv(double *data, double *respns, int n, int m, const int isign, double *ans)
{
	int i,no2;
	DP mag2,tmp;

	double *temp;
	temp = new double[n];
	temp[0]=respns[0];
	for (i=1;i<(m+1)/2;i++) {
		temp[i]=respns[i];
		temp[n-i]=respns[m-i];
	}
	for (i=(m+1)/2;i<n-(m-1)/2;i++)
		temp[i]=0.0;
	for (i=0;i<n;i++)
		ans[i]=data[i];
	realft(ans,n,1);
	realft(temp,n,1);
	no2=n>>1;
	if (isign == 1) {
		for (i=2;i<n;i+=2) {
			tmp=ans[i];
			ans[i]=(ans[i]*temp[i]-ans[i+1]*temp[i+1])/no2;
			ans[i+1]=(ans[i+1]*temp[i]+tmp*temp[i+1])/no2;
		}
		ans[0]=ans[0]*temp[0]/no2;
		ans[1]=ans[1]*temp[1]/no2;
	} else if (isign == -1) {
		for (i=2;i<n;i+=2) {
			if ((mag2=SQR(temp[i])+SQR(temp[i+1])) == 0.0)
				nrerror("Deconvolving at response zero in convlv");
			tmp=ans[i];
			ans[i]=(ans[i]*temp[i]+ans[i+1]*temp[i+1])/mag2/no2;
			ans[i+1]=(ans[i+1]*temp[i]-tmp*temp[i+1])/mag2/no2;
		}
		if (temp[0] == 0.0 || temp[1] == 0.0)
			nrerror("Deconvolving at response zero in convlv");
		ans[0]=ans[0]/temp[0]/no2;
		ans[1]=ans[1]/temp[1]/no2;
	} else nrerror("No meaning for isign in convlv");
	realft(ans,n,-1);
	delete [] temp;
}


//
//  Correlation using FFT
//  Computes the correlation of two real data sets data1[0..n-1] and data2[0..n-1]
//  (including any user-supplied zero pading). n MUST be an integer power of two.
//  The answer is returned in ans[0..n-1] stored in wrap-around order, i.e., 
//  correlations at increasingly negative lags are in ans[n-1] on down to ans[n/2],
//  while correlations at increasingly positive lags are in ans[0] (zero lag) on up
//  to ans[n/2-1]. Sign convention of this routine: if data1 lags data2, i.e.,
//  is shifted to the right of it, then ans will show a peak at positive lags.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void correl(double *data1, double *data2, int n, double *ans)
{
	int no2,i;
	DP tmp;

	double *temp;
	temp = new double[n];
	for (i=0;i<n;i++) {
		ans[i]=data1[i];
		temp[i]=data2[i];
	}
	realft(ans,n,1);
	realft(temp,n,1);
	no2=n>>1;
	for (i=2;i<n;i+=2) {
		tmp=ans[i];
		ans[i]=(ans[i]*temp[i]+ans[i+1]*temp[i+1])/no2;
		ans[i+1]=(ans[i+1]*temp[i]-tmp*temp[i+1])/no2;
	}
	ans[0]=ans[0]*temp[0]/no2;
	ans[1]=ans[1]*temp[1]/no2;
	realft(ans,n,-1);
	delete [] temp;
}


//
//  Linear Prediction
//  Given a real vector of data[0..n-1], this routine returns m linear prediction coefficients
//  as d[0..m-1], and returns the mean square discrepancy as xms.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void memcof(double *data, double &xms, double *d, int n, int m)
{
	int k,j,i;
	DP p=0.0;

	vector<double> wk1(n), wk2(n), wkm(m);
	for (j=0;j<n;j++) p += SQR(data[j]);
	xms=p/n;
	wk1[0]=data[0];
	wk2[n-2]=data[n-1];
	for (j=1;j<n-1;j++) {
		wk1[j]=data[j];
		wk2[j-1]=data[j];
	}
	for (k=0;k<m;k++) {
		DP num=0.0,denom=0.0;
		for (j=0;j<(n-k-1);j++) {
			num += (wk1[j]*wk2[j]);
			denom += (SQR(wk1[j])+SQR(wk2[j]));
		}
		d[k]=2.0*num/denom;
		xms *= (1.0-SQR(d[k]));
		for (i=0;i<k;i++)
			d[i]=wkm[i]-d[k]*wkm[k-1-i];
		if (k == m-1)
			return;
		for (i=0;i<=k;i++) wkm[i]=d[i];
		for (j=0;j<(n-k-2);j++) {
			wk1[j] -= (wkm[k]*wk2[j]);
			wk2[j]=wk2[j+1]-wkm[k]*wk1[j+1];
		}
	}
	nrerror("never get here in memcof.");
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Number-Theoretic Algorithms
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Euclid's Algorithm
//  Given two nonnegative integers a and b, find the greatest common divisors (gcd).
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
long Euclid( long a, long b )
{
	if ( b == 0 ) return a;
	else          return Euclid( b, a % b );
}


// 
//  Extended Euclid's Algorithm
//  Given two nonnegative integers a and b, find the greatest common divisors
//  as well as the integer x and y where gcd(a, b) = a*x + b*y.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//    
long Extended_Euclid( long a, long b, long &x, long & y)   
{   
	long t, d;   
	if ( b == 0 ) 
	{
		x=1;  
		y=0;
		return a;
	}   
	d = Extended_Euclid( b, a % b, x, y );   
    t = x;   
    x = y;   
    y = t - a / b * y;   
    return d;   
}   


//
//  Modular Linear Equations
//  Given nonnegative integers a, b, and n, solve the modular linear equation ax = b (mod n).
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Modular_Linear_Equation_Solver( long a, long b, long n )
{
	long d, x0, i, x, y;
	d = Extended_Euclid( a, n, x, y );
	if ( b % d == 0 )
	{
		x0 = ( x * ( b / d ) ) % n;
		for ( i = 0 ; i < d ; i++ )
			cout << "The " << i + 1 << "th solution is: " << ( x0 + i * ( n / d ) ) % n << endl;
	}
	else
		cout << "No solutions\n";
}


//
//  Modular Exponentiation
//  Given nonnegative integers a, b, and n, compute a^b mod n.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
long Modular_Exponentiation( long a, long b, long n )   
{   
	long c, d, bb, nbits, i;   

	c = 0;
	d = 1;
	bb = b;
	nbits = 0;
	
	while ( bb > 0 )
	{
		bb >>= 1;
		nbits++;
	}

	for ( i = nbits ; i > 0 ; i-- )
	{
		c *= 2;
		d = ( d * d ) % n;
		if ( ( b & ( 1 << ( i - 1 ) ) ) > 0 )
			d = ( d * a ) % n;
	}
	return d;
}  


//
//  Primality Test
//  Given an integer number n, test if the number n is a prime.
//  (Brute-Force Approach)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
bool Primality_Test( const long n )
{
	long i, limit;
	limit = (long) sqrt( (double) n );
	for ( i = 2 ; i <= limit ; i++ )
	{
		if ( n % i == 0 ) return false;
	}
	return true;
}


//
//  Miller-Rabin Primality Test
//  Given an input number n and a randomly chosen base value s, test if
//  the number n is a prime.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
bool Witness( long a, long n )
{
	long t, u, x0, x1, i;
	u = n - 1;
	t = 0;
	while ( u % 2 == 0 )
	{
		u /= 2;
		t++;
	}
	x0 = Modular_Exponentiation( a, u, n );
	for ( i = 0 ; i < t ; i++ )
	{
		x1 = ( x0 * x0 ) % n;
		if ( x1 == 1 && x0 != 1 && x0 != n - 1 )
			return true;
		x0 = x1;
	}
	if ( x0 != 1 )	
		return true;
	return false;
}
   

bool Miller_Rabin( long n, long s )   
{   
	long j, a;
	int seed = 30000;   
	for ( j = 0 ; j < s ; j++ )   
	{     
		a = (long)( ran1( seed ) * (double)( n - 2 ) + 1.0 );   
        if ( Witness( a, n ) ) 
			return false;   
    }   
    return true;         
}   


//
//  Eratosthenes Sieve Method
//  Generate all prime between 2 & the input number n (n > 2).
//  This is a very simple yet elegant method to generate primes.
//
//  Reference: 名題精選百則使用C語言 (冼鏡光)
//
void Sieve( unsigned long n )
{
	unsigned long i, k, m, count, prime;
	bool *sieve;

	m = ( n - 3 ) / 2;
	sieve = new bool[m];
	for ( i = 0 ; i <= m ; i++ )
		sieve[i] = true;

	count = 1;
	for ( i = 0 ; i <= m ; i++ )
	{
		if ( sieve[i] )
		{
			prime = i + i + 3;
			count++;
			for ( k = prime + i ; k <= m ; k += prime )
				sieve[k] = false;
		}
	}

	cout << setw(8) << "2";
	for ( i = 0, k = 2 ; i <= m ; i++ )
	{
		if ( sieve[i] )
		{
			if ( k > 9 )
			{
				cout << endl;
				k = 1;
			}
			cout << setw(8) << 2 * i + 3;
			k++;
		}
	}
	cout << "\nThere are " << count << " primes in total.\n";

	delete [] sieve;
}


//
//  Factorization
//  Given an integer n, this program finds all prime factors by using traditional division method.
//
//  Reference: 名題精選百則使用C語言 (冼鏡光)
//
#define SAVE_FACTOR( fact, exp ) { if(exp > 0) factors[count] = fact, exps[count++] = i; }

void Factorization( unsigned long n )
{
	unsigned long factors[100], exps[100];
	unsigned long work;
	int i, k, count = 0;

	for ( i = 0, work = n ; ( work & 0x01UL ) == 0 && work > 1 ; work >>= 1, i++ )
		;                                          // extract divisor 2
		SAVE_FACTOR( 2, i );                       // save it and its exp.   

	for ( k = 3 ; k <= work ; k+=2 ) 
	{
		for ( i = 0 ; work % k == 0 && work > 1 ; work /= k, i++ )  // for k = 3,5,7,9
			;		                               // extract divisor k 
			SAVE_FACTOR( k, i );                   // save it and its exp. 
	}
	cout << n << " = ";                            // display result          
	for ( i = 0 ; i < count ; i++ )
		cout << factors[i] << "(" << exps[i] << ")";
	cout << endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  String Processing & Matching
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Infix to Postfix conversion
//
//  Input:  One legal infix expression, which contains no blanks 
//          on one line of input.
//  Output: The given infix expression and the equivalent 
//          postfix expression.
//
//  Assumptions:
//     1. The input string is a legal infix expression, which can contain parentheses.
//     2. Only the following operators +, -, *, /, ^ are allowed.
//     3. Every character that is not an operator or a parenthesis is a legal operand.
//     4. The class Stack for a stack of characters is available.
//
int priority_op( char op )
{
	enum precedance { LOWEST, MID, HIGH, HIGHEST };
	if ( op == '^' )              return HIGHEST;
	if ( op == '*' || op == '/' ) return HIGH;
	if ( op == '+' || op == '-' ) return MID;
	return LOWEST;
}

void Infix_to_Postfix( char *in, char *out )
{
	int i, j, k;
	char ch, top;
	bool done;

	k = (int) strlen( in );
	Stack<char> S( k );
	strcpy_s( out, 10, "" );
	
	j = 0;
	for ( i = 0 ; i < k ; i++ )
	{
		ch = in[i];
		switch ( ch )
		{
		// Process Operators
		// Pops operators from stack S whose precedence is >= present operator ch and
		// appends them to the end of the postfix expression.
		case '+':  
		case '-':
		case '*':
		case '/':
		case '^':
			done = false;
			while ( !S.IsEmpty() && !done )
			{
				top = S.Pop();
				if ( ( top != '(' ) && ( priority_op( top ) >= priority_op( ch ) ) )
				{
					out[j++] = top;
				}
				else
				{
					S.Push( top );
					done = true;
				}
			}
			S.Push( ch );
			break;
	
		// Process Openparen - Push onto the stack 
		case '(':
			S.Push( ch );
			break;

		// Process close paren - Pop operators, appending 
		// them to the output string, until a matching 
		// open paren is found
		case ')': 
           while ( ( top = S.Pop() ) != '(' ) // Pop down to the matching open parenthesis
				out[j++] = top;
			break;

		default: // Operand
			out[j++] = ch;
			break;
		} // end switch
	} // end for

	while ( !S.IsEmpty() )
		out[j++] = S.Pop();
	out[j] = NULL;
}


//
//  Naive String Matching 
//  (Brute-Force Approach)
//  
//  Input:  Text[0..n-1] & Pattern[0..m-1]
//  Output: Valid shifts (matched string)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Naive_String_Matching( char *text, char *pattern )
{
	int  n, m, s, i;
	bool found, flag;
	
	n = (int) strlen( text );
	m = (int) strlen( pattern );

	found = false;
	for ( s = 0 ; s <= n - m ; s++ )
	{
		flag = true;
		for ( i = 0 ; i < m ; i++ )
		{
			if ( pattern[i] != text[s+i] )
			{
				flag = false;
				break;
			}
		}

		if ( flag ) 
		{
			cout << "Pattern occurs with shift " << s << endl;
			found = true;
		}
	}

	if ( !found ) cout << "Pattern not found!" << endl;
}


//
//  Knuth-Morris Pratt (KMP) String Matching 
//  
//  Input:  Text[0..n-1] & Pattern[0..m-1]
//  Output: Valid shifts (matched string)
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void KMP_String_Matching( char *text, char *pattern )
{
	int  n, m, k, q, i;
	int  *pi;
	bool found;

	n = (int) strlen( text );
	m = (int) strlen( pattern );
	pi = new int[m];
	found = false;

	// Compute Prefix Function 
	pi[0] = 0;
	k = 0;
	for ( q = 1 ; q < m ; q++ )
	{
		while ( k > 0 && pattern[k] != pattern[q] )
			k = pi[k-1];
		if ( pattern[k] == pattern[q] )
			k++;
		pi[q] = k;
	}

	// Matching
	q = 0;                                         // Number of characters matched
	for ( i = 0 ; i < n ; i++ )                    // Scan the text from left to right   
	{
		while ( q > 0 && pattern[q] != text[i] )
			q = pi[q-1];
		if ( pattern[q] == text[i] )
			q++;
		if ( q == m )
		{
			cout << "Pattern occurs with shift " << i-m+1 << endl;
			q = pi[q-1];
			found = true;
		}
	}
  	if ( !found ) cout << "Pattern not found!" << endl;  
	delete [] pi;
}


//
//  Boyer-Moore (BM) String Matching 
//  
//  Input:  Text[0..n-1] & Pattern[0..m-1]
//  Output: Valid shifts (matched string)
//
//  Reference: Introduction to the Design & Analysis of Algorithms (Levitin)
//
void BM_String_Matching( char *text, char *pattern )
{
	int  n, m, i, j, k;
	int  table[128];   // ASCII
	bool found;

	n = (int) strlen( text );
	m = (int) strlen( pattern );
	found = false;

	// Generate shift table
	for ( i = 0 ; i < 128 ; i++ )
		table[i] = m;

	for ( j = 0 ; j <= m - 2 ; j++ )
		table[pattern[j]] = m - 1 - j;

	// Matching
	i = m - 1;
	while ( i <= n-1 )
	{
		k = 0;
		while ( k <= m - 1 && pattern[m-1-k] == text[i-k] )
			k++;

		if ( k == m )
		{
			cout << "Pattern occurs with shift " << i - m + 1 << endl;
			i++;
			found = true;
		}
		else
			i += table[text[i]];
	}

  	if ( !found ) cout << "Pattern not found!" << endl;  
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Computational Geometry
//
////////////////////////////////////////////////////////////////////////////////////////////////////

struct Point { 
	double x;	
	double y; 
};


//
//  Determining whether two line segments intesect
//  Given two segments p1 p2 & p3 p4, determine if the two line segments intersect.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
double Direction( Point pi, Point pj, Point pk )
{
	double x1, x2, y1, y2;
	x1 = pk.x - pi.x;
	y1 = pk.y - pi.y;
	x2 = pj.x - pi.x;
	y2 = pj.y = pi.y;
	return ( x1 * y2 - x2 * y1 );
}


bool On_Segment( Point pi, Point pj, Point pk )
{
	if ( min( pi.x, pj.x ) <= pk.x && pk.x <= max( pi.x, pj.x ) &&
	     min( pi.y, pj.y ) <= pk.y && pk.y <= max( pi.y, pj.y ) )
	   return true;
	else
		return false;
}

bool Segments_Intersect( Point p1, Point p2, Point p3, Point p4 )
{
	double d1, d2, d3, d4;
	d1 = Direction( p3, p4, p1 );
	d2 = Direction( p3, p4, p2 );
	d3 = Direction( p1, p2, p3 );
	d4 = Direction( p1, p2, p4 );
	if ( ( ( d1 > 0 && d2 < 0 ) || ( d1 < 0 && d2 > 0 ) ) ||
	     ( ( d3 > 0 && d4 < 0 ) || ( d3 < 0 && d4 > 0 ) ) )
	   return true;
	else if ( d1 == 0 && On_Segment( p3, p4, p1 ) )
		return true;
	else if ( d2 == 0 && On_Segment( p3, p4, p2 ) )
		return true;
	else if ( d3 == 0 && On_Segment( p1, p2, p3 ) )
		return true;
	else if ( d4 == 0 && On_Segment( p1, p2, p4 ) )
		return true;
	else
		return false;
}


//
//  Polar Angle
//  Determine the polar angle of point p1 with respect to p0 in degrees.
//
double Polar_Angle( Point p0, Point p1 )
{
	const double Div180PI = 57.295779513082320876798154814105000;
	double x, r, angle;
	x = p1.x - p0.x;
	r = sqrt( ( p1.x - p0.x ) * ( p1.x - p0.x ) +
		      ( p1.y - p0.y ) * ( p1.y - p0.y ) );
	if ( p1.y >= p0.y )
		angle = acos( x / r ) * Div180PI;
	else
		angle = 360.0 - acos( x / r ) * Div180PI;
	return angle;
}


//
//  Jarvis's March (Gift Wrapping)
//  Given a set Q of n points, where n >= 3, find the convex hull.
//  On the output, the convex hull is returned as the first m points
//  of the set Q, where m is the number of vertices in the convex hull.
//
//  Reference: Introduction to Algorithms, 3rd Edition (Cormen)
//
void Jarvis_March( Point *Q, int n, int &m )
{
	int    i, j, min_indx;
	double xmin, ymin;
	double angle_min, angle_thresh;

	min_indx = 0;                 // Find the point with the minimum y-coordinate
	ymin = Q[0].y;
	for ( i = 1 ; i < n ; i++ )                   
	{
		if ( Q[i].y < ymin )              
		{
			ymin = Q[i].y;
			min_indx = i;
		}
	}
	xmin = Q[min_indx].x;         // Leftmost in case of a tie
	for ( i = 0 ; i < n ; i++ )
	{
		if ( Q[i].y == ymin && Q[i].x < xmin )  
			min_indx = i;
	}

	angle_thresh = 0.0;           
	for ( m = 0 ; m < n - 1 ; m++ )     // Find the convex hull
	{
		swap( Q[m].x, Q[min_indx].x );    
		swap( Q[m].y, Q[min_indx].y );

		min_indx  = n;
		angle_min = 360.0;
		for ( j = m + 1 ; j < n ; j++ )
		{
			if ( Polar_Angle(Q[m], Q[j]) >= angle_thresh && Polar_Angle(Q[m], Q[j]) < angle_min )
			{
				angle_min = Polar_Angle(Q[m], Q[j]);
				min_indx = j;
			}
			if ( m >= 2 && Polar_Angle(Q[0], Q[m]) >= angle_thresh && Polar_Angle(Q[0], Q[m]) < angle_min )
			{
				min_indx = n;
			}
		}
		if ( min_indx == n )          
		{
			m++;
			break;
		}	
		else
			angle_thresh = angle_min;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Approximation Algorithms
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Simulated Annealing
//  Solve the Traveling Salesman Problem (TSP) using simulated annealing.
//  The algorithm finds the shortest round-trip path to ncities whose coordinates are in
//  the array x[0..ncity-1], y[0..ncity-1]. The array iorder[0..ncity-1] specifies the
//  order in which the cities are visited. On input, the elements of iorder may be set to any
//  permutation of the numbers 0 to ncity-1. This routine will return the best alternative path
//  it can find.
//
//  Reference: Numerical Recipes in C++ (Press)
//
inline double alen( const double a, const double b, const double c, const double d )
{
	return sqrt( ( b - a ) * ( b - a ) + ( d - c ) * ( d - c ) );
}


double trncst( double *x, double *y, int *iorder, int *n, int ncity )
{
	int j,ii;
	DP de;
	double xx[6],yy[6];

	n[3] = ( n[2] + 1) % ncity;
	n[4] = ( n[0] + ncity - 1) % ncity;
	n[5] = ( n[1] + 1 ) % ncity;
	for ( j = 0 ; j < 6 ; j++ ) 
	{
		ii = iorder[n[j]];
		xx[j] = x[ii];
		yy[j] = y[ii];
	}
	de = -alen( xx[1], xx[5], yy[1], yy[5]);
	de -= alen( xx[0], xx[4], yy[0], yy[4]);
	de -= alen( xx[2], xx[3], yy[2], yy[3]);
	de += alen( xx[0], xx[2], yy[0], yy[2]);
	de += alen( xx[1], xx[3], yy[1], yy[3]);
	de += alen( xx[4], xx[5], yy[4], yy[5]);
	return de;
}


bool metrop( const double de, const double t )
{
	static int gljdum=1;

	return de < 0.0 || ran3( gljdum ) < exp( -de / t );
}


void trnspt( int *iorder, int *n, int ncity )
{
	int m1, m2, m3, nn, j, jj;
	int *jorder;
	jorder = new int[ncity];
	m1 = ( n[1] - n[0] + ncity ) % ncity;
	m2 = ( n[4] - n[3] + ncity ) % ncity;
	m3 = ( n[2] - n[5] + ncity ) % ncity;
	nn = 0;
	for ( j = 0 ; j <= m1 ; j++ ) 
	{
		jj = ( j + n[0] ) % ncity;
		jorder[nn++] = iorder[jj];
	}
	for ( j = 0 ; j <= m2 ; j++ ) 
	{
		jj = ( j + n[3] ) % ncity;
		jorder[nn++] = iorder[jj];
	}
	for ( j = 0 ; j <= m3 ; j++ ) 
	{
		jj = ( j + n[5] ) % ncity;
		jorder[nn++] = iorder[jj];
	}
	for ( j = 0 ; j < ncity ; j++ )
		iorder[j]=jorder[j];

	delete [] jorder;
}


double revcst( double *x, double *y, int *iorder, int *n, int ncity )
{
	int j,ii;
	DP de;
	double xx[4],yy[4];

	n[2] = ( n[0] + ncity - 1 ) % ncity;
	n[3] = ( n[1] + 1 ) % ncity;
	for ( j = 0 ; j < 4 ; j++ ) 
	{
		ii = iorder[n[j]];
		xx[j] = x[ii];
		yy[j] = y[ii];
	}
	de = -alen( xx[0], xx[2], yy[0], yy[2]);
	de -= alen( xx[1], xx[3], yy[1], yy[3]);
	de += alen( xx[0], xx[3], yy[0], yy[3]);
	de += alen( xx[1], xx[2], yy[1], yy[2]);
	return de;
}


void anneal( double *x, double *y, int *iorder, int ncity )
{
	const DP TFACTR = 0.9;
	bool ans;
	int i, i1, i2, idec, idum, j, k, nn, nover, nlimit, nsucc;
	static int n[6];
	unsigned long iseed;
	DP path, de, t;

	nover = 100 * ncity;
	nlimit = 10 * ncity;
	path = 0.0;
	t = 0.5;
	for ( i = 0 ; i < ncity - 1 ; i++ ) 
	{
		i1 = iorder[i];
		i2 = iorder[i+1];
		path += alen( x[i1], x[i2], y[i1], y[i2] );
	}
	i1 = iorder[ncity-1];
	i2 = iorder[0];
	path += alen( x[i1], x[i2], y[i1], y[i2] );
	idum = -1;
	iseed = 111;
	cout << fixed << setprecision(6);
	for ( j = 0 ; j < 100 ; j++ ) 
	{
		nsucc=0;
		for ( k = 0 ; k < nover ; k++ ) 
		{
			do {
				n[0] = int( ncity * ran3( idum ) );
				n[1] = int( ( ncity - 1 ) * ran3( idum ) );
				if ( n[1] >= n[0] ) ++n[1];
				nn = ( n[0] - n[1] + ncity - 1 ) % ncity;
			} while ( nn < 2 );
			idec = irbit1( iseed );
			if ( idec == 0 ) 
			{
				n[2] = n[1] + int( abs( nn - 1 ) * ran3( idum ) ) + 1;
				n[2] %= ncity;
				de = trncst( x, y, iorder, n, ncity );
				ans = metrop( de,t );
				if ( ans ) 
				{
					++nsucc;
					path += de;
					trnspt( iorder, n, ncity );
				}
			} 
			else 
			{
				de = revcst( x, y, iorder, n, ncity );
				ans = metrop( de, t );
				if ( ans ) 
				{
					++nsucc;
					path += de;
					reverse( iorder, n );
				}
			}
			if ( nsucc >= nlimit ) break;
		}

		cout << endl << "T = " << setw(12) << t;
		cout << "	 Path Length = " << setw(12) << path << endl;
		cout << "Successful Moves: " << nsucc << endl;
		
		t *= TFACTR;
		if ( nsucc == 0 ) return;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Integration of Functions
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Integration of Functions (Trapezoidal Rule)
//  This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
//  input as the function to be integrated between limits a and b, also input. When called with
//  n = 1, the routine returns the crudest estimate of integral_a to b f(x) dx. Subsequent calls 
//  with n = 2, 3, ...(in that sequential order) will improve the accuracy by adding 2^n-2 additional
//  interior points.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double trapzd(double func(const double), const double a, const double b, const int n)
{
	DP x,tnm,sum,del;
	static DP s;
	int it,j;

	if(n == 1) {
		return (s=0.5*(b-a)*(func(a)+func(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}


//
//  Integration of Functions (Improved Trapezoidal Rule)
//  Returns the integral of the function func from a to b. The constants EPS can be set to the
//  desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
//  number of steps. Integration is performed by the trapezoidal rule.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double qtrap(double func(const double), const double a, const double b)
{
	const int JMAX=20;
	const DP EPS=1.0e-10;
	int j;
	DP s,olds=0.0;

	for (j=0;j<JMAX;j++) {
		s=trapzd(func,a,b,j+1);
		if (j > 5)
			if (fabs(s-olds) < EPS*fabs(olds) ||
				(s == 0.0 && olds == 0.0)) return s;
		olds=s;
	}
	cout << "Too many steps in routine qtrap" << endl;
	return 0.0;
}


//
//  Integration of Functions (Simpson's Rule)
//  Returns the integral of the function func from a to b. The constants EPS can be set to the
//  desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
//  number of steps. Integration is performed by Simpson's rule.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double qsimp(double func(const double), const double a, const double b)
{
	const int JMAX=20;
	const DP EPS=1.0e-10;
	int j;
	DP s, st, ost=0.0,os=0.0;

	for (j=0;j<JMAX;j++) {
		st=trapzd(func,a,b,j+1);
		s=(4.0*st-ost)/3.0;
		if (j > 5)
			if (fabs(s-os) < EPS*fabs(os) ||
				(s == 0.0 && os == 0.0)) return s;
		os=s;
		ost=st;
	}
	cout << "Too many steps in routine qsimp" << endl;
	return 0.0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Interpolation and Extrapolation
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Polynomial Interpolation
//  Given arrays xa[0..n-1] and ya[0..n-1], and given a value x, this routine returns a value
//  y, and an error estimate dy. If P(x) is the polynomial of degree n - 1 such that 
//  P(xa_i) = ya_i, i = 0..n-1,then the returned value y = P(x).
//
//  Reference: Numerical Recipes in C++ (Press)
//
void polint(double *xa, double *ya, int n, const double x, double &y, double &dy)
{
	int i,m,ns=0;
	DP den,dif,dift,ho,hp,w;

	vector<double> c(n), d(n);
	dif=fabs(x-xa[0]);
	for (i=0;i<n;i++) {                 
		if ((dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ((den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
	}
}


//
//  Rational Function Interpolation
//  Given arrays xa[0..n-1] and ya[0..n-1], and given a value x, this routine returns a value
//  y, and an error estimate dy. The value returned is that of the diagonal rational function,
//  evaluated at x, that passes through the n points (xa_i, ya_i), i=0..n-1.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void ratint(double *xa, double *ya, int n, const double x, double &y, double &dy)
{
	const DP TINY=1.0e-25;
	int m,i,ns=0;
	DP w,t,hh,h,dd;

	vector<double> c(n), d(n);
	hh=fabs(x-xa[0]);
	for (i=0;i<n;i++) {
		h=fabs(x-xa[i]);
		if (h == 0.0) {
			y=ya[i];
			dy=0.0;
			return;
		} else if (h < hh) {
			ns=i;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]+TINY;
	}
	y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<n-m;i++) {
			w=c[i+1]-d[i];
			h=xa[i+m]-x;
			t=(xa[i]-x)*d[i]/h;
			dd=t-c[i+1];
			if (dd == 0.0) nrerror("Error in routine ratint");
			dd=w/dd;
			d[i]=c[i+1]*dd;
			c[i]=t*dd;
		}
		y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
	}
}


//
//  Cubic Spline
//  Given arrays x[0..n-1]and y[0..n-1] containing a tabulated function, i.e., yi = f(xi), with
//  x0 < x1 <...< x_n-1 and given values yp1 and ypn for the first derivative of the interpolating
//  function at points 0 and n-1, respectively, this routine returns an array y2[0..n-1] that
//  contains the second derivatives of the interpolating function at the tabulated points xi. If yp1
//  and/or ypn are equal to 1 x 10^30 or larger, the routine is signaled to set the corresponding
//  boundary condition for a natural spline, with zero second derivative on that boundary.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void spline(double *x, double *y, int n, const double yp1, const double ypn, double *y2)
{
	int i,k;
	DP p,qn,sig,un;

	vector<double> u(n-1);
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}


//
//  Cubic Spline Interpolation
//  Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function (with the xa_i's in
//  order), and given the array y2a[0..n-1], which is the output from spline above, and given
//  a value of x, this routine returns a cubic-spline interpolated value y.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void splint(double *xa, double *ya, double *y2a, int n, const double x, double &y)
{
	int k;
	DP h,b,a;

	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
		+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


//
//  Polynomial Interpolation Coefficients
//  Given arrays xa[0..n-1] and ya[0..n-1] containing a tabulated function ya_i = f(xa_i),
//  this routine returns an array of coefficients cof[0..n-1] such that ya_i = sum cof_j xa_i^j
//
//  Reference: Numerical Recipes in C++ (Press)
//
void polcof(double *xa, double *ya, int n, double *cof)
{
	int k,j,i;
	DP xmin,dy;

	vector<double> x(n), y(n);
	for (j=0;j<n;j++) {
		x[j]=xa[j];
		y[j]=ya[j];
	}
	for (j=0;j<n;j++) {
		double *x_t, *y_t;
		x_t = new double[n-j];
		y_t = new double[n-j];
		for (k=0;k<n-j;k++) {
			x_t[k]=x[k];
			y_t[k]=y[k];
		}
		polint(x_t,y_t,n,0.0,cof[j],dy);
		xmin=1.0e38;
		k = -1;
		for (i=0;i<n-j;i++) {
			if (fabs(x[i]) < xmin) {
				xmin=fabs(x[i]);
				k=i;
			}
			if (x[i] != 0.0)
				y[i]=(y[i]-cof[j])/x[i];
		}
		for (i=k+1;i<n-j;i++) {
			y[i-1]=y[i];
			x[i-1]=x[i];
		}
		delete [] x_t;
		delete [] y_t;
	}
}


//
//  Polynomial Interpolation in Two Dimensions
//  Given arrays x1a[0..m-1] n x2a[0..n-1] of independent variables, and a submatrix of
//  function values ya[0..m-1][0..n-1], tabulated at the grid points defined by x1a and x2a;
//  and given values x1 and x2 of the independent variables; this routine returns an interpolated
//  function value y, and an accuracy indication dy (based only on the interpolation in the x1
//  direction, however).
//
//  Reference: Numerical Recipes in C++ (Press)
//
void polint2(double *x1a, double *x2a, double **ya, int m, int n, double x1, double x2, double &y, double &dy)
{
	int j,k;
	double *ymtmp,*ya_t;
	
	ymtmp = new double[m];
	ya_t = new double[n];
	for(j=0;j<m;j++) {
		for(k=0;k<n;k++) ya_t[k]=ya[j][k];
		polint(x2a,ya_t,n,x2,ymtmp[j],dy);
	}
	polint(x1a,ymtmp,m,x1,y,dy);
	delete [] ymtmp;
	delete [] ya_t;
}


//
//  Bicubic Interpolation
//  Given arrays y[0..3], y1[0..3], y2[0..3], and y12[0..3], containing the function, gradients,
//  and cross derivative at the four grid points of a rectangular grid cell (numbered 
//  counterclockwise from the lower left), and given d1 and d2, the length of the grid cell in
//  the 1- and 2-directions, this routine returns the table c[0..3][0..3] that is used by routine
//  bcuint for bicubic interpolation
//
//  Reference: Numerical Recipes in C++ (Press)
//
void bcucof(double *y, double *y1, double *y2, double *y12, const double d1, const double d2, 
			double **c)
{
	static int wt_d[16][16]=
		{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
		-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
		 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
		 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
		 0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
		 0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
		-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
		 9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
		-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
		 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
		-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
		 4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};
	int l,k,j,i;
	DP xx,d1d2;
	vector<double> cl(16), x(16);

	d1d2=d1*d2;
	for (i=0;i<4;i++) {
		x[i]=y[i];
		x[i+4]=y1[i]*d1;
		x[i+8]=y2[i]*d2;
		x[i+12]=y12[i]*d1d2;
	}
	for (i=0;i<16;i++) {
		xx=0.0;
		for (k=0;k<16;k++) xx += wt_d[i][k]*x[k];
		cl[i]=xx;
	}
	l=0;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) c[i][j]=cl[l++];
}


//
//  Bicubic Interpolation
//  Bicubic interpolation within a grid square. Input quantities are y, y1, y2, y12 (as described
//  in bcucof); x1l and x1u, the lower and upper coordinates of the grid square in the 1-direction;
//  x2l and x2u likewise for the 2-direction; and x1, x2, the coordinates of the desired point for
//  the interpolation. The interpolated function value is returned as ansy, and the interpolated
//  gradient values as ansy1 and ansy2. This routine calls bcucof.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void bcuint(double *y, double *y1, double *y2, double *y12, 
			const double x1l, const double x1u, const double x2l, const double x2u,
			const double x1, const double x2, double &ansy, double &ansy1, double &ansy2)
{
	int i;
	DP t,u,d1,d2;
	double **c;
	c = Matrix<double> (4, 4);

	d1=x1u-x1l;
	d2=x2u-x2l;
	bcucof(y,y1,y2,y12,d1,d2,c);
	if (x1u == x1l || x2u == x2l)
		nrerror("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	ansy=ansy2=ansy1=0.0;
	for (i=3;i>=0;i--) {
		ansy=t*ansy+((c[i][3]*u+c[i][2])*u+c[i][1])*u+c[i][0];
		ansy2=t*ansy2+(3.0*c[i][3]*u+2.0*c[i][2])*u+c[i][1];
		ansy1=u*ansy1+(3.0*c[3][i]*t+2.0*c[2][i])*t+c[1][i];
	}
	ansy1 /= d1;
	ansy2 /= d2;
	FreeMatrix(c, 4, 4);
}


// 
//  Bicubic Spline
//  Given an m by n tabulatd function ya[0..m-1][0..n-1], and tabulated independent variables
//  x2a[0..n-1], this routine cnstructs one-dimensional natural cubic splines of the rows
//  of ya and returns the second-derivatives in the array y2a[0..m-1][0..n-1]. (The array
//  x1a[0..m-1] is included in the argument list merely for consistency with routine
//  splin2.)
//
//  Reference: Numerical Recipes in C++ (Press)
//
void splie2(double *x1a, double *x2a, double **ya, double **y2a, int m, int n)
{
	int j,k;

	double *ya_t, *y2a_t;
	ya_t  = new double[n];
	y2a_t = new double[n];
	for (j=0;j<m;j++) {
		for (k=0;k<n;k++) ya_t[k]=ya[j][k];
		spline(x2a,ya_t,n,1.0e30,1.0e30,y2a_t);
		for (k=0;k<n;k++) y2a[j][k]=y2a_t[k];
	}
	delete [] ya_t;
	delete [] y2a_t;
}


//
//  Bicubic Spline Interpolation
//  Given x1a, x2a, ya, m, n as described in splie2 and y2a as produced by that routine; and
//  given a desired interpolating point x1, x2; this routine returns an interpolated function 
//  value y by bicubic spline interpolation.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void splin2(double *x1a, double *x2a, double **ya, double **y2a, int m, int n,
			const double x1, const double x2, double &y)
{
	int j,k;

	double *ya_t, *y2a_t, *yytmp, *ytmp;
	ya_t  = new double[n];
	y2a_t = new double[n];
	yytmp = new double[m];
	ytmp  = new double[m];
	for (j=0;j<m;j++) {
		for (k=0;k<n;k++) {
			ya_t[k]=ya[j][k];
			y2a_t[k]=y2a[j][k];
		}
		splint(x2a,ya_t,y2a_t,n,x2,yytmp[j]);
	}
	spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
	splint(x1a,yytmp,ytmp,m,x1,y);
	
	delete [] ya_t;
	delete [] y2a_t;
	delete [] yytmp;
	delete [] ytmp;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Special Functions
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Gamma Function
//  Returns the value ln[Gamma(xx)] for xx > 0.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double gammln( double xx )
{
	int j;
	DP x,y,tmp,ser;
	static const DP cof[6]={76.18009172947146,-86.50532032941677,
		 24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
		-0.5395239384953e-5};

	y = x = xx;
	tmp = x + 5.5;
	tmp -= ( x + 0.5 ) * log( tmp );
	ser = 1.000000000190015;
	for ( j = 0 ; j < 6 ; j++ ) ser += cof[j] / ++y;
	return -tmp + log( 2.5066282746310005 * ser / x );
}


//
//  Factorial
//  Returns the value n! as a double-precision number
//
//  Reference: Numerical Recipes in C++ (Press)
//
double factrl( const int n )
{
	static int ntop = 4;
	static DP a[33] = { 1.0, 1.0, 2.0, 6.0, 24.0 };
	int j;

	if ( n < 0 ) nrerror( "Negative factorial in routine factrl" );
	if ( n > 32 ) return exp( gammln( n + 1.0 ) );
	while ( ntop < n ) 
	{
		j = ntop++;
		a[ntop] = a[j] * ntop;
	}
	return a[n];
}


//
//  Natural Log of Factorial
//  Returns ln(n!)
//
//  Reference: Numerical Recipes in C++ (Press)
//
double factln( const int n )
{
	static DP a[101];

	if ( n < 0 ) nrerror( "Negative factorial in routine factln" );
	if ( n <= 1 ) return 0.0;
	if ( n <= 100 )
		return ( a[n] != 0.0 ? a[n] : ( a[n] = gammln( n + 1.0 ) ) );
	else return gammln( n + 1.0 );
}


//
//  Binomial Coefficient
//  Returns the Binomial coefficient C(n, k) as a double-precision number
//
//  Reference: Numerical Recipes in C++ (Press)
//
double bico( const int n, const int k )
{
	return floor( 0.5 + exp( factln( n ) - factln( k ) - factln( n - k ) ) );
}


//
//  Beta Function
//  Returns the value of the beta function B(z, w)
//
//  Reference: Numerical Recipes in C++ (Press)
//
double beta( double z, double w )
{
	return exp( gammln( z ) * gammln( w ) - gammln( z + w ) );
}


//
//  Incomplete Gamma Function
//
//  Reference: Numerical Recipes in C++ (Press)
//
void gser( double &gamser, const double a, const double x, double &gln )
{
	const int ITMAX = 100;
	const DP EPS = numeric_limits<DP>::epsilon();
	int n;
	DP sum, del, ap;

	gln = gammln( a );
	if ( x <= 0.0 ) 
	{
		if ( x < 0.0 ) nrerror( "x less than 0 in routine gser" );
		gamser = 0.0;
		return;
	} 
	else 
	{
		ap = a;
		del = sum = 1.0 / a;
		for ( n = 0 ; n < ITMAX ; n++ ) 
		{
			++ap;
			del *= x / ap;
			sum += del;
			if ( fabs( del ) < fabs( sum ) * EPS ) 
			{
				gamser = sum * exp(-x + a * log( x ) - gln );
				return;
			}
		}
		nrerror( "a too large, ITMAX too small in routine gser" );
		return;
	}
}


//
//  Incomplete Gamma Function
//  Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction 
//  representation as gammcf. Also returns ln[Gamma(a)] as gln.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void gcf( double &gammcf, const double a, const double x, double &gln )
{
	const int ITMAX = 100;
	const DP EPS = numeric_limits<DP>::epsilon();
	const DP FPMIN = numeric_limits<DP>::min() / EPS;
	int i;
	DP an, b, c, d, del, h;

	gln = gammln( a );
	b = x + 1.0 - a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h = d;
	for ( i = 1 ; i <= ITMAX ; i++ ) 
	{
		an = -i * ( i - a );
		b += 2.0;
		d = an * d + b;
		if ( fabs( d ) < FPMIN ) d = FPMIN;
		c = b + an / c;
		if ( fabs( c ) < FPMIN ) c = FPMIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if ( fabs( del - 1.0 ) <= EPS ) break;
	}
	if ( i > ITMAX ) nrerror( "a too large, ITMAX too small in gcf" );
	gammcf = exp( -x + a * log( x ) - gln ) * h;
}


//
//  Incomplete Gamma Function
//  Returns the incomplete gamma function P(a, x)
//
//  Reference: Numerical Recipes in C++ (Press)
//
double gammp( double a, double x )
{
	DP gamser, gammcf, gln;

	if ( x < 0.0 || a <= 0.0 )
		nrerror( "Invalid arguments in routine gammp" );
	if ( x < a + 1.0 ) 
	{
		gser( gamser, a, x, gln );
		return gamser;
	} 
	else 
	{
		gcf( gammcf, a, x, gln );
		return 1.0 - gammcf;
	}
}


//
//  Incomplete Gamma Function
//  Returns the incomplete gamma function Q(a, x) = 1 - P(a, x)
//
//  Reference: Numerical Recipes in C++ (Press)
//
double gammq( double a, double x )
{
	DP gamser, gammcf, gln;

	if ( x < 0.0 || a <= 0.0 )
		nrerror( "Invalid arguments in routine gammq" );
	if ( x < a + 1.0 ) 
	{
		gser( gamser, a, x, gln );
		return 1.0 - gamser;
	} 
	else 
	{
		gcf( gammcf, a, x, gln );
		return gammcf;
	}
}


//
//  Error Function
//  Returns the error function erf(x)
//
//  Reference: Numerical Recipes in C++ (Press)
//
double erff(const double x)
{
	return x < 0.0 ? -gammp( 0.5, x * x ) : gammp( 0.5, x * x );
}


//
//  Complementary Error Function
//  Returns the complementary error function erfc(x)
//
//  Reference: Numerical Recipes in C++ (Press)
//
double erffc( const double x )
{
	return x < 0.0 ? 1.0 + gammp( 0.5, x * x ) : gammq( 0.5, x * x );
}


//
//  Complementary Error Function
//  Returns the complementary error function erfc(x) with fractional error everwhere less than
//  1.2 x 10^-7.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double erfcc( const double x )
{
	DP t, z, ans;

	z = fabs( x );
	t = 1.0 / ( 1.0 + 0.5 * z );
	ans = t * exp( -z * z - 1.26551223 + t * ( 1.00002368 + t * ( 0.37409196 + t * ( 0.09678418 +
		  t * ( -0.18628806 + t * ( 0.27886807 + t * ( -1.13520398 + t * ( 1.48851587 +
		  t * ( -0.82215223 + t * 0.17087277 ) ) ) ) ) ) ) ) );
	return ( x >= 0.0 ? ans : 2.0 - ans );
}


//
//  Incomplete Beta Function
//  Evalutes continued fraction for incomplete beta function by modified Lentz's method.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double betacf( double a, double b, double x )
{
	const int MAXIT = 100;
	const DP EPS = numeric_limits<DP>::epsilon();
	const DP FPMIN = numeric_limits<DP>::min()/EPS;
	int m, m2;
	DP aa, c, d, del, h, qab, qam, qap;

	qab = a + b;
	qap = a + 1.0;
	qam = a - 1.0;
	c = 1.0;
	d = 1.0 - qab * x / qap;
	if ( fabs( d ) < FPMIN ) d = FPMIN;
	d = 1.0 / d;
	h = d;
	for ( m = 1 ; m <= MAXIT ; m++ ) 
	{
		m2 = 2 * m;
		aa = m * ( b - m ) * x / ( ( qam + m2 ) * ( a + m2 ) );
		d = 1.0 + aa * d;
		if ( fabs( d ) < FPMIN ) d = FPMIN;
		c = 1.0 + aa / c;
		if ( fabs( c ) < FPMIN ) c = FPMIN;
		d = 1.0 / d;
		h *= d * c;
		aa = -( a + m ) * ( qab + m ) * x / ( ( a + m2 ) * ( qap + m2 ) );
		d = 1.0 + aa * d;
		if ( fabs( d ) < FPMIN ) d = FPMIN;
		c = 1.0 + aa / c;
		if ( fabs( c ) < FPMIN ) c = FPMIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if ( fabs( del - 1.0 ) <= EPS ) break;
	}
	if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}


//
//  Incomplete Beta Function
//  Returns the incomplete beta function Ix(a, b)
//
//  Reference: Numerical Recipes in C++ (Press)
//
double betai( double a, double b, double x )
{
	DP bt;

	if ( x < 0.0 || x > 1.0 ) nrerror( "Bad x in routine betai" );
	if ( x == 0.0 || x == 1.0 ) bt = 0.0;
	else
		bt = exp( gammln( a + b ) - gammln( a ) - gammln( b ) + a * log( x ) + b * log( 1.0 - x ) );
	if ( x < ( a + 1.0 ) / ( a + b + 2.0 ) )
		return bt * betacf( a, b, x ) / a;
	else
		return 1.0 - bt * betacf( b, a, 1.0 - x ) / b;
}


//
//  Bessel Function 
//  Returns the Bessel function J0(x) for any real x.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double bessj0( const double x )
{
	DP ax, z, xx, y, ans, ans1, ans2;

	if ( ( ax = fabs( x ) ) < 8.0 ) 
	{
		y = x * x;
		ans1 = 57568490574.0 + y * ( -13362590354.0 + y * ( 651619640.7
			   + y * ( -11214424.18 + y * ( 77392.33017 + y * ( -184.9052456 ) ) ) ) );
		ans2 = 57568490411.0 + y * ( 1029532985.0 + y * ( 9494680.718
			   + y * ( 59272.64853 + y * ( 267.8532712 + y * 1.0 ) ) ) );
		ans = ans1 / ans2;
	} 
	else 
	{
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 0.785398164;
		ans1 = 1.0 + y * ( -0.1098628627e-2 + y * ( 0.2734510407e-4
			   + y * ( -0.2073370639e-5 + y * 0.2093887211e-6 ) ) );
		ans2 = -0.1562499995e-1 + y * ( 0.1430488765e-3
			   + y * ( -0.6911147651e-5 + y * ( 0.7621095161e-6
			   - y * 0.934945152e-7 ) ) );
		ans = sqrt( 0.636619772 / ax ) * ( cos( xx ) * ans1 - z * sin( xx ) * ans2 );
	}
	return ans;
}


//
//  Bessel Function
//  Returns the Bessel function Y0(x) for positive x.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double bessy0( const double x )
{
	DP z, xx, y, ans, ans1, ans2;

	if ( x < 8.0 ) 
	{
		y = x * x;
		ans1 = -2957821389.0 + y * ( 7062834065.0 + y * ( -512359803.6
			   + y * ( 10879881.29 + y * ( -86327.92757 + y * 228.4622733 ) ) ) );
		ans2 = 40076544269.0 + y * ( 745249964.8 + y * ( 7189466.438
			   + y * ( 47447.26470 + y * ( 226.1030244 + y * 1.0 ) ) ) );
		ans = ( ans1 / ans2 ) + 0.636619772 * bessj0( x ) * log( x );
	} 
	else 
	{
		z = 8.0 / x;
		y = z * z;
		xx = x - 0.785398164;
		ans1 = 1.0 + y * ( -0.1098628627e-2 + y * ( 0.2734510407e-4
			   + y * ( -0.2073370639e-5 + y * 0.2093887211e-6 ) ) );
		ans2 = -0.1562499995e-1 + y * ( 0.1430488765e-3
			   + y * ( -0.6911147651e-5 + y * ( 0.7621095161e-6
			   + y * ( -0.934945152e-7 ) ) ) );
		ans = sqrt( 0.636619772 / x ) * ( sin( xx ) * ans1 + z * cos( xx ) * ans2 );
	}
	return ans;
}


//
//  Bessel Function
//  Returns the Bessel function J1(x) for any real x.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double bessj1( const double x )
{
	DP ax, z, xx, y, ans, ans1, ans2;

	if ( ( ax = fabs( x ) ) < 8.0 ) 
	{
		y = x * x;
		ans1 = x * ( 72362614232.0 + y * ( -7895059235.0 + y * ( 242396853.1
			   + y * ( -2972611.439 + y * ( 15704.48260 + y * ( -30.16036606 ) ) ) ) ) );
		ans2 = 144725228442.0 + y * ( 2300535178.0 + y * ( 18583304.74
			   + y * ( 99447.43394 + y * ( 376.9991397 + y * 1.0 ) ) ) );
		ans = ans1 / ans2;
	} 
	else 
	{
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 2.356194491;
		ans1 = 1.0 + y * ( 0.183105e-2 + y * ( -0.3516396496e-4
			   + y * ( 0.2457520174e-5 + y * ( -0.240337019e-6 ) ) ) );
		ans2 = 0.04687499995 + y * ( -0.2002690873e-3
			   + y * ( 0.8449199096e-5 + y * ( -0.88228987e-6
			   + y * 0.105787412e-6 ) ) );
		ans = sqrt( 0.636619772 / ax ) * ( cos( xx ) * ans1 - z * sin( xx ) * ans2 );
		if ( x < 0.0 ) ans = -ans;
	}
	return ans;
}


//
//  Bessel Function
//  Returns the Bessel function Y1(x) for any positive x.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double bessy1(const double x)
{
	DP z, xx, y, ans, ans1, ans2;

	if ( x < 8.0 ) 
	{
		y = x * x;
		ans1 = x * ( -0.4900604943e13 + y * ( 0.1275274390e13
			   + y * (-0.5153438139e11 + y * ( 0.7349264551e9
			   + y * (-0.4237922726e7 + y * 0.8511937935e4 ) ) ) ) );
		ans2 = 0.2499580570e14 + y * ( 0.4244419664e12
			   + y * ( 0.3733650367e10 + y * ( 0.2245904002e8
			   + y * ( 0.1020426050e6 + y * ( 0.3549632885e3 + y ) ) ) ) );
		ans = ( ans1 / ans2 ) + 0.636619772 * ( bessj1( x ) * log( x ) - 1.0 / x );
	} 
	else 
	{
		z = 8.0 / x;
		y = z * z;
		xx = x - 2.356194491;
		ans1 = 1.0 + y * ( 0.183105e-2 + y * ( -0.3516396496e-4
			   + y * ( 0.2457520174e-5 + y * ( -0.240337019e-6 ) ) ) );
		ans2 = 0.04687499995 + y * ( -0.2002690873e-3
			   + y * ( 0.8449199096e-5 + y * ( -0.88228987e-6
			   + y * 0.105787412e-6 ) ) );
		ans = sqrt( 0.636619772 / x ) * ( sin( xx ) * ans1 + z * cos( xx ) * ans2 );
	}
	return ans;
}


//
//  Bessel Function
//  Returns the Bessel function Jn(x) for any real x and n >= 2
//
//  Reference: Numerical Recipes in C++ (Press)
//
double bessj( const int n, const double x )
{
	const DP ACC = 160.0;
	const int IEXP = numeric_limits<DP>::max_exponent / 2;
	bool jsum;
	int j, k, m;
	DP ax, bj, bjm, bjp, dum, sum, tox, ans;

	if ( n < 2 ) nrerror( "Index n less than 2 in bessj" );
	ax = fabs( x );
	if ( ax * ax <= 8.0 * numeric_limits<DP>::min() ) 
		return 0.0;
	else if ( ax > DP( n ) ) 
	{
		tox = 2.0 / ax;
		bjm = bessj0( ax );
		bj = bessj1( ax );
		for ( j = 1 ; j < n ; j++ ) 
		{
			bjp = j * tox * bj - bjm;
			bjm = bj;
			bj = bjp;
		}
		ans = bj;
	} 
	else 
	{
		tox = 2.0 / ax;
		m = 2 * ( ( n + int( sqrt( ACC * n ) ) ) / 2 );
		jsum = false;
		bjp = ans = sum = 0.0;
		bj = 1.0;
		for ( j = m ; j > 0 ; j-- ) 
		{
			bjm = j * tox * bj - bjp;
			bjp = bj;
			bj = bjm;
			dum = frexp( bj, &k );
			if ( k > IEXP ) 
			{
				bj = ldexp( bj,-IEXP );
				bjp = ldexp( bjp,-IEXP );
				ans = ldexp( ans,-IEXP );
				sum = ldexp( sum,-IEXP );
			}
			if ( jsum ) sum += bj;
			jsum = !jsum;
			if ( j == n ) ans = bjp;
		}
		sum = 2.0 * sum - bj;
		ans /= sum;
	}
	return x < 0.0 && ( n & 1 ) ? -ans : ans;
}


//
//  Bessel Function
//  Returns the Bessel function Yn(x) for positive x and n >= 2
//
//  Reference: Numerical Recipes in C++ (Press)
//
double bessy( const int n, const double x )
{
	int j;
	DP by, bym, byp, tox;

	if ( n < 2 ) nrerror( "Index n less than 2 in bessy" );
	tox = 2.0 / x;
	by = bessy1( x );
	bym = bessy0( x );
	for ( j = 1 ; j < n ; j++ ) 
	{
		byp = j * tox * by - bym;
		bym = by;
		by = byp;
	}
	return by;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Root Finding and Nonlinear Sets of Equations
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Bracketing
//  Given a function func and an initial guessed range x1 to x2, the routine expands the range
//  geometrically until a root is bracketed by the returned values x1 and x2 (in which case zbrac
//  return true) or until the range becomes unacceptably large (in which case zbrac returns 
//  false).
//
//  Reference: Numerical Recipes in C++ (Press)
//
bool zbrac(double func(const double), double &x1, double &x2)
{
	const int NTRY=50;
	const DP FACTOR=1.6;
	int j;
	DP f1,f2;

	if (x1 == x2) nrerror("Bad initial range in zbrac");
	f1=func(x1);
	f2=func(x2);
	for (j=0;j<NTRY;j++) {
		if (f1*f2 < 0.0) return true;
		if (fabs(f1) < fabs(f2))
			f1=func(x1 += FACTOR*(x1-x2));
		else
			f2=func(x2 += FACTOR*(x2-x1));
	}
	return false;
}


//
//  Bracketing
//  Given a function fx defined on the interval from x1-x2 subdivide the interval into n equally
//  spaced segments, and search for zero crossings of the function. The arrays xb1[0..nb-1] and
//  xb2[0..nb-1] will be filled sequentially with any bracketing pairs that are found, and must be
//  provided with a size nb that is sufficient to hold the maximum number of roots sought. nroot
//  will be set to the number of bracketing pairs actually found.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void zbrak(double fx(const double), const double x1, const double x2, const int n,
	                 double *xb1, double *xb2, int nb, int &nroot)
{
	int i;
	DP x,fp,fc,dx;

	nroot=0;
	dx=(x2-x1)/n;
	fp=fx(x=x1);
	for (i=0;i<n;i++) {
		fc=fx(x += dx);
		if (fc*fp <= 0.0) {
			xb1[nroot]=x-dx;
			xb2[nroot++]=x;
			if(nroot == nb) return;
		}
		fp=fc;
	}
}


//
//  Bisection Method
//  Using bisection, find the root of a function func known to lie between x1 and x2. The root,
//  returned as rtbis, will be refined until its accuracy is +-xacc.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double rtbis(double func(const double), const double x1, const double x2, const double xacc)
{
	const int JMAX=40;
	int j;
	DP dx,f,fmid,xmid,rtb;

	f=func(x1);
	fmid=func(x2);
	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=0;j<JMAX;j++) {
		fmid=func(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}


//
//  False Position Method
//  Using the false position method, find the root of a function func known to lie
//  between x1 and x2. The root, returned as rtflsp, is refined until its accuracy 
//  is +-xacc.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double rtflsp(double func(const double), const double x1, const double x2, const double xacc)
{
	const int MAXIT=30;
	int j;
	DP fl,fh,xl,xh,dx,del,f,rtf;

	fl=func(x1);
	fh=func(x2);
	if (fl*fh > 0.0) nrerror("Root must be bracketed in rtflsp");
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xl=x2;
		xh=x1;
		swap(fl,fh);
	}
	dx=xh-xl;
	for (j=0;j<MAXIT;j++) {
		rtf=xl+dx*fl/(fl-fh);
		f=func(rtf);
		if (f < 0.0) {
			del=xl-rtf;
			xl=rtf;
			fl=f;
		} else {
			del=xh-rtf;
			xh=rtf;
			fh=f;
		}
		dx=xh-xl;
		if (fabs(del) < xacc || f == 0.0) return rtf;
	}
	nrerror("Maximum number of iterations exceeded in rtflsp");
	return 0.0;
}


//
//  Secant Method
//  Using the secant method, find the root of a function func thought to lie betwen
//  x1 and x2. The root, returned as rtsec, is refined until its accuracy is +-xacc.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double rtsec(double func(const double), const double x1, const double x2, const double xacc)
{
	const int MAXIT=30;
	int j;
	DP fl,f,dx,xl,rts;

	fl=func(x1);
	f=func(x2);
	if (fabs(fl) < fabs(f)) {
		rts=x1;
		xl=x2;
		swap(fl,f);
	} else {
		xl=x1;
		rts=x2;
	}
	for (j=0;j<MAXIT;j++) {
		dx=(xl-rts)*f/(f-fl);
		xl=rts;
		fl=f;
		rts += dx;
		f=func(rts);
		if (fabs(dx) < xacc || f == 0.0) return rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsec");
	return 0.0;
}


//
//  Ridders' Method
//  Using Ridders' method, return the root of a function func known to lie between
//  x1 and x2. The root, returned as zriddr, will be refined to an approximate
//  accuracy xacc.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double zriddr(double func(const double), const double x1, const double x2, const double xacc)
{
	const int MAXIT=60;
	const DP UNUSED=-1.11e30;
	int j;
	DP ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

	fl=func(x1);
	fh=func(x2);
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		xl=x1;
		xh=x2;
		ans=UNUSED;
		for (j=0;j<MAXIT;j++) {
			xm=0.5*(xl+xh);
			fm=func(xm);
			s=sqrt(fm*fm-fl*fh);
			if (s == 0.0) return ans;
			xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
			if (fabs(xnew-ans) <= xacc) return ans;
			ans=xnew;
			fnew=func(ans);
			if (fnew == 0.0) return ans;
			if (SIGN(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (SIGN(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (SIGN(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else nrerror("never get here.");
			if (fabs(xh-xl) <= xacc) return ans;
		}
		nrerror("zriddr exceed maximum iterations");
	}
	else {
		if (fl == 0.0) return x1;
		if (fh == 0.0) return x2;
		nrerror("root must be bracketed in zriddr.");
	}
	return 0.0;
}


//
//  Brent's Method
//  Using Brent's method, find the root of a function func known to lie between
//  x1 and x2. The root, returned as zbrent, will be refined until its accuracy
//  is tol.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double zbrent(double func(const double), const double x1, const double x2, const double tol)
{
	const int ITMAX=100;
	const DP EPS=numeric_limits<DP>::epsilon();
	int iter;
	DP a=x1,b=x2,c=x2,d,e,min1,min2;
	DP fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		nrerror("Root must be bracketed in zbrent");
	fc=fb;
	for (iter=0;iter<ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
			fb=func(b);
	}
	nrerror("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}


//
//  Newton-Raphson Method
//  Using the Newton-Raphson method, find the root of a function known to lie in the 
//  interval [x1, x2]. The root rtnewt will be refined until its accuracy is +-xacc. 
//  funcd is a user-supplied routine that returns both the function value and the 
//  first derivative of the function at the point x.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double rtnewt(void funcd(const double, double &, double &), const double x1, const double x2,
			  const double xacc)
{
	const int JMAX=20;
	int j;
	DP df,dx,f,rtn;

	rtn=0.5*(x1+x2);
	for (j=0;j<JMAX;j++) {
		funcd(rtn,f,df);
		dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			nrerror("Jumped out of brackets in rtnewt");
		if (fabs(dx) < xacc) return rtn;
	}
	nrerror("Maximum number of iterations exceeded in rtnewt");
	return 0.0;
}


//
//  Newton-Raphson Method & Bisection
//  Using a combination of Newton-Raphson and bisection, find the root of a function 
//  bracketed between x1 and x2. The root, returned as the function value rtsafe, will
//  be refined until its accuracy is known within +-xacc. funcd is a user-supplied 
//  routine that returns both the function value and the first derivative of the function.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double rtsafe(void funcd(const double, double &, double &), const double x1, const double x2,
		      const double xacc)
{
	const int MAXIT=100;
	int j;
	DP df,dx,dxold,f,fh,fl,temp,xh,xl,rts;

	funcd(x1,fl,df);
	funcd(x2,fh,df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		nrerror("Root must be bracketed in rtsafe");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	funcd(rts,f,df);
	for (j=0;j<MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		funcd(rts,f,df);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}


//
//  Laguerre's Method
//  Given the m+1 complex coefficients a[0..m] of the polynomial a[0] + a[1]x + ... + a[m]x^m, 
//  and given a complex value x, this routine improves x by Laguerre's method until it 
//  converges, within the achievable roundoff limit, to a root of the given polynomial.
//  The number of iterations is returned as its.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void laguer(complex<double> *a, int m, complex<double> &x, int &its)
{
	const int MR=8,MT=10,MAXIT=MT*MR;
	const DP EPS=numeric_limits<DP>::epsilon();
	static const DP frac[MR+1]=
		{0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
	int iter,j;
	DP abx,abp,abm,err;
	complex<DP> dx,x1,b,d,f,g,h,sq,gp,gm,g2;

	for (iter=1;iter<=MAXIT;iter++) {
		its=iter;
		b=a[m];
		err=abs(b);
		d=f=0.0;
		abx=abs(x);
		for (j=m-1;j>=0;j--) {
			f=x*f+d;
			d=x*d+b;
			b=x*b+a[j];
			err=abs(b)+abx*err;
		}
		err *= EPS;
		if (abs(b) <= err) return;
		g=d/b;
		g2=g*g;
		h=g2-2.0*f/b;
		sq=sqrt(DP(m-1)*(DP(m)*h-g2));
		gp=g+sq;
		gm=g-sq;
		abp=abs(gp);
		abm=abs(gm);
		if (abp < abm) gp=gm;
		dx=MAX(abp,abm) > 0.0 ? DP(m)/gp : polar(1+abx,DP(iter));
		x1=x-dx;
		if (x == x1) return;
		if (iter % MT != 0) x=x1;
		else x -= frac[iter/MT]*dx;
	}
	nrerror("too many iterations in laguer");
	return;
}


//
//  Laguerre's Method for Multiple Roots
//  Given the m+1 complx coefficients a[0..m] of the polynomial a[0] + a[1]x + ... + a[m]x^m, 
//  this routin successively calls laguer and find all m complx roots in roots[0..m-1]. 
//  The boolean variable polish should be input as true if polishing (also by Laguerre's 
//  method) is desired, false if the roots will be subsequently polished by other means.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void zroots( complex<double> *a, int m, complex<double> *roots, const bool &polish )
{
	const DP EPS = 1.0e-14;
	int i,its,j,jj;
	complex<DP> x,b,c;

	complex<DP> *ad;
	ad = new complex<double> [m + 1];
	for ( j = 0 ; j <= m ; j++ ) ad[j] = a[j];
	for ( j = m - 1 ; j >= 0 ; j-- ) {
		x = 0.0;
		complex<DP> *ad_v;
		ad_v = new complex<double> [j + 2];
		for ( jj = 0 ; jj < j + 2 ; jj++) ad_v[jj] = ad[jj];
		laguer( ad_v, j + 2, x, its );
		if ( fabs( imag( x ) ) <= 2.0 * EPS * fabs( real( x ) ) )
			x = complex<DP>( real( x ), 0.0 );
		roots[j] = x;
		b = ad[j + 1];
		for ( jj = j ; jj >= 0 ; jj--) {
			c = ad[jj];
			ad[jj] = b;
			b = x * b + c;
		}
		delete [] ad_v;
	}
	if ( polish )
		for ( j = 0 ; j < m ; j++ )
			laguer( a, m, roots[j], its );
	for ( j = 1 ; j < m ; j++ ) {
		x = roots[j];
		for ( i = j - 1 ; i >= 0 ; i-- ) {
			if ( real( roots[i] ) <= real( x ) ) break;
			roots[i + 1] = roots[i];
		}
		roots[i + 1] = x;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Statistical Description of Data
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//  
//  Moment
//  Given an array of data[0..n-1], this routine returns its mean, ave, average deviation adev,
//  standard deviation sdev, variance var, skewness skew, and kurtosis curt.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void moment( double *data, int n, double &ave, double &adev, double &sdev, double &var, double &skew, double &curt )
{
	int j;
	DP ep = 0.0, s, p;

	if ( n <= 1 ) nrerror("n must be at least 2 in moment");
	s = 0.0;
	for ( j = 0 ; j < n ; j++ ) s += data[j];
	ave = s / n;
	adev = var = skew = curt = 0.0;
	for ( j = 0 ; j < n ; j++ ) 
	{
		adev += fabs( s = data[j] - ave );
		ep += s;
		var += ( p = s * s );
		skew += ( p *= s );
		curt += ( p *= s );
	}
	adev /= n;
	var = ( var - ep * ep / n ) / ( n - 1 );
	sdev = sqrt( var );
	if ( var != 0.0 ) 
	{
		skew /= ( n * var * sdev );
		curt = curt / ( n * var * var ) - 3.0;
	} else 
		nrerror("No skew/kurtosis when variance = 0 (in moment)");
}


//
//  Average (Mean) and Variance
//  Given an array data[0..n-1], returns its mean as ave and its variance as var
//
//  Reference: Numerical Recipes in C++ (Press)
//
void avevar( double *data, int n, double &ave, double &var )
{
	DP s,ep;
	int j;

	ave = 0.0;
	for ( j = 0 ; j < n ; j++ ) 
		ave += data[j];
	ave /= n;
	
	var = ep = 0.0;
	for ( j = 0 ; j < n ; j++ ) 
	{
		s = data[j] - ave;
		ep += s;
		var += s * s;
	}
	var = ( var - ep * ep / n ) / ( n - 1 );
}


//  
//  Student's t-test
//  Given the array data1[0..n1-1] and data2[0..n2-1], return student's t as t, and its
//  significance as prob, small value of prob indicating that the arrays have significantly different
//  means. The data arrays are assumed to be drawn from populations with the same true variance.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void ttest( double *data1, double *data2, int n1, int n2, double &t, double &prob )
{
	DP var1, var2, svar, df, ave1, ave2;

	avevar( data1, n1, ave1, var1 );
	avevar( data2, n2, ave2, var2 );
	df = n1 + n2 - 2;
	svar = ( ( n1 - 1 ) * var1 + ( n2 - 1 ) * var2 ) / df;
	t = ( ave1 - ave2 ) / sqrt( svar * ( 1.0 / n1 + 1.0 / n2 ) );
	prob = betai( 0.5 * df, 0.5, df / ( df + t * t ) );
}


//  
//  Student's t-test
//  Given the array data1[0..n1-1] and data2[0..n2-1], this routine return student's t as 
//  t, and its significance as prob, small value of prob indicating that the arrays have significantly 
//  different means. The data arrays are allowed to be drawn from populations with unequal variances.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void tutest( double *data1, double *data2, int n1, int n2, double &t, double &prob )
{
	DP var1, var2, df, ave1, ave2;

	avevar( data1, n1, ave1, var1 );
	avevar( data2, n2, ave2, var2 );
	t = ( ave1 - ave2 ) / sqrt( var1 / n1 + var2 / n2 );
	df = SQR( var1 / n1 + var2 / n2 ) / ( SQR( var1 / n1 ) / ( n1 - 1 ) + SQR( var2 / n2 ) / ( n2 - 1 ) );
	prob = betai( 0.5 * df, 0.5, df / ( df + SQR( t ) ) );
}


//  
//  Student's t-test
//  Given the paired arrays data1[0..n-1] and data2[0..n-1], this routine return student's 
//  t for paired data as t, and its significance as prob, small value of prob indicating a significant
//  different of means. 
//
//  Reference: Numerical Recipes in C++ (Press)
//
void tptest( double *data1, double *data2, int n, double &t, double &prob )
{
	int j;
	DP var1, var2, ave1, ave2, sd, df, cov = 0.0;

	avevar( data1, n, ave1, var1 );
	avevar( data2, n, ave2, var2 );
	for ( j = 0 ; j < n ; j++ )
		cov += ( data1[j] - ave1 ) * ( data2[j] - ave2 );
	cov /= df = n - 1;
	sd = sqrt( ( var1 + var2 - 2.0 * cov ) / n );
	t = ( ave1 - ave2 ) / sd;
	prob = betai( 0.5 * df, 0.5, df / ( df + t * t ) );
}


//  
//  F-test
//  Given the arrays data1[0..n1-1] and data2[0..n2-1], this routine return the value of f,
//  and its significance as prob. Small values of prob indicate that the two arrays have significantly
//  different variances. 
//
//  Reference: Numerical Recipes in C++ (Press)
//
void ftest( double *data1, double *&data2, int n1, int n2, double &f, double &prob )
{
	DP var1, var2, ave1, ave2, df1, df2;

	avevar( data1, n1, ave1, var1 );
	avevar( data2, n2, ave2, var2 );
	if ( var1 > var2 ) 
	{
		f = var1 / var2;
		df1 = n1 - 1;
		df2 = n2 - 1;
	} 
	else 
	{
		f = var2 / var1;
		df1 = n2 - 1;
		df2 = n1 - 1;
	}
	prob = 2.0 * betai( 0.5 * df2, 0.5 * df1, df2 / ( df2 + df1 * f ) );
	if ( prob > 1.0 ) prob = 2.0 - prob;
}


//  
//  Chi-Square Test
//  Given the array bins[0..nbins-1] containting the observed numbers of events, and an array
//  ebins[0..nbins-1] containing the expected numbers of events, and given the number of 
//  constraints knstrn (normally one), this routine returns (trivially) the number of degrees of
//  freedom df, and (nontrivially) the chi-square chsq and the significance prob. A small value of
//  prob indicates a significant difference between the distributions bins and ebins. Note that
//  bins and ebins are both double arrays, although bins will normally contain integer values.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void chsone( double *bins, double *ebins, int nbins, const int knstrn, double &df, double &chsq, double &prob )
{
	int j;
	DP temp;

	df = nbins - knstrn;
	chsq = 0.0;
	for ( j = 0 ; j < nbins ; j++ ) 
	{
		if ( ebins[j] <= 0.0 ) nrerror( "Bad expected number in chsone" );
		temp = bins[j] - ebins[j];
		chsq += temp * temp / ebins[j];
	}
	prob = gammq( 0.5 * df, 0.5 * chsq );
}


//  
//  Chi-Square Test
//  Given the arrays bins1[0..nbins-1] and bins2[0..nbins-1], containting two sets of 
//  binned data, and given the number of constraints knstrn (normally 1 or 0), this routine
//  returns the number of degrees of freedom df, the chi-square chsq, and the significance prob.
//  A small value of prob indicates a significant difference between the distributions bins1 and
//  bins2. Note that bins1 and bins2 are both double arrays, although they will normally
//  contain integer values.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void chstwo( double *bins1, double *bins2, int nbins, const int knstrn, double &df, double &chsq, double &prob )
{
	int j;
	DP temp;

	df = nbins - knstrn;
	chsq = 0.0;
	for ( j = 0 ; j < nbins ; j++ )
	{
		if ( bins1[j] == 0.0 && bins2[j] == 0.0 )
		{
			--df;
		}
		else 
		{
			temp = bins1[j] - bins2[j];
			chsq += temp * temp / ( bins1[j] + bins2[j] );
		}
	}
	prob = gammq( 0.5 * df, 0.5 * chsq );
}


//
//  Kolmogorov-Smirnov Test
//  Given an array data1[0..n1-1], and an array data2[0..n2-1], this routine returns the
//  K-S statistic d, and the significance level prob for the null hypothesis that the data sets are
//  drawn from the same distribution. Small value of prob show that the cumulative distribution
//  function of data1 is significantly different from that of data2. The array data1 and data2
//  are modified by being sorted into ascending order.
//
//  Reference: Numerical Recipes in C++ (Press)
//
double probks( const double alam )
{
	const DP EPS1 = 1.0e-6, EPS2 = 1.0e-16;
	int j;
	DP a2, fac = 2.0, sum = 0.0, term, termbf = 0.0;

	a2 = -2.0 * alam * alam;
	for ( j = 1 ; j <= 100 ; j++ ) 
	{
		term = fac * exp( a2 * j * j );
		sum += term;
		if ( fabs( term ) <= EPS1 * termbf || fabs( term ) <= EPS2 * sum ) return sum;
		fac = -fac;
		termbf = fabs( term );
	}
	return 1.0;
}


void kstwo( double *data1, double *data2, int n1, int n2, double &d, double &prob )
{
	int j1 = 0, j2 =0;
	double d1, d2, dt, en1, en2, en, fn1 = 0.0, fn2 = 0.0;

	QuickSort( data1, n1 );
	QuickSort( data2, n2 );
	en1 = n1;
	en2 = n2;
	d = 0.0;
	while ( j1 < n1 && j2 < n2 ) 
	{
		if ( ( d1 = data1[j1] ) <= ( d2 = data2[j2] ) ) fn1 = j1++ / en1;
		if ( d2 <= d1 ) fn2 = j2++ / en2;
		if ( ( dt = fabs( fn2 - fn1 ) ) > d ) d = dt;
	}
	en = sqrt( en1 * en2 / ( en1 + en2 ) );
	prob=probks( ( en + 0.12 + 0.11 / en ) * d );
}


//
//  Contingency Table
//  Given a two-dimensional contingency table in the form of an array nn[0..ni-1][0..nj-1]
//  of integers, this routine returns the chi-square chisq, the number of degrees of freedom df,
//  the significance level prob (small values indicating a significant association), and two measures
//  of association, Cramer's V (cramrv) and the contingency coefficient C (ccc)
//
//  Reference: Numerical Recipes in C++ (Press)
//
void cntab1( int **nn, int ni, int nj, double &chisq, double &df, double &prob, double &cramrv, double &ccc )
{
	const DP TINY = 1.0e-30;
	int i, j, nnj, nni, minij;
	DP sum = 0.0, expctd, temp;

	vector<double> sumi( ni ), sumj( nj );
	nni = ni;
	nnj = nj;
	for ( i = 0 ; i < ni ; i++ ) 
	{
		sumi[i] = 0.0;
		for ( j = 0 ; j < nj ; j++ ) 
		{
			sumi[i] += nn[i][j];
			sum += nn[i][j];
		}
		if ( sumi[i] == 0.0 ) --nni;
	}

	for ( j = 0 ; j < nj ; j++ ) 
	{
		sumj[j] = 0.0;
		for ( i = 0 ; i < ni ; i++ ) sumj[j] += nn[i][j];
		if ( sumj[j] == 0.0 ) --nnj;
	}
	df = nni * nnj - nni - nnj + 1;

	chisq = 0.0;
	for ( i = 0 ; i < ni ; i++ ) 
	{
		for ( j = 0 ; j < nj ; j++ ) 
		{
			expctd = sumj[j] * sumi[i] / sum;
			temp = nn[i][j] - expctd;
			chisq += temp * temp / ( expctd + TINY );
		}
	}
	prob = gammq( 0.5 * df, 0.5 * chisq );
	minij = nni < nnj ? nni - 1 : nnj - 1;
	cramrv = sqrt( chisq / ( sum * minij ) );
	ccc = sqrt( chisq / ( chisq + sum ) );
}


//
//  Contingency Table
//  Given a two-dimensional contingency table in the form of an integer array nn[i][j], where i 
//  labels the x variable and ranges from 0 to ni-1, j labels the y variable and ranges from 0 to
//  nj-1, this routine returns the entropy h of the whole table, the entropy hx of the x distribution,
//  the entropy hy of the y distribution, the entropy hygx of y given x, the entropy hxgy of x
//  given y, the dependency uygx of y on x, the dependency uxgy of x on y and the symmetrical
//  dependency uxy.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void cntab2( int **nn, int ni, int nj, double &h, double &hx, double &hy, double &hygx, double &hxgy, double &uygx, double &uxgy, double &uxy )
{
	const DP TINY = 1.0e-30;
	int i, j;
	DP sum = 0.0, p;

	vector<double> sumi( ni ), sumj( nj );
	for ( i = 0 ; i < ni ; i++ ) 
	{
		sumi[i] = 0.0;
		for ( j = 0 ; j < nj ; j++ ) 
		{
			sumi[i] += nn[i][j];
			sum += nn[i][j];
		}
	}

	for ( j = 0 ; j < nj ; j++ ) 
	{
		sumj[j]=0.0;
		for ( i = 0 ; i < ni ; i++ )
			sumj[j] += nn[i][j];
	}

	hx = 0.0;
	for ( i = 0 ; i < ni ; i++ )
	{
		if ( sumi[i] != 0.0 ) 
		{
			p = sumi[i] / sum;
			hx -= p * log( p );
		}
	}

	hy = 0.0;
	for (j = 0 ; j < nj ; j++ )
	{
		if ( sumj[j] != 0.0 ) 
		{
			p = sumj[j] / sum;
			hy -= p * log( p );
		}
	}

	h = 0.0;
	for ( i = 0 ; i < ni ; i++ )
	{
		for ( j = 0 ; j < nj ; j++ )
		{
			if ( nn[i][j] != 0 ) 
			{
				p = nn[i][j] / sum;
				h -= p * log( p );
			}
		}
	}
	hygx = h - hx;
	hxgy = h - hy;
	uygx = ( hy - hygx ) / ( hy + TINY );
	uxgy = ( hx - hxgy ) / ( hx + TINY );
	uxy = 2.0 * ( hx + hy - h ) / ( hx + hy + TINY );
}


//
//  Linear Correlation Coefficient
//  Given two arrays x[0..n-1] and y[0..n-1], this routine computes their correlation coefficient
//  r (returned as r), the significance level at which the null hypothesis of zero correlation is 
//  disproved (prob whose small value indicates a significant correlation), and Fisher's z (returned
//  as z), whose value can be used in further statistical tests.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void pearsn( double *x, double *y, int n, double &r, double &prob, double &z )
{
	const DP TINY = 1.0e-20;
	int j;
	DP yt, xt, t, df;
	DP syy = 0.0, sxy = 0.0,sxx = 0.0,ay = 0.0, ax = 0.0;

	for ( j = 0 ; j < n ; j++ ) 
	{
		ax += x[j];
		ay += y[j];
	}
	ax /= n;
	ay /= n;

	for ( j = 0 ; j < n ; j++ ) 
	{
		xt = x[j] - ax;
		yt = y[j] - ay;
		sxx += xt * xt;
		syy += yt * yt;
		sxy += xt * yt;
	}

	r = sxy / ( sqrt( sxx * syy ) + TINY );
	z = 0.5 * log( ( 1.0 + r + TINY ) / ( 1.0 - r + TINY ) );
	df = n - 2;
	t = r * sqrt( df / ( ( 1.0 - r + TINY ) * ( 1.0 + r + TINY ) ) );
	prob = betai( 0.5 * df, 0.5, df / ( df + t * t ) );
}


//
//  Spearman Rank-Order Correlation Coefficient
//  Given two data arrays, data1[0..n-1] and data2[0..n-1], this routine returns their sum-
//  squared difference of ranks as D, the number of standard deviation by which D deviates
//  from its null-hypothesis expected value as zd, the two-sided significance level of this
//  deviation as probd, Spearman's rank correlation as rs, and the two-sided significance level
//  of its deviation from zero as probrs. The external routines crank (below) and quicksort
//  are used. A small value of either probd or probrs indicates a significant correlation (rs 
//  positive) or anticorrelation (rs negative).
//
//  Reference: Numerical Recipes in C++ (Press)
//
void crank( double *w, double &s, int n )
{
	int j = 1, ji, jt;
	DP t, rank;

	s = 0.0;
	while ( j < n ) 
	{
		if ( w[j] != w[j-1] ) 
		{
			w[j-1]=j;
			++j;
		} 
		else 
		{
			for ( jt = j + 1 ; jt <= n && w[jt-1] == w[j-1] ; jt++ );
			rank = 0.5 * ( j + jt - 1 );
			for ( ji = j ; ji <= ( jt - 1 ) ; ji++ )
				w[ji-1] = rank;
			t = jt - j;
			s += ( t * t * t - t );
			j = jt;
		}
	}
	if ( j == n ) w[n-1] = n;
}


void spear( double *data1, double *data2, int n, double &d, double &zd, double &probd, double &rs, double &probrs )
{
	int j;
	DP vard, t, sg, sf, fac, en3n, en, df, aved;

	double *wksp1, *wksp2;
	wksp1 = new double[n];
	wksp2 = new double[n];
	
	for ( j = 0 ; j < n ; j++ ) 
	{
		wksp1[j] = data1[j];
		wksp2[j] = data2[j];
	}
	QuickSort2( wksp1,wksp2,n );
	crank( wksp1,sf,n );
	QuickSort2( wksp2,wksp1,n );
	crank( wksp2,sg,n );
	d = 0.0;
	for ( j = 0 ; j < n ; j++ )
		d += SQR( wksp1[j] - wksp2[j] );
	en = n;
	en3n = en * en * en - en;
	aved = en3n / 6.0 - ( sf + sg ) / 12.0;
	fac = ( 1.0 - sf / en3n ) * ( 1.0 - sg / en3n );
	vard = ( ( en - 1.0 ) * en * en * SQR( en + 1.0 ) / 36.0 ) * fac;
	zd = ( d - aved ) / sqrt( vard );
	probd = erfcc( fabs( zd ) / 1.4142136 );
	rs = ( 1.0 - ( 6.0 / en3n ) * ( d + ( sf + sg ) / 12.0 ) ) / sqrt( fac );
	fac = ( rs + 1.0 ) * ( 1.0 - rs );
	if ( fac > 0.0 ) 
	{
		t = rs * sqrt( ( en - 2.0 ) / fac );
		df = en - 2.0;
		probrs = betai( 0.5 * df, 0.5, df / ( df + t * t ) );
	} else
		probrs = 0.0;
	
	delete [] wksp1;
	delete [] wksp2;
}


//
//  Kendall's Tau
//  Given data arrays data1[0..n-1] and data2[0..n-1], this program returns Kendall's tau as
//  tau, its number of standard deviations from zero as z, and its two-sided significance level
//  as prob. Small values of prob indicate a significant correlation (tau positive) or anti-
//  correlation (tau negative).
//
//  Reference: Numerical Recipes in C++ (Press)
//
void kendl1( double *data1, double *data2, int n, double &tau, double &z, double &prob )
{
	int is = 0, j, k, n2 = 0, n1 = 0;
	DP svar, aa, a2, a1;

	for ( j = 0 ; j < n - 1 ; j++ ) 
	{
		for ( k = j + 1 ; k < n ; k++ ) 
		{
			a1 = data1[j] - data1[k];
			a2 = data2[j] - data2[k];
			aa = a1 * a2;
			if ( aa != 0.0 ) 
			{
				++n1;
				++n2;
				aa > 0.0 ? ++is : --is;
			} 
			else 
			{
				if ( a1 != 0.0 ) ++n1;
				if ( a2 != 0.0 ) ++n2;
			}
		}
	}
	tau = is / ( sqrt( DP( n1 ) ) * sqrt( DP( n2 ) ) );
	svar = ( 4.0 * n + 10.0 ) / ( 9.0 * n * ( n - 1.0 ) );
	z = tau / sqrt( svar );
	prob = erfcc( fabs( z ) / 1.4142136 );
}


//
//  Kendall's Tau for Contingency Tables
//  Given a two-dimensional table tab[0..i-1][0..j-1], such that tab[k][l] contains
//  the number of evens falling in bin k of one variable and bin l of another, this program returns
//  Kendall's tau, its number of standard deviations from zero as z, and its two-sided significance
//  level as prob. Small values of prob indicate a significant correlation (tau positive) or
//  anticorrelation (tau negative) between the two variables. Although tab is a double array, it
//  will normally contain integral values.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void kendl2( int **tab, int i, int j, double &tau, double &z, double &prob )
{
	int k, l, nn, mm, m2, m1, lj, li, kj, ki;
	DP svar, s = 0.0, points, pairs, en2 = 0.0, en1 = 0.0;

	nn = i * j;
	points = tab[i-1][j-1];
	for ( k = 0 ; k <= nn - 2 ; k++ ) 
	{
		ki = ( k / j );
		kj = k - j * ki;
		points += tab[ki][kj];
		for ( l = k + 1 ; l <= nn - 1 ; l++ ) 
		{
			li = l / j;
			lj = l - j * li;
			mm = ( m1 = li - ki ) * ( m2 = lj - kj );
			pairs = tab[ki][kj] * tab[li][lj];
			if ( mm != 0 ) 
			{
				en1 += pairs;
				en2 += pairs;
				s += ( mm > 0 ? pairs : -pairs );
			} 
			else 
			{
				if ( m1 != 0 ) en1 += pairs;
				if ( m2 != 0 ) en2 += pairs;
			}
		}
	}
	tau = s / sqrt( en1 * en2 );
	svar = ( 4.0 * points + 10.0 ) / ( 9.0 * points * ( points - 1.0 ) );
	z = tau / sqrt( svar );
	prob = erfcc( fabs( z ) / 1.4142136 );
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Modeling of Data
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Fitting Data to a Straight Line
//  Given a set of data points x[0..ndata-1], y[0..ndata-1] with individual standard deviation
//  sig[0..ndata-1], fit them to a straight line y = a + bx by minimizing chi-square. Returned
//  are a, b and their respective probable uncertainties siga and sigb, the chi-square chi2, and
//  the goodness-of-fit probability q (that the fit would have chi-square this large or larger). If
//  mwt=false on input, then the standard deviations are assumed to be unavailable: q is returned as
//  1.0 and the normalization of chi2 is to unit standard deviation on all points.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void fit( double *x, double *y, int ndata, double *sig, const bool mwt, 
		  double &a, double &b, double &siga, double &sigb, double &chi2, double &q)
{
	int i;
	DP wt, t, sxoss, sx = 0.0, sy = 0.0, st2 = 0.0, ss, sigdat;

	b = 0.0;
	if ( mwt ) 
	{
		ss = 0.0;
		for ( i = 0 ; i < ndata ; i++ ) 
		{
			wt = 1.0 / SQR( sig[i] );
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	} 
	else 
	{
		for ( i = 0 ; i < ndata ; i++) 
		{
			sx += x[i];
			sy += y[i];
		}
		ss = ndata;
	}
	sxoss = sx / ss;
	
	if ( mwt ) 
	{
		for ( i = 0 ; i < ndata ; i++ ) 
		{
			t = ( x[i] - sxoss ) / sig[i];
			st2 += t * t;
			b += t * y[i] / sig[i];
		}
	} 
	else 
	{
		for ( i = 0 ; i < ndata ; i++ ) 
		{
			t = x[i] - sxoss;
			st2 += t * t;
			b += t * y[i];
		}
	}
	b /= st2;
	a = ( sy - sx * b ) / ss;
	siga = sqrt( (1.0 + sx * sx / ( ss * st2 ) ) / ss );
	sigb = sqrt( 1.0 / st2 );
	chi2 = 0.0;
	q = 1.0;
	if ( !mwt ) 
	{
		for ( i = 0 ; i < ndata ; i++ )
			chi2 += SQR( y[i] - a - b * x[i] );
		sigdat = sqrt( chi2 / ( ndata - 2 ) );
		siga *= sigdat;
		sigb *= sigdat;
	} 
	else 
	{
		for ( i = 0 ; i < ndata ; i++ )
			chi2 += SQR( ( y[i] - a - b * x[i] ) / sig[i] );
		if ( ndata > 2 ) q = gammq( 0.5 * ( ndata - 2 ), 0.5 * chi2 );
	}
}


// 
//  Fitting Data to a Straight Line
//  Given a set of data points x[0..ndata-1], y[0..ndata-1], fit them to a straight line 
//  y = a + bx by minimizing chi-square. This routine is a simplified version of fit.
//
//  Reference: Numerical Recipes in C++ (Press)
//
void fitline( double *x, double *y, int ndata, double &a, double &b )
{
	int i;
	double t, sxoss, sx = 0.0, sy = 0.0, st2 = 0.0, ss;

	b = 0.0;
	for ( i = 0 ; i < ndata ; i++ ) 
	{
		sx += x[i];
		sy += y[i];
	}
	ss = ndata;
	sxoss = sx / ss;
	for ( i = 0 ; i < ndata ; i++ ) 
	{
		t = x[i] - sxoss;
		st2 += t * t;
		b += t * y[i];
	}
	b /= st2;
	a = ( sy - sx * b ) / ss;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Games
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Magic Square
//  Given an order n, generate a n x n magic square A[0..n-1][0..n-1].
//  Fill 1..n^2 to a n x n matrix such that all summations of rows, columns, and two 
//  diagonals are equal.
//
//  There are three kinds of magic squares depending on n
//  (1) Odd          n = 2m + 1
//  (2) Singly-Even  n = 2(2m + 1)
//  (3) Doubly-Even  n = 4m
//
//  Reference: 名題精選百則使用C語言(冼鏡光)
//
void Magic_Square( int **A, int n )
{
	int i, j, k;

	if ( n % 2 != 0 )            // Odd
	{
		for ( i = 0 ; i < n ; i++ )        
		{
			for ( j = 0 ; j < n ; j++ )
			{
				A[i][j] = 0;
			}
		}

		i = 1;                            // i: row index  
		j = n / 2 - 1;                    // j: column index 
		for ( k = 1 ; k <= n * n ; k++ )  // k: number  
		{
			i--;                          // Start from middle
			j++;
			if ( i < 0 && j == n )        // Upper right 
			{
				i = 1;	
				j = n - 1;
			}
			else if ( i < 0 )             // Move up  
			{
				i = n - 1;
			}
			else if ( j == n )            // Move right 
			{
				j = 0;
			}
			if ( A[i][j] != 0 )           // Already occupied 
			{
				i += 2;
				j--;
			}
			A[i][j] = k;                  // Put number
		}
	}
	
	if ( n % 2 == 0 )               
	{
		int m = n / 2;
		if ( m % 2 == 1 )                 // Singly-Even
		{
			i = 0;
			j = m / 2;
			for ( k = 1 ; k <= m * m ; k++ )
			{
				A[i][j] = k;                // Put numbers in A
				A[i+m][j+m] = k + m * m;    // Put numbers in B
				A[i][j+m] = k + 2 * m * m;  // Put numbers in C
				A[i+m][j] = k + 3 * m * m;  // Put numbers in D
				if ( k % m == 0 )
					i++;
				else
				{
					i = (i == 0) ? m - 1 : i - 1;
					j = (j == m - 1) ? 0 : j + 1;
				}
			}
			
			int width = n / 4;
			int width1 = width - 1;
			for ( i = 0 ; i < n / 2 ; i++ )
			{
				if ( i != width )
				{
					for ( j = 0 ; j < width ; j++ )
						swap( A[i][j], A[n/2+i][j] );

					for ( j = 0 ; j < width1 ; j++ )
						swap( A[i][n-1-j], A[n/2+i][n-1-j] );
				}
				else
				{
					for ( j = 1 ; j <= width ; j++ )
						swap( A[width][j], A[n/2+width][j] );

					for ( j = 0 ; j < width1 ; j++ )
						swap( A[width][n-1-j], A[n/2+width][n-1-j] );
				}
			}
		}
		else                     // Doubly-Even
		{
			int marker = -1;
			for ( i = 0 ; i < n ; i++ )
				for ( j = 0 ; j < n ; j++ )
					A[i][j] = 0;

			for ( i = 0 ; i < n / 2 ; i++, marker = -marker )
				for ( j = 0 ; j < n / 2 ; j++, marker = -marker )
					A[i][j] = A[i][n-1-j] = marker;

			k = 1;                       // upward counter
			int inv_k = n * n;           // downward counter
			for ( i = 0 ; i < n / 2 ; i++ )
			{			
				for ( j = 0 ; j < n ; j++ )
				{
					if ( A[i][j] != -1 )  // marked
					{
						A[i][j] = k++;
						A[n-1-i][n-1-j] = inv_k--;
					}
					else                 // unmarked
					{ 
						A[i][j] = inv_k--;
						A[n-1-i][n-1-j] = k++;
					}
				}
			} 
		}
	}
}


//
//  N-Queen
//  Given a board size A[0..n-1][0..n-1] n > 3, compute one solution on the 
//  N queen's problem using a formula of legal positions. 
//
//  Remark: This algorithm is one of the coolest I have ever seen.
//
//  Reference: 名題精選百則使用C語言(冼鏡光)
//
void N_Queen( char **A, int n )
{
	int i, j;
	for ( i = 0 ; i < n ; i++ )
		for ( j = 0 ; j < n ; j++ )
			A[i][j] = '.';

	int m = n % 6;
	if ( m == 0 || m == 4 )  // Class 1
	{
		for ( i = 1 ; i <= n / 2 ; i++ )
			A[2*i-1][i-1] = A[2*i-2][n/2+i-1] = 'Q';
	}
	if ( m == 1 || m == 5 )  // Class 2
	{
		for(i=1;i<=(n-1)/2;i++)
			A[2*i-1][i-1] = 'Q';
		for(i=1;i<=(n+1)/2;i++)
			A[2*i-2][(n-1)/2+i-1] = 'Q';
	}
	if ( m == 3 )            // Class 3 
	{
		for ( i = 1 ; i <= ( n - 3 ) / 2 ; i++ )
			A[2*i+1][i-1] = A[2*i+2][(n-1)/2+i-1] = 'Q';
		A[0][n-2] = A[1][(n-1)/2-1] = A[2][n-1] = 'Q';
	}
	if ( m == 2 )            // Class 4
	{
		if ( n > 8 )
		{
			for ( i = 1 ; i <= 3 ; i++ )
				A[2*i-2][n/2-3+i] = 'Q';
			A[1][n-1] = A[3][0] = A[5][n-2] = 'Q';
			for ( i = 1 ; i <= n / 2 - 3 ; i++ )
				A[2*i+4][i] = A[2*i+5][n/2+i] = 'Q';
		}
		else
		{
			for ( i = 1 ; i <= 4 ; i++ )
				A[2*i-1][i-1] = 'Q';
			A[0][5] = A[2][4] = 'Q';
			A[4][7] = A[6][6] = 'Q';
		}
	}
}


//
//  Knight Tour
//  Given a n x n chess board A[0..n-1][0..n-1] and a starting position,
//  find a knight tour path through each square on the chess board 
//  exactly once. The algorithm uses the backtracking technique without
//  recursion. 
//
//  Warning: For large n (n >= 8), it may take very long time to solve.
//
//  Reference: 名題精選百則使用C語言(冼鏡光)
//
void Knight_Tour( int **A, int n, int start_x, int start_y )
{
	int x, y, new_x, new_y, d, count;
	int offset_x[] = { 2, 1, -1, -2, -2, -1,  1,  2 };
	int offset_y[] = { 1, 2,  2,  1, -1, -2, -2, -1 };
	bool found;	

	Stack<int> path_x( n * n );            // Stack for path x 
	Stack<int> path_y( n * n );            // Stack for path y 
	Stack<int> direction( n * n );         // Stack for direction of (x, y)    

	for ( x = 0 ; x < n ; x++ )            // Initialize the chess board
		for ( y = 0 ; y < n ; y++ )
			A[x][y] = 0;

	path_x.Push( start_x );                // Starting position 
	path_y.Push( start_y );
	direction.Push(0);

	count = 1;                             // Count the number of visits 
	A[start_x][start_y] = count;
	while ( count != n * n )
	{
		found = false;
		d = direction.Pop();
		direction.Push( d );
		while ( d < 8 )                    // Check 8 possible directions
		{
			x = path_x.Top();              // Get current position & compute new position
			y = path_y.Top();
			new_x = x + offset_x[d];
			new_y = y + offset_y[d];

			if ( ( new_x >= 0 && new_x < n && new_y >= 0 && new_y < n ) && ( A[new_x][new_y] == 0 ) )
			{
				path_x.Push( new_x );    // Found a new position available
				path_y.Push( new_y );
				direction.Push( 0 );
				found = true;
				count++;
				A[new_x][new_y] = count;
				break;
			}
			else
			{
				d = direction.Pop();
				d++;
				direction.Push(d);
			}
		}

		if ( !found )                   // All 8 positions checked, but no path exists
		{
			if ( !path_x.IsEmpty() )      
			{
				x = path_x.Pop();       // Backtrack
				y = path_y.Pop();
				d = direction.Pop();
				A[x][y] = 0;
				count--;

				d = direction.Pop();    // Continue to check next direction
				d++;
				direction.Push( d );
			}
			else
			{
				cout << "The knight tour is not found...\n";
				return;
			}  
		}  
	}
}


//
//  Maze 
//  Given a maze[0..nr-1][0..nc-1], find a path to solve the maze problem.
//  The method uses stacks to maintain the (x, y) coordinates being visited.
//
//  Input:
//     maze: integer array of nr x nc
//           0 = path  1 = wall
//     nr: #rows  nc: #cols
//
//     (start_x, start_y): start coordinate 
//     (end_x, end_y):     end coordinate
//     Both start & end coordinates must be 0 (path)
//   
//  Output:
//     The maze is marked with the path.
//     Print the path (solution) for the maze.
//  
#define MAZE_VISIT  2

void Maze( int **maze, int nr, int nc, int start_x, int start_y, int end_x, int end_y ) 
{
	int i, j, k;

	Stack<int> x( nr * nc );
	Stack<int> y( nr * nc );
	i = start_x;  
	j = start_y;
	while ( i != end_x || j != end_y )
	{
		if ( maze[i][j] == 0 )
		{
			x.Push( i );    
			y.Push( j );
			maze[i][j] = MAZE_VISIT;
		}
		if ( j + 1 < nc && maze[i][j+1] == 0 )       // Move Right
			j++; 
		else if ( i + 1 < nr && maze[i+1][j] == 0 )  // Move Down
			i++;
		else if ( j - 1 >= 0 && maze[i][j-1] == 0 )  // Move Left
			j--;
		else if ( i - 1 >= 0 && maze[i-1][j] == 0 )  // Move Up
			i--;
		else                                         // Dead End 
		{
			if ( x.IsEmpty() )                     
			{
				cout << "There is no path for the maze...\n";
				return;
			}
			else                                    // Backtrack
			{
				i = x.Pop();
				j = y.Pop();
				
				if ( ( j + 1 < nc && maze[i][j+1] == 0 ) || ( i + 1 < nr && maze[i+1][j] == 0 ) ||
				     ( j - 1 >= 0 && maze[i][j-1] == 0 ) || ( i - 1 >= 0 && maze[i-1][j] == 0 ) )
				{
					x.Push( i );
					y.Push( j );
				}  
			}
		}
	}

	// Push end coordinate
	x.Push( i ); 
	y.Push( j );
	for ( i = 0 ; i < nr ; i++ )
	{
		for ( j = 0 ; j < nc ; j++ )
		{
			if ( maze[i][j] == MAZE_VISIT )
				maze[i][j] = 0;
		}
	}

	// Output the path for the maze
	int *xx, *yy;
	k = x.NumKeys();
	xx = new int[k];
	yy = new int[k];
	
	i = k - 1;
	while ( !x.IsEmpty() )
	{
		xx[i] = x.Pop();
		yy[i] = y.Pop();
		i--;
	}
	
	cout << "The path for the maze is:" << endl;
	for ( i = 0 ; i < k ; i++ )
	{
		maze[xx[i]][yy[i]] = MAZE_VISIT;
		cout << "(" << xx[i] << "," << yy[i] << ")" << endl;
	}
	delete [] xx;
	delete [] yy;
}     