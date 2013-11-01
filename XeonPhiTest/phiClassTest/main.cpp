#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>

#include "TestClass.h"

using namespace std;

int main()
{
	//Timer start
	struct timeval start, mid, end;
	long mtime, seconds, useconds;
	gettimeofday( &start, NULL );
	clock_t startTime = clock();

	for ( unsigned int index = 0; index < SIZE; index++ )
	{
		TestClass timeWaster( index );
		timeWaster.Go();
	}

	//Timer middle
	gettimeofday( &mid, NULL );
	clock_t midTime = clock();
	cout << "Scalar done" << endl;

	TestClass inputArray[ SIZE ];
	for ( unsigned int index = 0; index < SIZE; index++ )
	{
		TestClass timeWaster( index );
		inputArray[ index ] = timeWaster;
	}

	#pragma offload target ( mic ) in ( inputArray )
	{
		#pragma omp parallel
		for ( unsigned int index = 0; index < SIZE; index++ )
		{
			inputArray[ index ].Go();
		}
	}

	//Timer end
	gettimeofday( &end, NULL );
	cout << "Scalar seconds: " << (float)( midTime - startTime ) / CLOCKS_PER_SEC << endl;
	cout << "Vector seconds: " << (float)( clock() - midTime ) / CLOCKS_PER_SEC << endl;
	seconds  = mid.tv_sec  - start.tv_sec;
	useconds = mid.tv_usec - start.tv_usec;
	mtime = ( seconds * 1000 + useconds / 1000 ) + 0.5;
	cout << "Scalar: time elapsed " << mtime << " milliseconds" << endl;
	seconds  = end.tv_sec  - mid.tv_sec;
	useconds = end.tv_usec - mid.tv_usec;
	mtime = ( seconds * 1000 + useconds / 1000 ) + 0.5;
	cout << "Vector: time elapsed " << mtime << " milliseconds" << endl;
}
