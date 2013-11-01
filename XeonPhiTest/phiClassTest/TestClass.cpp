#include "TestClass.h"

TestClass::TestClass()
{
}

TestClass::TestClass( int Index )
{
	//data.resize( ARRAYSIZE );
	//otherData.resize( ARRAYSIZE );
	for ( int i = 0; i < ARRAYSIZE; i++ )
	{
		data[ i ] = (double)( i + Index );
		otherData[ i ] = (double)( i * Index );
	}
}

TestClass::~TestClass()
{
}

void TestClass::Go()
{
	for ( int i = 0; i < ARRAYSIZE; i++ )
	{
		data[ i ] *= otherData[ i ];
	}
}
