#ifndef TEST_CLASS
#define TEST_CLASS

#pragma offload_attribute ( push, target( mic ) )

#include <vector>

#define SIZE 50000
#define ARRAYSIZE 50

class TestClass
{
	public:
		TestClass();
		TestClass( int Index );
		~TestClass();

		void Go();

	private:
		double data[ ARRAYSIZE ];
		double otherData[ ARRAYSIZE ];
//		std::vector< double > data, otherData;
};

#pragma offload_attribute ( pop )

#endif
