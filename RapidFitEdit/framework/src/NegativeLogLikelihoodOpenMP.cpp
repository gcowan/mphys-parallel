/**
        @class NegativeLogLikelihoodOpenMP

        A fit function with evaulate methods for an NLL calculation

        @author Benjamin M Wynne bwynne@cern.ch
        Edited by mark
	@date 2009-10-02
*/

//	RapidFit Headers
#include "NegativeLogLikelihoodOpenMP.h"
#include "ClassLookUp.h"
#include "IPDF.h"

//	System Headers
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <omp.h>

//Default constructor
NegativeLogLikelihoodOpenMP::NegativeLogLikelihoodOpenMP() : FitFunction()
{
	Name="NegativeLogLikelihoodOpenMP";
}

//Destructor
NegativeLogLikelihoodOpenMP::~NegativeLogLikelihoodOpenMP()
{
}

//Return the negative log likelihood for a PDF/DataSet result
double NegativeLogLikelihoodOpenMP::EvaluateDataSet( IPDF * TestPDF, IDataSet * TestDataSet, int number )
{
	(void)number;
	//Initialise the integral caching
	//ResultIntegrator->UpdateIntegralCache( TestDataSet->GetBoundary() );
    omp_set_num_threads(Threads);
    
	double total = 0.0;
	double integral = 0.0;
	double weight = 1.0;
	double value = 0.0;
	DataPoint* temporaryDataPoint=NULL;
	//bool flag = false;
    int numberOfPoints = TestDataSet->GetDataNumber();
    int numThreads;
    double * threadTotals;
    int threadNum;
    int lowerBound,upperBound;
    double pointValue;
    #pragma omp parallel default(none) shared(threadTotals,TestDataSet,TestPDF,numThreads,numberOfPoints) private(value,temporaryDataPoint,integral,weight,threadNum,lowerBound,upperBound,pointValue) 

    {
        
        
        //Getting the number of threads using only the master thread
        #pragma omp master
        {
            numThreads = omp_get_num_threads();     
            //setting up an array to store totals for each thread
            threadTotals = new double[numThreads];
            for(int i=0; i<numThreads; i++){
                threadTotals[i] = 0;
            }

        }
        //Adding this barrier to make sure lower and upper bounds are not calculated before master thread sets numThreads
        #pragma omp barrier
        //Getting thread number
        threadNum = omp_get_thread_num();
        
        //Now splitting up the loop
        lowerBound = (numberOfPoints/numThreads)*threadNum;
        //If last thread iterate upto end of loop
        if(threadNum == numThreads-1){
            upperBound = numberOfPoints;
        }

        else{
            upperBound = (numberOfPoints/numThreads)*(threadNum+1);
        }


        for(int i=lowerBound; i<upperBound; i++){
	
		    temporaryDataPoint = TestDataSet->GetDataPoint(i);
		    value = TestPDF->Evaluate(temporaryDataPoint);

		
		    //Find out the integral
		    integral = TestPDF->Integral( temporaryDataPoint, TestDataSet->GetBoundary() );
		

		    //Get the weight for this DataPoint (event)
		    weight = 1.0;
		    if( useWeights )
		    {   
			    weight = temporaryDataPoint->GetEventWeight();
		    }

		    pointValue= log( value / integral );

		    if( useWeights ) pointValue *= weight;
		    if( useWeights && weightsSquared ) pointValue *= weight;
		    threadTotals[threadNum]+=pointValue;
       }
    
	}
    //Now summing up all the totals
    for(int i=0; i<numThreads; i++){
        total += threadTotals[i];
    }
    delete threadTotals;
	if( false ) cerr << "PDF evaluates to " << value << endl;

	//Return negative log likelihood
	return -total;
}

//Return the up value for error calculations
double NegativeLogLikelihoodOpenMP::UpErrorValue( int Sigma )
{
	if ( Sigma == 1 )
	{
		return 0.5;
	}
	else if ( Sigma == 2 )
	{
		return 2.0;
	}
	else if ( Sigma == 3 )
	{
		return 4.5;
	}
	else
	{
		cerr << "I don't know UP for NLL sigma > 3" << endl;
		exit(1);
	}
}
