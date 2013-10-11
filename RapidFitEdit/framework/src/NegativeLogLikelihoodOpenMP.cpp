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

	//Loop over all data points
	double total = 0.0;
	double integral = 0.0;
	double weight = 1.0;
	double value = 0.0;
	DataPoint* temporaryDataPoint=NULL;
	//bool flag = false;
    ofstream myfile;
    myfile.open("time.txt");
    double time = omp_get_wtime();
    int numberOfPoints = TestDataSet->GetDataNumber();
    int numThreads;
    double * threadTotals;
    int threadNum;
    #pragma omp parallel default(none) shared(threadTotals,TestDataSet,TestPDF,numThreads,numberOfPoints) private(value,temporaryDataPoint,integral,weight,threadNum) 

    {
        //Getting the number of threads
        #pragma omp master
        {
            numThreads = omp_get_num_threads();     
            //setting up an array to store totals for each thread
            threadTotals = new double[numThreads];

        }

        //Getting thread number
        threadNum = omp_get_thread_num();
        
        //Now splitting up the loop
        int lowerBound = (numberOfPoints/numThreads)*threadNum;
        int upperBound = (numberOfPoints/numThreads)*(threadNum+1);

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

		    double pointValue= log( value / integral );

		    if( useWeights ) pointValue *= weight;
		    if( useWeights && weightsSquared ) pointValue *= weight;
		    threadTotals[threadNum]+=pointValue;
       }
       
       //Making sure we iterate through all the data points
       if(threadNum==numThreads-1 && upperBound<numberOfPoints){
           for(int i=upperBound; i<numberOfPoints; i++){
               
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

		    double pointValue= log( value / integral );

		    if( useWeights ) pointValue *= weight;
		    if( useWeights && weightsSquared ) pointValue *= weight;
		    threadTotals[threadNum]+=pointValue;
           }
       }

	}
    //Now summing up all the totals
    for(int i=0; i<numThreads; i++){
        total += threadTotals[i];
    }
    delete threadTotals;
    time = omp_get_wtime()-time;
    myfile<<"The time it took to calculate the negative log likelyhood was: "<<time<<std::endl;
    myfile.close();
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
