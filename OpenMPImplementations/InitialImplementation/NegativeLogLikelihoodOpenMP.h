/**
        @class NegativeLogLikelihoodOpenMP

        A fit function with evaulate methods for an NLL calculation

        @author Benjamin M Wynne bwynne@cern.ch
        edited by mark
	@date 2009-10-02
*/

#pragma once
#ifndef NEGATIVE_LOG_LIKELIHOODOPENMP_H
#define NEGATIVE_LOG_LIKELIHOODOPENMP_H

//	RapidFit Headers
#include "FitFunction.h"

class NegativeLogLikelihoodOpenMP : public FitFunction
{
	public:
		NegativeLogLikelihoodOpenMP();
		~NegativeLogLikelihoodOpenMP();

		virtual double UpErrorValue(int);

	protected:
		virtual double EvaluateDataSet( IPDF*, IDataSet*, int );
};
#endif

