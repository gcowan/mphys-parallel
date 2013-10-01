//Class for a gaussian function with peak at x=0

#ifndef GAUSSIANFUNC_H
#define GAUSSIANFUNC_H
class gaussianFunc{
    private:
        double sigma;

    public:
        //Constructor
        gaussianFunc(double sigmaVal);
        //Method to evaluate the function for a certain value of x
        double evaluate(double x);
        //Method to get the integral of the function
        double integral();
        //get and set methods 
        void setSigma(double newSigma);
        double getSigma();

  
};
#endif
