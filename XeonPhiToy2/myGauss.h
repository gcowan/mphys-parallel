#ifndef MYGAUSS_H
#define MYGAUSS_H

#pragma offload_attribute (push,target(mic))
#define SIZE 1
class myGauss{
    private:
        double sigma;
        double oneOverTwoSigSq;
        double sqrtTwoPi;
    public:
        myGauss(double sigma);
        myGauss();
        ~myGauss();
        double getSigma();
        void  setSigma(double newSigma);
        double evaluate(double xValue);
        double normValue();

};

#pragma offload_attribute(pop)
#endif
