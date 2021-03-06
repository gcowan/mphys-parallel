
#include "DPJpsiKSpinTwo.hh"
#include "DPBWResonanceShape.hh"
#include "DPBarrierFactor.hh"
#include "DPBarrierL0.hh"
#include "DPBarrierL1.hh"
#include "DPBarrierL2.hh"
#include "DPBarrierL3.hh"
#include "DPHelpers.hh"

#include <iostream>

DPJpsiKSpinTwo::DPJpsiKSpinTwo(int fLB, int fLR, double fmB, double mR, 
                               double gammaR, double m1, double m2, 
                               double RB, double RR, double fmJpsi):
A0(1,0),
Aplus(1,0),
Aminus(1,0)
{
  this->LB=fLB;
  this->LR=fLR;
  this->mB=fmB;
  this->mR=mR;
  this->m1=m1;
  this->m2=m2;
  mJpsi=fmJpsi;
  massShape = new DPBWResonanceShape(mR, gammaR, LR, m1, m2, RR);
  switch (LB)
  {
    case 0: barrierB=new DPBarrierL0(RB);
            break;
    case 1: barrierB=new DPBarrierL1(RB);
            break;
    case 2: barrierB=new DPBarrierL2(RB);
            break;
    case 3: barrierB=new DPBarrierL3(RB);
            break;
    default: std::cout<<"WARNING DPJpsiKSpinTwo (LB): Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrierB=new DPBarrierL0(RB);
             break;
  }
  switch (LR)
  {
    case 0: barrierR=new DPBarrierL0(RR);
            break;
    case 1: barrierR=new DPBarrierL1(RR);
            break;
    case 2: barrierR=new DPBarrierL2(RR);
            break;
    case 3: barrierR=new DPBarrierL3(RR);
            break;
    default: std::cout<<"WARNING DPJpsiKSpinTwo (LR): Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrierR=new DPBarrierL0(RR);
             break;
  }
}

DPJpsiKSpinTwo::~DPJpsiKSpinTwo()
{
  if ( massShape ) 
  {
    delete massShape;
  }
  if ( barrierR )
  {
    delete barrierR;
  }
  if ( barrierB )
  {
    delete barrierB;
  }
}

TComplex DPJpsiKSpinTwo::amplitude(double m23, double cosTheta1, 
                                   double cosTheta2, double phi, int twoLambda,
                                   int twoLambdaPsi)
{
  TComplex result(0,0);

  if ( twoLambda != 2 && twoLambda != -2 ) // Not right option
  {
    return result;
  }
  
  double pB = DPHelpers::daughterMomentum(this->mB, this->mJpsi, m23);
  double pR = DPHelpers::daughterMomentum(m23, this->m1, this->m2);

  double orbitalFactor = TMath::Power(pB/this->mB, this->LB)*
                         TMath::Power(pR/m23, this->LR);

  double barrierFactor = barrierB->barrier( DPHelpers::daughterMomentum(this->mB,
                                            this->mJpsi, this->mR), pB)*
                 barrierR->barrier(DPHelpers::daughterMomentum(this->mR,
                                   this->m1, this->m2), pR);

  TComplex massFactor = this->massShape->massShape(m23);

  // Angular part
  TComplex angular(0,0);
  angular+=Aminus*wigner1.function(cosTheta1,-1,twoLambda/2)*
           wigner2.function(cosTheta2,-1,0)*TComplex::Exp(+1.0*TComplex::I()*phi);
  // A0 has flat phi dependence, to save CPU do not multiply by complicated  1
  angular+=A0*wigner1.function(cosTheta1,0,twoLambda/2)*wigner2.function(cosTheta2,0,0);
  angular+=Aplus*wigner1.function(cosTheta1,1,twoLambda/2)*wigner2.function(cosTheta2,1,0)*
           TComplex::Exp(-1.0*TComplex::I()*phi);

  result = massFactor*barrierFactor*orbitalFactor*angular;

  return result;
}

void DPJpsiKSpinTwo::setHelicityAmplitudes(double magA0, double magAplus,
                     double magAminus, double phaseA0, double phaseAplus,
                     double phaseAminus)
{
  A0=TComplex(magA0*TMath::Cos(phaseA0),magA0*TMath::Sin(phaseA0));
  Aplus=TComplex(magAplus*TMath::Cos(phaseAplus),magAplus*TMath::Sin(phaseAplus));
  Aminus=TComplex(magAminus*TMath::Cos(phaseAminus),magAminus*TMath::Sin(phaseAminus));

}

void DPJpsiKSpinTwo::setResonanceParameters(double mass, double sigma)
{
        massShape->setResonanceParameters( mass, sigma );
}
