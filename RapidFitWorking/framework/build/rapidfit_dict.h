/********************************************************************
* framework/build/rapidfit_dict.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error framework/build/rapidfit_dict.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtablerapidfit_dict();
extern void G__cpp_setup_inheritancerapidfit_dict();
extern void G__cpp_setup_typetablerapidfit_dict();
extern void G__cpp_setup_memvarrapidfit_dict();
extern void G__cpp_setup_globalrapidfit_dict();
extern void G__cpp_setup_memfuncrapidfit_dict();
extern void G__cpp_setup_funcrapidfit_dict();
extern void G__set_cpp_environmentrapidfit_dict();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "framework/include/AcceptReject.h"
#include "framework/include/AngularAcceptance.h"
#include "framework/include/BasePDF.h"
#include "framework/include/BasePDF_Framework.h"
#include "framework/include/BasePDF_MCCaching.h"
#include "framework/include/Blinder.h"
#include "framework/include/ClassLookUp.h"
#include "framework/include/CombinedMistagCalib.h"
#include "framework/include/ComponentPlotter.h"
#include "framework/include/ComponentPlotter_config.h"
#include "framework/include/ComponentRef.h"
#include "framework/include/ConstraintFunction.h"
#include "framework/include/CorrectedCovariance.h"
#include "framework/include/DataPoint.h"
#include "framework/include/DataSetConfiguration.h"
#include "framework/include/DebugClass.h"
#include "framework/include/DoubleFixedResModel.h"
#include "framework/include/DoubleResolutionModel.h"
#include "framework/include/EdStyle.h"
#include "framework/include/ExternalConstraint.h"
#include "framework/include/FitAssembler.h"
#include "framework/include/FitFunction.h"
#include "framework/include/FitFunctionConfiguration.h"
#include "framework/include/FitResult.h"
#include "framework/include/FitResultVector.h"
#include "framework/include/FixedResolutionModel.h"
#include "framework/include/Foam.h"
#include "framework/include/FoamIntegrator.h"
#include "framework/include/FumiliFunction.h"
#include "framework/include/FumiliWrapper.h"
#include "framework/include/FunctionContour.h"
#include "framework/include/GoodnessOfFit.h"
#include "framework/include/IConstraint.h"
#include "framework/include/IDataGenerator.h"
#include "framework/include/IDataSet.h"
#include "framework/include/IFitFunction.h"
#include "framework/include/IMinimiser.h"
#include "framework/include/IMistagCalib.h"
#include "framework/include/InputParsing.h"
#include "framework/include/IntegratorFunction.h"
#include "framework/include/IPDF.h"
#include "framework/include/IPDF_Framework.h"
#include "framework/include/IPDF_MCCaching.h"
#include "framework/include/IPDF_NormalisationCaching.h"
#include "framework/include/IPrecalculator.h"
#include "framework/include/IResolutionModel.h"
#include "framework/include/IStudy.h"
#include "framework/include/JackKnife.h"
#include "framework/include/JPsiPhiDataGenerator.h"
#include "framework/include/main.h"
#include "framework/include/MakeFoam.h"
#include "framework/include/Mathematics.h"
#include "framework/include/MCStudy.h"
#include "framework/include/MemoryDataSet.h"
#include "framework/include/MinimiserConfiguration.h"
#include "framework/include/Minuit2Function.h"
#include "framework/include/Minuit2Wrapper.h"
#include "framework/include/MinuitWrapper.h"
#include "framework/include/MultiThreadedFunctions.h"
#include "framework/include/NegativeLogLikelihood.h"
#include "framework/include/NegativeLogLikelihoodNumerical.h"
#include "framework/include/NegativeLogLikelihoodThreaded.h"
#include "framework/include/NegativeLogLikelihoodThreadedNew.h"
#include "framework/include/NormalisedSumPDF.h"
#include "framework/include/Observable.h"
#include "framework/include/ObservableContinuousConstraint.h"
#include "framework/include/ObservableDiscreteConstraint.h"
#include "framework/include/ObservableRef.h"
#include "framework/include/OutputConfiguration.h"
#include "framework/include/ParameterSet.h"
#include "framework/include/ParseCommandLine.h"
#include "framework/include/PDFConfigurator.h"
#include "framework/include/PDFWithData.h"
#include "framework/include/PerEventAngularAcceptance.h"
#include "framework/include/PerEventResModel.h"
#include "framework/include/PhaseSpaceBoundary.h"
#include "framework/include/PhysicsBottle.h"
#include "framework/include/PhysicsParameter.h"
#include "framework/include/PrecalculatorConfig.h"
#include "framework/include/ProdPDF.h"
#include "framework/include/PseudoObservable.h"
#include "framework/include/RapidFitConfiguration.h"
#include "framework/include/RapidFitIntegrator.h"
#include "framework/include/RapidFitIntegratorConfig.h"
#include "framework/include/RapidFitMatrix.h"
#include "framework/include/RapidRun.h"
#include "framework/include/ResolutionModel.h"
#include "framework/include/ResultFormatter.h"
#include "framework/include/ResultParameter.h"
#include "framework/include/ResultParameterSet.h"
#include "framework/include/ScanParam.h"
#include "framework/include/ScanStudies.h"
#include "framework/include/SimpleMistagCalib.h"
#include "framework/include/SlicedAcceptance.h"
#include "framework/include/StatisticsFunctions.h"
#include "framework/include/StringProcessing.h"
#include "framework/include/SumPDF.h"
#include "framework/include/SWeightPrecalculator.h"
#include "framework/include/Threading.h"
#include "framework/include/ThreadingConfig.h"
#include "framework/include/ToyStudy.h"
#include "framework/include/TripleFixedResModel.h"
#include "framework/include/VectoredFeldmanCousins.h"
#include "framework/include/XMLConfigReader.h"
#include "framework/include/XMLTag.h"
#include "pdfs/include/Bd2JpsiKstar_sWave.h"
#include "pdfs/include/Bd2JpsiKstar_sWave_Fs.h"
#include "pdfs/include/Bd2JpsiKstar_sWave_Fs_withAcc.h"
#include "pdfs/include/Bd2JpsiKstar_sWave_Fscopy.h"
#include "pdfs/include/Bd2JpsiKstar_sWave_KpiBins.h"
#include "pdfs/include/Bd2JpsiKstar_sWave_rTerms.h"
#include "pdfs/include/Bd2JpsiKstar_withTimeRes_withAverageAngAcc.h"
#include "pdfs/include/Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms.h"
#include "pdfs/include/BoxPDF.h"
#include "pdfs/include/Bs2DsPi.h"
#include "pdfs/include/Bs2DsPi_acc.h"
#include "pdfs/include/Bs2DsPi_lowmassbkg.h"
#include "pdfs/include/Bs2DsPi_lowmassbkg_updated.h"
#include "pdfs/include/Bs2DsPi_mistagParameter.h"
#include "pdfs/include/Bs2DsPiBkg_withTimeRes.h"
#include "pdfs/include/Bs2DsPiMassSignal.h"
#include "pdfs/include/Bs2Jpsifzero_Signal_v5.h"
#include "pdfs/include/Bs2Jpsifzero_Signal_v5_forcombined.h"
#include "pdfs/include/Bs2Jpsifzero_Signal_v5a.h"
#include "pdfs/include/Bs2Jpsifzero_Signal_v6.h"
#include "pdfs/include/Bs2Jpsifzero_SignalAlt_BaseClass_dev.h"
#include "pdfs/include/Bs2Jpsifzero_SignalAlt_MO_dev.h"
#include "pdfs/include/Bs2Jpsifzero_SignalAlt_MP_dev.h"
#include "pdfs/include/Bs2JpsiPhi_Angluar_Terms.h"
#include "pdfs/include/Bs2JpsiPhi_Signal_v5.h"
#include "pdfs/include/Bs2JpsiPhi_Signal_v5_old.h"
#include "pdfs/include/Bs2JpsiPhi_Signal_v6.h"
#include "pdfs/include/Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4.h"
#include "pdfs/include/Bs2JpsiPhi_SignalAlt_MO_1angle_v4.h"
#include "pdfs/include/Bs2JpsiPhi_SignalAlt_MO_v4.h"
#include "pdfs/include/Bs2JpsiPhiLongLivedBkg.h"
#include "pdfs/include/Bs2JpsiPhiLongLivedBkg_II.h"
#include "pdfs/include/Bs2JpsiPhiLongLivedBkg_withTimeRes.h"
#include "pdfs/include/Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist.h"
#include "pdfs/include/Bs2JpsiPhiMassBkg.h"
#include "pdfs/include/Bs2JpsiPhiMassBkgLL.h"
#include "pdfs/include/Bs2JpsiPhiMassSignal.h"
#include "pdfs/include/Bs2JpsiPhiPromptBkg_tripleGaussian.h"
#include "pdfs/include/Bs2JpsiPhiPromptBkg_withTimeRes.h"
#include "pdfs/include/Bs2JpsiPhiPromptBkg_withTimeResDouble.h"
#include "pdfs/include/Bs2PhiPhi.h"
#include "pdfs/include/BsMass.h"
#include "pdfs/include/CrystalBall.h"
#include "pdfs/include/DoubleExponential.h"
#include "pdfs/include/DPBackground.h"
#include "pdfs/include/DPHistoBackground.h"
#include "pdfs/include/DPTotalAmplitudePDF.h"
#include "pdfs/include/DPTotalAmplitudePDF_withAcc.h"
#include "pdfs/include/DPTotalAmplitudePDF_withAcc_withBkg.h"
#include "pdfs/include/Exponential.h"
#include "pdfs/include/ExponentialWithDeltaGamma.h"
#include "pdfs/include/FlatPDF.h"
#include "pdfs/include/GammaDistribution.h"
#include "pdfs/include/Landau.h"
#include "pdfs/include/LandauGauss.h"
#include "pdfs/include/LogNormalDistribution.h"
#include "pdfs/include/LongLivedBkg.h"
#include "pdfs/include/LongLivedBkg_3Dangular.h"
#include "pdfs/include/MistagDistribution.h"
#include "pdfs/include/Novosibirsk.h"
#include "pdfs/include/OptimisedDoubleGauss.h"
#include "pdfs/include/OptimisedGauss.h"
#include "pdfs/include/PerEventErrorHistogram.h"
#include "pdfs/include/PolyPDF.h"
#include "pdfs/include/RapidFit_Pdf_Exponential.h"
#include "pdfs/include/SimpleDoubleGauss.h"
#include "pdfs/include/SimpleGauss.h"
#include "pdfs/include/SimpleGauss2D.h"
#include "pdfs/include/SimpleGauss3D.h"
#include "pdfs/include/SingleGauss.h"
#include "pdfs/include/Sinusoid.h"
#include "pdfs/include/StudentT.h"
#include "pdfs/include/TemplatePDF.h"
#include "pdfs/include/WrongPVAssocBkg.h"
#include "pdfs/dalitz/include/CalculateAngles.hh"
#include "pdfs/dalitz/include/DPBarrierFactor.hh"
#include "pdfs/dalitz/include/DPBarrierL0.hh"
#include "pdfs/dalitz/include/DPBarrierL1.hh"
#include "pdfs/dalitz/include/DPBarrierL2.hh"
#include "pdfs/dalitz/include/DPBarrierL3.hh"
#include "pdfs/dalitz/include/DPBarrierL4.hh"
#include "pdfs/dalitz/include/DPBarrierL5.hh"
#include "pdfs/dalitz/include/DPBWResonanceShape.hh"
#include "pdfs/dalitz/include/DPComponent.hh"
#include "pdfs/dalitz/include/DPGLassShape.hh"
#include "pdfs/dalitz/include/DPHelpers.hh"
#include "pdfs/dalitz/include/DPJpsiKaon.hh"
#include "pdfs/dalitz/include/DPLassShape.hh"
#include "pdfs/dalitz/include/DPMassShape.hh"
#include "pdfs/dalitz/include/DPNonresonant.hh"
#include "pdfs/dalitz/include/DPTotalAmplitude.hh"
#include "pdfs/dalitz/include/DPWignerFunction.hh"
#include "pdfs/dalitz/include/DPWignerFunctionGeneral.hh"
#include "pdfs/dalitz/include/DPWignerFunctionJ0.hh"
#include "pdfs/dalitz/include/DPWignerFunctionJ1.hh"
#include "pdfs/dalitz/include/DPWignerFunctionJ1over2.hh"
#include "pdfs/dalitz/include/DPWignerFunctionJ2.hh"
#include "pdfs/dalitz/include/DPWignerFunctionJ3over2.hh"
#include "pdfs/dalitz/include/DPZplusK.hh"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__rapidfit_dictLN_TClass;
extern G__linked_taginfo G__rapidfit_dictLN_TBuffer;
extern G__linked_taginfo G__rapidfit_dictLN_TMemberInspector;
extern G__linked_taginfo G__rapidfit_dictLN_TObject;
extern G__linked_taginfo G__rapidfit_dictLN_TList;
extern G__linked_taginfo G__rapidfit_dictLN_RapidRun;
extern G__linked_taginfo G__rapidfit_dictLN_auto_ptrlETListgR;

/* STUB derived class for protected member access */
