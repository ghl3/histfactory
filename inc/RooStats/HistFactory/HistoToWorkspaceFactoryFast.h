// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, Akira Shibata
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_HISTOTOWORKSPACEFACTORYFAST
#define ROOSTATS_HISTOTOWORKSPACEFACTORYFAST

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>

#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <TObject.h>
#include <TH1.h>
#include <TDirectory.h>

//#include "RooStats/HistFactory/EstimateSummary.h"
//#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/Systematics.h"
class ParamHistFunc;


namespace RooStats{
  namespace HistFactory{

    // Forward Declarations FTW
    class Measurement;
    class Channel;
    class Sample;


    class HistoToWorkspaceFactoryFast: public TObject {
    
    public:
      typedef std::map<std::string, double> param_map;
      HistoToWorkspaceFactoryFast(  std::string, std::string , std::vector<std::string> , double =200, double =20, int =0, int =6, TFile* =NULL, param_map = param_map() );
      HistoToWorkspaceFactoryFast(  RooStats::HistFactory::Measurement& Meas );
      static void ConfigureWorkspaceForMeasurement( const std::string& ModelName, RooWorkspace* ws_single, Measurement& measurement );
    
      HistoToWorkspaceFactoryFast();
      virtual ~HistoToWorkspaceFactoryFast();
    
      RooWorkspace* MakeSingleChannelModel( Measurement& measurement, Channel& channel );
      static RooWorkspace* MakeCombinedModel( Measurement& measurement );
    
      void SetFunctionsToPreprocess(std::vector<std::string> lines){ fPreprocessFunctions = lines; }
    
      //ES//
      /*
	void AddEfficiencyTerms(RooWorkspace* proto, std::string prefix, std::string interpName,
	std::map<std::string,std::pair<double,double> > systMap,
	std::vector<std::string>& likelihoodTermNames, std::vector<std::string>& totSystTermNames);
      */			    
      void AddEfficiencyTerms(RooWorkspace* proto, std::string prefix, std::string interpName,
			      std::vector<OverallSys>& systList, 			 
			      std::vector<std::string>& likelihoodTermNames, std::vector<std::string>& totSystTermNames);

      //ES// std::string AddNormFactor(RooWorkspace *, std::string & , std::string & , EstimateSummary & , bool );
      std::string AddNormFactor(RooWorkspace* proto, std::string& channel, std::string& sigmaEpsilon, Sample& sample, bool doRatio);

      void AddMultiVarGaussConstraint(RooWorkspace* proto, std::string prefix,int lowBin, int highBin, std::vector<std::string>& likelihoodTermNames);
    
      void AddPoissonTerms(RooWorkspace* proto, std::string prefix, std::string obsPrefix, std::string expPrefix, int lowBin, int highBin,
			   std::vector<std::string>& likelihoodTermNames);
    
      //void Combine_old();
    
      RooWorkspace *  MakeCombinedModel(std::vector<std::string>, std::vector<RooWorkspace*>);
    
      //void Combine_ratio(std::vector<std::string> , std::vector<RooWorkspace*>);
    
      void Customize(RooWorkspace* proto, const char* pdfNameChar, std::map<std::string,std::string> renameMap);
    
      static void EditSyst(RooWorkspace* proto, const char* pdfNameChar, 
			   std::map<std::string,double> gammaSyst, std::map<std::string,double> uniformSyst, std::map<std::string,double> logNormSyst, std::map<std::string,double> noSyst);
    
      //void FormatFrameForLikelihood(RooPlot* frame, std::string XTitle=std::string("#sigma / #sigma_{SM}"), std::string YTitle=std::string("-log likelihood"));
    
    
      //ES// void LinInterpWithConstraint(RooWorkspace* proto, TH1* nominal, std::vector<TH1*> lowHist, std::vector<TH1*> highHist,
      //ES// 				 std::vector<std::string> sourceName, std::string prefix, std::string productPrefix, std::string systTerm,
      //ES// 				 int lowBin, int highBin, std::vector<std::string>& likelihoodTermNames);
      void LinInterpWithConstraint(RooWorkspace* proto, TH1* nominal, std::vector<HistoSys>,  std::string prefix, std::string productPrefix, std::string systTerm,
				   int lowBin, int highBin, std::vector<std::string>& likelihoodTermNames);
    
      TDirectory* Makedirs( TDirectory* file, std::vector<std::string> names );
    
      //ES// RooWorkspace* MakeSingleChannelModel(std::vector<RooStats::HistFactory::EstimateSummary> summary, std::vector<std::string> systToFix, bool doRatio=false);
      // GHL: Renaming to "MakeSingleChannelWorkspace"
      RooWorkspace* MakeSingleChannelWorkspace(Measurement& measurement, Channel& channel);
    
      void  MakeTotalExpected(RooWorkspace* proto, std::string totName, std::string /**/, std::string /**/,
			      int lowBin, int highBin, std::vector<std::string>& syst_x_expectedPrefixNames,
			      std::vector<std::string>& normByNames);
    
      TDirectory* Mkdir( TDirectory * file, std::string name );
    
      RooDataSet* MergeDataSets(RooWorkspace* combined,
				std::vector<RooWorkspace*> wspace_vec, 
				std::vector<std::string> channel_names, 
				std::string dataSetName,
				RooArgList obsList,
				RooCategory* channelCat);

      static void PrintCovarianceMatrix(RooFitResult* result, RooArgSet* params, std::string filename);
      void ProcessExpectedHisto(TH1* hist,RooWorkspace* proto, std::string prefix, std::string productPrefix, std::string systTerm, double low, double high, int lowBin, int highBin);
      void SetObsToExpected(RooWorkspace* proto, std::string obsPrefix, std::string expPrefix, int lowBin, int highBin);
      //void FitModel(const std::string& FileNamePrefix, RooWorkspace *, std::string, std::string, TFile*, FILE*);
      //std::string FilePrefixStr(std::string);
    
      TH1* MakeScaledUncertaintyHist( const std::string& Name, std::vector< std::pair<TH1*,TH1*> > HistVec );
      TH1* MakeAbsolUncertaintyHist( const std::string& Name, const TH1* Hist );
      RooArgList createStatConstraintTerms( RooWorkspace* proto, std::vector<std::string>& constraintTerms, ParamHistFunc& paramHist, TH1* uncertHist, 
					    Constraint::Type type, Double_t minSigma );
    
      inline void SetObsNameVec(const std::vector<std::string>& obsNameVec) { fObsNameVec = obsNameVec; }
      inline void SetObsName(const std::string& obsName) { fObsNameVec.clear(); fObsNameVec.push_back(obsName); fObsName = obsName; }
      inline void AddObsName(const std::string& obsName) { fObsNameVec.push_back(obsName); }

      void ConfigureHistFactoryDataset(RooDataSet* obsData, TH1* nominal, RooWorkspace* proto,
				       std::vector<std::string> obsNameVec);
    
    
      //string fFileNamePrefix;
      //string fRowTitle;
      std::vector<std::string> fSystToFix;
      std::map<std::string, double> fParamValues;
      double fNomLumi, fLumiError;
      int fLowBin, fHighBin;    
      //std::stringstream fResultsPrefixStr;
      //TFile * fOut_f;
      // FILE * pFile;
    
    private:
    
      void GuessObsNameVec(TH1* hist);
    
      std::vector<std::string> fObsNameVec;
      std::string fObsName;
      std::vector<std::string> fPreprocessFunctions;
    
      ClassDef(RooStats::HistFactory::HistoToWorkspaceFactoryFast,3)
    };
  
  }
}

#endif
