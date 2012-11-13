// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, Akira Shibata
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//_________________________________________________
/*
  BEGIN_HTML
  <p>
  This is a package that creates a RooFit probability density function from ROOT histograms 
  of expected distributions and histograms that represent the +/- 1 sigma variations 
  from systematic effects. The resulting probability density function can then be used
  with any of the statistical tools provided within RooStats, such as the profile 
  likelihood ratio, Feldman-Cousins, etc.  In this version, the model is directly
  fed to a likelihodo ratio test, but it needs to be further factorized.</p>

  <p>
  The user needs to provide histograms (in picobarns per bin) and configure the job
  with XML.  The configuration XML is defined in the file config/Config.dtd, but essentially
  it is organized as follows (see config/Combination.xml and config/ee.xml for examples)</p>

  <ul>
  <li> - a top level 'Combination' that is composed of:</li>
  <ul>
  <li>- several 'Channels' (eg. ee, emu, mumu), which are composed of:</li>
  <ul>
  <li>- several 'Samples' (eg. signal, bkg1, bkg2, ...), each of which has:</li>
  <ul>
  <li> - a name</li>
  <li> - if the sample is normalized by theory (eg N = L*sigma) or not (eg. data driven)</li>
  <li> - a nominal expectation histogram</li>
  <li> - a named 'Normalization Factor' (which can be fixed or allowed to float in a fit)</li>
  <li> - several 'Overall Systematics' in normalization with:</li>
  <ul>
  <li> - a name</li>
  <li> - +/- 1 sigma variations (eg. 1.05 and 0.95 for a 5% uncertainty)</li>
  </ul>
  <li>- several 'Histogram Systematics' in shape with:</li>
  <ul>
  <li>- a name (which can be shared with the OverallSyst if correlated)</li>
  <li>- +/- 1 sigma variational histograms</li>
  </ul>
  </ul>
  </ul>
  <li>- several 'Measurements' (corresponding to a full fit of the model) each of which specifies</li>
  <ul>
  <li>- a name for this fit to be used in tables and files</li>
  <ul>
  <li>      - what is the luminosity associated to the measurement in picobarns</li>
  <li>      - which bins of the histogram should be used</li>
  <li>      - what is the relative uncertainty on the luminosity </li>
  <li>      - what is (are) the parameter(s) of interest that will be measured</li>
  <li>      - which parameters should be fixed/floating (eg. nuisance parameters)</li>
  </ul>
  </ul>
  </ul>
  END_HTML
*/
//


// from std
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

// from root
#include "TFile.h"
#include "TH1F.h"
#include "TDOMParser.h"
#include "TXMLAttr.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"

// from roofit
#include "RooStats/ModelConfig.h"

// from this package
#include "Helper.h"
//#include "RooStats/HistFactory/ConfigParser.h"
#include "RooStats/HistFactory/EstimateSummary.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"
#include "RooStats/HistFactory/HistFactoryException.h"

#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

using namespace RooFit;
//using namespace RooStats;
//using namespace HistFactory;

//using namespace std;



RooWorkspace* RooStats::HistFactory::MakeModelAndMeasurementFast( RooStats::HistFactory::Measurement& measurement ) {

  // This will be returned
  RooWorkspace* ws = NULL;

  try {

    std::cout << "Making Model and Measurements (Fast) for measurement: " << measurement.GetName() << std::endl;

    double lumiError = measurement.GetLumi()*measurement.GetLumiRelErr();

    std::cout << "using lumi = " << measurement.GetLumi() << " and lumiError = " << lumiError
	 << " including bins between " << measurement.GetBinLow() << " and " << measurement.GetBinHigh() << std::endl;
    std::cout << "fixing the following parameters:"  << std::endl;

    for(std::vector<std::string>::iterator itr=measurement.GetConstantParams().begin(); itr!=measurement.GetConstantParams().end(); ++itr){
      std::cout << "   " << *itr << std::endl;
    }

    std::string rowTitle = measurement.GetName();
    
    std::vector<RooWorkspace*> channel_workspaces;
    std::vector<std::string>        channel_names;

    // This holds the TGraphs that are created during the fit
    std::string outputFileName = measurement.GetOutputFilePrefix() + "_" + measurement.GetName() + ".root";
    std::cout << "Creating the output file: " << outputFileName << std::endl;
    TFile* outFile = new TFile(outputFileName.c_str(), "recreate");

    // This holds the table of fitted values and errors
    std::string tableFileName = measurement.GetOutputFilePrefix() + "_results.table";
    std::cout << "Creating the table file: " << tableFileName << std::endl;
    FILE*  tableFile =  fopen( tableFileName.c_str(), "a"); 

    std::cout << "Creating the HistoToWorkspaceFactoryFast factory" << std::endl;

    HistoToWorkspaceFactoryFast factory( measurement );

    std::cout << "Setting preprocess functions" << std::endl;

    // Make the factory, and do some preprocessing
    // HistoToWorkspaceFactoryFast factory(measurement, rowTitle, outFile);
    factory.SetFunctionsToPreprocess( measurement.GetPreprocessFunctions() );
  
    // for results tables
    fprintf(tableFile, " %s &", rowTitle.c_str() );
  
    // First: Loop to make the individual channels

    for( unsigned int chanItr = 0; chanItr < measurement.GetChannels().size(); ++chanItr ) {
    
      HistFactory::Channel& channel = measurement.GetChannels().at( chanItr );

      if( ! channel.CheckHistograms() ) {
	std::cout << "MakeModelAndMeasurementsFast: Channel: " << channel.GetName()
		  << " has uninitialized histogram pointers" << std::endl;
	throw hf_exc();
	return NULL;
      }

      std::string ch_name = channel.GetName();
      channel_names.push_back(ch_name);

      std::cout << "Starting to process channel: " << ch_name << std::endl;

      RooWorkspace* ws_single = factory.MakeSingleChannelModel( measurement, channel );

      channel_workspaces.push_back(ws_single);

      // Get the Paramater of Interest as a RooRealVar
      RooRealVar* poi = (RooRealVar*) ws_single->var( (measurement.GetPOI()).c_str() );
      
      // Make the output
      std::string ChannelFileName = measurement.GetOutputFilePrefix() + "_" 
	+ ch_name + "_" + rowTitle + "_model.root";
      ws_single->writeToFile( ChannelFileName.c_str() );
    
      // Now, write the measurement to the file
      // Make a new measurement for only this channel
      RooStats::HistFactory::Measurement meas_chan( measurement );
      meas_chan.GetChannels().clear();
      meas_chan.GetChannels().push_back( channel );
      std::cout << "Opening File to hold channel: " << ChannelFileName << std::endl;
      TFile* chanFile = TFile::Open( ChannelFileName.c_str(), "UPDATE" );
      std::cout << "About to write channel measurement to file" << std::endl;
      meas_chan.writeToFile( chanFile );
      std::cout << "Successfully wrote channel to file" << std::endl;
      chanFile->Close();

      // do fit unless exportOnly requested
      if(! measurement.GetExportOnly()){
	if(!poi){
	  std::cout << "Can't do fit for: " << measurement.GetName() 
		    << ", no parameter of interest" << std::endl;
	} else {
	  if(ws_single->data("obsData")){
	    FitModelAndPlot(measurement.GetName(), measurement.GetOutputFilePrefix(), ws_single, 
			    ch_name, "obsData",    outFile, tableFile);
	  } else {
	    FitModelAndPlot(measurement.GetName(), measurement.GetOutputFilePrefix(), ws_single, 
			    ch_name, "asimovData", outFile, tableFile);
	  }
	}
      }

      fprintf(tableFile, " & " );
    } // End loop over channels
  
    /***
	Second: Make the combined model:
	If you want output histograms in root format, create and pass it to the combine routine.
	"combine" : will do the individual cross-section measurements plus combination	
    ***/
  
    // Use HistFactory to combine the individual channel workspaces
    ws = factory.MakeCombinedModel(channel_names, channel_workspaces);

    // Configure that workspace
    HistoToWorkspaceFactoryFast::ConfigureWorkspaceForMeasurement( "simPdf", ws, measurement );

    // Get the Parameter of interest as a RooRealVar
    RooRealVar* poi = (RooRealVar*) ws->var( (measurement.GetPOI()).c_str() );

    std::string CombinedFileName = measurement.GetOutputFilePrefix() + "_combined_"
      + rowTitle + "_model.root";
    std::cout << "Writing combined workspace to file: " << CombinedFileName << std::endl;
    ws->writeToFile( CombinedFileName.c_str() );
    std::cout << "Writing combined measurement to file: " << CombinedFileName << std::endl;
    TFile* combFile = TFile::Open( CombinedFileName.c_str(), "UPDATE" );
    if( combFile == NULL ) {
      std::cout << "Error: Failed to open file " << CombinedFileName << std::endl;
      throw hf_exc();
    }
    measurement.writeToFile( combFile );
    combFile->Close();
    
    // Fit the combined model
    if(! measurement.GetExportOnly()){
      if(!poi){
	std::cout << "Can't do fit for: " << measurement.GetName() 
		  << ", no parameter of interest" << std::endl;
      } 
      else {
	if(ws->data("obsData")){
	  FitModelAndPlot(measurement.GetName(), measurement.GetOutputFilePrefix(), ws,"combined", 
			  "obsData",    outFile, tableFile);
	} 
	else {
	  FitModelAndPlot(measurement.GetName(), measurement.GetOutputFilePrefix(), ws,"combined", 
			  "asimovData", outFile, tableFile);
	}
      }
    }
  
    fprintf(tableFile, " \\\\ \n");

    outFile->Close();
    delete outFile;

    fclose( tableFile );

  }
  catch(std::exception& e)
    {
      std::cout << e.what() << std::endl;
      throw hf_exc();
      return NULL;
    }

  return ws;


}


  ///////////////////////////////////////////////
void RooStats::HistFactory::FitModelAndPlot(const std::string& MeasurementName, 
					    const std::string& FileNamePrefix, 
					    RooWorkspace * combined, std::string channel, 
					    std::string data_name, 
					    TFile* outFile, FILE* tableFile  ) {
    
    ModelConfig* combined_config = (ModelConfig *) combined->obj("ModelConfig");
    if(!combined_config){
      std::cout << "Error: no ModelConfig found in Measurement: "
		<< MeasurementName <<  std::endl;
      throw hf_exc();
    }

    RooAbsData* simData = combined->data(data_name.c_str());
    if(!simData){
      std::cout << "Error: Failed to get dataset: " << data_name
		<< " in measurement: " << MeasurementName << std::endl;
      throw hf_exc();
    }
    
    const RooArgSet* POIs = combined_config->GetParametersOfInterest();
    if(!POIs) {
      std::cout << "Not Fitting Model for measurement: " << MeasurementName
		<< ", no poi found" << std::endl;
      // Should I throw an exception here?
      return;
    }

    RooAbsPdf* model = combined_config->GetPdf();
    if( model==NULL ) {
      std::cout << "Error: Failed to find pdf in ModelConfig: " << combined_config->GetName()
		<< std::endl;
      throw hf_exc();
    }

    // Save a Snapshot
    RooArgSet PoiPlusNuisance( *combined_config->GetNuisanceParameters() );
    PoiPlusNuisance.add( *combined_config->GetParametersOfInterest() );
    combined->saveSnapshot("InitialValues", PoiPlusNuisance);

    ///////////////////////////////////////
    //Do combined fit
    std::cout << "\n\n---------------" << std::endl;
    std::cout << "---------------- Doing "<< channel << " Fit" << std::endl;
    std::cout << "---------------\n\n" << std::endl;
    model->fitTo(*simData, Minos(kTRUE), PrintLevel(1));
    //    PrintCovarianceMatrix(result, allParams, "results/"+FilePrefixStr(channel)+"_corrMatrix.table" );

    // ONE COULD ADD HistFactoryNavigation OUTPUT HERE

    // Print the output table file 
    // (Would like to update this...)
    if( outFile != NULL ) {

      //
      // We only print info for the FIRST POI,
      // but the rest are still saved in the ModelConfig
      // 

      RooRealVar* poi = NULL; // (RooRealVar*) POIs->first();
      // for results tables
      TIterator* params_itr=POIs->createIterator();
      TObject* params_obj=0;
      while( (params_obj=params_itr->Next()) ) {
	poi = (RooRealVar*) params_obj;
	std::cout << "printing results for " << poi->GetName() 
		  << " at " << poi->getVal()<< " high " 
		  << poi->getErrorLo() << " low " 
		  << poi->getErrorHi() << std::endl;
      }

      fprintf(tableFile, " %.4f / %.4f  ", poi->getErrorLo(), poi->getErrorHi());

      RooAbsReal* nll = model->createNLL(*simData);
      RooAbsReal* profile = nll->createProfile(*poi);
      if( profile==NULL ) {
	std::cout << "Error: Failed to make ProfileLikelihood for: " << poi->GetName() 
		  << " using model: " << model->GetName()
		  << " and data: " << simData->GetName()
		  << std::endl;
	throw hf_exc();
      }

      RooPlot* frame = poi->frame();
      if( frame == NULL ) {
	std::cout << "Error: Failed to create RooPlot frame for: " << poi->GetName() << std::endl;
	throw hf_exc();
      }

      // Draw the likelihood curve
      FormatFrameForLikelihood(frame);
      TCanvas* c1 = new TCanvas( channel.c_str(), "",800,600);
      nll->plotOn(frame, ShiftToZero(), LineColor(kRed), LineStyle(kDashed));
      profile->plotOn(frame);
      frame->SetMinimum(0);
      frame->SetMaximum(2.);
      frame->Draw();
      c1->SaveAs( (FileNamePrefix+"_"+channel+"_"+MeasurementName+"_profileLR.eps").c_str() );
      delete c1;

      // Save to the output file
      TDirectory* channel_dir = outFile->mkdir(channel.c_str());
      if( channel_dir == NULL ) {
	std::cout << "Error: Failed to make channel directory: " << channel << std::endl;
	throw hf_exc();
      }
      TDirectory* summary_dir = channel_dir->mkdir("Summary");
      if( summary_dir == NULL ) {
	std::cout << "Error: Failed to make Summary directory for channel: " 
		  << channel << std::endl;
	throw hf_exc();
      }
      summary_dir->cd();
    
      RooCurve* curve=frame->getCurve();
      Int_t curve_N=curve->GetN();
      Double_t* curve_x=curve->GetX();
      delete frame;

      Double_t * x_arr = new Double_t[curve_N];
      Double_t * y_arr_nll = new Double_t[curve_N];
      
      for(int i=0; i<curve_N; i++){
	double f=curve_x[i];
	poi->setVal(f);
	x_arr[i]=f;
	y_arr_nll[i]=nll->getVal();
      }

      TGraph* g = new TGraph(curve_N, x_arr, y_arr_nll);
      g->SetName( (FileNamePrefix +"_nll").c_str() );
      g->Write(); 
      delete g;
      delete [] x_arr;
      delete [] y_arr_nll;

    } // End: if( outFile==NULL )

    // Finally, restore the initial values
    combined->loadSnapshot("InitialValues");

  }


void RooStats::HistFactory::FitModel(RooWorkspace * combined, std::string data_name ) {

    std::cout << "In Fit Model" << std::endl;
    ModelConfig * combined_config = (ModelConfig *) combined->obj("ModelConfig");
    if(!combined_config){
      std::cout << "no model config " << "ModelConfig" << " exiting" << std::endl;
      return;
    }
    
    RooAbsData* simData = combined->data(data_name.c_str());
    if(!simData){
      std::cout << "no data " << data_name << " exiting" << std::endl;
      return;
    }

    const RooArgSet * POIs=combined_config->GetParametersOfInterest();
    if(!POIs){
      std::cout << "no poi " << data_name << " exiting" << std::endl;
      return;
    }

    RooAbsPdf* model=combined_config->GetPdf();
    model->fitTo(*simData, Minos(kTRUE), PrintLevel(1));
    
  }


void RooStats::HistFactory::FormatFrameForLikelihood(RooPlot* frame, std::string /*XTitle*/, 
						     std::string YTitle){

    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(255);
    gStyle->SetTitleFillColor(255);
    gStyle->SetFrameFillColor(0);  
    gStyle->SetStatColor(255);
    
    RooAbsRealLValue* var = frame->getPlotVar();
    double xmin = var->getMin();
    double xmax = var->getMax();
    
    frame->SetTitle("");
    //      frame->GetXaxis()->SetTitle(XTitle.c_str());
    frame->GetXaxis()->SetTitle(var->GetTitle());
    frame->GetYaxis()->SetTitle(YTitle.c_str());
    frame->SetMaximum(2.);
    frame->SetMinimum(0.);
    TLine * line = new TLine(xmin,.5,xmax,.5);
    line->SetLineColor(kGreen);
    TLine * line90 = new TLine(xmin,2.71/2.,xmax,2.71/2.);
    line90->SetLineColor(kGreen);
    TLine * line95 = new TLine(xmin,3.84/2.,xmax,3.84/2.);
    line95->SetLineColor(kGreen);
    frame->addObject(line);
    frame->addObject(line90);
    frame->addObject(line95);
}


