



#include "RooStats/HistFactory/HistFactoryNavigation.h"



namespace RooStats {
  namespace HistFactory {

    HistFactoryNavigation::HistFactoryNavigation(RooAbsPdf* model) {
      fModel = model;
    }



    void HistFactoryNavigation::_GetNodes(std::vector< std::vector< TH1* > >& channelHistListVec, std::vector< TH1* >& channelDataVec,
					  ModelConfig* mc, RooAbsData* data) {
      
      // Get the pdf from the ModelConfig
      RooAbsPdf* modelPdf = mc->GetPdf();

      // Create vectors to hold the channel pdf's
      // as well as the set of observables for each channel
      std::map< RooAbsPdf* >  channelPdfVec;
      std::map< RooArgSet* >  channelObservVec;
      
      // Check if it is a simultaneous pdf or not
      // (if it's an individual channel, it won't be, if it's
      // combined, it's simultaneous)
      // Fill the channel vectors based on the structure
      // (Obviously, if it's not simultaneous, there will be
      // only one entry in the vector for the single channel)
      if(strcmp(modelPdf->ClassName(),"RooSimultaneous")==0){

	// If so, get a list of the component pdf's:
	RooSimultaneous* simPdf = (RooSimultaneous*) modelPdf;
	RooCategory* channelCat = (RooCategory*) (&simPdf->indexCat());
	
	// Iterate over the categories and get the
	// pdf and observables for each category
	TIterator* iter = channelCat->typeIterator() ;
	RooCatType* tt = NULL;
        while((tt=(RooCatType*) iter->Next())) {
	  std::string ChannelName = tt->GetName();
	  RooAbsPdf* pdftmp = simPdf->getPdf(ChannelName.c_str()) ;
	  RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;
	  channelPdfVec.push_back( pdftmp );
	  channelObservVec.push_back( obstmp );
	  fChannelNameVec.push_back( ChannelName );
	}

      } else { 
	RooArgSet* obstmp = modelPdf->getObservables(*mc->GetObservables()) ;	
	channelPdfVec.push_back( modelPdf );
	channelObservVec.push_back( obstmp );
	channelNameVec.push_back( modelPdf->GetName() );
      } 


      // Okay, now we have a list of the pdf's 
      // and the observable list per channel
      // Now loop over the channel pdf's:
      // and find their RooRealSumPdf's
      std::map< std::string, RooRealSumPdf* > channelSumNodeVec;

      for( unsigned int i = 0; i < channelPdfVec.size(); ++i ) {

	RooAbsPdf* pdf   = channelPdfVec.at(i);
	std::string Name = channelNameVec.at(i);

	// Loop over the pdf's components
	RooArgSet* components = pdf->getComponents();
	TIterator* argItr = components->createIterator();
	RooAbsArg* arg = NULL;
	while( (arg=(RooAbsArg*)argItr->Next()) ) {
	  std::string ClassName = arg->ClassName();
	  if( ClassName == "RooRealSumPdf" ) {
	    channelSumNodeVec.push_back( (RooRealSumPdf*) arg );
	    std::cout << "Found RooRealSumPdf: " << arg << " " << arg->GetName() << std::endl;
	    break;
	  }
	}

      }
      
      // Okay, now we have all necessary
      // nodes filled for each channel.

      for( unsigned int i = 0; i < channelPdfVec.size(); ++i ) {

	std::string    Channel  = channelNameVec.at(i);
	RooRealSumPdf* sumPdf   = (RooRealSumPdf*) channelSumNodeVec.at(i);
	RooArgList* observables = (RooArgList*) channelObservVec.at(i);
      
	std::cout << "Looping over sum pdf: " << Channel << " " << sumPdf << std::endl;
	
	// Print the nodes:
	RooArgList nodes = sumPdf->funcList();
	nodes.Print();

	// Create a vector to store the sample
	// histograms for this channel
	std::vector< TH1* > histList;

	// Loop over the sample nodes in this
	// channel's RooRealSumPdf
	TIterator* sampleItr = nodes.createIterator();
	RooAbsArg* sample;
	while( (sample=(RooAbsArg*)sampleItr->Next()) ) {

	  if( observables->getSize() == 0 ) {
	    std::cout << "Error: No observables found in channel: " << Channel << std::endl;
	  }
	  if( observables->getSize() > 1 ) {
	    std::cout << "Error: This function for now only works with 1-d histograms" << std::endl;
	  }


	  RooRealVar* xVar = (RooRealVar*) observables->at(0);

	  // First, get the Data histogram
	  // for this channel
	  TH1* dataHist = new TH1F( (Channel+"_data").c_str(), "", 
				    xVar->numBins(), xVar->getMin(), xVar->getMax() );
	  //dataHist->Sumw2();
	  data->fillHistogram( dataHist, RooArgSet( *xVar ) );
	  // Fix the histogram's bin errors:
	  for( int j_bin = 0; j_bin < dataHist->GetNbinsX(); ++j_bin ) {
	    dataHist->SetBinError( j_bin + 1, TMath::Sqrt( dataHist->GetBinContent( j_bin + 1 ) ) );
	  }
	  channelDataVec.push_back( dataHist );

	  // Get a list of the RooAbsReal nodes:
	  RooAbsReal* func = (RooAbsReal*) sample;

	  // Do a bit of work to get the name of each sample
	  std::string SampleName = sample->GetName();
	  if( SampleName.find("L_x_") != std::string::npos ) {
	    size_t index = SampleName.find("L_x_");
	    SampleName.replace( index, 4, "" );
	  }
	  if( SampleName.find(Channel.c_str()) != std::string::npos ) {
	    size_t index = SampleName.find(Channel.c_str());
	    SampleName = SampleName.substr(0, index-1);
	  }

	  TH1* sampleHist = MakeHistFromRooFunction( func, xVar, SampleName /*sample->GetName()*/ ); 
	  sampleHist->SetTitle( SampleName.c_str() );
	  //	  sampleHist->SetName( SampleName.c_str() );
	  histList.push_back( sampleHist );
	  std::cout << "Adding histogram for sample: " << SampleName /*sample->GetName()*/ << std::endl;

	}

	channelHistListVec.push_back( histList );

	// Okay, now we have a list of histograms
	// representing the samples for this channel.

      }


    }



    TH1* HistFactoryNavigation::MakeHistFromRooFunction( RooAbsReal* func, RooRealVar* var, std::string name ) {

      // Turn a RooAbsReal* into a TH1* based 
      // on a template histogram.  
      // Loop over the bins of the Template,
      // find the bin centers, 
      // Scan the input Var over those bin centers,
      // and use the value of the function
      // to make the new histogram


      // Make the new histogram
      // Cone and empty the template
      //      TH1* hist = (TH1*) histTemplate.Clone( name.c_str() );

      TH1* hist = (TH1*) var->createHistogram( name.c_str(), RooFit::Binning(var->getBinning()) );
      hist->Reset();

      std::cout << "Created histogram: " << hist->GetName() 
		<< " with binning: " 
		<< "n=" << hist->GetNbinsX()
		<< " [" << hist->GetXaxis()->GetXmin() 
		<< ", " << hist->GetXaxis()->GetXmax() 
		<< "] " << std::endl;

      int nBins = hist->GetNbinsX();

      //std::vector< double > binCenters;

      for( int i = 0; i < nBins; ++i ) {

	double binCenter = hist->GetBinCenter(i+1);
	var->setVal( binCenter );

	double val = func->getVal();

	// double binWidth = var->getBinWidth( i );
	// hist->Fill( binCenter, val*binWidth );

	hist->Fill( binCenter, val );
	
	std::cout << "Filling: " << func->GetName() << " "
		  << " (" << binCenter << "," << val << ")"
		  << std::endl;

      }

      return hist;
      

    }



  }
}




