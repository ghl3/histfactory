

#include "RooStats/HistFactory/HistFactoryNavigation.h"
#include "RooStats/HistFactory/HistFactoryException.h"

namespace RooStats {
  namespace HistFactory {

    HistFactoryNavigation::HistFactoryNavigation(ModelConfig* mc) {

      // Save the model pointer
      fModel = mc->GetPdf();
      fObservables = (RooArgSet*) mc->GetObservables();
      
      // Initialize the rest of the members
      _GetNodes(fModel, fObservables);

    }


    RooAbsReal* HistFactoryNavigation::SampleFunction(const std::string& channel, const std::string& sample){

      std::map< std::string, std::map< std::string, RooAbsReal*> >::iterator channel_itr;
      channel_itr = fChannelSampleFunctionMap.find(channel);
      if( channel_itr==fChannelSampleFunctionMap.end() ){
	std::cout << "Error: Channel: " << channel << " not found in Navigation" << std::endl;
	throw hf_exc();
      }

      std::map< std::string, RooAbsReal*>& SampleMap = channel_itr->second;
      std::map< std::string, RooAbsReal*>::iterator sample_itr;
      sample_itr = SampleMap.find(sample);
      if( sample_itr==SampleMap.end() ){
	std::cout << "Error: Sample: " << sample << " not found in Navigation" << std::endl;
	throw hf_exc();
      }
      
      return sample_itr->second;

    }

    RooArgSet* HistFactoryNavigation::GetObservableSet(const std::string& channel) {

      std::map< std::string, RooArgSet*>::iterator channel_itr;
      channel_itr = fChannelObservMap.find(channel);
      if( channel_itr==fChannelObservMap.end() ){
	std::cout << "Error: Channel: " << channel << " not found in Navigation" << std::endl;
	throw hf_exc();
      }
      
      return channel_itr->second;

    }


    TH1* HistFactoryNavigation::GetSampleHist(const std::string& channel, const std::string& sample) {

      // Convert the ArgSet to an ArgList
      // This is temporary!
      RooArgList observable_list( *GetObservableSet(channel) );//observable_list( *(fChannelObservMap[channel]) );
      RooRealVar* observable = (RooRealVar*) observable_list.at(0);
      
      std::string HistName = channel + "_" + sample + "_hist";

      RooAbsReal* sample_function = SampleFunction(channel, sample); //fChannelSampleFunctionMap[channel][sample];

      return MakeHistFromRooFunction( sample_function, observable, HistName );
				     
    }


    void HistFactoryNavigation::_GetNodes(RooAbsPdf* modelPdf, const RooArgSet* observables) {
      
      // Get the pdf from the ModelConfig
      //RooAbsPdf* modelPdf = mc->GetPdf();
      //RooArgSet* observables = mc->GetObservables();

      // Create vectors to hold the channel pdf's
      // as well as the set of observables for each channel
      //std::map< std::string, RooAbsPdf* >  channelPdfMap;
      //std::map< std::string, RooArgSet* >  channelObservMap;
      
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
	  fChannelNameVec.push_back( ChannelName );
	  RooAbsPdf* pdftmp = simPdf->getPdf(ChannelName.c_str()) ;
	  RooArgSet* obstmp = pdftmp->getObservables(*observables) ;
	  fChannelPdfMap[ChannelName] = pdftmp;
	  fChannelObservMap[ChannelName] =  obstmp;
	}

      } else { 
	RooArgSet* obstmp = modelPdf->getObservables(*observables) ;	
	std::string ChannelName = modelPdf->GetName();
	fChannelNameVec.push_back(ChannelName);
	fChannelPdfMap[ChannelName] = modelPdf;
	fChannelObservMap[ChannelName] = obstmp;

      }

      // Okay, now we have maps of the pdfs 
      // and the observable list per channel
      // We then loop over the channel pdfs:
      // and find their RooRealSumPdfs
      // std::map< std::string, RooRealSumPdf* > channelSumNodeMap;

      for( unsigned int i = 0; i < fChannelNameVec.size(); ++i ) {

	std::string ChannelName = fChannelNameVec.at(i);
	RooAbsPdf* pdf = fChannelPdfMap[ChannelName];
	//std::string Name = fChannelNameMap[ChannelName];

	// Loop over the pdf's components and find
	// the (one) that is a RooRealSumPdf
	// Based on the mode, we assume that node is 
	// the "unconstrained" pdf node for that channel
	RooArgSet* components = pdf->getComponents();
	TIterator* argItr = components->createIterator();
	RooAbsArg* arg = NULL;
	while( (arg=(RooAbsArg*)argItr->Next()) ) {
	  std::string ClassName = arg->ClassName();
	  if( ClassName == "RooRealSumPdf" ) {
	    fChannelSumNodeMap[ChannelName] = (RooRealSumPdf*) arg;
	    break;
	  }
	}
      }
      
      // Okay, now we have all necessary
      // nodes filled for each channel.
      for( unsigned int i = 0; i < fChannelNameVec.size(); ++i ) {

	std::string ChannelName = fChannelNameVec.at(i);
	RooRealSumPdf* sumPdf = (RooRealSumPdf*) fChannelSumNodeMap[ChannelName];
	
	// We now take the RooRealSumPdf and loop over
	// its component functions.  The RooRealSumPdf turns
	// a list of functions (expected events or bin heights
	// per sample) and turns it into a pdf.
	// Therefore, we loop over it to find the expected
	// height for the various samples

	// First, create a map to store the function nodes
	// for each sample in this channel
	std::map< std::string, RooAbsReal*> sampleFunctionMap;

	// Loop over the sample nodes in this
	// channel's RooRealSumPdf
	RooArgList nodes = sumPdf->funcList();
	TIterator* sampleItr = nodes.createIterator();
	RooAbsArg* sample;
	while( (sample=(RooAbsArg*)sampleItr->Next()) ) {

	  // Cast this node as a function
	  RooAbsReal* func = (RooAbsReal*) sample;

	  // Do a bit of work to get the name of each sample
	  std::string SampleName = sample->GetName();
	  if( SampleName.find("L_x_") != std::string::npos ) {
	    size_t index = SampleName.find("L_x_");
	    SampleName.replace( index, 4, "" );
	  }
	  if( SampleName.find(ChannelName.c_str()) != std::string::npos ) {
	    size_t index = SampleName.find(ChannelName.c_str());
	    SampleName = SampleName.substr(0, index-1);
	  }

	  // And simply save this node into our map
	  sampleFunctionMap[SampleName] = func;

	}

	fChannelSampleFunctionMap[ChannelName] = sampleFunctionMap;

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

    // A simple wrapper to use a ModelConfig
    void HistFactoryNavigation::_GetNodes(ModelConfig* mc) {
      RooAbsPdf* modelPdf = mc->GetPdf();
      const RooArgSet* observables = mc->GetObservables();
      _GetNodes(modelPdf, observables);
    }

  }
}




