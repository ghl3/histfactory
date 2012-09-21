
#include <map>

#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TMath.h"

#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooRealSumPdf.h"
#include "RooStats/HistFactory/Measurement.h"

#include "RooStats/ModelConfig.h"


namespace RooStats {
  namespace HistFactory { 

    class HistFactoryNavigation {
  
    public:

      // Initialze based on an already-created HistFactory Model
      HistFactoryNavigation(ModelConfig* mc);

      // Return the total hist (fitted) in the channel
      // This is harder, implement later
      //TH1* GetChannelHist(const std::string& channel); 

      // Get the RooAbsReal function for a given sample in a given channel
      RooAbsReal* SampleFunction(const std::string& channel, const std::string& sample);

      // Get the set of observables for a given channel
      RooArgSet* GetObservableSet(const std::string& channel);

      // The (current) histogram for that sample 
      // This includes all parameters and interpolation
      TH1* GetSampleHist(const std::string& channel, const std::string& sample, const std::string& name="");  

      // Get the total channel histogram for this channel
      // This is found by summing over the sample histograms
      TH1* GetChannelHist(const std::string& channel, const std::string& name="");

      // The value of the ith bin for the total in that channel
      // double GetBinValue(int bin, const std::string& channel);  
      // The value of the ith bin for that sample and channel 
      // double GetBinValue(int bin, const std::string& sample, const std::string& channel);  

      // Get the nuisance parameters for this channell
      // TO BE IMPLEMENTED
      // RooArgSet* GetNuisanceParameters(const std::string& channel);

      // Should pretty print all channels and the current values 
      void PrintState();
      // Should pretty print this and the current values
      void PrintState(const std::string& channel);

      // Print a "HistFactory style" RooDataSet in a readable way
      static void PrintDataSet(RooDataSet* data);
  
    private:

      // Internal Methods
      
      // To be used in the constructor, and results should be cached
      // Given a ModelConfig* and a RooAbsData* (optionally),
      // return vectors of the channel names, TH1*'s for the channels, 
      // and TH1's for the data
      // (For this method, we don't want to return histograms, we want
      // pointers for the pdf objects)
      void _GetNodes(ModelConfig* mc);
      void _GetNodes(RooAbsPdf* model, const RooArgSet* observables);
      
      void PrettyPrintHistogram(TH1* hist);

      // Make a histogram from a funciton
      // Edit so it can take a RooArgSet of parameters
      TH1* MakeHistFromRooFunction( RooAbsReal* func, RooArgList vars, std::string name="Hist" );

      std::map< std::string, RooAbsReal*> GetSampleFunctionMap(const std::string& channel);
      
      // The HistFactory Pdf Pointer
      RooAbsPdf* fModel;

      // The observables
      RooArgSet* fObservables;

      // The list of channels
      std::vector<std::string> fChannelNameVec;

      // Map of channel names to their full pdf's
      std::map< std::string, RooAbsPdf* > fChannelPdfMap;  

      // Map of channel names to pdf without constraint
      std::map< std::string, RooAbsPdf* > fChannelSumNodeMap;  
      
      // Map of channel names to their set of ovservables
      std::map< std::string, RooArgSet*> fChannelObservMap;

      // Map of Map of Channel, Sample names to Function Nodes
      // Used by doing: fChannelSampleFunctionMap["MyChannel"]["MySample"]
      std::map< std::string, std::map< std::string, RooAbsReal*> > fChannelSampleFunctionMap;

    };

  }
}
