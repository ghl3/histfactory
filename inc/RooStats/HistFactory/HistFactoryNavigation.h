
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
      HistFactoryNavigation(RooAbsPdf* model);

      // Return the total hist (fitted) in the channel
      TH1* GetChannelHist(const std::string& channel); 
  
      //  The fitted histogram for that sample
      TH1* GetSampleHist(const std::string& sample, const std::string& channel);  

      // The value of the ith bin for the total in that channel
      double GetBinValue(int bin, const std::string& channel);  
      // The value of the ith bin for that sample and channel 
      double GetBinValue(int bin, const std::string& sample, const std::string& channel);  

      // Get the observables for this channel
      RooArgSet* GetChannelObservables(const std::string& channel);
      // Get the nuisance parameters for this channell
      RooArgSet* GetNuisanceParameters(const std::string& channel);


      // <-- Should pretty print all channels and the current values 
      void PrintState(); 
      // <-- Should pretty print this and the current values
      void PrintState(const std::string& channel); 
  
  
  
  
    private:

      // Internal Methods

      
      // To be used in the constructor, and results should be cached
      // Given a ModelConfig* and a RooAbsData* (optionally),
      // return vectors of the channel names, TH1*'s for the channels, 
      // and TH1's for the data
      // (For this method, we don't want to return histograms, we want
      // pointers for the pdf objects)
      void _GetNodes(ModelConfig* mc);
      void _GetNodes(RooAbsPdf* model, RooArgSet* observables);


      TH1* MakeHistFromRooFunction( RooAbsReal* func, RooRealVar* var, std::string name="Hist" );
      

      // The HistFactory Pdf Pointer
      RooAbsPdf* fModel;

      // The list of channels
      std::vector<std::string> fChannelNameVec;

      // Map of channel names to their full pdf's
      std::map< std::string, RooAbsPdf* > fChannelPdfMap;  

      // Map of channel names to pdf without constraint
      std::map< std::string, RooAbsPdf* > fChannelSumNodeMap;  
      
      // Map of channel names to their set of ovservables
      std::map< std::string, RooArgAet*> fChannelObservMap;
      

    };

  }
}
