

// Blah
y
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

  // The HistFactory Pdf Pointer
  RooAbsPdf* fModel;


};
