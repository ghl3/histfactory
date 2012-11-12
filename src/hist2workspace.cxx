
#include <string>
#include <exception>
#include <vector>

//void topDriver(string input); // in MakeModelAndMeasurements
//void fastDriver(string input); // in MakeModelAndMeasurementsFast

//#include "RooStats/HistFactory/MakeModelAndMeasurements.h"
#include "RooStats/HistFactory/ConfigParser.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

//_____________________________batch only_____________________
#ifndef __CINT__

namespace RooStats {
  namespace HistFactory {
    void fastDriver(std::string input) {

      // Create the initial list of measurements and channels
      std::vector< HistFactory::Measurement > measurement_list;
      // std::vector< HistFactory::Channel >     channel_list;

      // Fill them using the XML parser
      HistFactory::ConfigParser xmlParser;
      measurement_list = xmlParser.GetMeasurementsFromXML( input );
	
      // At this point, we have all the information we need
      // from the xml files.

      // We will make the measurements 1-by-1
      for(unsigned int i = 0; i < measurement_list.size(); ++i) {
	HistFactory::Measurement measurement = measurement_list.at(i);
	measurement.CollectHistograms();
	MakeModelAndMeasurementFast( measurement );
      }

      return;

    }
  } // namespace RooStats
} // namespace HistFactory


int main(int argc, char** argv) {

  if( !(argc>1) ) {
    std::cerr << "need input file" << std::endl;
    exit(1);
  }
  
  if(argc==2){
    std::string input(argv[1]);
    try {
      RooStats::HistFactory::fastDriver(input);
    }
    catch (std::string str) {
      std::cerr << "caught exception: " << str << std::endl ;
    }
    catch( const std::exception& e ) {
      std::cerr << "Caught Exception: " << e.what() << std::endl;
    }
    catch(...) {
      exit(1);
    }
  }
  
  if(argc==3){
    std::string flag(argv[1]);
    std::string input(argv[2]);

    if(flag=="-standard_form") {
      try {
	RooStats::HistFactory::fastDriver(input);
      }
      catch (std::string str) {
	std::cerr << "caught exception: " << str << std::endl ;
      }
      catch( const std::exception& e ) {
	std::cerr << "Caught Exception: " << e.what() << std::endl;
      }
      catch(...) {
	exit(1);
      }
    }
      
    else if(flag=="-number_counting_form") {
      std::cout << "ERROR: 'number_counting_form' is now depricated." << std::endl;
      /*
      try {
	//topDriver(input);
      }
      catch (std::string str) {
	std::cerr << "caught exception: " << str << std::endl ;
      }
      catch( const std::exception& e ) {
	std::cerr << "Caught Exception: " << e.what() << std::endl;
	}
      */
      return 255;
    }
    
    else {
      std::cerr << "Unrecognized flag.  " << std::endl;
      return 255;
    }
  }

  return 0;

}

#endif
