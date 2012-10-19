
#ifndef HISTFACTORY_ASIMOV_H
#define HISTFACTORY_ASIMOV_H

#include <string>
#include <vector>
#include <map>

#include "RooWorkspace.h"

namespace RooStats{
  namespace HistFactory {

    class Asimov {

    public:
      
      Asimov() {;}
      Asimov(std::string Name) : fName(Name) {;}

      void ConfigureWorkspace( RooWorkspace* );

      std::string GetName() { return fName; }

      void SetFixedParam(const std::string& param, bool constant=true) { fParamsToFix[param] = constant; }
      void SetParamValue(const std::string& param, double value) { fParamValsToSet[param] = value; }
      
      std::map< std::string, bool >& GetParamsToFix() { return fParamsToFix; }
      std::map< std::string, double >& GetParamsToSet() { return fParamValsToSet; }

    protected:

      std::string fName;

      std::map<std::string, bool> fParamsToFix;
      std::map< std::string, double > fParamValsToSet;

    };


  } // namespace HistFactory
} // namespace RooStats

#endif
