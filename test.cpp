#include <armadillo>
#include <libconfig.h++>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>

using namespace arma;
using namespace std;
using namespace libconfig;

int main(int argc, char **argv){
  Config cfg;
  string file_name;
  file_name = argv[1];
  // Read the file. If there is an error, report it and exit.
  // cout << file_name << endl;
  
  try
  {
    cfg.readFile(file_name);
  }
  catch(const FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    return(EXIT_FAILURE);
  }
  catch(const ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    return(EXIT_FAILURE);
  }

  // Get the store name.
  try
  {
    string name = cfg.lookup("name");
    cout << "Store name: " << name << endl << endl;
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No 'name' setting in configuration file." << endl;
  }



}