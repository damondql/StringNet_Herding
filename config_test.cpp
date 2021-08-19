#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;

int main(int argc, char **argv) {
    Config cfg;
    // string file_name;
    // file_name = argv[1];
    try
    {
    cfg.readFile("example.cfg");
    }
    catch(const FileIOException &fioex)
    {
    std::cerr << "File I/O error" << std::endl;
    return(EXIT_FAILURE);
    }
    catch(const ParseException &pex)
    {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                << " - " << pex.getError() << std::endl;
    return(EXIT_FAILURE);
    }
    
    double a;
    cfg.lookupValue("bg",a);
    cout << "a: " << a << endl;

    double bg;
    cfg.lookupValue("bg",bg);
    cout << "bg: " << bg << endl;
}