
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;

// This example demonstrates the handling of parsing errors in
// 'invalid.cfg'.

int main(int argc, char **argv)
{
  Config cfg;

  try
  {
    cfg.readFile("/home/dasclab/catkin_ws_pxdrones/src/stringnet_herding/src/config/origin_config.cfg");
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

  double dt = 0;
  cfg.lookupValue("dt",dt);
  cout << "dt: " <<dt << endl;

  double Rs;
  cfg.lookupValue("Rs", Rs);
  cout << "Rs: " <<Rs << endl;

  double rho_P = 0;
  cfg.lookupValue("rho_P", rho_P);
  cout << "rho_P: " << rho_P <<endl;

  double wW = 0;
  cfg.lookupValue("wW", wW);
  cout<< "wW: " << wW << endl;

  double lW = 0;
  cfg.lookupValue("lW", lW);
  cout << "lW: " << lW << endl;

  double rho_A = 0;
  cfg.lookupValue("rho_A", rho_A);
  cout << "rho_A: " << rho_A << endl;

  double rho_c_A = 0;
  cfg.lookupValue("rho_c_A", rho_c_A);
  cout << "rho_c_A: " << rho_c_A << endl;

  double C_d = 0;
  cfg.lookupValue("C_d", C_d);
  cout << "C_d: " << C_d << endl;

  double vma = 0;
  cfg.lookupValue("vma", vma);
  cout << "vma: " << vma << endl;

  double rDmin = 0;
  cfg.lookupValue("rDmin", rDmin);
  cout << "rDmin: " << rDmin << endl;

  double rho_D = 0;
  cfg.lookupValue("rho_D", rho_D);
  cout << "rho_D: " << rho_D << endl;

  double rho_sn_max = 0;
  cfg.lookupValue("rho_sn_max", rho_sn_max);
  cout << "rho_sn_max: " << rho_sn_max << endl;

  int largeP = 0;
  cfg.lookupValue("largeP", largeP);
  cout << "largeP: " << largeP << endl;

  double kr1 = 0;
  cfg.lookupValue("kr1", kr1);
  cout << "kr1: " << kr1 << endl;

  double kr2 = 0;
  cfg.lookupValue("kr2", kr2);
  cout << "kr2: " << kr2 << endl;

  double kv1 = 0;
  cfg.lookupValue("kv1", kv1);
  cout << "kv1: " << kv1 << endl;

  double kv2 = 0;
  cfg.lookupValue("kv2", kv2);
  cout << "kv2: " << kv2 << endl << endl;

  double kAFr = 0;
  cfg.lookupValue("kAFr", kAFr);
  cout << "kAFr: " << kAFr << endl;

  double kAOr = 0;
  cfg.lookupValue("kAOr", kAOr);
  cout << "kAOr: " << kAOr << endl;

  double kAOr2 = 0;
  cfg.lookupValue("kAOr2", kAOr2);
  cout << "kAOr2: " << kAOr2 << endl;

  double kADr = 0;
  cfg.lookupValue("kADr", kADr);
  cout << "kADr: " << kADr << endl;

  double kADv = 0;
  cfg.lookupValue("kADv", kADv);
  cout << "kADv: " << kADv << endl;

  double kAFv = 0;
  cfg.lookupValue("kAFv", kAFv);
  cout << "kAFv: " << kAFv << endl;

  double alphaAFv = 0;
  cfg.lookupValue("alphaAFv", alphaAFv);
  cout << "alphaAFv: " << alphaAFv << endl;

  double kAOv = 0;
  cfg.lookupValue("kAOv", kAOv);
  cout << "kAOv: " << kAOv << endl;

  double kAOv2 = 0;
  cfg.lookupValue("kAOv2", kAOv2);
  cout << "kAOv2: " << kAOv2 << endl;

  double alphaAOv = 0;
  cfg.lookupValue("alphaAOv", alphaAOv);
  cout << "alphaAOv: " << alphaAOv << endl;

  double alphaADv = 0;
  cfg.lookupValue("alphaADv", alphaADv);
  cout << "alphaADv: " << alphaADv << endl;

  double kAPr = 0;
  cfg.lookupValue("kAPr", kAPr);
  cout << "kAPr: " << kAPr << endl;

  double kAPv = 0;
  cfg.lookupValue("kAPv", kAPv);
  cout << "kAPv: " << kAPv << endl << endl;

  double kDOr = 0;
  cfg.lookupValue("kDOr", kDOr);
  cout << "kDOr: " << kDOr << endl;

  double kDDr = 0;
  cfg.lookupValue("kDDr", kDDr);
  cout << "kDDr: " << kDDr << endl;

  double kDOv1 = 0;
  cfg.lookupValue("kDOv1", kDOv1);
  cout << "kDOv1: " << kDOv1 << endl;

  double kDOv2 = 0;
  cfg.lookupValue("kDOv2", kDOv2);
  cout << "kDOv2: " << kDOv2 << endl;

  double kDFr = 0;
  cfg.lookupValue("kDFr", kDFr);
  cout << "kDFr: " << kDFr << endl;

  double kDFv = 0;
  cfg.lookupValue("kDFv", kDFv);
  cout << "kDFv: " << kDFv << endl;

  double alphaDFv = 0;
  cfg.lookupValue("alphaDFv", alphaDFv);
  cout << "alphaDFv: " << alphaDFv << endl;

  double alphaDOv = 0;
  cfg.lookupValue("alphaDOv", alphaDOv);
  cout << "alphaDOv: " << alphaDOv << endl;

  double kDDv = 0;
  cfg.lookupValue("kDDv", kDDv);
  cout << "kDDv: " << kDDv << endl;

  double alphaDDv = 0;
  cfg.lookupValue("alphaDDv", alphaDDv);
  cout << "alphaDDv: " << alphaDDv << endl;

  double kDRr = 0;
  cfg.lookupValue("kDRr", kDRr);
  cout << "kDRr: " << kDRr << endl;

  double kDRv = 0;
  cfg.lookupValue("kDRv", kDRv);
  cout << "kDRv: " << kDRv << endl;

  double kDFphi = 0;
  cfg.lookupValue("kDFphi", kDFphi);
  cout << "kDFphi: " << kDFphi << endl;

  double kDFphid = 0;
  cfg.lookupValue("kDFphid", kDFphid);
  cout << "kDFphid: " << kDFphid << endl;

  double kDFphir = 0;
  cfg.lookupValue("kDFphir", kDFphir);
  cout << "kDFphir: " << kDFphir << endl;

  double kDFphiv = 0;
  cfg.lookupValue("kDFphiv", kDFphiv);
  cout << "kDFphiv: " << kDFphiv << endl;

  double kDDesr = 0;
  cfg.lookupValue("kDDesr", kDDesr);
  cout << "kDDesr: " << kDDesr << endl;

  double kDDesv = 0;
  cfg.lookupValue("kDDesv", kDDesv);
  cout << "kDDesv: " << kDDesv << endl << endl;

  double bd = 0;
  cfg.lookupValue("bd", bd);
  cout << "bd: " << bd << endl;

  double GO = 0;
  cfg.lookupValue("GO", GO);
  cout << "GO: " << GO << endl;

  double E_OF1 = 0;
  cfg.lookupValue("E_OF1", E_OF1);
  cout << "E_OF1: " << E_OF1 << endl;

  double E_OF2 = 0;
  cfg.lookupValue("E_OF2", E_OF2);
  cout << "E_OF2: " << E_OF2 << endl;

  return(EXIT_SUCCESS);
}


