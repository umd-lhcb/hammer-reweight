// Author: Yipeng Sun, Alex Fernez
// License: BSD 2-clause
// Last Change: Thu Mar 23, 2023 at 02:12 AM -0400

// Alex note to self: here the filename "default" just refers to the fact that the nominal
// parameter values for the D** decays have been reset to the paper values (ie. not our shifted
// values that we chose because the fit seemed to want to move these FF params far from the paper
// values); the corresponding set nominal values and corresponding variations are thus different
// from the reweighter script in this folder without "default" in the filename

#include <algorithm>
#include <array>
#include <exception>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <stdio.h>

#include <TMath.h>
#include <TString.h>
#include <ROOT/RDataFrame.hxx>

#include <Hammer/Hammer.hh>
#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>
#include <Hammer/Process.hh>

#include <cxxopts.hpp>

#include "const.h"
#include "utils_general.h"
#include "utils_ham.h"

using namespace std;
using ROOT::RDataFrame;
using ROOT::RDF::RNode;
using ROOT::VecOps::RVec;

#define DEBUG_OUT cout << __LINE__ << endl;

///////////////////
// Configurables //
///////////////////

//#define DEBUG_CLI
//#define FORCE_MOMENTUM_CONSERVATION_LEPTONIC
#define RADIATIVE_CORRECTION
#define SOFT_PHOTON_THRESH 0.1

typedef map<vector<Int_t>, unsigned long> DecayFreq;

void setInputFF(Hammer::Hammer& ham, TString run) {
  if (run == "run1") {
    ham.setFFInputScheme({
        {"BD", "ISGW2"},  // 12573010
        {"BD*", "ISGW2"}  // 11574010
        // for run 1, B -> D*MuNu is modelled w/ HQET, which is not implemented
        // in HAMMER
    });
  } else if (run == "run2") {
    ham.setFFInputScheme({{"BD", "CLN_1"},
                          {"BD*", "CLN_1"},
                          {"BD**0*", "ISGW2"},
                          {"BD**1", "ISGW2"},
                          {"BD**1*", "ISGW2"},
                          {"BD**2*", "ISGW2"},
                          {"BsDs**1", "ISGW2"},
                          {"BsDs**2*", "ISGW2"}});

    // 12573001, 12573012
    // clang-format off
    // HQET2(hqetrho2, hqetv1_1, indelta): 1.131 1.035 0.38
    ham.setOptions("BtoDCLN_1: {RhoSq: 1.131, Delta: 0.38, G1: 1.035}");  // HQET3
    // 11574011, 11574021
    // HQET2(hqetrho2, hqetha1_1, hqetr1_1, hqetr2_1, hqetr0_1): 1.122 0.908 1.270 0.852 1.15
    ham.setOptions("BtoD*CLN_1: {RhoSq: 1.122, F1: 0.908, R1: 1.270, R2: 0.852, R0: 1.15}");  // HQET2 11874440, rest of D**
    // ISGW2, which has no configurable parameter
    // clang-format on
  }
}

map<string, complex<double>> specializedWC = {
    {"SM", 1},     {"S_qLlL", 0}, {"S_qRlL", 0}, {"V_qLlL", 0},
    {"V_qRlL", 0}, {"T_qLlL", 0}, {"S_qLlR", 0}, {"S_qRlR", 0},
    {"V_qLlR", 0}, {"V_qRlR", 0}, {"T_qRlR", 0}};

// clang-format off
// +, -, +, -, ...
vector<vector<string>> BtoDVars = {
  // shifting 1-th param in + direction....
  {
    "{ap: [0.01564266061259705, -0.034768835742855116, -0.09341578728622298, 0.0]}",
    "{a0: [0.07921507727285786, -0.20200574935567273, -0.32999651965104015, 0.0]}",
  },
  // shifting 1-th param in - direction....
  {
    "{ap: [0.015677339387402952, -0.03363116425714489, -0.08658421271377702, 0.0]}",
    "{a0: [0.07948492272714215, -0.20799425064432725, -0.13000348034895984, 0.0]}",
  },
  // shifting 2-th param in + direction....
  {
    "{ap: [0.0156957728134992, -0.03274798737619618, -0.11173093154013228, 0.0]}",
    "{a0: [0.07952054730874657, -0.20479549469640895, -0.22925983421113325, 0.0]}",
  },
  // shifting 2-th param in - direction....
  {
    "{ap: [0.015624227186500802, -0.035652012623803825, -0.06826906845986772, 0.0]}",
    "{a0: [0.07917945269125344, -0.20520450530359102, -0.23074016578886677, 0.0]}",
  },
  // shifting 3-th param in + direction....
  {
    "{ap: [0.015686808500749508, -0.031712849151909125, -0.08969198838997502, 0.0]}",
    "{a0: [0.07939582669134693, -0.19132668876796238, -0.22961524687243018, 0.0]}",
  },
  // shifting 3-th param in - direction....
  {
    "{ap: [0.015633191499250493, -0.03668715084809088, -0.09030801161002497, 0.0]}",
    "{a0: [0.07930417330865308, -0.2186733112320376, -0.23038475312756984, 0.0]}",
  },
  // shifting 4-th param in + direction....
  {
    "{ap: [0.015564244420542984, -0.03420228477497694, -0.09000030258603699, 0.0]}",
    "{a0: [0.07887140238985467, -0.20499939148119928, -0.22999994184165481, 0.0]}",
  },
  // shifting 4-th param in - direction....
  {
    "{ap: [0.015755755579457017, -0.03419771522502306, -0.089999697413963, 0.0]}",
    "{a0: [0.07982859761014534, -0.2050006085188007, -0.2300000581583452, 0.0]}",
  },
  // shifting 5-th param in + direction....
  {
    "{ap: [0.015685136124387915, -0.03519607262617201, -0.09006434570610047, 0.0]}",
    "{a0: [0.07914343789044714, -0.2048177906662987, -0.22998668416129328, 0.0]}",
  },
  // shifting 5-th param in - direction....
  {
    "{ap: [0.015634863875612085, -0.03320392737382799, -0.08993565429389952, 0.0]}",
    "{a0: [0.07955656210955286, -0.20518220933370127, -0.23001331583870674, 0.0]}",
  }
};

vector<vector<string>> BtoDstVars = {
  // shifting 1-th param in + direction....
  {
    "{avec: [0.0012794288639999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 1-th param in - direction....
  {
    "{avec: [0.0012021219839999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 2-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.004483799039999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 2-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.006880312319999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 3-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0166209792]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 3-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0320823552]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 4-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.0004912852223999999, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 4-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.0004742777088, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 5-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 0.000247382016, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 5-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, -8.5037568e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 6-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0069576192]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 6-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, -0.0015461375999999997]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 7-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [5.7207091199999995e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 7-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [-9.276825599999997e-06, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 8-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0033241958399999995, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 8-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0013142169599999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 9-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.017780582399999998]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 9-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.054887884799999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 10-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.0021066124799999997, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 10-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.001998382848, -0.00776934144, 2.7057407999999998e-05]}",
  },
  // shifting 11-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00614589696, 2.7057407999999998e-05]}",
  },
  // shifting 11-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00939278592, 2.7057407999999998e-05]}",
  },
  // shifting 12-th param in + direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, 0.034737846527999994]}",
  },
  // shifting 12-th param in - direction....
  {
    "{avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}",
    "{bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}",
    "{cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}",
    "{dvec: [0.002052497664, -0.00776934144, -0.034683731711999996]}",
  }
};

vector<vector<string>> BtoD0starVars = {
  // shifting 1-th param in + direction....
  {
    "{ztp: 3.0}",
    "{zeta1: 0.6}",
  },
  // shifting 1-th param in - direction....
  {
    "{ztp: -2.5999999999999996}",
    "{zeta1: 0.6}",
  },
  // shifting 2-th param in + direction....
  {
    "{ztp: 0.2}",
    "{zeta1: 1.2}",
  },
  // shifting 2-th param in - direction....
  {
    "{ztp: 0.2}",
    "{zeta1: 0.0}",
  }
};

vector<vector<string>> BtoD1Vars = {
  // shifting 1-th param in + direction....
  {
    "{tp: -1.2000000000000002}",
    "{tau1: -0.5}",
    "{tau2: 2.9}",
  },
  // shifting 1-th param in - direction....
  {
    "{tp: -2.0}",
    "{tau1: -0.5}",
    "{tau2: 2.9}",
  },
  // shifting 2-th param in + direction....
  {
    "{tp: -1.6}",
    "{tau1: 0.09999999999999998}",
    "{tau2: 2.9}",
  },
  // shifting 2-th param in - direction....
  {
    "{tp: -1.6}",
    "{tau1: -1.1}",
    "{tau2: 2.9}",
  },
  // shifting 3-th param in + direction....
  {
    "{tp: -1.6}",
    "{tau1: -0.5}",
    "{tau2: 5.699999999999999}",
  },
  // shifting 3-th param in - direction....
  {
    "{tp: -1.6}",
    "{tau1: -0.5}",
    "{tau2: 0.10000000000000009}",
  },
};

vector<vector<string>> BtoD1starVars = {
  // shifting 1-th param in + direction....
  {
    "{ztp: 3.0}",
    "{zeta1: 0.6}",
  },
  // shifting 1-th param in - direction....
  {
    "{ztp: -2.5999999999999996}",
    "{zeta1: 0.6}",
  },
  // shifting 2-th param in + direction....
  {
    "{ztp: 0.2}",
    "{zeta1: 1.2}",
  },
  // shifting 2-th param in - direction....
  {
    "{ztp: 0.2}",
    "{zeta1: 0.0}",
  },
};

vector<vector<string>> BtoD2starVars = {
  // shifting 1-th param in + direction....
  {
    "{tp: -1.2000000000000002}",
    "{tau1: -0.5}",
    "{tau2: 2.9}",
  },
  // shifting 1-th param in - direction....
  {
    "{tp: -2.0}",
    "{tau1: -0.5}",
    "{tau2: 2.9}",
  },
  // shifting 2-th param in + direction....
  {
    "{tp: -1.6}",
    "{tau1: 0.09999999999999998}",
    "{tau2: 2.9}",
  },
  // shifting 2-th param in - direction....
  {
    "{tp: -1.6}",
    "{tau1: -1.1}",
    "{tau2: 2.9}",
  },
  // shifting 3-th param in + direction....
  {
    "{tp: -1.6}",
    "{tau1: -0.5}",
    "{tau2: 5.699999999999999}",
  },
  // shifting 3-th param in - direction....
  {
    "{tp: -1.6}",
    "{tau1: -0.5}",
    "{tau2: 0.10000000000000009}",
  },
};

// Various FF config
map<string, vector<vector<string>>> ffVarSpecs = {
  {"BD", BtoDVars},
  {"BD*", BtoDstVars},
  {"BD**0*", BtoD0starVars},
  {"BD**1", BtoD1Vars},
  {"BD**1*", BtoD1starVars},
  {"BD**2*", BtoD2starVars},
  {"BsDs**1", BtoD1Vars},
  {"BsDs**2*", BtoD2starVars}
};

map<string, string> ffSchemeByDecay = {
  {"BD", "BGL"},
  {"BD*", "BGL"},
  {"BD**0*", "BLR"},
  {"BD**1", "BLR"},
  {"BD**1*", "BLR"},
  {"BD**2*", "BLR"},
  {"BsDs**1", "BLR"},
  {"BsDs**2*", "BLR"}
};

void setBtoDBGLDefault(Hammer::Hammer& ham, const string scheme) {
  ham.setOptions(scheme + ": {ChiT: 0.0005131}");
  ham.setOptions(scheme + ": {ChiL: 0.006332}");
  ham.setOptions(scheme + ": {BcStatesp: [6.329, 6.92, 7.02]}");
  ham.setOptions(scheme + ": {BcStates0: [6.716, 7.121]}");
  ham.setOptions(scheme + ": {ap: [0.01566, -0.0342, -0.09, 0.0]}");
  ham.setOptions(scheme + ": {a0: [0.07935, -0.205, -0.23, 0.0]}");
}

void setBtoDstarBGLDefault(Hammer::Hammer& ham, const string scheme) {
  ham.setOptions(scheme + ": {Vcb: 0.0384}");
  ham.setOptions(scheme + ": {Chim: 0.0003894}");
  ham.setOptions(scheme + ": {Chip: 0.0005131}");
  ham.setOptions(scheme + ": {ChimL: 0.019421}");
  ham.setOptions(scheme + ": {BcStatesf: [6.739, 6.75, 7.145, 7.15]}");
  ham.setOptions(scheme + ": {BcStatesg: [6.329, 6.92, 7.02]}");
  ham.setOptions(scheme + ": {BcStatesP1: [6.275, 6.842, 7.25]}");
  ham.setOptions(scheme + ": {avec: [0.0012407754239999998, -0.005682055679999999, -0.0243516672]}");
  ham.setOptions(scheme + ": {bvec: [0.00048278146559999996, 8.117222399999999e-05, 0.0027057408]}");
  ham.setOptions(scheme + ": {cvec: [2.39651328e-05, 0.0023192063999999996, -0.036334233599999995]}");
  ham.setOptions(scheme + ": {dvec: [0.002052497664, -0.00776934144, 2.7057407999999998e-05]}");
  // Alex note to self: gen_ham_params outputs setOptions(... abcderr: [...10 values...]) too,
  // but I don't think Hammer needs the errors (they weren't set in previous iterations of Hammer
  // reweighting, either)
}

void setBtoD0starBLRDefault(Hammer::Hammer& ham, const string scheme) {
  ham.setOptions(scheme + ": {as: 0.26}");
  ham.setOptions(scheme + ": {mb: 4.71}");
  ham.setOptions(scheme + ": {mc: 1.31}");
  ham.setOptions(scheme + ": {zt1: 0.7}");
  ham.setOptions(scheme + ": {ztp: 0.2}");
  ham.setOptions(scheme + ": {zeta1: 0.6}");
  ham.setOptions(scheme + ": {chi1: 0.0}");
  ham.setOptions(scheme + ": {chi2: 0.0}");
  ham.setOptions(scheme + ": {laB: 0.4}");
  ham.setOptions(scheme + ": {laS: 0.76}");
};

void setBtoD1BLRDefault(Hammer::Hammer& ham, const string scheme) {
  ham.setOptions(scheme + ": {as: 0.26}");
  ham.setOptions(scheme + ": {mb: 4.71}");
  ham.setOptions(scheme + ": {mc: 1.31}");
  ham.setOptions(scheme + ": {t1: 0.7}");
  ham.setOptions(scheme + ": {tp: -1.6}");
  ham.setOptions(scheme + ": {tau1: -0.5}");
  ham.setOptions(scheme + ": {tau2: 2.9}");
  ham.setOptions(scheme + ": {eta1: 0.0}");
  ham.setOptions(scheme + ": {eta2: 0.0}");
  ham.setOptions(scheme + ": {eta3: 0.0}");
  ham.setOptions(scheme + ": {laB: 0.4}");
  ham.setOptions(scheme + ": {laP: 0.8}");
};

void setBtoD1starBLRDefault(Hammer::Hammer& ham, const string scheme) {
  ham.setOptions(scheme + ": {as: 0.26}");
  ham.setOptions(scheme + ": {mb: 4.71}");
  ham.setOptions(scheme + ": {mc: 1.31}");
  ham.setOptions(scheme + ": {zt1: 0.7}");
  ham.setOptions(scheme + ": {ztp: 0.2}");
  ham.setOptions(scheme + ": {zeta1: 0.6}");
  ham.setOptions(scheme + ": {chi1: 0.0}");
  ham.setOptions(scheme + ": {chi2: 0.0}");
  ham.setOptions(scheme + ": {laB: 0.4}");
  ham.setOptions(scheme + ": {laS: 0.76}");
};

void setBtoD2starBLRDefault(Hammer::Hammer& ham, const string scheme) {
  ham.setOptions(scheme + ": {as: 0.26}");
  ham.setOptions(scheme + ": {mb: 4.71}");
  ham.setOptions(scheme + ": {mc: 1.31}");
  ham.setOptions(scheme + ": {t1: 0.7}");
  ham.setOptions(scheme + ": {tp: -1.6}");
  ham.setOptions(scheme + ": {tau1: -0.5}");
  ham.setOptions(scheme + ": {tau2: 2.9}");
  ham.setOptions(scheme + ": {eta1: 0.0}");
  ham.setOptions(scheme + ": {eta2: 0.0}");
  ham.setOptions(scheme + ": {eta3: 0.0}");
  ham.setOptions(scheme + ": {laB: 0.4}");
  ham.setOptions(scheme + ": {laP: 0.8}");
};

map<string, function<void(Hammer::Hammer&, const string)>>
ffSchemeDefaultsByDecay = {
  {"BD", setBtoDBGLDefault},
  {"BD*", setBtoDstarBGLDefault},
  {"BD**0*", setBtoD0starBLRDefault},
  {"BD**1", setBtoD1BLRDefault},
  {"BD**1*", setBtoD1starBLRDefault},
  {"BD**2*", setBtoD2starBLRDefault},
  {"BsDs**1", setBtoD1BLRDefault},
  {"BsDs**2*", setBtoD2starBLRDefault}
};
// clang-format on

const int numOfFFVar = 24;

string decayDescr(const string decay) {
  // the pattern seems to be you just need to insert the word "to"
  // I'm going to assume the fact (that should be fine for our analysis) that the
  // decay string only contains 'D' once (where the name of the daughter begins)
  auto daughterPos = decay.find("D");
  auto bMeson = decay.substr(0,daughterPos);
  auto daughter = decay.substr(daughterPos);
  return bMeson + "to" + daughter;
}

vector<string> setOutputFF(Hammer::Hammer& ham) {
  vector<string> hamFFSchemes{"OutputFF"};

  ham.addFFScheme("OutputFF", ffSchemeByDecay);
  // Set defaults
  for (auto [decay, ffName] : ffSchemeByDecay) {
    auto fullDescr = decayDescr(decay) + ffName;
    cout << "Decay: " << decay << "; default FF: " << fullDescr << endl;
    ffSchemeDefaultsByDecay[decay](ham, fullDescr);
  }

  // Set variations
  for (int i = 1; i <= numOfFFVar; i++) {
    cout << "Configuring FF scheme: OutputFFVar" + to_string(i) << endl;
    map<string, string> schemes{};
    for (auto const& [decay, vars] : ffVarSpecs) {
      if (i <= vars.size()) {
        // Need to reweight the decay in this HAMMER scheme
        auto ffName    = ffSchemeByDecay[decay] + "_" + to_string(i);
        auto descr     = decayDescr(decay);
        schemes[decay] = ffName;
        // Configure the FF scheme defaults for this decay
        cout << "  Variation for decay: " << decay
             << "; with FF: " << descr + ffName << endl;
        ffSchemeDefaultsByDecay[decay](ham, descr + ffName);
        for (auto const& shift : vars[i - 1])
          ham.setOptions(descr + ffName + ": " +
                         shift);  // configure FF variations
      }
    }
    ham.addFFScheme("OutputFFVar" + to_string(i), schemes);
  }

  return hamFFSchemes;
}

void setDecays(Hammer::Hammer& ham) {
  ham.includeDecay("BDTauNu");
  ham.includeDecay("BDMuNu");

  ham.includeDecay("BD*TauNu");
  ham.includeDecay("BD*MuNu");

  ham.includeDecay("BD**0*TauNu");
  ham.includeDecay("BD**1TauNu");
  ham.includeDecay("BD**1*TauNu");
  ham.includeDecay("BD**2*TauNu");

  ham.includeDecay("BD**0*MuNu");
  ham.includeDecay("BD**1MuNu");
  ham.includeDecay("BD**1*MuNu");
  ham.includeDecay("BD**2*MuNu");

  ham.includeDecay("BsDs**1MuNu");
  ham.includeDecay("BsDs**2*MuNu");
}

/////////////
// Helpers //
/////////////

// For muting hammer logging
// NOTE: Hammer always uses 'stdout' for logging
int muteStdout() {
  fflush(stdout);
  int  fd      = dup(STDOUT_FILENO);
  auto devNull = freopen("/dev/null", "w", stdout);
  return fd;
}

void restoreStdout(int fd) {
  fflush(stdout);
  dup2(fd, fileno(stdout));
  close(fd);
}

/////////////
// Filters //
/////////////

bool truthMatchOk(double q2True, bool isTauDecay, int bMesonId, int dau1Id,
                  int dau2Id, int muID) {
  double q2Min = 100 * 100;
  if (isTauDecay) q2Min = 1700 * 1700;

  // we require that there's ONE and ONLY ONE D meson
  return findIn(LEGAL_B_MESON_IDS, TMath::Abs(bMesonId)) && q2True > q2Min &&
         isDMeson(TMath::Abs(dau1Id)) && !isDMeson(TMath::Abs(dau2Id)) &&
         TMath::Abs(muID) == 13;
}

/////////////////////////
// Branch name helpers //
/////////////////////////

vector<string> getBTrueP(string particle, string bId) {
  return setBrPrefix(particle, {"TRUEP_E", "TRUEP_X", "TRUEP_Y", "TRUEP_Z"},
                     {bId});
}

vector<string> getDauTrueP(string particle, string dau) {
  return setBrPrefix(particle + "_TrueHadron_" + dau,
                     {"PE", "PX", "PY", "PZ", "ID"});
}

////////////////////
// FSR correction //
////////////////////

typedef tuple<double, double, double, double, int> HamPartCtn;

bool isSoftPhoton(Hammer::Particle part) {
  if (part.p().E() < SOFT_PHOTON_THRESH) return true;
  return false;
}

bool isSoftPhoton(HamPartCtn part) {
  auto pe = get<0>(part);
  if (pe < SOFT_PHOTON_THRESH) return true;
  return false;
}

vector<HamPartCtn> buildPhotonVec(RVec<float>& arrPe, RVec<float>& arrPx,
                                  RVec<float>& arrPy, RVec<float>& arrPz,
                                  RVec<float>& arrId, int size) {
  vector<HamPartCtn> result{};
  for (auto idx = 0; idx != size; idx++) {
    auto part = buildPartVec(arrPe[idx], arrPx[idx], arrPy[idx], arrPz[idx],
                             static_cast<int>(arrId[idx]));
    if (!isSoftPhoton(part)) result.emplace_back(part);
  }
  return result;
}

string printP(double pe, double px, double py, double pz) {
  stringstream buffer;
  buffer << pe << "," << px << "," << py << "," << pz;
  return buffer.str();
}

string printP(Hammer::Particle part) {
  stringstream buffer;
  auto         p = part.p();
  buffer << p.E() << "," << p.px() << "," << p.py() << "," << p.pz();
  return buffer.str();
}

string addRadiativePhotons(Hammer::Process& proc, Hammer::ParticleIndices& idx,
                           int refMomId, vector<HamPartCtn>& photons) {
  stringstream buffer;

  for (const auto& p : photons) {
    auto [pe, px, py, pz, photonMomId] = p;
    if (photonMomId == refMomId || photonMomId == -refMomId) {
      auto partPhoton = buildHamPart(pe, px, py, pz, 22);
      idx.emplace_back(proc.addParticle(partPhoton));
      buffer << "  Adding photon: " << printP(pe, px, py, pz) << " to "
             << refMomId << endl;
    }
  }

  return buffer.str();
}

/////////////////
// Reweighting //
/////////////////

pair<RNode, vector<string>> prepAuxOutput(RNode df, string bMesonName) {
  auto outputBrs = vector<string>{};

  // truth variables
  df = df.Define("q2_true", bMesonName + "_True_Q2 / 1000 / 1000");
  outputBrs.emplace_back("q2_true");

  df = df.Define("is_tau", bMesonName + "_True_IsTauDecay");
  outputBrs.emplace_back("is_tau");

  auto kinematicSuffix = {"_PE", "_PX", "_PY", "_PZ"};
  for (int i = 0; i < 2; i++) {
    auto partName = "d_meson" + to_string(i + 1);

    df = df.Define(partName + "_true_id",
                   bMesonName + "_TrueHadron_D" + to_string(i) + "_ID");
    outputBrs.emplace_back(partName + "_true_id");

    auto sIdx         = to_string(i);
    auto kinematicBrs = vector<string>{};
    for (const auto& suf : kinematicSuffix) {
      kinematicBrs.emplace_back("TrueHadron_D" + sIdx + suf);
    }

    df = df.Define(partName + "_true_m", invM,
                   setBrPrefix(bMesonName, kinematicBrs));
    outputBrs.emplace_back(partName + "_true_m");
  }

  // truth-matching ok
  df = df.Define("ham_tm_ok", truthMatchOk,
                 setBrPrefix(bMesonName,
                             {"True_Q2", "True_IsTauDecay", "TRUEID",
                              "TrueHadron_D0_ID", "TrueHadron_D1_ID"},
                             {"mu_TRUEID"}));
  outputBrs.emplace_back("ham_tm_ok");

  return {df, outputBrs};
}

// NOTE: Refine branch names if input ntuples have different tree structrue!
pair<RNode, vector<string>> prepHamInput(RNode df, string bMesonName) {
  auto outputBrs = vector<string>{};

  // B meson
  df = df.Define("b_id_fixed", bIdFix,
                 setBrPrefix(bMesonName, {"TRUEID", "TrueHadron_D0_ID"}));
  df = df.Define("part_B", buildPartVec, getBTrueP(bMesonName, "b_id_fixed"));

  // D meson (can be a D**, D*, or D)
  df = df.Define("part_D", buildPartVec, getDauTrueP(bMesonName, "D0"));

  // D meson daughters
  df =
      df.Define("part_D_dau0", buildPartVec, getDauTrueP(bMesonName, "D0_GD0"));
  df =
      df.Define("part_D_dau1", buildPartVec, getDauTrueP(bMesonName, "D0_GD1"));
  df =
      df.Define("part_D_dau2", buildPartVec, getDauTrueP(bMesonName, "D0_GD2"));

#ifdef RADIATIVE_CORRECTION
  // the radiative photons container
  df = df.Define("part_photon_arr", buildPhotonVec,
                 setBrPrefix(bMesonName, {"MCTrue_gamma_E", "MCTrue_gamma_PX",
                                          "MCTrue_gamma_PY", "MCTrue_gamma_PZ",
                                          "MCTrue_gamma_mother_ID",
                                          "MCTrue_gamma_ArrayLength"}));

  // Tau/Mu, Nu_Tau/Nu_Mu associated w/ B -> D decay
  df = df.Define("part_Tau_id", tauIdFix, {"mu_TRUEID"});
  df = df.Define(
      "part_Tau", buildPartVec,
      setBrPrefix(bMesonName,
                  {"TrueTau_PE", "TrueTau_PX", "TrueTau_PY", "TrueTau_PZ"},
                  {"part_Tau_id"}));
  df = df.Define("part_Mu_id", muIdFix, {"mu_TRUEID"});
  df = df.Define(
      "part_Mu", buildPartVec,
      setBrPrefix(bMesonName,
                  {"TrueMu_PE", "TrueMu_PX", "TrueMu_PY", "TrueMu_PZ"},
                  {"part_Mu_id"}));
  df = df.Define("part_L",
                 [](HamPartCtn pTau, HamPartCtn pMu, bool isTau) {
                   if (isTau) return pTau;
                   return pMu;
                 },
                 {"part_Tau", "part_Mu", bMesonName + "_True_IsTauDecay"});

  df = df.Define("part_NuL_id", nuIdFix,
                 {"mu_TRUEID", bMesonName + "_True_IsTauDecay"});
  df = df.Define("part_NuL", buildPartVec,
                 setBrPrefix(bMesonName,
                             {"TrueNeutrino_PE", "TrueNeutrino_PX",
                              "TrueNeutrino_PY", "TrueNeutrino_PZ"},
                             {"part_NuL_id"}));
  // for leptons produced by primary B -> D decays, ALWAYS use part_L and
  // part_NuL

  // TauNu, MuNu associated w/ Tau -> Mu decay
  df = df.Define("part_NuTau_id", tauNuTauIdFix, {"mu_TRUEID"});
  df = df.Define("part_NuTau", buildPartVec,
                 setBrPrefix(bMesonName,
                             {"TrueTauNuTau_PE", "TrueTauNuTau_PX",
                              "TrueTauNuTau_PY", "TrueTauNuTau_PZ"},
                             {"part_NuTau_id"}));
  df = df.Define("part_NuMu_id", tauNuMuIdFix, {"mu_TRUEID"});
  df = df.Define("part_NuMu", buildPartVec,
                 setBrPrefix(bMesonName,
                             {"TrueTauNuMu_PE", "TrueTauNuMu_PX",
                              "TrueTauNuMu_PY", "TrueTauNuMu_PZ"},
                             {"part_NuMu_id"}));
#else
  // this is just a place holder in case we don't add radiative photons back!
  df = df.Define("part_photon_arr", []() { return vector<HamPartCtn>{}; }, {});
#endif

  return {df, outputBrs};
}

auto reweightWrapper(Hammer::Hammer& ham, unsigned long& numOfEvt,
                     unsigned long& numOfEvtOk, vector<string>& schemes) {
  return [&](bool truthMatchOk, bool isTau, HamPartCtn pB, HamPartCtn pD,
             HamPartCtn pDDau0, HamPartCtn pDDau1, HamPartCtn pDDau2,
             HamPartCtn pL, HamPartCtn pNuL, HamPartCtn pMu, HamPartCtn pNuMu,
             HamPartCtn pNuTau, vector<HamPartCtn> pPhotons) {
    bool                      hamOk = true;
    double                    wtFF  = 1.0;
    array<double, numOfFFVar> wtFFVars;
    fill_n(wtFFVars.begin(), numOfFFVar, 1.0);  // all default to 1.0

    string debugMsg = "====\n";
    numOfEvt += 1;
    if (!truthMatchOk)
      // FIXME: it's ugly! but not sure RDataFrame frame can deduce type if we
      // use anything fancy.
      return tuple<bool, double, double, double, double, double, double, double,
                   double, double, double, double, double, double, double,
                   double, double, double, double, double, double, double,
                   double, double, double, double>{
          false,        wtFF,         wtFFVars[0],  wtFFVars[1],  wtFFVars[2],
          wtFFVars[3],  wtFFVars[4],  wtFFVars[5],  wtFFVars[6],  wtFFVars[7],
          wtFFVars[8],  wtFFVars[9],  wtFFVars[10], wtFFVars[11], wtFFVars[12],
          wtFFVars[13], wtFFVars[14], wtFFVars[15], wtFFVars[16], wtFFVars[17],
          wtFFVars[18], wtFFVars[19], wtFFVars[20], wtFFVars[21], wtFFVars[22],
          wtFFVars[23]};

    Hammer::Process proc;
    auto            partB   = buildHamPart(pB);
    auto            partD   = buildHamPart(pD);
    auto            partL   = buildHamPart(pL);
    auto            partNuL = buildHamPart(pNuL);
    auto particles = vector<Hammer::Particle>{partB, partD, partL, partNuL};

    // add B meson
    auto partBIdx = proc.addParticle(partB);

    debugMsg += "  B meson 4-mom: " + printP(partB) + '\n';
    debugMsg += "  D meson 4-mom: " + printP(partD) + '\n';
    debugMsg += "  primary charged lepton 4-mom: " + printP(partL) + '\n';
    debugMsg += "  primary neutrino 4-mom: " + printP(partNuL) + '\n';

    // add direct B daughters
    auto                    partDIdx   = proc.addParticle(partD);
    auto                    partLIdx   = proc.addParticle(partL);
    auto                    partNuLIdx = proc.addParticle(partNuL);
    Hammer::ParticleIndices partBDauIdx{partDIdx, partLIdx, partNuLIdx};
#ifdef RADIATIVE_CORRECTION
    debugMsg += addRadiativePhotons(proc, partBDauIdx, partB.pdgId(), pPhotons);
#endif
    proc.addVertex(partBIdx, partBDauIdx);

    // in case of a D*, add its daughters as well
    Hammer::ParticleIndices partDDauIdx{};
#ifdef RADIATIVE_CORRECTION
    debugMsg += addRadiativePhotons(proc, partDDauIdx, partD.pdgId(), pPhotons);
#endif

    bool isDst = isDstMeson(get<4>(pD));
    if (isDst) {
      auto partDDaus = {buildHamPart(pDDau0), buildHamPart(pDDau1),
                        buildHamPart(pDDau2)};
      for (const auto& p : partDDaus) {
        // don't add placeholder particle
        // don't add photons as it will be a duplicate
        if (p.pdgId() == 0 || p.pdgId() == 22) continue;
        partDDauIdx.emplace_back(proc.addParticle(p));
        particles.emplace_back(p);
        debugMsg += "  D daughters: " + printP(p) + "\n";
      }
    }
    if (partDDauIdx.size()) proc.addVertex(partDIdx, partDDauIdx);

    // in case of a Tau, add its daughters
    Hammer::ParticleIndices partLDauIdx{};
#ifdef RADIATIVE_CORRECTION
    debugMsg += addRadiativePhotons(proc, partLDauIdx, partL.pdgId(), pPhotons);
#endif

    if (isTau) {
      auto partLDaus = {buildHamPart(pMu), buildHamPart(pNuMu),
                        buildHamPart(pNuTau)};
      for (const auto& p : partLDaus) {
        partLDauIdx.emplace_back(proc.addParticle(p));
        particles.emplace_back(p);
        debugMsg += "  secondary leptons: " + printP(p) + "\n";
      }
    }
    if (partLDauIdx.size()) proc.addVertex(partLIdx, partLDauIdx);

    // make sure invariant mass is not negative
    auto partInvMOk = vector<bool>{};
    for (const auto& p : particles) {
      if (p.p().mass() >= 0)
        partInvMOk.emplace_back(true);
      else
        partInvMOk.emplace_back(false);
    }
    bool allPartsOk =
        find(partInvMOk.begin(), partInvMOk.end(), false) == partInvMOk.end();
    if (!allPartsOk) {
      cout << "  WARN: Bad kinematics for candidate: " << numOfEvt << endl;
      hamOk = false;
    }

#ifdef DEBUG_CLI
    cout << debugMsg;
#endif

    // add the whole decay chain to hammer and see if it likes it
    if (hamOk) {
      ham.initEvent();
      int procId = 0;
      try {
        procId = ham.addProcess(proc);
      } catch (const exception& e) {
        cout << "  WARN: HAMMER doesn't add process properly: " << numOfEvt
             << endl;
        cout << e.what() << endl;
        hamOk = false;
      }
      if (procId == 0) hamOk = false;
    }

    // compute FF weight
    auto nominalFFScheme = schemes[0];
    if (hamOk) {
      try {
        ham.processEvent();
        wtFF = ham.getWeight(nominalFFScheme);
      } catch (const exception& e) {
        cout << "  WARN: HAMMER doesn't like candidate for reweighting: "
             << numOfEvt << endl;
        cout << e.what() << endl;
        hamOk = false;
      }

      if (isnan(wtFF) || isinf(wtFF)) hamOk = false;

      if (hamOk) {
        numOfEvtOk += 1;
        // compute various FF variation weights here
        // failure's allowed because not all FF schemes require 20 variations!
        auto fd = muteStdout();
        try {
          for (int i = 0; i < numOfFFVar; i++) {
            auto ffName = "OutputFFVar" + to_string(i + 1);
            wtFFVars[i] = ham.getWeight(ffName);
          }
        } catch (const exception& e) {
        }
        restoreStdout(fd);
      }
    }

#ifdef DEBUG_CLI
    if (hamOk) cout << "  FF weight: " << wtFF << endl;
#endif

    // FIXME: it's ugly! but not sure RDataFrame frame can deduce type if we
    // use anything fancy.
    return tuple<bool, double, double, double, double, double, double, double,
                 double, double, double, double, double, double, double, double,
                 double, double, double, double, double, double, double, double,
                 double, double>{
        hamOk,        wtFF,         wtFFVars[0],  wtFFVars[1],  wtFFVars[2],
        wtFFVars[3],  wtFFVars[4],  wtFFVars[5],  wtFFVars[6],  wtFFVars[7],
        wtFFVars[8],  wtFFVars[9],  wtFFVars[10], wtFFVars[11], wtFFVars[12],
        wtFFVars[13], wtFFVars[14], wtFFVars[15], wtFFVars[16], wtFFVars[17],
        wtFFVars[18], wtFFVars[19], wtFFVars[20], wtFFVars[21], wtFFVars[22],
        wtFFVars[23]};
  };
}

//////////
// Main //
//////////

int main(int argc, char** argv) {
  cxxopts::Options argOpts("ReweightRDX", "Reweight RDX FF w/ HAMMER.");

  // clang-format off
  argOpts.add_options()
    // positional
    ("ntpIn", "specify input ntuple.", cxxopts::value<string>())
    ("ntpOut", "specify output ntuple.", cxxopts::value<string>())
    ("extra", "unused.", cxxopts::value<vector<string>>())
    // keyword
    ("h,help", "print help.")
    ("t,trees", "specify tree name.",
     cxxopts::value<vector<string>>()
     ->default_value("TupleBminus/DecayTree,TupleB0/DecayTree"))
    ("b,bMesons", "specify B meson name.",
     cxxopts::value<vector<string>>()->default_value("b,b0"))
    ("r,run", "specify run.", cxxopts::value<string>()->default_value("run2"))
  ;
  // setup positional argument
  argOpts.parse_positional({"ntpIn", "ntpOut", "extra"});
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  auto ntpIn   = parsedArgs["ntpIn"].as<string>();
  auto ntpOut  = parsedArgs["ntpOut"].as<string>();
  auto trees   = parsedArgs["trees"].as<vector<string>>();
  auto bMesons = parsedArgs["bMesons"].as<vector<string>>();
  auto run     = parsedArgs["run"].as<string>();

  Hammer::Hammer ham{};

  setDecays(ham);
  setInputFF(ham, run);
  auto ffSchemes = setOutputFF(ham);

  ham.setUnits("MeV");
  ham.setOptions("ProcessCalc: {CheckForNaNs: true}");
  ham.initRun();

  // only use SM Wilson coefficients
  ham.specializeWCInWeights("BtoCTauNu", specializedWC);
  ham.specializeWCInWeights("BtoCMuNu", specializedWC);

  // output option
  auto writeOpts  = ROOT::RDF::RSnapshotOptions{};
  writeOpts.fMode = "UPDATE";

  for (int idx = 0; idx != trees.size(); idx++) {
    RNode          df     = static_cast<RNode>(RDataFrame(trees[idx], ntpIn));
    auto           bMeson = bMesons[idx];
    vector<string> outputBrs{"runNumber", "eventNumber"};
    unsigned long  numOfEvt   = 0;
    unsigned long  numOfEvtOk = 0;

    cout << "Handling " << trees[idx] << " with B meson name " << bMeson
         << endl;

    // prepare aux output branches like q2_true
    auto [dfAux, outputBrsAux] = prepAuxOutput(df, bMeson);
    df                         = dfAux;
    for (const auto& br : outputBrsAux) outputBrs.emplace_back(br);

    // prepare HAMMER particles
    tie(df, ignore) = prepHamInput(df, bMeson);

    // reweight FF
    auto reweight = reweightWrapper(ham, numOfEvt, numOfEvtOk, ffSchemes);
    df            = df.Define("ff_result", reweight,
                   {"ham_tm_ok", "is_tau", "part_B", "part_D", "part_D_dau0",
                    "part_D_dau1", "part_D_dau2", "part_L", "part_NuL",
                    "part_Mu", "part_NuMu", "part_NuTau", "part_photon_arr"});
    df            = df.Define("ham_ok", "get<0>(ff_result)");
    df            = df.Define("wff_dstnocorr", "get<1>(ff_result)");
    for (int i = 0; i < numOfFFVar; i++) {
      auto outputBrName = "wff_dstnocorr_var" + to_string(i + 1);
      df = df.Define(outputBrName, "get<" + to_string(i + 2) + ">(ff_result)");
      outputBrs.emplace_back(outputBrName);
    }
    outputBrs.emplace_back("wff_dstnocorr");

    df.Snapshot(trees[idx], ntpOut, outputBrs, writeOpts);

    cout << "Total number of candidates: " << numOfEvt << endl;
    cout << "Hammer reweighted candidates: " << numOfEvtOk << endl;
    cout << "Reweighted fraction: "
         << static_cast<float>(numOfEvtOk) / static_cast<float>(numOfEvt)
         << endl;
  }
}
