// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Wed Nov 02, 2022 at 06:28 AM -0400

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
                          {"BD**2*", "ISGW2"}});

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
    "{avec: [0.0013093595787969515, -0.0052989262672254065, -0.015913932066091246]}",
    "{bvec: [0.0005197681308535841, 0.00016603505416844483, -0.0006321111943059936]}",
    "{cvec: [-5.177382160552352e-06, 0.0039643426962044336]}",
    "{dvec: [0.002193046250379089, -0.007623983661755219]}",
  },
  // shifting 1-th param in - direction....
  {
    "{avec: [0.0013558152412030484, -0.006899052532774594, -0.03421474793390876]}",
    "{bvec: [0.0005178955451464161, 0.00014309180583155521, 0.002303067194305994]}",
    "{cvec: [1.770955216055235e-05, 0.002552385703795567]}",
    "{dvec: [0.002209922809620911, -0.008667837338244782]}",
  },
  // shifting 2-th param in + direction....
  {
    "{avec: [0.0013397362480110139, -0.006521632793600443, -0.02421161899889002]}",
    "{bvec: [0.0005189331126016106, 2.9282075140068343e-05, 0.0056238390883810215]}",
    "{cvec: [2.935628992630371e-06, 0.003093501118943172]}",
    "{dvec: [0.00221141314116393, -0.008757525011895096]}",
  },
  // shifting 2-th param in - direction....
  {
    "{avec: [0.001325438571988986, -0.005676346006399558, -0.025917061001109982]}",
    "{bvec: [0.0005187305633983896, 0.0002798447848599317, -0.003952883088381022]}",
    "{cvec: [9.596541007369628e-06, 0.0034232272810568285]}",
    "{dvec: [0.00219155591883607, -0.007534295988104905]}",
  },
  // shifting 3-th param in + direction....
  {
    "{avec: [0.0013144346207163281, -0.005254689981076678, -0.0250185679052699]}",
    "{bvec: [0.0005203025768633124, 0.00013216463996638218, 0.0007760568769642677]}",
    "{cvec: [2.0998023917081746e-05, 0.00207355363197959]}",
    "{dvec: [0.002218941226886135, -0.00880679120989766]}",
  },
  // shifting 3-th param in - direction....
  {
    "{avec: [0.0013507401992836718, -0.006943288818923322, -0.025110112094730103]}",
    "{bvec: [0.0005173610991366878, 0.00017696222003361786, 0.0008948991230357324]}",
    "{cvec: [-8.465853917081748e-06, 0.0044431747680204105]}",
    "{dvec: [0.002184027833113865, -0.007485029790102343]}",
  },
  // shifting 4-th param in + direction....
  {
    "{avec: [0.001323811332292626, -0.00573634484126035, -0.02510240054335281]}",
    "{bvec: [0.0005190844466766717, 0.00018295337353295966, 0.0010424221684713713]}",
    "{cvec: [2.196485973692794e-05, 0.0026878612381293792]}",
    "{dvec: [0.0021847908345044085, -0.0066818895558519685]}",
  },
  // shifting 4-th param in - direction....
  {
    "{avec: [0.001341363487707374, -0.006461633958739651, -0.02502627945664719]}",
    "{bvec: [0.0005185792293233285, 0.00012617348646704038, 0.0006285338315286289]}",
    "{cvec: [-9.43268973692794e-06, 0.003828867161870621]}",
    "{dvec: [0.0022181782254955917, -0.009609931444148033]}",
  },
  // shifting 5-th param in + direction....
  {
    "{avec: [0.001327500295987637, -0.005996563846581542, -0.025076318272573946]}",
    "{bvec: [0.0005177666045874265, 0.00023488155209841835, 0.0008509946079322073]}",
    "{cvec: [1.4376101549438943e-06, 0.003329470782126069]}",
    "{dvec: [0.0021908139443399567, -0.008147735240384731]}",
  },
  // shifting 5-th param in - direction....
  {
    "{avec: [0.0013376745240123629, -0.0062014149534184585, -0.025052361727426056]}",
    "{bvec: [0.0005198970714125737, 7.42453079015817e-05, 0.0008199613920677929]}",
    "{cvec: [1.1094559845056106e-05, 0.0031872576178739315]}",
    "{dvec: [0.0022121551156600435, -0.00814408575961527]}",
  },
  // shifting 6-th param in + direction....
  {
    "{avec: [0.0013319207323129788, -0.006040683488732809, -0.025071829661523392]}",
    "{bvec: [0.0005183670734954005, 3.909704317488329e-05, 0.0008407401360212947]}",
    "{cvec: [-1.6555060704176955e-05, 0.0032994547146866014]}",
    "{dvec: [0.002192088448421568, -0.008142906670424709]}",
  },
  // shifting 6-th param in - direction....
  {
    "{avec: [0.0013332540876870211, -0.006157295311267191, -0.02505685033847661]}",
    "{bvec: [0.0005192966025045997, 0.00027002981682511674, 0.0008302158639787055]}",
    "{cvec: [2.9087230704176955e-05, 0.003217273685313399]}",
    "{dvec: [0.002210880611578432, -0.008148914329575293]}",
  },
  // shifting 7-th param in + direction....
  {
    "{avec: [0.0013522950438472476, -0.0060937603087213845, -0.02506493006423695]}",
    "{bvec: [0.0005260461199218985, 0.0001531369470211833, 0.0008361024288432554]}",
    "{cvec: [1.1212341603788612e-05, 0.003262183032498195]}",
    "{dvec: [0.0022554689492709874, -0.00814511416298112]}",
  },
  // shifting 7-th param in - direction....
  {
    "{avec: [0.0013128797761527524, -0.006104218491278616, -0.02506374993576305]}",
    "{bvec: [0.0005116175560781017, 0.00015598991297881675, 0.0008348535711567447]}",
    "{cvec: [1.3198283962113878e-06, 0.003254545367501805]}",
    "{dvec: [0.002147500110729013, -0.008146706837018881]}",
  },
  // shifting 8-th param in + direction....
  {
    "{avec: [0.0013181982077841504, -0.006100303569829522, -0.02506423265469129]}",
    "{bvec: [0.0005198266992863658, 0.0001561459019956034, 0.0008353656699881132]}",
    "{cvec: [-9.138585693877095e-06, 0.003257469951532203]}",
    "{dvec: [0.002208249015602, -0.00814578957470808]}",
  },
  // shifting 8-th param in - direction....
  {
    "{avec: [0.0013469766122158496, -0.0060976752301704785, -0.025064447345308712]}",
    "{bvec: [0.0005178369767136344, 0.00015298095800439664, 0.0008355903300118869]}",
    "{cvec: [2.1670755693877094e-05, 0.0032592584484677973]}",
    "{dvec: [0.0021947200443980003, -0.008146031425291921]}",
  },
  // shifting 9-th param in + direction....
  {
    "{avec: [0.0013213816393571078, -0.006098125355004741, -0.02506448307732861]}",
    "{bvec: [0.0005156160808193653, 0.0001529212683065104, 0.000835591898588393]}",
    "{cvec: [1.7663582719550113e-05, 0.003259342036434862]}",
    "{dvec: [0.002204761239756867, -0.00814588292849938]}",
  },
  // shifting 9-th param in - direction....
  {
    "{avec: [0.0013437931806428922, -0.006099853444995259, -0.02506419692267139]}",
    "{bvec: [0.0005220475951806349, 0.00015620559169348964, 0.0008353641014116071]}",
    "{cvec: [-5.131412719550114e-06, 0.003257386363565138]}",
    "{dvec: [0.002198207820243133, -0.008145938071500622]}",
  },
  // shifting 10-th param in + direction....
  {
    "{avec: [0.0013333256099200574, -0.006099027247184637, -0.025064332250233773]}",
    "{bvec: [0.0005138848807118044, 0.0001546573726872871, 0.0008354723294001544]}",
    "{cvec: [5.481032602412766e-06, 0.0032583117728193444]}",
    "{dvec: [0.002201957965286372, -0.008145903278427446]}",
  },
  // shifting 10-th param in - direction....
  {
    "{avec: [0.0013318492100799426, -0.006098951552815363, -0.02506434774976623]}",
    "{bvec: [0.0005237787952881958, 0.00015446948731271293, 0.0008354836705998458]}",
    "{cvec: [7.051137397587234e-06, 0.003258416627180656]}",
    "{dvec: [0.002201011094713628, -0.008145917721572556]}",
  }
};

vector<vector<string>> BtoD0starVars = {
  // shifting 1-th param in + direction....
  {
    "{ztp: 1.3199999999999998}",
    "{zeta1: 1.98}",
  },
  // shifting 1-th param in - direction....
  {
    "{ztp: -4.279999999999999}",
    "{zeta1: 1.98}",
  },
  // shifting 2-th param in + direction....
  {
    "{ztp: -1.48}",
    "{zeta1: 2.58}",
  },
  // shifting 2-th param in - direction....
  {
    "{ztp: -1.48}",
    "{zeta1: 1.38}",
  }
};

vector<vector<string>> BtoD1Vars = {
  // shifting 1-th param in + direction....
  {
    "{tp: -0.4}",
    "{tau1: 1.2999999999999998}",
    "{tau2: -0.7399999999999998}",
  },
  // shifting 1-th param in - direction....
  {
    "{tp: -1.2000000000000002}",
    "{tau1: 1.2999999999999998}",
    "{tau2: -0.7399999999999998}",
  },
  // shifting 2-th param in + direction....
  {
    "{tp: -0.8}",
    "{tau1: 1.9}",
    "{tau2: -0.7399999999999998}",
  },
  // shifting 2-th param in - direction....
  {
    "{tp: -0.8}",
    "{tau1: 0.6999999999999998}",
    "{tau2: -0.7399999999999998}",
  },
  // shifting 3-th param in + direction....
  {
    "{tp: -0.8}",
    "{tau1: 1.2999999999999998}",
    "{tau2: 2.06}",
  },
  // shifting 3-th param in - direction....
  {
    "{tp: -0.8}",
    "{tau1: 1.2999999999999998}",
    "{tau2: -3.5399999999999996}",
  }
};

vector<vector<string>> BtoD1starVars = {
  // shifting 1-th param in + direction....
  {
    "{ztp: 1.3199999999999998}",
    "{zeta1: 1.98}",
  },
  // shifting 1-th param in - direction....
  {
    "{ztp: -4.279999999999999}",
    "{zeta1: 1.98}",
  },
  // shifting 2-th param in + direction....
  {
    "{ztp: -1.48}",
    "{zeta1: 2.58}",
  },
  // shifting 2-th param in - direction....
  {
    "{ztp: -1.48}",
    "{zeta1: 1.38}",
  }
};

vector<vector<string>> BtoD2starVars = {
  // shifting 1-th param in + direction....
  {
    "{tp: -0.4}",
    "{tau1: 1.2999999999999998}",
    "{tau2: -0.7399999999999998}",
  },
  // shifting 1-th param in - direction....
  {
    "{tp: -1.2000000000000002}",
    "{tau1: 1.2999999999999998}",
    "{tau2: -0.7399999999999998}",
  },
  // shifting 2-th param in + direction....
  {
    "{tp: -0.8}",
    "{tau1: 1.9}",
    "{tau2: -0.7399999999999998}",
  },
  // shifting 2-th param in - direction....
  {
    "{tp: -0.8}",
    "{tau1: 0.6999999999999998}",
    "{tau2: -0.7399999999999998}",
  },
  // shifting 3-th param in + direction....
  {
    "{tp: -0.8}",
    "{tau1: 1.2999999999999998}",
    "{tau2: 2.06}",
  },
  // shifting 3-th param in - direction....
  {
    "{tp: -0.8}",
    "{tau1: 1.2999999999999998}",
    "{tau2: -3.5399999999999996}",
  }
};

// Various FF config
map<string, vector<vector<string>>> ffVarSpecs = {
  {"BD", BtoDVars},
  {"BD*", BtoDstVars},
  {"BD**0*", BtoD0starVars},
  {"BD**1", BtoD1Vars},
  {"BD**1*", BtoD1starVars},
  {"BD**2*", BtoD2starVars},
};

map<string, string> ffSchemeByDecay = {
  {"BD", "BGL"},
  {"BD*", "BGL"},
  {"BD**0*", "BLR"},
  {"BD**1", "BLR"},
  {"BD**1*", "BLR"},
  {"BD**2*", "BLR"},
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
  ham.setOptions(scheme + ": {Vcb: 0.0415}");
  ham.setOptions(scheme + ": {Chim: 0.0003068}");
  ham.setOptions(scheme + ": {Chip: 0.000528}");
  ham.setOptions(scheme + ": {ChimL: 0.002466}");
  ham.setOptions(scheme + ": {BcStatesf: [6.73, 6.736, 7.135, 7.142]}");
  ham.setOptions(scheme + ": {BcStatesg: [6.337, 6.899, 7.012, 7.28]}");
  ham.setOptions(scheme + ": {BcStatesP1: [6.275, 6.842, 7.25]}");
  ham.setOptions(scheme + ": {avec: [0.00133258741, -0.0060989894, -0.02506434]}");
  ham.setOptions(scheme + ": {bvec: [0.0005188318380000001, 0.00015456343000000002, 0.0008354780000000001]}");
  ham.setOptions(scheme + ": {cvec: [6.266085e-06, 0.0032583642]}");
  ham.setOptions(scheme + ": {dvec: [0.00220148453, -0.0081459105]}");
  // Alex note to self: gen_ham_params outputs setOptions(... abcderr: [...10 values...]) too,
  // but I don't think Hammer needs the errors (they weren't set in previous iterations of Hammer
  // reweighting, either)
}

void setBtoD0starBLRDefault(Hammer::Hammer& ham, const string scheme) {
  ham.setOptions(scheme + ": {as: 0.26}");
  ham.setOptions(scheme + ": {mb: 4.71}");
  ham.setOptions(scheme + ": {mc: 1.31}");
  ham.setOptions(scheme + ": {zt1: 0.7}");
  ham.setOptions(scheme + ": {ztp: -1.48}");
  ham.setOptions(scheme + ": {zeta1: 1.98}");
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
  ham.setOptions(scheme + ": {tp: -0.8}");
  ham.setOptions(scheme + ": {tau1: 1.2999999999999998}");
  ham.setOptions(scheme + ": {tau2: -0.7399999999999998}");
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
  ham.setOptions(scheme + ": {ztp: -1.48}");
  ham.setOptions(scheme + ": {zeta1: 1.98}");
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
  ham.setOptions(scheme + ": {tp: -0.8}");
  ham.setOptions(scheme + ": {tau1: 1.2999999999999998}");
  ham.setOptions(scheme + ": {tau2: -0.7399999999999998}");
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
  {"BD**2*", setBtoD2starBLRDefault}
};
// clang-format on

const int numOfFFVar = 20;

string decayDescr(const string decay) {
  auto daughter = decay.substr(1);
  return "Bto" + daughter;
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
                   double, double, double, double, double, double, double>{
          false,        wtFF,         wtFFVars[0],  wtFFVars[1],  wtFFVars[2],
          wtFFVars[3],  wtFFVars[4],  wtFFVars[5],  wtFFVars[6],  wtFFVars[7],
          wtFFVars[8],  wtFFVars[9],  wtFFVars[10], wtFFVars[11], wtFFVars[12],
          wtFFVars[13], wtFFVars[14], wtFFVars[15], wtFFVars[16], wtFFVars[17],
          wtFFVars[18], wtFFVars[19]};

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
                 double, double, double, double, double, double>{
        hamOk,        wtFF,         wtFFVars[0],  wtFFVars[1],  wtFFVars[2],
        wtFFVars[3],  wtFFVars[4],  wtFFVars[5],  wtFFVars[6],  wtFFVars[7],
        wtFFVars[8],  wtFFVars[9],  wtFFVars[10], wtFFVars[11], wtFFVars[12],
        wtFFVars[13], wtFFVars[14], wtFFVars[15], wtFFVars[16], wtFFVars[17],
        wtFFVars[18], wtFFVars[19]};
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
    df            = df.Define("wff", "get<1>(ff_result)");
    for (int i = 0; i < numOfFFVar; i++) {
      auto outputBrName = "wff_var" + to_string(i + 1);
      df = df.Define(outputBrName, "get<" + to_string(i + 2) + ">(ff_result)");
      outputBrs.emplace_back(outputBrName);
    }
    outputBrs.emplace_back("ham_ok");
    outputBrs.emplace_back("wff");

    df.Snapshot(trees[idx], ntpOut, outputBrs, writeOpts);

    cout << "Total number of candidates: " << numOfEvt << endl;
    cout << "Hammer reweighted candidates: " << numOfEvtOk << endl;
    cout << "Reweighted fraction: "
         << static_cast<float>(numOfEvtOk) / static_cast<float>(numOfEvt)
         << endl;
  }
}
