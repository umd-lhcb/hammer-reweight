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
  // Note: these variations are all manually set to 10sigma instead of 1sigma
  // shifting 1-th param in + direction....
  {
    "{avec: [0.0012465300969462544, -0.006078741080525021, -0.027022644719643378]}",
    "{bvec: [0.00048156025250519207, 1.9109299873591896e-06, 0.00506796892822629]}",
    "{cvec: [3.355953586166967e-06, 0.002083037477675795, -0.02204522052306778]}",
    "{dvec: [0.0020582655123657667, -0.014622323576421525, 0.3470464727296727]}",
  },
  // shifting 1-th param in - direction....
  {
    "{avec: [0.001235020751053745, -0.005285370279474977, -0.02168068968035662]}",
    "{bvec: [0.00048400267869480785, 0.00016043351801264078, 0.0003435126717737098]}",
    "{cvec: [4.4574312013833035e-05, 0.002555375322324204, -0.05062324667693221]}",
    "{dvec: [0.002046729815634233, -0.000916359303578476, -0.34699235791367267]}",
  },
  // shifting 2-th param in + direction....
  {
    "{avec: [0.001164237594974296, -0.0026653850340459186, -0.04633165839455134]}",
    "{bvec: [0.00048505335943473967, -1.735253300333417e-05, 0.002325477279678842]}",
    "{cvec: [-2.1779232994205136e-05, -0.005533233191349908, 0.1484452829218463]}",
    "{dvec: [0.0020339935574234572, -0.00982515574159766, -0.007790606662095515]}",
  },
  // shifting 2-th param in - direction....
  {
    "{avec: [0.0013173132530257035, -0.00869872632595408, -0.0023716760054486585]}",
    "{bvec: [0.00048050957176526025, 0.00017969698100333414, 0.0030860043203211583]}",
    "{cvec: [6.970949859420514e-05, 0.010171645991349907, -0.22111375012184628]}",
    "{dvec: [0.0020710017705765425, -0.005713527138402339, 0.007844721478095515]}",
  },
  // shifting 3-th param in + direction....
  {
    "{avec: [0.0014937800485939367, -0.014790260007525405, -0.09839611436453814]}",
    "{bvec: [0.0004692282904880357, 0.00041544409117401447, 0.004262515284802589]}",
    "{cvec: [9.011383666413758e-05, 0.0019066083870788136, -0.0450451444754436]}",
    "{dvec: [0.0020714194151520552, -0.01010790006508856, -0.00025156835522024593]}",
  },
  // shifting 3-th param in - direction....
  {
    "{avec: [0.0009877707994060628, 0.003426148647525407, 0.04969277996453815]}",
    "{bvec: [0.0004963346407119641, -0.0002530996431740145, 0.001148966315197411]}",
    "{cvec: [-4.218357106413758e-05, 0.0027318044129211855, -0.027623322724556386]}",
    "{dvec: [0.0020335759128479445, -0.005430782814911439, 0.0003056831712202459]}",
  },
  // shifting 4-th param in + direction....
  {
    "{avec: [0.001369170146042264, -0.01230600055107656, -0.022543062861824722]}",
    "{bvec: [0.0004783648808389333, -0.0008636463452013094, 0.04508860307725271]}",
    "{cvec: [-6.0607917454812375e-05, 0.0038488052109782065, -0.03593271723285831]}",
    "{dvec: [0.002101459341855659, -0.012862188209125897, -0.00037139218316628895]}",
  },
  // shifting 4-th param in - direction....
  {
    "{avec: [0.0011123807019577355, 0.0009418891910765621, -0.026160271538175277]}",
    "{bvec: [0.0004871980503610666, 0.0010259907932013094, -0.039677121477252705]}",
    "{cvec: [0.00010853818305481237, 0.0007896075890217924, -0.03673574996714168]}",
    "{dvec: [0.0020035359861443405, -0.0026764946708741034, 0.0004255069991662889]}",
  },
  // shifting 5-th param in + direction....
  {
    "{avec: [0.0012196205558262715, -0.0052523382682886114, -0.024793214060702617]}",
    "{bvec: [0.0004793723178580663, 0.00044236264531985784, 0.0044529584647468486]}",
    "{cvec: [9.99709987125693e-05, 0.0013156807256480539, -0.036272666561984045]}",
    "{dvec: [0.0018314591531018646, 0.005668526224921964, 0.00027450204937167936]}",
  },
  // shifting 5-th param in - direction....
  {
    "{avec: [0.001261930292173728, -0.006111773091711387, -0.023910120339297382]}",
    "{bvec: [0.0004861906133419336, -0.00028001819731985787, 0.0009585231352531522]}",
    "{cvec: [-5.204073311256931e-05, 0.0033227320743519453, -0.036395800638015945]}",
    "{dvec: [0.002273536174898135, -0.021207209104921965, -0.00022038723337167938]}",
  },
  // shifting 6-th param in + direction....
  {
    "{avec: [0.0012814138164973018, -0.008267343026409062, -0.024139009449046952]}",
    "{bvec: [0.0004794221071018516, -0.0003744019336643319, 0.002138039907403325]}",
    "{cvec: [-0.0002497820085051978, 0.008283293510078829, -0.036007641895928866]}",
    "{dvec: [0.001932584608310705, -0.0071501850939849135, 3.232260675255132e-05]}",
  },
  // shifting 6-th param in - direction....
  {
    "{avec: [0.0012001370315026978, -0.0030967683335909354, -0.024564324950953047]}",
    "{bvec: [0.0004861408240981483, 0.0005367463816643319, 0.0032734416925966755]}",
    "{cvec: [0.00029771227410519774, -0.003644880710078829, -0.036660825304071123]}",
    "{dvec: [0.0021724107196892948, -0.008388497786015086, 2.1792209247448674e-05]}",
  },
  // shifting 7-th param in + direction....
  {
    "{avec: [0.0012843702142687969, -0.006412653041316254, -0.02426888037804516]}",
    "{bvec: [0.0004912622442030115, -0.0010792733850466924, 0.0025823443016507385]}",
    "{cvec: [-4.289238173878899e-05, 0.001892036455619719, -0.03633097247280449]}",
    "{dvec: [0.0020932378153919453, -0.007726833683069959, 2.784348983304206e-05]}",
  },
  // shifting 7-th param in - direction....
  {
    "{avec: [0.0011971806337312027, -0.004951458318683744, -0.02443445402195484]}",
    "{bvec: [0.0004743006869969884, 0.0012416178330466924, 0.0028291372983492618]}",
    "{cvec: [9.0822647338789e-05, 0.00274637634438028, -0.036337494727195496]}",
    "{dvec: [0.0020117575126080544, -0.007811849196930041, 2.6271326166957937e-05]}",
  },
  // shifting 8-th param in + direction....
  {
    "{avec: [0.0014180123942665649, -0.005679643681207845, -0.02435154108227566]}",
    "{bvec: [0.0005499882615659136, 9.74964954688785e-05, 0.002706034825926937]}",
    "{cvec: [5.749844471477067e-05, 0.0023304544770001152, -0.03633354666453006]}",
    "{dvec: [0.0025221654709425127, -0.0077612251640925755, 2.7193996188637165e-05]}",
  },
  // shifting 8-th param in - direction....
  {
    "{avec: [0.0010635384537334347, -0.005684467678792153, -0.02435179331772434]}",
    "{bvec: [0.00041557466963408623, 6.484795253112148e-05, 0.0027054467740730633]}",
    "{cvec: [-9.56817911477067e-06, 0.002307958322999884, -0.03633492053546993]}",
    "{dvec: [0.001582829857057487, -0.007777457715907424, 2.692081981136283e-05]}",
  },
  // shifting 9-th param in + direction....
  {
    "{avec: [0.0011849957135970663, -0.005974991668828165, -0.024315485343657205]}",
    "{bvec: [0.0005056941533157017, 0.0003198652591399119, 0.0026680212185089757]}",
    "{cvec: [-7.230744495556924e-05, 0.0022016338683811198, -0.03633012659798523]}",
    "{dvec: [0.0020731556917226524, -0.007768293450027011, 2.7078831778645917e-05]}",
  },
  // shifting 9-th param in - direction....
  {
    "{avec: [0.0012965551344029332, -0.005389119691171833, -0.024387849056342795]}",
    "{bvec: [0.00045986877788429815, -0.00015752081113991193, 0.0027434603814910246]}",
    "{cvec: [0.00012023771055556925, 0.0024367789316188794, -0.03633834060201476]}",
    "{dvec: [0.0020318396362773473, -0.007770389429972989, 2.703598422135408e-05]}",
  },
  // shifting 10-th param in + direction....
  {
    "{avec: [0.0010843791398429444, -0.005667196887414759, -0.024354078616807308]}",
    "{bvec: [0.00046690408606668424, 6.459274908984185e-05, 0.0027080342611609798]}",
    "{cvec: [2.6334724289269157e-05, 0.00232700335213522, -0.036334482774126486]}",
    "{dvec: [0.0021139160274426254, -0.007768418630867613, 2.7071881786327663e-05]}",
  },
  // shifting 10-th param in - direction....
  {
    "{avec: [0.0013971717081570551, -0.005696914472585239, -0.02434925578319269]}",
    "{bvec: [0.0004986588451333157, 9.775169891015812e-05, 0.0027034473388390205]}",
    "{cvec: [2.159554131073084e-05, 0.002311409447864779, -0.0363339844258735]}",
    "{dvec: [0.0019910793005573743, -0.007770264249132386, 2.7042934213672333e-05]}",
  },
  // shifting 11-th param in + direction....
  {
    "{avec: [0.0012346853058009702, -0.005685677566489024, -0.024351239731133584]}",
    "{bvec: [0.000521997411348442, 8.206954100809469e-05, 0.0027052592996004163]}",
    "{cvec: [4.981979169892542e-05, 0.0023187949863201483, -0.0363341398420916]}",
    "{dvec: [0.0020473372545730516, -0.007769434970123697, 2.7055913157374695e-05]}",
  },
  // shifting 11-th param in - direction....
  {
    "{avec: [0.0012468655421990294, -0.005678433793510974, -0.024352094668866415]}",
    "{bvec: [0.0004435655198515579, 8.027490699190528e-05, 0.002706222300399584]}",
    "{cvec: [-1.889526098925419e-06, 0.002319617813679851, -0.03633432735790839]}",
    "{dvec: [0.002057658073426948, -0.007769247909876303, 2.70589028426253e-05]}",
  },
  // shifting 12-th param in + direction....
  {
    "{avec: [0.0012414250753808162, -0.005685699220444434, -0.024351209123777216]}",
    "{bvec: [0.0004730355493736248, 8.299681595606507e-05, 0.002705248325108172]}",
    "{cvec: [3.83162945171957e-05, 0.002318349403817377, -0.03633415212052943]}",
    "{dvec: [0.002052598312588332, -0.007769340223737234, 2.705742918688994e-05]}",
  },
  // shifting 12-th param in - direction....
  {
    "{avec: [0.0012401257726191834, -0.005678412139555564, -0.024352125276222783]}",
    "{bvec: [0.0004925273818263752, 7.93476320439349e-05, 0.002706233274891828]}",
    "{cvec: [9.613971082804295e-06, 0.002320063396182622, -0.03633431507947056]}",
    "{dvec: [0.0020523970154116677, -0.007769342656262766, 2.7057386813110055e-05]}",
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
    df            = df.Define("wff_dst10sig", "get<1>(ff_result)");
    for (int i = 0; i < numOfFFVar; i++) {
      auto outputBrName = "wff_dst10sig_var" + to_string(i + 1);
      df = df.Define(outputBrName, "get<" + to_string(i + 2) + ">(ff_result)");
      outputBrs.emplace_back(outputBrName);
    }
    outputBrs.emplace_back("wff_dst10sig");

    df.Snapshot(trees[idx], ntpOut, outputBrs, writeOpts);

    cout << "Total number of candidates: " << numOfEvt << endl;
    cout << "Hammer reweighted candidates: " << numOfEvtOk << endl;
    cout << "Reweighted fraction: "
         << static_cast<float>(numOfEvtOk) / static_cast<float>(numOfEvt)
         << endl;
  }
}
