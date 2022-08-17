// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Wed Aug 17, 2022 at 05:56 AM -0400

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
    "{ap: [0.015642660612597052, -0.034035685234997386, -0.0898070271651518, 0.0]}",
    "{a0: [0.0858147593040351, -0.3050702334404309, -0.22752112980378628, 0.0]}"
  },
  //shifting 1-th param in - direction....
  {
    "{ap: [0.01567733938740295, -0.03436431476500262, -0.09019297283484819, 0.0]}",
    "{a0: [0.07288524069596491, -0.10492976655956911, -0.23247887019621374, 0.0]}"
  },
  // shifting 2-th param in + direction....
  {
    "{ap: [0.015536159303514138, -0.03274798737619622, -0.08610235496526895, 0.0]}",
    "{a0: [0.07939786111920759, -0.2055198293137311, -0.2513856759832376, 0.0]}"
  },
  // shifting 2-th param in - direction....
  {
    "{ap: [0.015783840696485862, -0.03565201262380378, -0.09389764503473104, 0.0]}",
    "{a0: [0.07930213888079242, -0.20448017068626886, -0.20861432401676247, 0.0]}"
  },
  // shifting 3-th param in + direction....
  {
    "{ap: [0.015185466182351984, -0.04806686174041869, -0.08969198838997502, 0.0]}",
    "{a0: [0.07255753045478437, -0.20504393049104058, -0.23088155901344717, 0.0]}"
  },
  // shifting 3-th param in - direction....
  {
    "{ap: [0.016134533817648015, -0.020333138259581315, -0.09030801161002497, 0.0]}",
    "{a0: [0.08614246954521564, -0.2049560695089594, -0.22911844098655285, 0.0]}"
  },
  // shifting 4-th param in + direction....
  {
    "{ap: [0.01566286514976162, -0.034199101149958354, -0.08990582046864679, 0.0]}",
    "{a0: [0.07936645416843903, -0.20499939148119928, -0.22998280566971113, 0.0]}"
  },
  // shifting 4-th param in - direction....
  {
    "{ap: [0.01565713485023838, -0.03420089885004165, -0.09009417953135321, 0.0]}",
    "{a0: [0.07933354583156098, -0.2050006085188007, -0.2300171943302889, 0.0]}"
  },
  //shifting 5-th param in + direction....
  {
    "{ap: [0.014646019937572174, -0.034165525575076655, -0.08997191655772308, 0.0]}",
    "{a0: [0.07430176773341135, -0.2049993836925806, -0.22998668416129328, 0.0]}"
  },
  // shifting 5-th param in - direction....
  {
    "{ap: [0.016673980062427828, -0.03423447442492335, -0.09002808344227692, 0.0]}",
    "{a0: [0.08439823226658866, -0.20500061630741936, -0.23001331583870674, 0.0]}"
  }
};

vector<vector<string>> BtoDstVars = {
  // shifting 1-th param in + direction....
  {
    "{avec: [0.001309359578796952, -0.006085426132949664, -0.025170335831386732]}",
    "{bvec: [0.00046842382043563353, -0.00016195225293225718, 0.000790434370119119]}",
    "{cvec: [0.003152772380981239, -0.0027771913155404675]}",
    "{dvec: [-0.0040444465715125735, -0.0067898966427725775]}"
  },
  // shifting 1-th param in - direction....
  {
    "{avec: [0.001355815241203048, -0.006112552667050336, -0.02495834416861327]}",
    "{bvec: [0.0005692398555643667, 0.0004710791129322572, 0.0008805216298808812]}",
    "{cvec: [-0.003140240210981239, 0.009293919715540468]}",
    "{dvec: [0.008447415631512573, -0.009501924357227423]}"
  },
  // shifting 2-th param in + direction....
  {
    "{avec: [0.001754279448766413, -0.0065216327936004435, -0.02246589784587869]}",
    "{bvec: [0.0016167032260642866, 0.0035135110015019643, 0.0029118310116814116]}",
    "{cvec: [0.00044630540683350854, 0.0029678260484947024]}",
    "{dvec: [0.0024553262270693298, -0.008182553823460335]}"
  },
  // shifting 2-th param in - direction....
  {
    "{avec: [0.000910895371233587, -0.005676346006399557, -0.027662782154121313]}",
    "{bvec: [-0.0005790395500642864, -0.0032043841415019645, -0.0012408750116814113]}",
    "{cvec: [-0.00043377323683350857, 0.003548902351505298]}",
    "{dvec: [0.0019476428329306705, -0.008109267176539666]}"
  },
  // shifting 3-th param in + direction....
  {
    "{avec: [0.0028996813693635847, -0.005821918831159753, -0.0250185679052699]}",
    "{bvec: [0.00048139245157527037, 2.6927621190305163e-05, 0.0007488148973385366]}",
    "{cvec: [-9.868135740501089e-06, 0.00326607533281066]}",
    "{dvec: [0.0021878267263853394, -0.008143472502158787]}"
  },
  // shifting 3-th param in - direction....
  {
    "{avec: [-0.00023450654936358478, -0.0063760599688402475, -0.025110112094730103]}",
    "{bvec: [0.0005562712244247298, 0.0002821992388096949, 0.0009221411026614636]}",
    "{cvec: [2.240030574050109e-05, 0.0032506530671893405]}",
    "{dvec: [0.002215142333614661, -0.008148348497841214]}"
  },
  // shifting 4-th param in + direction....
  {
    "{avec: [0.001332750419363132, -0.006098955947371931, -0.02506284486010554]}",
    "{bvec: [0.0005190844466766718, 0.0001430243978124457, 0.0008300109725351035]}",
    "{cvec: [0.0002068001106268701, 0.0033310155864089584]}",
    "{dvec: [0.0018894234909381357, -0.009727993036208988]}",
  },
  // shifting 4-th param in - direction....
  {
    "{avec: [0.0013324244006368679, -0.006099022852628069, -0.025065835139894462]}",
    "{bvec: [0.0005185792293233284, 0.00016610246218755435, 0.0008409450274648966]}",
    "{cvec: [-0.00019426794062687007, 0.003185712813591042]}",
    "{dvec: [0.0025135455690618643, -0.006563827963791013]}"
  },
  // shifting 5-th param in + direction....
  {
    "{avec: [0.0013327717846240666, -0.006102809647810181, -0.025066442065359646]}",
    "{bvec: [0.0005214526744740508, 0.00023488155209841835, 0.0007100919344854683]}",
    "{cvec: [2.605614158483439e-06, 0.003269032440650809]}",
    "{dvec: [0.0021867733464291986, -0.008143136992461518]}"
  },
  // shifting 5-th param in - direction....
  {
    "{avec: [0.0013324030353759333, -0.006095169152189819, -0.025062237934640356]}",
    "{bvec: [0.0005162110015259494, 7.424530790158168e-05, 0.0009608640655145318]}",
    "{cvec: [9.926555841516562e-06, 0.0032476959593491913]}",
    "{dvec: [0.0022161957135708016, -0.008148684007538484]}"
  },
  // shifting 6-th param in + direction....
  {
    "{avec: [0.0013108660484635445, -0.005964527795849014, -0.025069475335398918]}",
    "{bvec: [0.0005364246349187906, 0.0001688524711814243, 0.0008407401360212947]}",
    "{cvec: [7.741654312295005e-06, 0.0032576668370395798]}",
    "{dvec: [0.002202424158845132, -0.008146064670639655]}"
  },
  // shifting 6-th param in - direction....
  {
    "{avec: [0.0013543087715364554, -0.006233451004150986, -0.025059204664601084]}",
    "{bvec: [0.0005012390410812096, 0.00014027438881857574, 0.0008302158639787055]}",
    "{cvec: [4.790515687704995e-06, 0.0032590615629604206]}",
    "{dvec: [0.002200544901154868, -0.008145756329360347]}"
  },
  // shifting 7-th param in + direction....
  {
    "{avec: [0.001332515735692203, -0.0060990289765959944, -0.025063801220319485]}",
    "{bvec: [0.0005193966069284976, 0.00015268177666130997, 0.0008258205874372128]}",
    "{cvec: [1.1212341603788608e-05, 0.0032178936756208686]}",
    "{dvec: [0.0022412742473197095, -0.00815494272313176]}"
  },
  // shifting 7-th param in - direction....
  {
    "{avec: [0.001332659084307797, -0.006098949823404006, -0.025064878779680517]}",
    "{bvec: [0.0005182670690715026, 0.00015644508333869007, 0.0008451354125627873]}",
    "{cvec: [1.319828396211392e-06, 0.003298834724379132]}",
    "{dvec: [0.0021616948126802907, -0.008136878276868242]}"
  },
  // shifting 8-th param in + direction....
  {
    "{avec: [0.0013342705139521996, -0.006099735113036324, -0.025080833537241616]}",
    "{bvec: [0.0005110195808831549, 0.0001651110149861937, 0.0008420967885866525]}",
    "{cvec: [7.719682599049204e-06, 0.003257469951532203]}",
    "{dvec: [0.0022027839248535833, -0.008146140096944605]}",
  },
  // shifting 8-th param in - direction....
  {
    "{avec: [0.0013309043060478003, -0.006098243686963677, -0.025047846462758386]}",
    "{bvec: [0.0005266440951168453, 0.00014401584501380635, 0.0008288592114133476]}",
    "{cvec: [4.812487400950795e-06, 0.0032592584484677973]}",
    "{dvec: [0.002200185135146417, -0.008145680903055397]}"
  },
  // shifting 9-th param in + direction....
  {
    "{avec: [0.0013325722709559543, -0.006098955604264764, -0.02506415712597228]}",
    "{bvec: [0.000518659811170497, 0.0001533723071420947, 0.0008343390391924923]}",
    "{cvec: [2.1729559912085935e-05, 0.003263454693060284]}",
    "{dvec: [0.002204761239756867, -0.008144350246025094]}"
  },
  // shifting 9-th param in - direction....
  {
    "{avec: [0.0013326025490440457, -0.0060990231957352365, -0.02506452287402772]}",
    "{bvec: [0.0005190038648295032, 0.00015575455285790534, 0.0008366169608075079]}",
    "{cvec: [-9.197389912085936e-06, 0.0032532737069397165]}",
    "{dvec: [0.002198207820243133, -0.008147470753974907]}"
  },
  // shifting 10-th param in + direction....
  {
    "{avec: [0.001332871541573003, -0.006099621108185075, -0.0250664407692947]}",
    "{bvec: [0.0005234096327192811, 0.00015450162342007755, 0.0008355884848210036]}",
    "{cvec: [6.335300106495867e-06, 0.003258391812615387]}",
    "{dvec: [0.0022014928961516023, -0.008145903278427446]}"
  },
  // shifting 10-th param in - direction....
  {
    "{avec: [0.001332303278426997, -0.006098357691814925, -0.025062239230705302]}",
    "{bvec: [0.0005142540432807191, 0.00015462523657992249, 0.0008353675151789965]}",
    "{cvec: [6.196869893504133e-06, 0.0032583365873846133]}",
    "{dvec: [0.002201476163848398, -0.008145917721572556]}"
  }
};

vector<vector<string>> BtoD0starVars = {
  // shifting 1-th param in + direction....
  {
    "{zt1: 1.098150778538033}",
    "{ztp: -2.5782005874353002}",
    "{zeta1: 1.1090232758998373}",
  },
  // shifting 1-th param in - direction....
  {
    "{zt1: 0.3018492214619671}",
    "{ztp: 2.9782005874353006}",
    "{zeta1: 0.09097672410016278}",
  },
  // shifting 2-th param in + direction....
  {
    "{zt1: 0.5963724470564943}",
    "{ztp: 0.18310737440727248}",
    "{zeta1: 0.5888576579281841}",
  },
  // shifting 2-th param in - direction....
  {
    "{zt1: 0.8036275529435056}",
    "{ztp: 0.21689262559272754}",
    "{zeta1: 0.6111423420718158}",
  },
  // shifting 3-th param in + direction....
  {
    "{zt1: 0.637119328388831}",
    "{ztp: 0.27680365896258735}",
    "{zeta1: 1.0683714283645744}",
  },
  // shifting 3-th param in - direction....
  {
    "{zt1: 0.7628806716111689}",
    "{ztp: 0.12319634103741266}",
    "{zeta1: 0.13162857163542552}",
  }
};

vector<vector<string>> BtoD1Vars = {
  // shifting 1-th param in + direction....
  {
    "{t1: 0.7694432788984408}",
    "{tp: -2.372337805334135}",
    "{tau1: -3.2416986114558015}",
    "{tau2: 2.7747263953487225}",
  },
  // shifting 1-th param in - direction....
  {
    "{t1: 0.6305567211015591}",
    "{tp: -0.8276621946658654}",
    "{tau1: 2.2416986114558015}",
    "{tau2: 3.0252736046512774}",
  },
  // shifting 2-th param in + direction....
  {
    "{t1: 0.6945480372695823}",
    "{tp: -1.2056731495789785}",
    "{tau1: -0.6152812019868762}",
    "{tau2: 2.988878858706443}",
  },
  // shifting 2-th param in - direction....
  {
    "{t1: 0.7054519627304177}",
    "{tp: -1.9943268504210216}",
    "{tau1: -0.38471879801312386}",
    "{tau2: 2.811121141293557}",
  },
  // shifting 3-th param in + direction....
  {
    "{t1: 0.708079553210943}",
    "{tp: -1.6088894147532562}",
    "{tau1: -0.49916538047982556}",
    "{tau2: 2.9410176316310555}",
  },
  // shifting 3-th param in - direction....
  {
    "{t1: 0.691920446789057}",
    "{tp: -1.591110585246744}",
    "{tau1: -0.5008346195201745}",
    "{tau2: 2.8589823683689444}",
  },
  // shifting 4-th param in + direction....
  {
    "{t1: 0.44568095475294084}",
    "{tp: -1.6153253350947394}",
    "{tau1: -0.5042655344258667}",
    "{tau2: 2.946860606539276}",
  },
  // shifting 4-th param in - direction....
  {
    "{t1: 0.9543190452470591}",
    "{tp: -1.5846746649052608}",
    "{tau1: -0.49573446557413325}",
    "{tau2: 2.8531393934607236}",
  }
};

vector<vector<string>> BtoD1starVars = {
    // shifting 1-th param in + direction....
  {
    "{zt1: 1.098150778538033}",
    "{ztp: -2.5782005874353002}",
    "{zeta1: 1.1090232758998373}",
  },
  // shifting 1-th param in - direction....
  {
    "{zt1: 0.3018492214619671}",
    "{ztp: 2.9782005874353006}",
    "{zeta1: 0.09097672410016278}",
  },
  // shifting 2-th param in + direction....
  {
    "{zt1: 0.5963724470564943}",
    "{ztp: 0.18310737440727248}",
    "{zeta1: 0.5888576579281841}",
  },
  // shifting 2-th param in - direction....
  {
    "{zt1: 0.8036275529435056}",
    "{ztp: 0.21689262559272754}",
    "{zeta1: 0.6111423420718158}",
  },
  // shifting 3-th param in + direction....
  {
    "{zt1: 0.637119328388831}",
    "{ztp: 0.27680365896258735}",
    "{zeta1: 1.0683714283645744}",
  },
  // shifting 3-th param in - direction....
  {
    "{zt1: 0.7628806716111689}",
    "{ztp: 0.12319634103741266}",
    "{zeta1: 0.13162857163542552}",
  }
};

vector<vector<string>> BtoD2starVars = {
    // shifting 1-th param in + direction....
  {
    "{t1: 0.7694432788984408}",
    "{tp: -2.372337805334135}",
    "{tau1: -3.2416986114558015}",
    "{tau2: 2.7747263953487225}",
  },
  // shifting 1-th param in - direction....
  {
    "{t1: 0.6305567211015591}",
    "{tp: -0.8276621946658654}",
    "{tau1: 2.2416986114558015}",
    "{tau2: 3.0252736046512774}",
  },
  // shifting 2-th param in + direction....
  {
    "{t1: 0.6945480372695823}",
    "{tp: -1.2056731495789785}",
    "{tau1: -0.6152812019868762}",
    "{tau2: 2.988878858706443}",
  },
  // shifting 2-th param in - direction....
  {
    "{t1: 0.7054519627304177}",
    "{tp: -1.9943268504210216}",
    "{tau1: -0.38471879801312386}",
    "{tau2: 2.811121141293557}",
  },
  // shifting 3-th param in + direction....
  {
    "{t1: 0.708079553210943}",
    "{tp: -1.6088894147532562}",
    "{tau1: -0.49916538047982556}",
    "{tau2: 2.9410176316310555}",
  },
  // shifting 3-th param in - direction....
  {
    "{t1: 0.691920446789057}",
    "{tp: -1.591110585246744}",
    "{tau1: -0.5008346195201745}",
    "{tau2: 2.8589823683689444}",
  },
  // shifting 4-th param in + direction....
  {
    "{t1: 0.44568095475294084}",
    "{tp: -1.6153253350947394}",
    "{tau1: -0.5042655344258667}",
    "{tau2: 2.946860606539276}",
  },
  // shifting 4-th param in - direction....
  {
    "{t1: 0.9543190452470591}",
    "{tp: -1.5846746649052608}",
    "{tau1: -0.49573446557413325}",
    "{tau2: 2.8531393934607236}",
  }
};

// Various FF config
map<string, vector<vector<string>>> ffVarSpecs = {
  {"BD", BtoDVars},
  {"BD*", BtoDstVars},
  {"BD**0", BtoD0starVars},
  {"BD**1", BtoD1Vars},
  {"BD**1*", BtoD1starVars},
  {"BD**2*", BtoD2starVars},
};

map<string, string> ffSchemeByDecay = {
  {"BD", "BGL"},
  {"BD*", "BGL"},
  {"BD**0", "BLR"},
  {"BD**1", "BLR"},
  {"BD**1*", "BLR"},
  {"BD**2", "BLR"},
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
  ham.setOptions(scheme + ": {laS: 0.8}");
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
  ham.setOptions(scheme + ": {laS: 0.8}");
};

map<string, function<void(Hammer::Hammer&, const string)>>
ffSchemeDefaultsByDecay = {
  {"BD", setBtoDBGLDefault},
  {"BD*", setBtoDstarBGLDefault},
  {"BD**0", setBtoD0starBLRDefault},
  {"BD**1", setBtoD1BLRDefault},
  {"BD**1*", setBtoD1starBLRDefault},
  {"BD**2", setBtoD2starBLRDefault}
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
