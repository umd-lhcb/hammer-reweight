// Author: Yipeng Sun, Alex Fernez
// License: BSD 2-clause

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
  // shifting in 1-th direction (+)...
  {
    "{ap: [0.012351197624241281, -0.027367962562844683, -0.07277729207579971, 0.0]}",
    "{a0: [0.062484046787995004, -0.1572277107072667, -0.2805601501901888, 0.0]}",
  },
  // shifting in 1-th direction (-)...
  {
    "{ap: [0.01896880237575872, -0.04103203743715532, -0.10722270792420029, 0.0]}",
    "{a0: [0.096215953212005, -0.25277228929273327, -0.1794398498098112, 0.0]}",
  },
  // shifting in 2-th direction (+)...
  {
    "{ap: [0.01591678690539162, -0.036804793342622426, -0.07049365753599796, 0.0]}",
    "{a0: [0.08068819570228931, -0.2123544213516238, -0.23067308360617908, 0.0]}",
  },
  // shifting in 2-th direction (-)...
  {
    "{ap: [0.01540321309460838, -0.03159520665737758, -0.10950634246400204, 0.0]}",
    "{a0: [0.0780118042977107, -0.19764557864837617, -0.22932691639382094, 0.0]}",
  },
  // shifting in 3-th direction (+)...
  {
    "{ap: [0.016629800510417916, -0.03439666118284376, -0.09000068878903926, 0.0]}",
    "{a0: [0.08411975918337586, -0.20489880083482004, -0.2299946563017477, 0.0]}",
  },
  // shifting in 3-th direction (-)...
  {
    "{ap: [0.014690199489582086, -0.03400333881715624, -0.08999931121096073, 0.0]}",
    "{a0: [0.07458024081662415, -0.20510119916517994, -0.23000534369825232, 0.0]}",
  },
  // shifting in 4-th direction (+)...
  {
    "{ap: [0.01587525531571869, -0.03322237162082708, -0.08993374264698996, 0.0]}",
    "{a0: [0.08074892773777949, -0.20516187356690724, -0.23001236118642798, 0.0]}",
  },
  // shifting in 4-th direction (-)...
  {
    "{ap: [0.015444744684281312, -0.035177628379172925, -0.09006625735301003, 0.0]}",
    "{a0: [0.07795107226222052, -0.20483812643309274, -0.22998763881357204, 0.0]}",
  }
};

vector<vector<string>> BtoDstVars = {
  // shifting in 1-th direction (+)...
  {
    "{avec: [0.001260431099035367, -0.005814457482657006, -0.024987372258585112]}",
    "{bvec: [0.0004900145218689943, 7.466869133956658e-05, 0.0029848016188998186]}",
    "{cvec: [2.2336792682128603e-05, 0.002339522831233444, -0.03566541618423161]}",
    "{dvec: [0.0020844011242765284, -0.008571735209323101, 0.03473567640725931]}",
  },
  // shifting in 1-th direction (-)...
  {
    "{avec: [0.0012211197489646326, -0.005549653877342992, -0.023715962141414887]}",
    "{bvec: [0.0004755484093310056, 8.76757566604334e-05, 0.0024266799811001817]}",
    "{cvec: [2.5593472917871395e-05, 0.002298889968766555, -0.03700305101576838]}",
    "{dvec: [0.0020205942037234714, -0.006966947670676899, -0.03468156159125931]}",
  },
  // shifting in 2-th direction (+)...
  {
    "{avec: [0.0015387628369786265, -0.007131957792233949, -0.035101031616549414]}",
    "{bvec: [0.0005970822644012318, 0.00010237454067166451, 0.0034572861427716907]}",
    "{cvec: [2.780849764268198e-05, 0.0021380081108407323, -0.02835811244570046]}",
    "{dvec: [0.002539299800309922, -0.009928456447903296, -0.0003855357904247226]}",
  },
  // shifting in 2-th direction (-)...
  {
    "{avec: [0.000942788011021373, -0.004232153567766049, -0.01360230278345059]}",
    "{bvec: [0.0003684806667987682, 5.9969907328335465e-05, 0.0019541954572283096]}",
    "{cvec: [2.0121767957318017e-05, 0.002500404689159267, -0.04431035475429953]}",
    "{dvec: [0.0015656955276900776, -0.005610226432096703, 0.00043965060642472254]}",
  },
  // shifting in 3-th direction (+)...
  {
    "{avec: [0.0011944443127184525, -0.005038648081395574, -0.024998540920203566]}",
    "{bvec: [0.0004678529136803674, 0.00017872786959570949, -0.0015031853253623566]}",
    "{cvec: [3.315476147439584e-05, 0.002133188698409768, -0.03648934835435856]}",
    "{dvec: [0.001983601899842271, -0.007033613662738017, 7.172198921323677e-05]}",
  },
  // shifting in 3-th direction (-)...
  {
    "{avec: [0.001287106535281547, -0.006325463278604424, -0.023704793479796434]}",
    "{bvec: [0.0004977100175196325, -1.638342159570951e-05, 0.006914666925362357]}",
    "{cvec: [1.477550412560416e-05, 0.0025052241015902314, -0.03617911884564143]}",
    "{dvec: [0.0021213934281577285, -0.008505069217261983, -1.7607173213236773e-05]}",
  },
  // shifting in 4-th direction (+)...
  {
    "{avec: [0.0010525834985406642, -0.005539398065100158, -0.024839002422353804]}",
    "{bvec: [0.00040007024225617046, 0.00010215297346911433, 0.00316670932085098]}",
    "{cvec: [2.9873722179234278e-05, 0.002089597326599368, -0.036454496301279]}",
    "{dvec: [0.001695750537909285, -0.0058208626148416615, 6.292224991778452e-05]}",
  },
  // shifting in 4-th direction (-)...
  {
    "{avec: [0.0014289673494593353, -0.00582471329489984, -0.023864331977646196]}",
    "{bvec: [0.0005654926889438294, 6.019147453088565e-05, 0.0022447722791490202]}",
    "{cvec: [1.805654342076572e-05, 0.002548815473400631, -0.03621397089872099]}",
    "{dvec: [0.0024092447900907147, -0.009717820265158338, -8.807433917784529e-06]}",
  },
  // shifting in 5-th direction (+)...
  {
    "{avec: [0.0012075850974850877, -0.005408495873000874, -0.024379032588291473]}",
    "{bvec: [0.0004702385858323742, 0.00011930190186244815, 0.0027596695898884447]}",
    "{cvec: [4.9439312756620866e-05, 0.001729530841348389, -0.03636946480320376]}",
    "{dvec: [0.0020141697196403383, -0.007891848135727951, 2.5423063900131175e-05]}",
  },
  // shifting in 5-th direction (-)...
  {
    "{avec: [0.0012739657505149119, -0.005955615486999124, -0.024324301811708526]}",
    "{bvec: [0.0004953243453676257, 4.3042546137551814e-05, 0.0026518120101115556]}",
    "{cvec: [-1.5090471566208716e-06, 0.0029088819586516102, -0.03629900239679623]}",
    "{dvec: [0.0020908256083596614, -0.007646834744272048, 2.869175209986882e-05]}",
  },
  // shifting in 6-th direction (+)...
  {
    "{avec: [0.0013941961600586292, -0.005752250078470359, -0.024335912951066713]}",
    "{bvec: [0.0005505690398239217, 0.0001129750440107342, 0.0027021039962866504]}",
    "{cvec: [3.481670698011398e-05, 0.00224601582532898, -0.036331843535330095]}",
    "{dvec: [0.002323926197539937, -0.0077008863683725864, 2.834272417602903e-05]}",
  },
  // shifting in 6-th direction (-)...
  {
    "{avec: [0.0010873546879413703, -0.005611861281529639, -0.024367421448933287]}",
    "{bvec: [0.0004149938913760782, 4.936940398926576e-05, 0.00270937760371335]}",
    "{cvec: [1.3113558619886017e-05, 0.002392396974671019, -0.036336623664669894]}",
    "{dvec: [0.0017810691304600624, -0.007837796511627413, 2.5772091823970966e-05]}",
  },
  // shifting in 7-th direction (+)...
  {
    "{avec: [0.0012454313752453725, -0.005613848189523626, -0.024359267901695007]}",
    "{bvec: [0.00048589371414307, 0.0001990769235573227, 0.002717452094888066]}",
    "{cvec: [3.1204932268956716e-05, 0.0023597717881083022, -0.03633446023379283]}",
    "{dvec: [0.0020649995839086663, -0.0077726602657066145, 2.6995779467860616e-05]}",
  },
  // shifting in 7-th direction (-)...
  {
    "{avec: [0.001236119472754627, -0.005750263170476372, -0.024344066498304992]}",
    "{bvec: [0.00047966921705692994, -3.673247555732273e-05, 0.002694029505111934]}",
    "{cvec: [1.672533333104328e-05, 0.002278641011891697, -0.03633400696620716]}",
    "{dvec: [0.0020399957440913334, -0.007766022614293385, 2.711903653213938e-05]}",
  },
  // shifting in 8-th direction (+)...
  {
    "{avec: [0.0012490993646979737, -0.005653192783471472, -0.024355228700926254]}",
    "{bvec: [0.00048165109815985875, 5.796526338073458e-05, 0.002709457132419868]}",
    "{cvec: [3.397776021702497e-05, 0.0023310263169516133, -0.036334623627886614]}",
    "{dvec: [0.0020612874833442064, -0.00776927100384874, 2.7058188885903403e-05]}",
  },
  // shifting in 8-th direction (-)...
  {
    "{avec: [0.0012324514833020259, -0.005710918576528526, -0.024348105699073745]}",
    "{bvec: [0.00048391183304014116, 0.00010437918461926539, 0.0027020244675801325]}",
    "{cvec: [1.3952505382975026e-05, 0.002307386483048386, -0.036333843572113375]}",
    "{dvec: [0.0020437078446557933, -0.0077694118761512596, 2.7056627114096593e-05]}",
  },
  // shifting in 9-th direction (+)...
  {
    "{avec: [0.0012246529983434185, -0.005681364543193438, -0.024351814454451295]}",
    "{bvec: [0.00048047004114849676, 8.023782795287083e-05, 0.002705870051600037]}",
    "{cvec: [2.3352662048822048e-05, 0.002319727168433616, -0.03633424482385493]}",
    "{dvec: [0.0020626132108799377, -0.007769190131293707, 2.705982770755252e-05]}",
  },
  // shifting in 9-th direction (-)...
  {
    "{avec: [0.001256897849656581, -0.0056827468168065606, -0.024351519945548704]}",
    "{bvec: [0.00048509289005150316, 8.210662004712915e-05, 0.002705611548399963]}",
    "{cvec: [2.457760355117795e-05, 0.002318685631566383, -0.03633422237614506]}",
    "{dvec: [0.002042382117120062, -0.0077694927487062925, 2.7054988292447474e-05]}",
  },
  // shifting in 10-th direction (+)...
  {
    "{avec: [0.0012396272967539704, -0.0056823690258229675, -0.024351632250586955]}",
    "{bvec: [0.0004871191595490027, 8.117865822067813e-05, 0.0027056998296810207]}",
    "{cvec: [2.6922346617834093e-05, 0.002319199324290693, -0.03633422467676812]}",
    "{dvec: [0.0020518615497555942, -0.007769352606630588, 2.705722714449867e-05]}",
  },
  // shifting in 10-th direction (-)...
  {
    "{avec: [0.001241923551246029, -0.005681742334177031, -0.024351702149413044]}",
    "{bvec: [0.0004784437716509972, 8.116578977932184e-05, 0.0027057817703189796]}",
    "{cvec: [2.1007918982165905e-05, 0.0023192134757093062, -0.03633424252323187]}",
    "{dvec: [0.0020531337782444055, -0.007769330273369412, 2.7057588855501326e-05]}",
  },
  // shifting in 11-th direction (+)...
  {
    "{avec: [0.0012407115735699984, -0.0056816925430860435, -0.02435171285388519]}",
    "{bvec: [0.00048376543120609697, 8.099043285636128e-05, 0.0027057898823278845]}",
    "{cvec: [2.2535039512706306e-05, 0.0023192918178923396, -0.03633424172055911]}",
    "{dvec: [0.0020524868475665435, -0.007769341577600906, 2.7057405639191926e-05]}",
  },
  // shifting in 11-th direction (-)...
  {
    "{avec: [0.0012408392744300012, -0.005682418816913955, -0.02435162154611481]}",
    "{bvec: [0.00048179749999390295, 8.13540151436387e-05, 0.002705691717672116]}",
    "{cvec: [2.5395226087293692e-05, 0.0023191209821076596, -0.03633422547944088]}",
    "{dvec: [0.0020525084804334562, -0.0077693413023990935, 2.705741036080807e-05]}",
  }
};

vector<vector<string>> BtoD0starVars = {
  // shifting in 1-th direction (+)...
  {
    "{ztp: 0.3392445647927047}",
    "{zeta1: 3.339839371663234}",
  },
  // shifting in 1-th direction (-)...
  {
    "{ztp: -3.299244564792705}",
    "{zeta1: 0.6201606283367658}",
  }
};

vector<vector<string>> BtoD1Vars = {
  // shifting in 1-th direction (+)...
  {
    "{tp: -1.3651051259161324}",
    "{tau1: 2.2486709629509813}",
    "{tau2: 1.5375085845529477}",
  },
  // shifting in 1-th direction (-)...
  {
    "{tp: -0.23489487408386767}",
    "{tau1: 0.3513290370490184}",
    "{tau2: -3.0175085845529472}",
  },
  // shifting in 2-th direction (+)...
  {
    "{tp: -0.40472797320903636}",
    "{tau1: 1.541751591399067}",
    "{tau2: -0.7426223683972742}",
  },
  // shifting in 2-th direction (-)...
  {
    "{tp: -1.1952720267909638}",
    "{tau1: 1.0582484086009327}",
    "{tau2: -0.7373776316027253}",
  }
};

vector<vector<string>> BtoD1starVars = {
  // shifting in 1-th direction (+)...
  {
    "{ztp: 0.3392445647927047}",
    "{zeta1: 3.339839371663234}",
  },
  // shifting in 1-th direction (-)...
  {
    "{ztp: -3.299244564792705}",
    "{zeta1: 0.6201606283367658}",
  }
};

vector<vector<string>> BtoD2starVars = {
  // shifting in 1-th direction (+)...
  {
    "{tp: -1.3651051259161324}",
    "{tau1: 2.2486709629509813}",
    "{tau2: 1.5375085845529477}",
  },
  // shifting in 1-th direction (-)...
  {
    "{tp: -0.23489487408386767}",
    "{tau1: 0.3513290370490184}",
    "{tau2: -3.0175085845529472}",
  },
  // shifting in 2-th direction (+)...
  {
    "{tp: -0.40472797320903636}",
    "{tau1: 1.541751591399067}",
    "{tau2: -0.7426223683972742}",
  },
  // shifting in 2-th direction (-)...
  {
    "{tp: -1.1952720267909638}",
    "{tau1: 1.0582484086009327}",
    "{tau2: -0.7373776316027253}",
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
    df            = df.Define("wff_norescale", "get<1>(ff_result)");
    for (int i = 0; i < numOfFFVar; i++) {
      auto outputBrName = "wff_norescale_var" + to_string(i + 1);
      df = df.Define(outputBrName, "get<" + to_string(i + 2) + ">(ff_result)");
      outputBrs.emplace_back(outputBrName);
    }
    outputBrs.emplace_back("wff_norescale");

    df.Snapshot(trees[idx], ntpOut, outputBrs, writeOpts);

    cout << "Total number of candidates: " << numOfEvt << endl;
    cout << "Hammer reweighted candidates: " << numOfEvtOk << endl;
    cout << "Reweighted fraction: "
         << static_cast<float>(numOfEvtOk) / static_cast<float>(numOfEvt)
         << endl;
  }
}
