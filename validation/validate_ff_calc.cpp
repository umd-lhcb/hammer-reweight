#include <ff_dstaunu.hpp>

#include <TF1.h>

using namespace std;

enum BMeson { Charged = 1, Neutral = 0 };

enum FFType { ISGW2 = 0, CLN = 1 };

// in GeV
const Double_t m_B0    = 5.2792;
const Double_t m_Bplus = 5.2792;

// q2 distributions with a particular FF parameterization
Double_t q2_cln_B0ToDstTauNu_pdf(Double_t x) {
  BToDstaunu ff_calc{};
  ff_calc.SetMasses(BMeson::Neutral);
  return ff_calc.Compute(x, FFType::CLN, m_B0);
};

int main() {
  auto plot_q2_cln_B0ToDstTauNu =
      TF1("plot_q2_cln_B0ToDstTauNu", "q2_cln_B0ToDstTauNu_pdf(x)", 3, 14);
}
