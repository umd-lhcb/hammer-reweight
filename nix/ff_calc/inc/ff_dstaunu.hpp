//------------------------------------------------------------------------------
// File and Version Information:
//      $Id: RateCalc.hh,v 1.2 2021/09/09 00:44:41 yipengsun Exp $
//
// Description:
//    Calculation of the B->DlNu rates based on CLN FF from hep-ph/1203.2654
//    for a lepton of mass ml.
//    The angular dependence is from Korner, Shuler (1990).
//    Normalization is calculated for each set of FF values.
//    The New Physics dependence is parameterized in terms of
//    gSR = -mb (tanBeta/mH)^2 from hep-ph/1203.2654
//
//    FromSP8ToThisModel(q2, ctl,  ctv, chi, isDgamma, lplus, ml)
//       Re-weights SP8 MC to the CLN parameterization.
//    SetMasses(isBm)
//       Sets the masses of the B and D* depending on the charge.
//    Rate(isCLN, ml)
//       Branching fraction.
//    Gamma_q2(q2, A1, V, A2, A0, ml)
//       q2 spectrum.
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//      Michael Mazur                             INFN Pisa
//
/// Revision History:
//      21/09/09 yipengsun -- Updated primary constants
//      21/08/26 yipengsun -- Updated default FF parameters
//      20/10/27 yipengsun -- Reformatted with clang-format
//      12/05/10 manuelf   -- Normalization validated with EvtGen, including
//                            Higgs
//                            Added the theta spectrum integrated over q2
//                            *NOT VALIDATED*
//      12/05/09 manuelf   -- Corrected a missing sin(thetaL)sin(thetaV) and
//                            added D*->Dgamma
//      12/05/03 manuelf   -- Added the NP dependence in terms of gSR
//      12/04/02 manuelf   -- Created off XSLBToDstrtaunu_CLN.cc by M. Mazur
//------------------------------------------------------------------------------

#ifndef BTODSTAUNU
#define BTODSTAUNU

#include <fstream>
#include <iostream>
#include "TMath.h"
#include "TString.h"

using namespace std;
using std::cout;
using std::endl;

class BToDstaunu {
 public:
  BToDstaunu(double rho2 = 1.122, double R1 = 1.270, double R2 = 0.852,
             double R0 = 1.15, double gSR = 0);
  ~BToDstaunu();

  // Primary constants from HAMMER v1.1, PDG and PRD 82 112007 (2010)
  static constexpr double mTau = 1.77686;
  static constexpr double mMu  = 0.10566;
  static constexpr double mE   = 0.000511;
  static constexpr double Vcb  = 0.0415;  // HAMMER
  static constexpr double F1   = 0.908;
  static constexpr double GF   = 0.000011664;
  static constexpr double PI   = 3.14159265;
  static constexpr double hbar = 6.582119e-25;
  static constexpr double mb_quark =
      4.18;  // [GeV] in the MSbar scheme, evaluated at mb.
  static constexpr double mc_quark = 1.27;  // from PDG 2020

  // Secondary constants
  double _mB, _mBSP8, _BLifeTime;
  double _mDs, _mDsSP8, _Dsmaxq2;
  double _rho2, _R1, _R2, _R0, _gSR;
  double _isBm;

  void   SetMasses(int isBm);
  void   ComputeCLN(double q2, double &A1, double &V, double &A2, double &A0);
  void   ComputeISGW2(double q2, double &A1, double &V, double &A2, double &A0);
  void   ComputeLinearQ2(double q2, double &A1, double &V, double &A2,
                         double &A0);
  void   HadronicAmp(double q2, double A1, double V, double A2, double A0,
                     double &H0, double &Ht, double &Hplus, double &Hminus);
  double EvtGetas(double massq, double massx);
  double EvtGetGammaji(double z);

  double Compute(double q2, double ctl, double ctv, double chi, int isDgamma,
                 bool lplus, int isCLN, double ml);
  double Normalization(double ml);
  double Gamma_q2Angular(double q2, double ctl, double ctv, double chi,
                         int isDgamma, bool lplus, double A1, double V,
                         double A2, double A0, double ml);
  double FromSP8ToThisModel(double q2, double ctl, double ctv, double chi,
                            int isDgamma, bool lplus, double ml);

  double Compute(double q2, int isCLN, double ml);
  double Compute(double q2, double ctl, int isCLN, double ml);
  double ComputetL(double thetaL, int isCLN, double ml);
  double Gamma_q2tL(double q2, double ctl, double A1, double V, double A2,
                    double A0, double ml);

  // Functions only used to calculate the Normalization once
  double IntRate(double minX, double maxX, int isCLN, double ctl, double ml,
                 int nPoints = 10000);
  double Gamma_q2(double q2, double A1, double V, double A2, double A0,
                  double ml);
  double Rate(int isCLN, double ml);
  double SubPoly(double Coef[3][7], int indexR);
  void   Polynomial(double ml);
  double SumPoly(double Coef[3][7]);
};

#endif
