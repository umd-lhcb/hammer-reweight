// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon May 02, 2022 at 05:34 PM -0400

#pragma once

#include <tuple>

#include <TString.h>

#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>

using std::tuple;

////////////////////////////
// HAMMER-related helpers //
////////////////////////////

TString printP(const Hammer::FourMomentum& p) {
  char tmp[80];
  sprintf(tmp, "%.2f, %.2f, %.2f, %.2f", p.E(), p.px(), p.py(), p.pz());
  return TString(tmp);
}

auto buildPartVec(double pe, double px, double py, double pz, int pid) {
  return tuple<double, double, double, double, int>{pe, px, py, pz, pid};
}

auto buildHamPart(double pe, double px, double py, double pz, int pid) {
  auto fourMom = Hammer::FourMomentum(pe, px, py, pz);
  auto partId  = static_cast<Hammer::PdgId>(pid);

  return Hammer::Particle(fourMom, partId);
}

auto buildHamPart(Hammer::FourMomentum fourMom, int pid) {
  auto partId = static_cast<Hammer::PdgId>(pid);
  return Hammer::Particle(fourMom, partId);
}

auto buildHamPart(tuple<double, double, double, double, int> pack) {
  auto [pe, px, py, pz, pid] = pack;
  return buildHamPart(pe, px, py, pz, pid);
}
