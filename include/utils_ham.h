// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Tue Apr 26, 2022 at 11:43 PM -0400

#pragma once

#include <TString.h>

#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>

////////////////////////////
// HAMMER-related helpers //
////////////////////////////

TString printP(const Hammer::FourMomentum& p) {
  char tmp[80];
  sprintf(tmp, "%.2f, %.2f, %.2f, %.2f", p.E(), p.px(), p.py(), p.pz());
  return TString(tmp);
}

auto particle(double pe, double px, double py, double pz, int pid) {
  auto fourMom = Hammer::FourMomentum(pe, px, py, pz);
  auto partId  = static_cast<Hammer::PdgId>(pid);

  return Hammer::Particle(fourMom, partId);
}

auto particle(Hammer::FourMomentum fourMom, int pid) {
  auto partId = static_cast<Hammer::PdgId>(pid);
  return Hammer::Particle(fourMom, partId);
}
