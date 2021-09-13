final: prev:

{
  ff_calc = prev.callPackage ./ff_calc { };
  hammer-reweight = prev.callPackage ./default.nix { };
}
