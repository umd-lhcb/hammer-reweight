self: super:

let
  pythonOverrides = {
    packageOverrides = self: super: {
      cymove = super.callPackage ./cymove {};
    };
  };
in

{
  python3 = super.python3.override pythonOverrides;

  hammer-phys = super.callPackage ./hammer-phys {};
  ff_calc = super.callPackage ./ff_calc {};
}
