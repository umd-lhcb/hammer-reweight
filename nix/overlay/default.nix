self: super:

let
  pythonOverrides = {
    packageOverrides = self: super: {
      cymove = super.callPackage ./cymove {};
    };
  };
in

{
  hammer-phys = super.callPackage ./hammer-phys {};
  python3 = super.python3.override pythonOverrides;
}
