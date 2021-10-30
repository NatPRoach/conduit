from src.package_version import ConduitVersionString

# Package

version       = ConduitVersionString
author        = "Nathan Roach"
description   = "De novo transcriptome assembler"
license       = "GPLv2"

# Dependencies

requires "hts >= 0.3.1", "nim >= 1.4.0"

srcDir = "src/"

before install:
  echo "Building poaV2"
  withDir "src/poaV2":
    exec "make clean"
  if existsEnv("CC"):
    withDir "src/poaV2":
      exec "make CC=$CC"
  else:
    withDir "src/poaV2":
      exec "make"

before build:
  echo "Building poaV2"
  withDir "src/poaV2":
    exec "make clean"
  if existsEnv("CC"):
    withDir "src/poaV2":
      exec "make CC=$CC"
  else:
    withDir "src/poaV2":
      exec "make"

bin = @["conduit", "conduitUtils", "conduit_clustering"]
