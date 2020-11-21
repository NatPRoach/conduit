# Package

version       = "0.1.0"
author        = "Nathan Roach"
description   = "De novo transcriptome assembler"
license       = "GPLv2"

# Dependencies

requires "hts >= 0.3.1", "nim >= 1.0.0"

srcDir = "src/"

before install:
  echo "Building poaV2"
  withDir "src/poaV2":
    exec "make"

bin = @["conduit", "conduitUtils", "conduit_clustering"]
