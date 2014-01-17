#!/bin/bash
# Example environmental variables to add for IDLUTILS/kCorrect/iSEDfit
# add this text to your .bashrc file

export IDL_LIBS=/global/data/products

# IDLUtils -- Helper Library
export IDLUTILS_DIR=$IDL_LIBS/idlutils
export IDL_PATH=$IDL_PATH:+$IDLUTILS_DIR/pro
export IDL_PATH=$IDL_PATH:+$IDLUTILS_DIR/goddard/pro

# RED -- Cosmological utilities
export RED_DIR=$HOME/idl/red
export IDL_PATH=$IDL_PATH:+$RED_DIR

# KCorrect
export KCORRECT_DIR=$IDL_LIBS/kcorrect
export PATH=$PATH:$KCORRECT_DIR/bin
export IDL_PATH=$IDL_PATH:+$KCORRECT_DIR/pro
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$KCORRECT_DIR/lib

# IMPRO / ISEDFIT
export IMPRO_DIR=$HOME/idl/impro
export IDL_PATH=$IDL_PATH:+$IMPRO_DIR

