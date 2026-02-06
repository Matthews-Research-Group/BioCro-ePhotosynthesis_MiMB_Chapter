#!/bin/bash

# this is using clang++ on Mac. Change it to match your computer's
# aliase for using c++20
alias myCppCompiler='clang++ -std=c++20'
# Get the directory of the script
SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname "$SCRIPT")
# Now go back once 
MAIN_DIR=$(dirname "$SCRIPTPATH")
ephoto_path="$MAIN_DIR/models/ePhotosynthesis_C"
sundial_path="/opt/homebrew/include"

# myephoto_single.exe will be needed for plotting
myCppCompiler -o myephoto_single.exe -I$sundial_path -I$ephoto_path/include -I$ephoto_path/build run_ephoto_EPS_singlerun.cpp -L$ephoto_path/build -lEPhotosynthesis -Wl,-rpath,$ephoto_path/build
