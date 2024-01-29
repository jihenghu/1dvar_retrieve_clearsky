#!/bin/sh

# Test case input data
COEF_FILENAME="/home/jihenghu/rttov13/rtcoef_rttov13/rttov13pred54L/rtcoef_gpm_1_gmi.dat"  # Location of this file is set below in $COEF_DIR
PROF_FILENAME="prof.dat"                        # Input profile(s), usually found in $TEST_DIR set below
NPROF=1                                         # Number of profiles defined in prof.dat
NLEVELS=51                                      # Number of profile levels
DO_SOLAR=0                                      # 0 = solar off / 1 = solar on

NCHAN=9                                        # Number of channels to simulate for each profile
CHAN_LIST=$(seq -s ' ' $NCHAN)                  # Space-separated channel-list

NTHREADS=1                                      # Number of threads to use (compile RTTOV with OpenMP to exploit this)

echo " "
echo " "
echo " Test K "
echo " "

echo  "Coef filename:      ${COEF_FILENAME}"
echo  "Input profile file: ${PROF_FILENAME}"
echo  "Number of profiles: ${NPROF}"
echo  "Number of levels:   ${NLEVELS}"
echo  "Do solar:           ${DO_SOLAR}"
echo  "Number of channels: ${NCHAN}"
echo  "Channel list:       ${CHAN_LIST}"
echo  "Number of threads:  ${NTHREADS}"

./example_k.exe << EOF
"${COEF_FILENAME}", Coefficient filename
${PROF_FILENAME}  , Input profile filename
${NPROF}          , Number of profiles
${NLEVELS}        , Number of levels
${DO_SOLAR}       , Turn solar radiation on/off
${NCHAN}          , Number of channels
${CHAN_LIST}      , Channel numbers
${NTHREADS}       , Number of threads
EOF

exit
