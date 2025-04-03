np=1

# generate observed data at original datum

cp redatum_obs_org.par redatum.par
time mpirun -n  $np ../../bin/main_redatum_DAW #> out1

mv dataobs_00001.bin dataobs_org_00001.bin

# compute DAWs

cp redatum_DAW.par redatum.par
time mpirun -n  $np ../../bin/main_redatum_DAW #> out1

# resample the DAW at new datum

cp redatum_Resample.par redatum.par
time mpirun -n  $np ../../bin/main_redatum_DAW #> out1

# generate observed data at new datum

cp redatum_obs_new.par redatum.par
time mpirun -n  $np ../../bin/main_redatum_DAW #> out1

cp -r dataobs_00001.bin dataobs_new_00001.bin
