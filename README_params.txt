Explanation of Parameter file params.dat

(int)     !nxb - number of collocation points in x
(int)     !nyb - in y
(int)     !nzb - in z  - note, as coded, x,y,z need to be the same
(double)  !mass - mass of black hole when used
(int)     !iLapse - initialize lapse, 0-->one  30-->PlaneWave(AWA)  8-->CosineBump  81-->CossquaredBump 
(int)     !eLapse - lapse evolution,  0-->do not evolve  1-->Harmonic
(int)     !iShift - initialize shift, 0-->zero 30-->PlaneWave (use 0 instead)  8-->CosineBump  81-->CossquaredBump
(int)     !eShift - shift evolution,
(int)     !igij - initialize metric,  0-->Flat  1-->Flat w/rand pert  30-->PlaneWave(AWA)  8-->CosineBump  81-->CossquaredBump  90-->trK distortion
(int)     !iKij - initialize K,       0-->Flat  1-->Zero w/rand pert  30-->PW  8-->CosineBump  81-->CossquaredBump  90-->trK distortion
(int)     !bcflag_lapse - unused
(int)     !bcflag_shift - unused 
(int)     !bcflag_gij   - unused
(int)     !bcflag_Kij   - unused
(int)     !start - from which timestep to start (nonzero means use restart dump)
(int)     !timesteps - number of timesteps before stopping
(int)     !outsteps - output file every (outsteps) number of timesteps
(int)     !dumpsteps - write restart dump (dumpsteps) number of timesteps
(double)  !dt - time interval
(double)  !xleft  - x minimum boundary 
(double)  !xright - x maximum boundary
(double)  !yleft  - y min   Note, these need to be -pi to pi for now
(double)  !yright - y max    unless equations are rewritten.  
(double)  !zleft  - z min  
(double)  !zright - z max  
(int      !iadm - 0-->SmarrYork equations (usually called "ADM")  1-->ADM
equations (true ADM)
