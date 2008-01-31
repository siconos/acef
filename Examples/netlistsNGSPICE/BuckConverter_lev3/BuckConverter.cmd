* lancement simu
.control
set width = 230
set nobreak
rusage everything
source BuckConverter_lev3.cir
run
echo Simulation BuckConverter_lev3 completed !
rusage everything
print col v(6) v(2) @L1[i] v(3) v(7) v(5) @Dpswitch[id] @Dnswitch[id] i(Vmespswitch) i(Vmesnswitch) > BuckConverter_lev3.dat
rusage everything
quit
.endc
.end
