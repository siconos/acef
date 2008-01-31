* lancement simu
.control
set width = 230
set nobreak
rusage everything
source BuckConverterModSah_abs.cir
run
echo Simulation BuckConverterModSah_abs completed !
rusage everything
print col v(6) v(2) @L1[i] v(3) v(7) v(5) @Dpswitch[id] @Dnswitch[id] i(Vmespswitch) i(Vmesnswitch) > res.txt
rusage everything
quit
.endc
.end
