print "file = $0 ," , "time step = $1 ns ,", "time begin = $2 ns ,", "time end = $3 ns ," , "file sampling = 1 out of $4 ," , "starting at $5"
system "head --lines=`echo 'scale=0;$3/$1' | bc -l` $0 | tail --lines=+`echo 'scale=0;$2/$1' | bc -l` | sed -n $5~$4p > temp.dat"

reset

xinc = ($3e-9 - $2e-9)/5.0
set xtics xinc
set grid xtics
set grid ytics
set lmargin 10
set rmargin 28
set tmargin 0
set bmargin 0

set format x ""
set xrange [$2e-9:$3e-9]

xper = 1
set key outside nobox samplen 2 spacing 1.0 title "NGSPICE results"

set multiplot

set size 1.0,0.3
set xlabel

set origin 0.,0.7
set ylabel "A"
set yrange[-0.2:1]
plot\
    "temp.dat" every xper u 2:5 t  "L current" w l,\
    "temp.dat" every xper u 2:10 t  "NMOS diode current" w l

set origin 0.,0.4
set yrange[-0.2:1]
plot\
    "temp.dat" every xper u 2:11 t "PMOS ISD current" w l,\
    "temp.dat" every xper u 2:12 t "NMOS ISD current" w l

set origin 0.,0.1
set format x "%g"
set xlabel "time in s"
set ylabel "V"
set yrange[-1:4]
plot\
    "temp.dat" every xper u 2:6 t "Output voltage" w l,\
    "temp.dat" every xper u 2:3 t "ramp voltage" w l,\
    "temp.dat" every xper u 2:4 t "MOS P drain potential" w l,\
    "temp.dat" every xper u 2:8 t "Error voltage" w l

set nomultiplot

#    "temp.dat" every xper u 1:6 t "Gate voltage" w l
