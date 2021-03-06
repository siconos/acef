print "file = $0 ," , "time step = $1 ns ,", "time begin = $2 ns ,", "time end = $3 ns ," , "file sampling = 1 out of $4"
system "head --lines=`echo 'scale=0;$3/$1' | bc -l` $0 | tail --lines=+`echo 'scale=0;$2/$1' | bc -l` | sed -n 1~$4p > temp.dat"

reset

xinc = ($3e-9 - $2e-9)/5.0
set grid xtics
set grid ytics
set xtics xinc
set format x ""
set lmargin 10
set rmargin 28
set tmargin 0
set bmargin 0
set xrange [$2e-9:$3e-9]

xper = 1
set key outside nobox samplen 2 spacing 1.0 title "NGSPICE results"

set multiplot

set size 1.0,0.3
set origin 0.,0.7
set xlabel
set ylabel "V"
set yrange[-1:6]
plot\
    "temp.dat" every xper u 2:3 t "ramp voltage" w l,\
    "temp.dat" every xper u 2:4 t "MOS P drain potential" w l,\
    "temp.dat" every xper u 2:6 t "Output voltage" w l,\
    "temp.dat" every xper u 2:7 t "Gate voltage" w l,\
    "temp.dat" every xper u 2:8 t "Verror voltage" w l

set origin 0.,0.4
set format x "%g"
set xlabel "time in s"
set ylabel "A"
set yrange[-0.1:0.7]
plot\
    "temp.dat" every xper u 2:5      t "L1 current" w l,\
    "temp.dat" every xper u 2:11     t "MOSP current" w l,\
    "temp.dat" every xper u 2:12     t "MOSN current" w l

set origin 0.,0.1
set format x "%g"
set xlabel "time in s"
set ylabel "A"
set yrange[-0.02:0.32]
plot\
    "temp.dat" every xper u 2:10     t "diode NMOS current" w l

set nomultiplot

#    "temp.dat" every xper u 1:8     t "diode PMOS current" w l
#    "temp.dat" u 1:10 t "Comparator output voltage SPICE" w l
#    "temp.dat" every xper u 2:9      t "diode PMOS current" w l,\
#    "temp.dat" every xper u 2:(-$11) t "ValimI current" w l,\
