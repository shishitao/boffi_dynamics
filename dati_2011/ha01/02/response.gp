reset
set term pdfcairo enhanced mono  dash  size 16cm,04cm font "TeXGyreSchola,6"
file='response.dat'

set xlabel 'Time/s' font ",7"
set xrange [0:10]
set xtics 1
set format x "%2.0f"
set format y "%2.1f"
set key top left
set grid

set out 'displacement.pdf'
set ylabel 'Displacement/mm' font ",7"
plot file u 1:($2*1000) w l t ''

set out 'velocity.pdf'
set ytics 5
set mytics 5
set mgrid
set ylabel 'Velocity/(mm/s)' font ",7"
plot file u 1:($4*1000) w l t '' 

set out 'force.pdf'
set ytics 2.5
set ylabel 'Force/kN' font ",7"
plot file u 1:($3/1000) w l t ''

set out 'acceleration.pdf'
set ylabel 'Acceleration/(m/s^2)' font ",5"
set ytics  0.10 font ",4"
set format y "%+4.1f"
set y2tics 0.25 font ",4"
set y2label "External Force  (kN)" font ",5"
plot file u 1:5 w l lw 3 t 'Acceleration',\
     file u 1:($6/1000.) axis x1y2 w l lt 1 t 'Load' 

set out
