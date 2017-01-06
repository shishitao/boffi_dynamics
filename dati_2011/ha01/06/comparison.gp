reset

set term pdfcairo enhanced color solid size 20cm,04cm font "TeXGyreSchola,6"
set out 'comparison.pdf'

set dummy a
set xrange [0:6]
set xlabel '{/Symbol a} = {/Symbol w}_o t' font ",6"
set ylabel '200 x_3({/Symbol a}) / L' font ",6"
set ytics 0.2
set mxtics 0.1
wf=        7.0;w1=+1.191108;     w2=+1.77818;   w3=+3.3330414
cf=-0.11692366;c1=-0.00070510474;c2=-0.02192027;c3=+0.2574789

set grid; set samples 1000; set key center bottom horizontal font ',04'

plot cf*sin(wf*a)+c1*sin(w1*a)+c2*sin(w2*a)+c3*sin(w3*a) \
     lw 9  lc rgb '#ff8080' t 'analytical', \
     cf*sin(wf*a)+c1*sin(w1*a)+c2*sin(w2*a)+c3*sin(w3*a) \
     lw 04  lc rgb '#ffffff' t '', \
     'num_res.dat' u 1:04 w l lw 1  lc 'black' t 'numerical'

set out
