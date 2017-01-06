reset
set key outside autotitle columnheader
set yrange [7:0] ; set ytics (1,2,3,4,5,6)
set term pdfcairo enh size 4in, 1.8in mono dashed font "mono"
set out 'r_n.pdf'
plot for [i = 2:7] 'ex1.dat' u 1:i w l lw 3, for [i=1:6] i lt 1
set out

