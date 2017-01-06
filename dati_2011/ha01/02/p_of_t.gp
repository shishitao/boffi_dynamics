reset
set term pdfcairo enhanced mono  solid size 16cm,04cm font "TeXGyreSchola,6"

# plot the forcing function for problem #2 in 1st homework 2010-11

set dummy t           ; set key bottom left
set samples 2000      ; set grid
set xrange [-1:7]     ; set yrange [-1.2:1.2]
set xlabel "Time (s)" ; set ylabel "Force (kN)"
set xtics 1           ; set ytics 1

set out 'p_of_t.pdf'
plot t<0?0:t<6?t/6*sin(pi*5.*t*t/6.):sin(pi*10.*t) t "p(t)"
set out ''
