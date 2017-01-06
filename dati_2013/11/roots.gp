reset
set term X11 enhanced
set xrange [0:5]
set yrange [-.35:.05]
set xtics 0.5, 1 format "%.1f{/Symbol p}"
set ytics -0.3, 0.1
set grid linestyle 1
set key right  bottom horizontal
set samples 2000
set term pdfcairo enh size 4in, 1in mono dashed
set out 'roots2.pdf'
plot cos(x*pi) lt 2 lw 3 t "cos({/Symbol b}L)",-1./cosh(x*pi) lt 3 lw 3 t "-1/cosh({/Symbol b}L)"
set out
!pdflatex lez11.tex 2>&1 > /dev/null &