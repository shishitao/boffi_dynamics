reset
unset xlabel
unset ylabel
set xrange [0:1]
set yrange [*:1.4]
set zeroaxis lt 1
set xtics 0,.25`
set ytics 0,.5
set key top center horizontal samplen 2
#set ylabel "{/Symbol f}_n"
#set xlabel "x/L"
coeff(b)=(cosh(b)+cos(b))/(sinh(b)+sin(b)) 
f1(b)=cosh(b)-cos(b) -coeff(b)*(sinh(b)-sin(b)) 
f(x)=(cosh(b*x)-cos(b*x)-coeff(b)*(sinh(b*x)-sin(b*x)))/f1(b)
set term pdfcairo enh size 4in, 1.2in # mono dashed
set out 'clamped_eigenfunctions.pdf'
plot b=1.8751, f(x) lw 3 lt 7 t "n=1",\
     b=4.6941, f(x) lw 3 lt 1 t "2",\
     b=7.8548, f(x) lw 3 t "3",\
     b=10.996, f(x) lw 3 t "4",\
     b=4.5*pi, f(x) lw 3 t "5"
set out
!pdflatex lez11.tex 2>&1 > /dev/null &