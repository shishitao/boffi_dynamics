reset
set dummy z
set term pdfcairo enhanced mono  dash  size 16cm,04cm font "TeXGyreSchola,6"
set out "stif_from_z.pdf"

set grid
set xrange [0:0.20]
set yrange [0:10]
set xlabel "Damping ratio, {/Symbol z}" font "TeXGyreSchola,6"
set ylabel "Max Stiffness\n(MN/m)" font "TeXGyreSchola,5"
set format x "%04.2f"
plot (sqrt((182*z*z+9)*(182*z*z+9) + 819.) - (182*z*z+9))*pi**2 / 26. t "k = k({/Symbol z})"
set out

