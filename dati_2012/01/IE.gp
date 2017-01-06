# ISOLATION EFFICIENCY
#
set xrange [0:2]; set yrange [0:50]
set grid
set dummy delta
set size 0.5, 0.5
set terminal postscript eps enhanced monochrome dashed defaultplex 'Helvetica' 12
set output "IE.eps"
set xlabel "{/Symbol D}_{st} [cm]"
set ylabel "Input frequency [Hz]"
set ytics 5; set xtics .2
freq(delta)=(113./(2*355.))*sqrt((981./delta)*(2.-IE)/(1.-IE))
plot IE=0.00, freq(delta) linewidth 2 title "IE=0.00", \
     IE=0.50, freq(delta) linewidth 3 title "IE=0.50", \
     IE=0.60, freq(delta) linewidth 2 title "IE=0.60", \
     IE=0.70, freq(delta) linewidth 3 title "IE=0.70", \
     IE=0.80, freq(delta) linewidth 2 title "IE=0.80", \
     IE=0.90, freq(delta) linewidth 3 title "IE=0.90", \
     IE=0.95, freq(delta) linewidth 2 title "IE=0.95", \
     IE=0.98, freq(delta) linewidth 3 title "IE=0.98", \
     IE=0.99, freq(delta) linewidth 2 title "IE=0.99"
set out
