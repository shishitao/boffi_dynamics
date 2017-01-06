reset
set term pdfcairo enhanced color solid size 16cm,11cm font "TeXGyreSchola,6"

# plot the peak force in problem 2, 1st homework 2010-11,
#  a) as a function of z for a fixed durationof the external force transient
#  b) as a color/contour map as a function of z and t_0

set format x  '%04.2f'
set cbtics (0, 1.25, 2., 3., 4., 5., 10.)
set clabel    '%04.2f'
set ytics 1

# first, plot the peak force response for a duration of the transient
# equal to 6.02 s (i'm lying in the output, saying it's for t_0=6.0)

set xlabel 'Damping ratio {/Symbol z}'
set ylabel 'p_{max}({/Symbol z},t_0=6.0)/kN'
set out 'f_of_z.pdf'
plot  "< awk 'NF==03&&$2==6.02{print $1,$3/1000}' maximum.dat" w l title 'p/kN'

# next, plot a contour+color map of the peak force using damping ratio
# in abscissa and external force transient duration in ordinates

set out 'map.pdf'

# we want a color map
set cbrange [0:10]
set pm3d map
set pm3d implicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean

# there is a nice predefined palette for color maps, this one is
# different...
set palette defined ( 0.000 '#ffffff', 1.188 '#fdf2f3', 1.250 'yellow',\
                      1.312 '#fdf1f1', 1.900 '#fdebec', 2.000 'cyan',\
                      2.100 '#fce8ea', 2.850 '#fce1e2', 3.000 'yellow',\
                      3.150 '#fbdddf', 3.800 '#fbd7d9', 4.000 'cyan',\
                      4.200 '#fad2d5', 4.750 '#facdcf', 5.000 'red',\
                      5.250 '#f9c7ca', 10.000 '#f5969b')

# we want also contour lines at the same values where the palette
# above is discontinous
set contour
set cntrparam levels discrete 1.25, 2.00, 3.00, 4.00, 5.00,10.00
# but not the labels on the contour lines
unset clabel

set yrange [3:8]
set ylabel 'Transient duration, t_0/s'
splot 'maximum.dat' using 1:2:($3/1000.) title 'p_{max}/kN'

# close the output
set out ''
