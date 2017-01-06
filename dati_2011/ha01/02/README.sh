#!/bin/sh
#
# this procedure is a shell script for a *nix system,
# to execute the procedure type "sh README.sh <RET>" 
# at a command prompt

echo 'This procedure recreates the graphics in this directory.
Press <Enter> to execute the next step of the procedure.

(you may want to change the font name in *.gp files)
'
echo -n '% gnuplot p_of_t.gp' ; read
pwd;gnuplot p_of_t.gp
echo -n 'Graphical file "p_of_t.pdf" created.
% gnuplot stif_from_z.gp' ; read
gnuplot stif_from_z.gp
echo -n 'Graphical file "stif_from_z.pdf" created.
% python integ.py > response.dat' ; read
python integ.py > response.dat
echo -n 'Data file "response.dat" created.
% gnuplot response.gp' ; read
gnuplot response.gp
echo -n 'Graphical files "displacement.pdf", "velocity.pdf",
        "acceleration.pdf" and "force.pdf" created.
% python maximum.py > maximum.dat' ; read
echo "This is going to take some time..."
time python maximum.py > maximum.dat
echo -n 'Data file "maximum.dat" created.
% gnuplot map.gp' ; read
gnuplot map.gp
echo 'Graphical file "map.pdf" created.
All files updated.'
