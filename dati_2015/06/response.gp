set term pdf mono dashed size 10cm, 5cm enh
set dummy t
set samples 2000
set xrange [0:30]
set key top center horizontal
ww1 = (7.0-sqrt(33.0))/4.0
ww2 = (7.0+sqrt(33.0))/4.0
w1  = sqrt(ww1)
w2  = sqrt(ww2)
#
p21 = 1.0
p11 = 2.0/(3.0-2.0*ww1)
p22 = 1.0
p12 = 2.0/(3.0-2.0*ww2)
#
m1  = 2.0*p11*p11+1.0
m2  = 2.0*p12*p12+1.0
#
p11 = p11/sqrt(m1)
p21 = p21/sqrt(m1)
p12 = p12/sqrt(m2)
p22 = p22/sqrt(m2)
# print p11,p12
# print p21,p22
b1(w) = w/w1; b2(w) = w/w2
dst1 = p21/ww1; dst2 = p22/ww2
# print dst1, dst2
q1(t)=dst1*(sin(w*t)-b1(w)*sin(w1*t))/(1-b1(w)*b1(w))
q2(t)=dst2*(sin(w*t)-b2(w)*sin(w2*t))/(1-b2(w)*b2(w))
x1(t)=p11*q1(t)+p12*q2(t)
x2(t)=p21*q1(t)+p22*q2(t)
w = 2.0
set xlabel "{/Symbol a} = {/Symbol w}_o t"
#
set out 'x_response.pdf'
set yrange [*:2.5]
set ylabel "x_i/{/Symbol D}_{st}"
plot x1(t) lw 3 t "x_1({/Symbol a})/{/Symbol D}_st",\
     x2(t) lw 3 t "x_2({/Symbol a})/{/Symbol D}_st"
#
set out 'q_response.pdf'
set yrange [*:2.5]
set ylabel "q_i/{/Symbol D}_{st}"
plot q1(t) lw 3 t "q_1({/Symbol a})/{/Symbol D}_st",\
     q2(t) lw 3 t "q_2({/Symbol a})/{/Symbol D}_st"
#
set out
