$include gamsdeclarations.gms
X1.fx = 4.317044786795817; 
X2.fx = 3.0; 
X3.fx = 2.1894282915266405; 
solve TEST using NLP minimizing Z; 
file fresa /'gamsoutput.txt'/ ;
put fresa ;
put z.l /;
put res.l /;
put TEST.modelstat /;
