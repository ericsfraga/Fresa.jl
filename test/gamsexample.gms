$include gamsdeclarations.gms
X1.fx = 4.0; 
X2.fx = 2.0; 
X3.fx = 2.0; 
solve TEST using NLP minimizing Z; 
file fresa /'gamsoutput.txt'/ ;
put fresa ;
put z.l /;
put res.l /;
put TEST.modelstat /;
