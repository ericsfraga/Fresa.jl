$include gamsdeclarations.gms
X1.fx = 2.387182340440622; 
X2.fx = 0.8810162763650371; 
X3.fx = 0.0; 
solve TEST using NLP minimizing Z; 
file fresa /'gamsoutput.txt'/ ;
put fresa ;
put z.l /;
put res.l /;
put TEST.modelstat /;
