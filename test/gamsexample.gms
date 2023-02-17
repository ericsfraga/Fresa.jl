$include gamsdeclarations.gms
X1.fx = 1.546940052510804; 
X2.fx = 0.7699931892418711; 
X3.fx = 0.005636823583293402; 
solve TEST using NLP minimizing Z; 
file fresa /'gamsoutput.txt'/ ;
put fresa ;
put z.l /;
put res.l /;
put TEST.modelstat /;
