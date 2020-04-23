* [[file:~/s/research/julia/Fresa.jl/fresa.org::gamsmodel][gamsmodel]]
$TITLE Test Problem 
$OFFDIGIT
$OFFSYMXREF 
$OFFSYMLIST 

VARIABLES X1, X2, X3, Z, res ; 
POSITIVE VARIABLES X1, X2, X3 ; 

EQUATIONS CON1, CON2, CON3, OBJ ;

CON1..  X2 - X3 =G= 0 ; 
CON2..  X1 - X3 =G= 0 ; 
CON3..  X1 - X2**2 + X1*X2 - 4 =E= res ;
OBJ..   Z =E= SQR(X1) + SQR(X2) + SQR(X3) ; 

* Upper bounds 
X1.UP = 5 ; 
X2.UP = 3 ; 
X3.UP = 3 ; 

* Initial point 
X1.L = 4 ; 
X2.L = 2 ; 
X3.L = 2 ; 

MODEL TEST / ALL / ; 

OPTION LIMROW = 0; 
OPTION LIMCOL = 0;
* gamsmodel ends here
