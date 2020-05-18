
* Minimal Regulator Identifier (MiReN) GAMS code
* Mohammad Mazharul Islam
* emailto: mislam3@huskers.unl.edu
* University of Nebraska-Lincoln, Lincoln, NE 68588
* GAMS code for identifying minimal regulatory network and global regulators
* in a gene expression network 
* from experimental gene expression data
* (network inference using the slack variables)
*
* single level formulation

OPTION decimals = 8
       sysout = off
       solprint = on
       resLim = 10000000
       iterLim = 10000000
       domlim = 1000
       limcol = 1000
       limrow = 1000
       optCA = 1e-9
       optCR = 1e-9
       work = 10000000
       dmpopt
       dmpSym
       dualcheck = 1
       eject
       measure
       memoryStat
       profile = 1
       threads = 2
       LP = cplex
       nlp = minos5
       mip = cplex
       minlp = baron;

Sets
ii   set of regulators
/
$include reg_list.txt
/

i   set of genes
/
$include gene_list.txt
/

t   time points
/
$include time_points.txt
/

t_sub(t) all time points but 1st
/
$include time_sub_points.txt
/

;

Parameters
X(i,t)  experimental expression values considering genes
/
$include gene_expression_levels.txt
/

TF(ii,t) experimental tf list
/
$include reg_expression_level.txt
/

time_step(t) time step from time t to time t+1
/
$include time_steps.txt
/
;

Variables
A(ii,i)  Interaction coefficient of transcription factor on gene

z_outer     outer objective 

z_single    single level objective

z   objective
z_dual  dual objective

p(ii,i)  linearizing variables 
q(ii,i)  linearizing variables 


lambda(i,t) dual variables associated with equality constraints 
delta(i)    dual variable associated with equality contraint on the first time point
;

Positive Variables
Dpos(i,t)     positive slack variable
Dneg(i,t)     negative slack variable

alpha(ii,i)
beta(ii,i)
;

Binary Variables
y(ii,i)  number of transcription factors acting on a gene
;

Scalar
total_tf
m
epsilon
negepsilon
total_edges
W
;

total_tf = 10;
m=1000;
W = 1000;
epsilon=0.0000001;
negepsilon=-0.0000001;
total_edges=9;

Equations
    outer
    
    primaldual
    
    single_objective
    
    obj
    main(i,t)
    main2(i)
    tf_on_gene_low(ii,i)
    tf_on_gene_high(ii,i)
    max_tf(i)
    
    dualobj
    dual1(ii,i)
    dual2
    dual3
    dual4
    dual5
    
    linearize11(ii,i)
    linearize12(ii,i)
    linearize13(ii,i)
    linearize14(ii,i)

    linearize21(ii,i)
    linearize22(ii,i)
    linearize23(ii,i)
    linearize24(ii,i)
;

***********************************Equations***********************************
* Outer level objective 
outer.. z_outer =e= sum(ii,sum(i,y(ii,i)));

max_tf(i)..sum(ii,y(ii,i)) =l= total_tf;

primaldual..    sum(i,sum(t,Dpos(i,t)+Dneg(i,t))) =e= - sum(i,sum(ii,(m*p(ii,i)) + (m*q(ii,i)))) + sum(i,sum(t$t_sub(t),lambda(i,t)*(X(i,t)-X(i,t-1))));

*single level primal objective
single_objective..  z_single =e= sum(i,sum(t,Dpos(i,t)+Dneg(i,t))) + sum(ii,sum(i,y(ii,i)));

****************Primal equations******************
obj..   z =e= sum(i,sum(t,Dpos(i,t)+Dneg(i,t)));
main(i,t)$t_sub(t)..    X(i,t) - X(i,t-1) - time_step(t)*sum(ii,A(ii,i)*TF(ii,t-1)) =e= Dpos(i,t) - Dneg(i,t);

*ignoring/including the main2 constraint (ignoring the first time point)
main2(i)..  X(i,'1')-X(i,'7')-time_step('1')*sum(ii,A(ii,i)*TF(ii,'7'))=e=Dpos(i,'1')-Dneg(i,'1');

tf_on_gene_low(ii,i)..   A(ii,i) =l= m*y(ii,i);
tf_on_gene_high(ii,i)..  A(ii,i) =g= -m*y(ii,i);


******************dual equations*****************
dualobj..   z_dual =e= - sum(i,sum(ii,(m*p(ii,i)) + (m*q(ii,i)))) + sum(i,sum(t$t_sub(t),lambda(i,t)*(X(i,t)-X(i,t-1))));
dual1(ii,i).. sum(t$t_sub(t),lambda(i,t)*time_step(t)*TF(ii,t)) + alpha(ii,i) - beta(ii,i) =l= 0;

dual2(i,t)..    lambda(i,t) =g= -1;
dual3(i,t)..    lambda(i,t) =l= 1;

dual4(i)..    delta(i) =g= -1;
dual5(i)..    delta(i) =l= 1;

* linearization constraints 
linearize11(ii,i)..   p(ii,i) =l= W*y(ii,i);
linearize12(ii,i)..   p(ii,i) =g= -W*y(ii,i);
linearize13(ii,i)..   p(ii,i) =l= alpha(ii,i) + W*(1-y(ii,i));
linearize14(ii,i)..   p(ii,i) =g= alpha(ii,i) - W*(1-y(ii,i));

linearize21(ii,i)..   q(ii,i) =l= W*y(ii,i);
linearize22(ii,i)..   q(ii,i) =g= -W*y(ii,i);
linearize23(ii,i)..   q(ii,i) =l= beta(ii,i) + W*(1-y(ii,i));
linearize24(ii,i)..   q(ii,i) =g= beta(ii,i) - W*(1-y(ii,i));

********************************************************************************

Model Single_level
/
single_objective

max_tf

main
main2
tf_on_gene_low
tf_on_gene_high
/
;

alias(i,i1);


file Gene_TF_single /MiReN_results_single_level.txt/;
PUT Gene_TF_single;

Single_level.optfile=1;
Solve Single_level using mip minimizing z_single

PUT "-----------------------------------------------------------------------"//;

PUT "Solver termination condition: ", Single_level.solveStat/; 
PUT "Model stat: ", Single_level.modelstat/;
PUT "ELAPSED TIME = ", Single_level.etSolve, " seconds"//;

PUT "objective = ", z_single.l:10:5//;

PUT "*************************Interaction coefficients***********************"//;

loop(i,
       PUT "gene: ";
       PUT i.tl:0:0/;
       loop(ii,
              PUT "A(",ii.tl:0:0,",",i.tl:0:0,") = ",A.l(ii,i):10:5,"    y(",ii.tl:0:0,",",i.tl:0:0,") = ",  y.l(ii,i)/; 
       );
       PUT /;
);

PUT /;
