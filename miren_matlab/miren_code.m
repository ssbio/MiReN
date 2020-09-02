% MiReN (Minimum Regulator Identifier) Matlab code
% Author: Mohammad mazharul Islam
% emailto: mislam3@huskers.unl.edu
% MiReN is a multi-objective Mixed-Integer Linear Programming (MILP)
% framework to infer minimal regulatory influences on 
% genes of interest from temporal transcriptomic data.

% Notes:
% Require matlab version 2019a or newer


clear all;

miren = optimproblem;

%% Reading inputs

% Time points
time_point = readcell('time_points.txt');

% time sub points
time_sub_point = readcell('time_sub_points.txt');

% time steps
d = load('time_steps.txt');
time_step = d(:,2);

% Gene expression
gene_exp = readtable('gene_expression_levels_df.txt');
gene_exp = gene_exp(2:end,2:end);
gene_exp = table2array(gene_exp);

% Regualtor expression
reg_exp = readtable('reg_expression_levels_df.txt');
reg_exp = reg_exp(2:end,2:end);
reg_exp = table2array(reg_exp);


%% Other parameters
% number of max regulatory connections to a gene
K = 10; % Max number of regulatory connection for a gene
lambda = 10; % Coeffcient for the L1-regularization (penalty) term in the objective function


%% Sets

II = length(reg_exp);
I = length(gene_exp);
T = length(time_point);
Tsub = length(time_sub_point);


%% Variable declaration

% y(ii,i) 
y = optimvar('y',[II,I],'Type','integer','LowerBound',0,'UpperBound',1);

% A(ii, i)
A = optimvar('A',[II,I],'Type','continuous','LowerBound',-1000,'UpperBound',1000);

% Dpos(i,t)
Dpos = optimvar('Dpos',[I,T],'Type','continuous','LowerBound',0,'UpperBound',100);

% Dneg(i,t)
Dneg = optimvar('Dneg',[I,T],'Type','continuous','LowerBound',0,'UpperBound',100);

% Qpos(ii,i)
Qpos = optimvar('Qpos',[II,I],'Type','continuous','LowerBound',0,'UpperBound',100);

% Qneg(ii,i)
Qneg = optimvar('Qneg',[II,I],'Type','continuous','LowerBound',0,'UpperBound',100);


%% Constraints

maxreg = optimconstr(I);
for i = 1:I
    maxreg(i) = sum(y(:,i)) == K;
end
% show(maxreg)

regression = optimconstr(I*T);
for i = 1:I
    for t = 2:T
        regression(i,t) = Dpos(i,t) - Dneg(i,t) - gene_exp(i,t) + gene_exp(i,t-1) + time_step(t)*sum(A(:,i).*reg_exp(:,t-1)) == 0;
    end
end
% show(regression)

LB = optimconstr(II*I);
for ii = 1:II 
    for i = 1:I
        LB(ii,i) = A(ii,i) + 1000*y(ii,i) >= 0;
    end
end
% show(LB)

UB = optimconstr(II*I);
for ii = 1:II 
    for i = 1:I
        UB(ii,i) = A(ii,i) - 1000*y(ii,i) <= 0;
    end
end
% show(UB)

Q = optimconstr(II*I);
for ii = 1:II 
    for i = 1:I
        Q(ii,i) = Qpos(ii,i) - Qneg(ii,i) - A(ii,i) == 0;
    end
end
% show(Q)


%% objective function

z = sum(sum(y)) + sum(sum(Dpos)) + sum(sum(Dneg)) + lambda*(sum(sum(Qpos)) + sum(sum(Qneg)));
miren.Objective = z;
% show(z)

%% Optimization model definition

miren.Constraints.max_reg = maxreg;
miren.Constraints.regressions = regression;
miren.Constraints.lowerbound = LB;
miren.Constraints.upperbound = UB;
miren.Constraints.Qslack = Q;

%% Solution

[sol,fval] = solve(miren);

