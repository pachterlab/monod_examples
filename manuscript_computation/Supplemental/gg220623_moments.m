function gg220525_moments
clc;clear;close all;
syms s b un um bet gam tau alf

assume(s,'positive')
assume(b,'positive')
assume(bet,'positive')
assume(gam,'positive')
assume(tau,'positive')
assume(alf,'positive')

f = bet/(bet-gam);

%%% use this for the null case
% xn = un;
% xm = um;

%%% use this for the Bernoulli case
% syms pn pm
% assume(pn,'positive')
% assume(pm,'positive')
% xn = un*pn;
% xm = um*pm;

%%% use this for the Poisson case
syms Ln Lm
assume(Ln,'positive')
assume(Lm,'positive')
xn = exp(un*Ln)-1;
xm = exp(um*Lm)-1;


U = xm*f*exp(-gam*s) + (xn-xm*f)*exp(-bet*s);
Uint = int(U,s,0,inf);
%%%%%%%%%%%%%%%%%%%
fprintf('\nbursty\n')
G1 = 1/(1-b*U)-1;

m1int = int(subs(subs(diff(G1,un),un,0),um,0));
m1 = simplify(subs(m1int,s,inf)-subs(m1int,s,0));

m2int = int(subs(subs(diff(G1,um),un,0),um,0));
m2 = simplify(subs(m2int,s,inf)-subs(m2int,s,0));

v1int = expand(int(subs(subs(diff(G1,un,2),un,0),um,0))) ;
v1 = simplify((subs(v1int,s,inf)-subs(v1int,s,0))/m1); %Poisson

v2int = expand(int(subs(subs(diff(G1,um,2),un,0),um,0)));
v2 = simplify((subs(v2int,s,inf)-subs(v2int,s,0))/m2); %Poisson

covint = expand(int(subs(subs(diff(diff(G1,um),un),un,0),um,0)));
cov_ = simplify((subs(covint,s,inf)-subs(covint,s,0))); %Poisson

fprintf('m1: %s\n',m1)
fprintf('m2: %s\n',m2)
fprintf('v1: %s\n',v1)
fprintf('v2: %s\n',v2)
fprintf('cov: %s\n',cov_)

% fprintf('')

%%%%%%%%%%%%%%%%%%%%
fprintf('\nconstitutive\n')

G2 = Uint;

m1 = subs(subs(diff(G2,un),un,0),um,0);
m2 = subs(subs(diff(G2,um),un,0),um,0);
v1 = subs(subs(diff(G2,un,2),un,0),um,0)/m1;
v2 = subs(subs(diff(G2,um,2),un,0),um,0)/m2;
cov = simplify(subs(subs(diff(diff(G2,um),un),un,0),um,0));

fprintf('m1: %s\n',m1)
fprintf('m2: %s\n',m2)
fprintf('v1: %s\n',v1)
fprintf('v2: %s\n',v2)
fprintf('cov: %s\n',cov)

%%%%%%%%%%%%%%%%%%%%

fprintf('\nextrinsic\n')

syms alf 
assume(alf,'positive');
G3 = -alf*log(1-Uint);

m1 = subs(subs(diff(G3,un),un,0),um,0);
m2 = subs(subs(diff(G3,um),un,0),um,0);
v1 = simplify(subs(subs(diff(G3,un,2),un,0),um,0)/m1);
v2 = simplify(subs(subs(diff(G3,um,2),un,0),um,0)/m2);
cov = simplify(subs(subs(diff(diff(G3,um),un),un,0),um,0));


fprintf('m1: %s\n',m1)
fprintf('m2: %s\n',m2)
fprintf('v1: %s\n',v1)
fprintf('v2: %s\n',v2)
fprintf('cov: %s\n',cov)

syms M1 V1
assume(M1,'positive');
assume(V1,'positive')
% eqn = m1^2/(m1*(1+v1) - m1)==R;
disp('MoM solution:')
disp(solve(M1==m1,V1==m1*(1+v1),'ReturnConditions',true))

%%%%%%%%%%%%%%%%%%%%

fprintf('\nCIR\n')

G4 = 1/2*(1-sqrt(1-4*b*U));
m1int = int(subs(subs(diff(G4,un),un,0),um,0));
m1 = simplify(subs(m1int,s,inf)-subs(m1int,s,0));

m2int = int(subs(subs(diff(G4,um),un,0),um,0));
m2 = simplify(subs(m2int,s,inf)-subs(m2int,s,0));

v1int = expand(int(subs(subs(diff(G4,un,2),un,0),um,0))) ;
v1 = simplify((subs(v1int,s,inf)-subs(v1int,s,0))/m1);

v2int = expand(int(subs(subs(diff(G4,um,2),un,0),um,0)));
v2 = simplify((subs(v2int,s,inf)-subs(v2int,s,0))/m2);

covint = expand(int(subs(subs(diff(diff(G4,um),un),un,0),um,0)));
cov_ = simplify((subs(covint,s,inf)-subs(covint,s,0))); %Poisson

fprintf('m1: %s\n',m1)
fprintf('m2: %s\n',m2)
fprintf('v1: %s\n',v1)
fprintf('v2: %s\n',v2)
fprintf('cov: %s\n',cov_)

%%%%%%%%%%%%%%%%%%%%

fprintf('\nDelayed splicing\n')

G5 = tau*b*xn/(1-b*xn) - 1/gam * log(1-b*xm);

m1 = subs(subs(diff(G5,un),un,0),um,0);
m2 = subs(subs(diff(G5,um),un,0),um,0);
v1 = simplify(subs(subs(diff(G5,un,2),un,0),um,0)/m1);
v2 = simplify(subs(subs(diff(G5,um,2),un,0),um,0)/m2);
cov = simplify(subs(subs(diff(diff(G5,um),un),un,0),um,0));


fprintf('m1: %s\n',m1)
fprintf('m2: %s\n',m2)
fprintf('v1: %s\n',v1)
fprintf('v2: %s\n',v2)
fprintf('cov: %s\n',cov)

syms M1 V1
assume(M1,'positive');
assume(V1,'positive')
% eqn = m1^2/(m1*(1+v1) - m1)==R;
disp('MoM solution:')
disp(solve(M1==m1,V1==m1*(1+v1),'ReturnConditions',true))

%%%%%%%%%%%%%%%%%%%%

fprintf('\nDelayed degradation\n')

Up = xm+(xn-xm)*exp(-bet*tau);
% G6 =  1/bet / (1-b*um) * log((b*Up-1)/(b*un-1));
G6 = tau*b*xm/(1-b*xm) - 1/bet * log(1-b*Up) + 1/bet / (1-b*xm) * log((b*Up-1)/(b*xn-1));
% G6 = exp(G6);
subs(subs(G6,un,0),um,0);
m1 = subs(subs(diff(G6,un),un,0),um,0);
m2 = subs(subs(diff(G6,um),un,0),um,0);
v1 = simplify(subs(subs(diff(G6,un,2),un,0),um,0)/m1);
v2 = simplify(subs(subs(diff(G6,um,2),un,0),um,0)/m2);
cov = simplify(subs(subs(diff(diff(G6,um),un),un,0),um,0));


fprintf('m1: %s\n',m1)
fprintf('m2: %s\n',m2)
fprintf('v1: %s\n',v1)
fprintf('v2: %s\n',v2)
fprintf('cov: %s\n',cov)

syms M1 V1
assume(M1,'positive');
assume(V1,'positive')
% eqn = m1^2/(m1*(1+v1) - m1)==R;
disp('MoM solution:')
disp(solve(M1==m1,V1==m1*(1+v1),'ReturnConditions',true))



return
