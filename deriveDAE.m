% Symbolic derivation of equations of motion for a Rocket 

clear

syms xG  yG   theta phi  phidot   real
syms vGx vGy  omega         real
syms aGx aGy  alpha         real
syms d m IG                 real  
syms N  T  g                real

i = [1 0 0]'; j=[0 1 0]'; k=cross(i,j);  % standard unit vectors
body_lambda = cos(theta)*i + sin(theta)*j;    % unit vector along   Rocket
body_normal = cross(k,body_lambda);                % unit vector perp  to sleigh  

fin_lambda = cos(theta+phi)*i + sin(theta+phi)*j;    % unit vector along wing   
fin_normal = -sin(theta+phi)*i + cos(theta+phi)*j;   % unit vector per to wing
rGC  = -d*body_lambda;       

vG   = vGx*i  +  vGy* j; 
aG   = aGx*i  +  aGy* j;
LMB = T*body_lambda + N*fin_normal -m*g*j -  m*(aGx*i + aGy*j);   

%AMB

%Torque of constraint force N
MG   = cross(rGC, N*fin_normal);
AMB=cross(d*body_lambda,aG)*m+ cross(d*body_lambda,-j*m*g) + IG*alpha*k;


%Constraint Equation. Need to express in terms of state variables
vC       = vG +  cross(omega*k,rGC);
% The constraint that the finn can't move perpendicular to itself
velconst = dot(fin_normal,vC);  

% We can only solve for the constraint by differentiating it.

dbydtvelconst =  diff(velconst,theta) * omega + ...
                 diff(velconst,vGx)   * aGx + ...
                 diff(velconst,vGy)   * aGy + ...
                 diff(velconst,omega) * alpha + ...
                 diff(velconst,phi)   * phidot;
dbydtvelconst = simplify(dbydtvelconst);  


% All of the stuff above gives 4 scalar equations in the variables of
% interest:  aGx, aGy, alpha, and N
eqn1 = dot(LMB, i);
eqn2 = dot(LMB, j); 
eqn3 = dot(AMB, k);
eqn4 = dbydtvelconst;   % Constraint equation in terms of accelerations

eqns = [eqn1, eqn2, eqn3, eqn4]; % Pack up the 4 equations in a vector
vars = [aGx, aGy, alpha, N];     % Pack up the 4 key unknowns in a vector

%TWO WAYS TO GET THESE EQUATIONS INTO A RHS FILE

%Method 1:  Use Solve (not good for systems with many DoF)
%[aGx, aGy, alpha, N] = solve(eqns,vars);
%    aGx   = simplify(aGx);
%    aGy   = simplify(aGy);
%    alpha = simplify(alpha);
%    N     = simplify(N);


%Method 2: Use equationsToMatrix (this is preferred for big problems)
[A,b] = equationsToMatrix(eqns,vars);  
    A  = simplify(A);  
    b  = simplify(b);

%Create two matlab functions that rhsrocket can call
matlabFunction(A,'file','Amatrix','Optimize',true);
matlabFunction(b,'file','bvector','Optimize',true)  ; % Eom are M*vars = b;
