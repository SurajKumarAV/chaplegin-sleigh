function zdot  = rhsdynamics(t,z,p)
% Right Hand Side calculation of ODES: zdot = f(t,z)

% Unpack state z
xG  = z(1); yG =z(2);  theta =z(3);
vGx = z(4); vGy=z(5);  omega =z(6);
phi = z(7);

% control here
u = controls(t,z,p);
phidot = u.phidot;
T = u.T;

% Unpack parameters
m = p.m; IG = p.IG;  d = p.d; g = p.g;

%rhs evaluation
xGdot     = vGx ;  % The first three ODEs are easy
yGdot     = vGy ;
thetadot  = omega ;
posdots   = [xGdot;yGdot;thetadot];  %store in one 3-element vector

A = Amatrix(IG,d,m,phi,theta);
b = bvector(T,d,g,m,omega,phi,phidot,theta,vGx,vGy);

% Solve for the acceleration terms and the constraint force
w = A\b;  

veldots = w(1:3);   
phidots = phidot;
%pack up
zdot= [posdots; veldots;phidots];  % 6-element right hand side vector
end

