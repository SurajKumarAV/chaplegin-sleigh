clear; clf
close all
deriveDAE  % Derives the equations of motion using deriveDAE.m

% ANIMATION PARAMETERS 
dur           = 10;                % time duration of simulation (world time)
speedupfactor = 1;                % for animating faster or slower
tspan         = linspace(0,dur);  % Only two numbers, duration in world time

% PROBLEM SETUP
% PARAMETERS
% Physical parameters 
d=1;        
m=1;        
IG=1;  
g = 0;
L = 2;  % This is only for graphics, its the drawn length of sleigh
p.d = d;   
p.m  = m;   
p.IG = IG;   % pack the parameters into struct p
p.g = g;

%Useful unit vectors for setting up initial conditions:
uniti  = [ 1 0 0]';  unitj = [ 0 1 0]';  unitk = [ 0 0 1]'; 
save = 0;

% INITIAL CONDITIONS

rG0     = [ 0 0 ]';         % initial position of G, x and y coordinates
vC0     = 1;                % initial speed of point C
theta0  = 0*pi/180;         % initial angle of sleigh
omega0  = 0;                % initial angular velocity of sleigh
phi0    = 0*pi/180;         % initial fin angle

% The things below are created; by the numbers above.
% The velocities of C and G are calculated to satisfy
% the constraint that vC is parallel to lambda
fin_lambda0 = cos(theta0+phi0)*uniti + sin(theta0+phi0)*unitj;
lambda0 = cos(theta0) * uniti + sin(theta0) * unitj;
vC0     = vC0 * fin_lambda0; % make vC0 a vector in the right direction
%Next, find velocity of G given velocity of C, theta0 and omega0
vG0     = vC0 + cross((omega0 * unitk), d*lambda0);
vG0     = vG0(1:2);  %  Only want a 2D vector, not 3D

% Put all the stuff above into a single 6-element column vector 
z0 = [rG0; theta0;   vG0; omega0; phi0];  %3 configs,3 velocities and fin angle

% SOLVE ODEs
small = 1e-4;  
options  = odeset('RelTol', small , 'AbsTol', small );
f        =  @(t,z)   rhsdynamics(t,z,p);
%SOLVE THE ODES !!!  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soln     = ode45(f,tspan, z0,options);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PLOTS and ANIMATION

%Make plot of the trajectory of point C
tspan = linspace(0,dur,100);  % Times at which we want the solution for plot
z     = deval(soln,tspan);    % Each column of z is a solution at
                              % the corresponding element of tspan
xG = z(1,:); yG = z(2,:);  theta = z(3,:); %configuration variables.

xC = xG - cos(theta) * d;
yC = yG - sin(theta) * d;

if (save)
    writerObj = VideoWriter('animation.avi');
    writerObj.FrameRate = 30;
    open(writerObj);
end

figure;
plot(xC,yC, '-b','LineWidth', 3)
title('Trajectory, curve shows motion of point C')
xlabel({'$x_C$'; 'Note rocket is always tangent to curve at C'}); 
ylabel('$y_C$')

xmax  = max(soln.y(1,:)); ymax  = max(soln.y(2,:)); %
xmin  = min(soln.y(1,:)); ymin  = min(soln.y(2,:));

axis equal  
axis([xmin-L, xmax+L, ymin-L, ymax+L]) 
hold on  
% ANIMATION
% Animation starts with the first frame
z = deval(soln,0);  % Evaluate solution at t = 0

%Unpack z at t=0
xG0     = z(1);
yG0     = z(2);
theta0  = z(3);
phi0    = z(7);

lambda = cos(theta0) * uniti + sin(theta0) * unitj;
rGC    = -d*lambda(1:2);      % pos of C relative to G
rGE    = (L-d)*lambda(1:2);   % pos of E relative to G

xC0    = xG0 + rGC(1);     % initial position of C
yC0    = yG0 + rGC(2);

xE0    = xG0 + rGE(1);     % initial position of C
yE0    = yG0 + rGE(2);

%Plot line segment from C to E (full length of sleigh)
x =  [xC0, xE0  ];
y =  [yC0, yE0  ];
xf0 = [xC0 + (d/4)*cos(theta0+phi0), xC0 - (d/4)*cos(theta0+phi0)];
yf0 = [yC0 + (d/4)*sin(theta0+phi0), yC0 - (d/4)*sin(theta0+phi0)];
h1 = plot(x,y,'LineWidth',10); % create plot object (struct)
h2 = text(xC0,yC0, 'C','FontSize', 15);  % Put a 'C' at the fin
h3 = text(xG0,yG0, 'G','FontSize', 15);  % Put a 'G' at the center of mass
h4 = plot(xf0,yf0,'LineWidth',2);
shg  

%Animate  Plot
tic     
t = 0;  
while t<dur    
    z = deval(soln,t);  % Evaluates ODE solution at time t
    xG = z(1); yG = z(2); theta = z(3); phi = z(7);
    lambda = cos(theta) * uniti + sin(theta) * unitj;
    rGC    = -d*lambda;
    xC     = xG + rGC(1);
    yC     = yG + rGC(2);
    rGE    = (L-d)*lambda;
    xE     = xG + rGE(1);
    yE     = yG + rGE(2);
    x   =  [xC, xE];
    y   =  [yC, yE];
    xf = [xC + (d/4)*cos(theta+phi), xC - (d/4)*cos(theta+phi)];
    yf = [yC + (d/4)*sin(theta+phi), yC - (d/4)*sin(theta+phi)];
    h1.XData = x;  h1.YData = y; 
    h4.XData = xf; h4.YData = yf;
    h2.Position = [xC-.05 yC 0]; 
    h3.Position = [xG-.05 yG 0];
    drawnow 
    if (save) 
        im = frame2im(getframe(gcf));
        writeVideo(writerObj,im);    
    end
    t = toc * speedupfactor; 
end
if (save)
    close(writerObj);
end
hold off


