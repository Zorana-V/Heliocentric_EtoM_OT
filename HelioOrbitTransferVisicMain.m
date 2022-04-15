clc; clear;

%-------------------------------------------------------%
% Helicenric Portion of Earth to Mars Orbit Transfer    %
%-------------------------------------------------------%
% Objective: Find a solution for the orbit transfer     %
%   between Earth and Mars during the heliocent         %
%   phase.                                              %
% Constraints:                                          %
%   Diff Eqns of Motion:                                %
%       dr/dt = v_r                                     %
%       dtheta/dt = (v_theta)/r                         %
%       dv_r/dt = (v_theta^2)/r + (T/m)sinBeta - mu/r^2 %
%       dv_theta/dt = -(v_r*v_theta)/r + (T/m)cosBeta   %
%       dm/dt = -(T/v_theta)                            %
%   Bounds:                                             %
%       (r(0),r(tf)) = (r0,rf)                          %
%       (theta(0),theta(tf)) = (theta0,free)            %
%       (v_r(0),v_r(tf)) = (v_r0,v_rf)                  %
%       (v_theta(0),v_theta(tf)) = (v_theta0,v_theta)   %
%       (m(0),m(tf)) = (m0,mf)                          %
% ------------------------------------------------------%

%-------------------------------------------------------%
% Defining Initial Conditions, Universal Constants, and %
% Boundary Conditions                                   %
%-------------------------------------------------------%

muSun = 1;            %gravity parameter of the Sun [km^3 * s^-1]
T               = 0.1405;
v_e             = 1.87583;          %escape velocity


% AUXDATA = Auxiliary Data [Structure]
auxdata.muSun   = muSun;
auxdata.v_e     = v_e;

% BOUNDS ON VARIABLES
t0              = 0;                %initial time
r0              = 1;                %initial radius (Earth Radius) [AU]
theta0          = 0;                %initial theta
v_r0            = 0;                %initial velocity in the r direction
v_theta0        = sqrt(muSun/r0);   %initial velocity in the theta direction
m0              = 1;                %initial mass
rf              = 1.5;              %final radius (Mars Radius) [AU]
v_rf            = 0;                %final velocity in the r direction
v_thetaf        = sqrt(muSun/rf);   %final velocity in the theta direction

tfmin       = t0;       %lower tolerance on final time
tfmax       = +10;      %upper tolerance on final time
rfmin       = r0;       %lower tolerance on final radius
rfmax       = +10;      %upper tolerance on final radius
v_rmin      = -10;      %lower tolerance on v_r
v_rmax      = +10;      %upper tolerance on v_r
v_thetamin  = -10;      %lower tolerance on v_theta
v_thetamax  = +10;      %upper tolerance on v_theta
thetamin    = theta0;   %lower tolerance on angle theta
thetamax    = +4*pi;    %upper tolerance on angle theta
mmin        = 0.1;      %lower tolerance on mass
mmax        = m0;       %upper tolerance on mass
betamin     = (-2*pi);  %lower tolerance on beta
betamax     = 2*pi;     %upper tolerance on beta
w_rmin      = -10;      %lower tolerance on w_r
w_rmax      = +10; 
w_thetamin  = -10;
w_thetamax  = +10;
Tmin        = 0;
Tmax        = T;

% BOUNDS [STRUCTURE]
bounds.phase.initialtime.lower      = t0;
bounds.phase.initialtime.upper      = t0;
bounds.phase.finaltime.lower        = tfmin;
bounds.phase.finaltime.upper        = tfmax;
bounds.phase.initialstate.lower     = [r0,theta0,v_r0,v_theta0,m0];
bounds.phase.initialstate.upper     = [r0,theta0,v_r0,v_theta0,m0];
bounds.phase.state.lower            = [rfmin,thetamin,v_rmin,v_thetamin,mmin];
bounds.phase.state.upper            = [rfmax,thetamax,v_rmax,v_thetamax,mmax];
bounds.phase.finalstate.lower       = [rf,thetamin,v_rf,v_thetaf,mmin];
bounds.phase.finalstate.upper       = [rf,thetamax,v_rf,v_thetaf,mmax];
%bounds.phase.control.lower          = [betamin];
%bounds.phase.control.upper          = [betamax];
bounds.phase.control.lower          = [w_rmin,w_thetamin,Tmin];
bounds.phase.control.upper          = [w_rmax,w_thetamax,Tmax];

bounds.phase.path.lower             = 1;
bounds.phase.path.upper             = 1;

% GUESS [STRUCTURE]
tGuess              = [t0; tfmax];
rGuess              = [r0; rf];
thetaGuess          = [theta0; thetamax];
v_rGuess            = [v_r0; v_rf];
v_thetaGuess        = [v_theta0; v_thetaf];
mGuess              = [m0; mmax];
%betaGuess           = [betamin; betamax];
w_rGuess            = [1; 1];
w_thetaGuess        = [0; 0];
TGuess              = [Tmax; Tmax];

guess.phase.time    = tGuess;
guess.phase.state   = [rGuess,thetaGuess,v_rGuess,v_thetaGuess,mGuess];
guess.phase.control = [w_rGuess,w_thetaGuess,TGuess];

% MESH [STRUCTURE]
numIntervals         = 10;
mesh.method          = 'hp-LiuRao-Legendre';
mesh.tolerance       = 1e-8;
mesh.maxiterations   = 2;
mesh.colpointsmin    = 3;
mesh.colpointsmax    = 10;
mesh.phase.colpoints = 3*ones(1,numIntervals);
mesh.phase.fraction  = ones(1,numIntervals)/numIntervals;

% SETUP [STRUCTURE]
setup.name                           = 'Heliocentric-Orbit-Transfer-Problem';
setup.functions.continuous           = @HelioOrbitTransferVisicContinuous;
setup.functions.endpoint             = @HelioOrbitTransferVisicEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'none';
setup.method                         = 'RPM-Differentiation';

% CALL GPOPS-II
output = gpops2(setup); %function to use gpops2 w Setup Structure

%finding vectors that represent the orbit
solution    = output.result.solution;
t           = solution.phase.time(:,1);
r           = solution.phase.state(:,1);
theta       = solution.phase.state(:,2);
v_r         = solution.phase.state(:,3);
v_theta     = solution.phase.state(:,4);
m           = solution.phase.state(:,5);
w_r         = solution.phase.control(:,1);
w_theta     = solution.phase.control(:,2);
T           = solution.phase.control(:,3);

beta        = atan2(w_theta,w_r);
iT          = find(T>0.01);
thetaiT     = theta(iT);
riT         = r(iT);

figure(1)
plot(t,r,t,theta,t,v_r,t,v_theta,t,T)

figure(2)
plot(t,w_r,t,w_theta)

figure(3)
plot(t,beta*180/pi)

figure(4)
polarplot(theta,r,thetaiT,riT,'-s')

figure(5)
plot(t,T)





