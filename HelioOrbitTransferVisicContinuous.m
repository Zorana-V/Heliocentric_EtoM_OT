function output = HelioOrbitTransferVisicContinuous(input)

%-------------------------------------------------------%
% Input:                                                %
%   auxdata - a structure containing the following:     %
%       time                                            %
%       state                                           %
%       control                                         %
%                                                       %
% Output:                                               %
%   output - a structure containing the following:      %
%       dynamics                                        %
%       path                                            %
%       integrand                                       %
%-------------------------------------------------------%

muSun     = input.auxdata.muSun; %calling the gparameter for the sun
v_e       = input.auxdata.v_e;

%calling the initial guesses for each of the following:
t           = input.phase.time(:,1);
r           = input.phase.state(:,1);
theta       = input.phase.state(:,2);
v_r         = input.phase.state(:,3);
v_theta     = input.phase.state(:,4);
m           = input.phase.state(:,5);
%beta        = input.phase.control(:,1);
w_r         = input.phase.control(:,1);
w_theta     = input.phase.control(:,2);
T           = input.phase.control(:,3);

%calculating the values of the diff. equations:
rdot            = v_r;
thetadot        = (v_theta./r);
v_rdot          = (((v_theta).^2)./r) + ((T./m).*w_theta) - (muSun./(r.^2));
v_thetadot      = (-(v_r.*v_theta)./r) + ((T./m).*w_r);
mdot            = (-T./v_e);

%updated guess as function is running
output(1).dynamics  = [rdot,thetadot,v_rdot,v_thetadot,mdot];
output(1).path      = w_r.^2 + w_theta.^2;