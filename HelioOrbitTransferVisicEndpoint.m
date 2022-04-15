function output = HelioOrbitTransferVisicEndpoint(input)

%-------------------------------------------------------%
% Input:                                                %
%   auxdata                                             %
%   input.phase - a structure containing:               %
%       initialtime                                     %
%       initialstate                                    %
%       finaltime                                       %
%       finalstate                                      %
%                                                       %
% Output:                                               %
%   output - a structure containing the following:      %
%       objective                                       %
%       event                                           %
%-------------------------------------------------------%

%computing the minimal time required
mf               = input.phase.finalstate(end,5);
output.objective = -mf;