function PB = PSB(x,y,uf,ub)
% PSF function represents backward laser field propagation equation
% 
% Inputs:
% - z : Position along the longitudinal distance of the cavity
% - pp: pump power
% - psf: Forward laser power 
% - psb: Backward laser power
%  
% Outputs: 
%   PB : Backward signal power as a function of distance z
%  
% Comments:
%   - Computation is done in the ideal case of a two level energy system
%
% References: 
%  
% Written by Justice Sompo, University of Johannesburg, South Africa    
close all
clc
% 
format longe
% Check inputs
if nargin <4
    error(message('MATLAB:PSB:NotEnoughInputs'));
end
% Rate equation Parameters
sigma12S = 2.33e-25;   % Laser absorption cross section
sigma21S = 2.64e-25;   % Laser emission cross section
sigma12P = 2e-25;      % Pump absorption cross section
sigma21P = 0;          % Pump emission cross section
% Lifetimes
t21 = 10e-3;           % Lifetime of the metastable energy level of Erbium
A21 = 1/t21;
% Concentrations
Ner = 1.2e25;          % Erbium ions concentration   
% Overlap factors and background losses
Gammap = 0.64;         % Overlap factor at pump wavelength
Gammas = 0.43;         % Overlap factor at signal wavelength
alphap = 0.005;        % Background loss at pump wavelength
alphas = 0.005;        % Background loss at laser wavelength
lambdap = 980e-9; % Pumb wavelenght
lambdas = 1550e-9; %signal wavelength 
% Other physical constants
r = 2.3e-6;           % Doped fibre core
h = 6.626e-34;        % Plank's constant
cel = 3e8;            % speed of light in free space
% Calculated parameters
A = pi*r^2;  
F_s = cel/lambdas;den1 = A*h*F_s;
F_p = cel/lambdap;  den2 = A*h*F_p;
% Fiber Bragg Grating Parameters
W12 = ((sigma12S*Gammas)*(uf+ub))/den1;
W21 = ((sigma21S*Gammas)*(uf+ub))/den1;
R12 = ((sigma12P*Gammap)*(y))/den2;
%Rate equations solving
N1 = Ner*(W21+A21)./(W12+R12+W21+A21);
N2 = Ner*(W12+R12)./(W21+A21+W12+R12);
PB= (-1)*(ub.*Gammas).*(sigma21S*N2-sigma12S*N1)+(alphas.*ub);
