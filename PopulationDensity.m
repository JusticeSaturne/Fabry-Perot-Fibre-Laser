function [N1,N2] = PopulationDensity(Pp,Psb,Psf)
% PopulationDensity function represents the system of linear equations for
% population densty calculation.
% 
% Inputs:
% -The function takes the pump power and initial guesses of output powers
%  
% Outputs: 
%   N1: Population density of ground level as a function of distance z
%   N2: Population density of metastable level as a function of
%   distance z
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
if nargin <3
    error(message('MATLAB:PopulationDensity:NotEnoughInputs'));
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
Gammap = 0.81;         % Overlap factor at pump wavelength
Gammas = 0.6;         % Overlap factor at signal wavelength
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
F_s = cel/lambdas;
den1 = A*h*F_s;
F_p = cel/lambdap;  
den2 = A*h*F_p;
% Fiber Bragg Grating Parameters
W12 = ((sigma12S*Gammas).*(Psb+Psf))/den1;
W21 = ((sigma21S*Gammas).*(Psb+Psf))/den1;
R12 = ((sigma12P*Gammap).*(Pp))/den2;
%Rate equations solving
N1 = Ner*(W21+A21)./(W12+R12+W21+A21);
N2 = Ner*(W12+R12)./(W21+A21+W12+R12);