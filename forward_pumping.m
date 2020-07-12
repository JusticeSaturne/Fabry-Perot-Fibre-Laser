% This script implements a relaxation algorithm to solve the propagation
% equations of a Fabry-Perot erbium doped fiber laser, in the forward pumping
% configuration.
% The cavity of the laser is divided into 100 sections and power in each 
% section is computed using a fourth-order Runge-Kutta shcheme.
%  Default values of parameters
%   - Pump power: 100 mW
%   - Doped fibre length: 10 m
%   - Number of sections: 100
%   - Output mirror reflectivity: 10%
%
% Comments:
%   The default values can be changed by the user to characterise the fiber
%   laser
%
% Rerences:
%
% Written by Justice Sompo, University of Johannesburg, South-Africa
clear all
close all
clc
format longe
tic
%
%------------------------COMPUTATION PARAMETERS----------------------------
FiberLength = 10;           % meters
n = 100;                    % Number of sections
step = (FiberLength-0)/n;   % step size for computation
R1 = 0.98;                  % Reflectivity of mirror 1
R2 = 0.1;                   % Reflectivity of mirror 2

%==========================================================================
%                     TRIP 1 : FORWARD PROPAGARTION
%==========================================================================
%----------------------------INITIAL CONDITIONS----------------------------
ppf = 60e-3;               % Forward launched power at z = 0 (100 mW)
z_initial = FiberLength;    % origin of longitudinal distance (z=L for backward displacement)
pp_initial = 15;            % Initial guess of the pump power at z = L 
psf_initial = 1e-3;         % initial guess of the forward propagating signal at z = L
psb_initial = 10e-3;        % initial guess of the backward propagating signal at z = L
z = z_initial;
pp = pp_initial;
psf = psf_initial;
psb = psb_initial;

%==========================================================================
%         WHILE LOOP TO COMPUTE LASER PARAMETERS
%==========================================================================
%
err1 = 10;                % Arbitrarily set value to enter the while loop
err2 = 10;
err3 = 10;
tol = 1e-4;
IterationNumber = 0;       % Initialise the counter

while (err1 >= tol && err2 >= tol || err3 >= tol) && IterationNumber <= 40
    
    %----------------------CREATE VECTORS FOR PLOTTING---------------------
    Z = z_initial;
    PP = pp_initial;
    PSF = psf_initial;
    PSB = psb_initial;
    %----------------Iterative calculations using RK4 method---------------
    p = n-1 ;                  % Integration from z = L to z = 0
    for w = p:-1:0             % Number of loops
        z = z-step;            % Doing Reverse Integration from "L" to Zero(0)
        
        k21 = step*feval('PPF',z,pp,psf,psb);
        f21 = step*feval('PSF',z,pp,psf,psb);
        b21 = step*feval('PSB',z,pp,psf,psb);
        
        k22 = step*feval('PPF',(z+0.5*step),(pp+0.5*k21), (psf+0.5*f21),(psb+0.5*b21));
        f22 = step*feval('PSF',z+0.5*step,pp+0.5*k21, psf+0.5*f21, psb+0.5*b21);
        b22 = step*feval('PSB',z+0.5*step,pp+0.5*k21, psf+0.5*f21, psb+0.5*b21);
        
        k23 = step*feval('PPF',(z+0.5*step),(pp+0.5*k22), (psf+0.5*f22), (psb+0.5*b22));
        f23 = step*feval('PSF',(z+0.5*step),(pp+0.5*k22), (psf+0.5*f22), (psb+0.5*b22));
        b23 = step*feval('PSB',(z+0.5*step),(pp+0.5*k22), (psf+0.5*f22), (psb+0.5*b22));
        
        k24 = step*feval('PPF',z+step,pp+k23,psf+f23,psb+b23);
        f24 = step*feval('PSF',z+step,pp+k23,psf+f23,psb+b23);
        b24 = step*feval('PSB',z+step,pp+k23,psf+f23,psb+b23);
        
        pp = pp-(1/6)*(k21+2*(k22+k23)+k24);      % Pump power evolution
        psf = psf -(1/6)*(f21+2*(f22+f23)+f24);   % Forward laser power evolution  
        psb = psb -(1/6)*(b21+2*(b22+b23)+b24);   % Backward laser power evolution
        
        Z = [Z,z];
        PP = [PP,pp];
        PSF = [PSF, psf];
        PSB = [PSB, psb];
    end
    %------------STORE THE VALUES OF PUMP, SIGNAL AND POSITION-------------
    z_FirstPass = z;               
    pp_FirstPass = pp;
    psf_FirstPass = psf;
    psb_FirstPass = psb;
    %----------------TEST CONVERGENCE AFTER THE FIRST PASS-----------------
    
    err1 = abs(psf-R1*psb);
    err2 = abs(pp-ppf);
    if err1 <=tol && err2 <=tol
        break;
    end
    
    %======================================================================
    %                     TRIP 2 : BACKWARD PROPAGARTION
    %======================================================================
    %-------------------APPLY BC FOR NEW INITIAL CONDITIONS----------------
    z = 0;
    pp = ppf;
    psf = psb_FirstPass*R1;
    psb = psb_FirstPass;
    %--------------------CREATE VECTORS FOR PLOTTING-----------------------
    Z = z;
    PP = pp;
    PSF = psf;
    PSB = psb;
    for w = 0:n-1           
        z = z+step;        % Integration is from z = 0 to z = L
        
        k31 = step*feval('PPF',z,pp,psf,psb);
        f31 = step*feval('PSF',z,pp,psf,psb);
        b31 = step*feval('PSB',z,pp,psf,psb);
        
        k32 = step*feval('PPF',(z+0.5*step),(pp+0.5*k31), (psf+0.5*31),(psb+0.5*b31));
        f32 = step*feval('PSF',z+0.5*step,pp+0.5*k31, psf+0.5*f31, psb+0.5*b31);
        b32 = step*feval('PSB',z+0.5*step,pp+0.5*k31, psf+0.5*f31, psb+0.5*b31);
        
        k33 = step*feval('PPF',(z+0.5*step),(pp+0.5*k32), (psf+0.5*f32), (psb+0.5*b32));
        f33 = step*feval('PSF',(z+0.5*step),(pp+0.5*k32), (psf+0.5*f32), (psb+0.5*b32));
        b33 = step*feval('PSB',(z+0.5*step),(pp+0.5*k32), (psf+0.5*f32), (psb+0.5*b32));
        
        k34 = step*feval('PPF',z+step,pp+k33,psf+f33,psb+b33);
        f34 = step*feval('PSF',z+step,pp+k33,psf+f33,psb+b33);
        b34 = step*feval('PSB',z+step,pp+k33,psf+f33,psb+b33);
        
        pp = pp+(1/6)*(k31+2*(k32+k33)+k34);      % Pump power evolution
        psf = psf +(1/6)*(f31+2*(f32+f33)+f34);   % Forward laser power evolution
        psb = psb +(1/6)*(b31+2*(b32+b33)+b34);   % Backward laser power evolution
        
        Z = [Z,z];
        PP = [PP,pp];
        PSF = [PSF, psf];
        PSB = [PSB, psb];
    end
    %--------------STORE THE VALUES OF PUMP, SIGNAL AND POSITION-----------
    z_SecondPass = z;
    pp_SecondPass = pp;
    psf_SecondPass = psf;
    psb_SecondPass = psb;  
    
    err3 = abs(psb-psf*R2); % Test convergence criteria
    z = z_SecondPass;
    pp = pp_SecondPass;
    psb = psf_SecondPass*R2;
    psf = psf_SecondPass;
    IterationNumber = (IterationNumber + 1) ;
end
%---------------------COMPUTE POPULATIONS DENSITY--------------------------
Pout = psf*(1-R2)
[N1,N2] = PopulationDensity(PP,PSB,PSF);
%==========================================================================
%                           PLOTTING
%==========================================================================
figure(1)
plot(Z,N1,'m','linewidth',2)
hold on
plot(Z,N2,'g','linewidth',2)
xlabel('Fiber length (m)') 
ylabel('Population density (ions/m^3)')
figure(2)
Pout = psf*(1-R2);
plot(Z,PSF*1000,'linewidth',2)
hold on
plot(Z,PSB*1000,'r','linewidth',2)
hold on
plot(Z,PP*1000,'m','linewidth',2)
xlabel('Fiber length (m)') 
ylabel('Power (mW))') 
toc;
