% File Explaining the NSA Code for the Fat/Water Toolbox
% Author: Emily Bice 
% Email: emily.bice@gmail.com
% Advisor: Angel Pineda
% Email: pineda@fullerton.edu
% Updated 2012-01-19

% We assume that a set of input and output parameters are in the
% MATLAB workspace.  In our case, we used 'test_hernando_110818'.


% Add to matlab path
imDataParams.TE = [-0.4 1.2 2.8]*1e-3;
imDataParams.TE = [1.5500    3.8200    6.0900]*1e-3;
imDataParams.FieldStrength = 1.5;
imDataParams.PrecessionIsClockwise = 1;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
% $$$ algoParams.species(2).frequency = [-3.40];
% $$$ algoParams.species(2).relAmps = [1.0];

% For the first part of the NSA analysis, we look at the NSA as a function
% of fat fraction.  We assume a two-peak model where the fat and water are
% aligned at the echo at pi/4 radians.

% First we generate the NSA normalized to be between 0 and N points for
% the field map and no R2* model:
% $$$ 
% $$$ FMK_NSA = fw_nsa_2pluspoint(imDataParams, algoParams);
% $$$ % For testing, we have plots of the resulting NSA.
% $$$ ha=axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
% $$$ text(0.5,1,'Theoretical NSA for a range of fat fractions','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)
% $$$ 
% $$$ pause %pressing any key will make MATLAB continue
% $$$ close
% $$$ 
% $$$ % By adding arguments we can compute the NSA for the model with unknown
% $$$ % Field Map, given in the structure simParams:
simParams = [];
simParams.fieldmap = 110;
% $$$ 
% $$$ IDEAL_NSA = fw_nsa_2pluspoint(imDataParams, algoParams, simParams);
% $$$ ha=axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
% $$$ text(0.5,1,'Theoretical NSA for a range of fat fractions; \psi = 110','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)
% $$$ 
% $$$ pause
% $$$ close 

IDEAL_R2_NSA_MC = fw_mcsim_2pluspoint(imDataParams, algoParams, simParams, 11, 200, 100);
ha=axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
text(0.5,1,'Monte Carlo NSA for a range of fat fractions; \psi = 110, R_2^* = 50','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)


% And for the model including a single R2*:

% $$$ simParams.r2star = 50;
% $$$ 
% $$$ IDEAL_R2_NSA = fw_nsa_2pluspoint(imDataParams, algoParams, simParams);
% $$$ ha=axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
% $$$ text(0.5,1,'Theoretical NSA for a range of fat fractions; \psi = 110, R_2^* = 50','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)
% $$$ 
% $$$ pause
% $$$ close
% $$$ 
% $$$ %We also have a sample code which carries out a Monte Carlo using least
% $$$ %squares estimation to test the NSA.  Note that for testing purposes we
% $$$ %used less fat/water ratios and a small number of realizations.  With these
% $$$ %parameters, the Monte Carlo took 10 seconds in my laptop.
% $$$ 
% $$$ IDEAL_R2_NSA_MC = fw_mcsim_2pluspoint(imDataParams, algoParams, simParams, 11, 200, 100);
% $$$ ha=axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
% $$$ text(0.5,1,'Monte Carlo NSA for a range of fat fractions; \psi = 110, R_2^* = 50','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)
% $$$ 
% $$$ pause
% $$$ close
% $$$ 
