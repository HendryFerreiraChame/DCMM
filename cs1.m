##########################################################################################################
#
#   <This program implements the Dynamic Computational Model of Motivation (DCMM)
#   which is a conectionist model of motivation based on SDT and HMIEM. The framework
#   consists on interconected Continuos Attractor Neural Netoworks (CANN)>
#
#   Copyright (C) <2018>  <Hendry Ferreira Chame>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, version 3.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#    
#   Programmed by: Hendry Ferreira Chame
#   Institution: Federal University of Rio Grande (FURG), Rio Grande, BRAZIL
#   E-mail: hendrychame@furg.br, hendryf@gmail.com
#     
#   Paper:  A dynamic computational model of motivation based on self-determination theory and CANN
#   Authors: Hendry Ferreira Chame, Fernanda Pinto Mota, Silvia Silva da Costa Botelho
#   Information Sciences, Elsevier, 2018. https://doi.org/10.1016/j.ins.2018.09.055
#
#   Description:  this program implenets Case Study 1 (CS1). 
#                 The plots generated are in the folder 'plots' 
#
##########################################################################################################
close all;
clear all;
clc;

display("#################################################################################################");
display("")
display("  <GNU Octave Implementation of DCMM>  Copyright (C) <2018>  <Hendry Ferreira Chame>          ");
display("  This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome    "); 
display("  to redistribute it under certain conditions; for details consult the file 'LICENSE.txt'.      ");
display("")
display("  For more details about the framework, please consult the following reference:                 "); 
display("")
display(" \"A dynamic computational model of motivation based on self-determination theory and CANN\"    ");  
display("  Hendry Ferreira Chame, Fernanda Pinto Mota, Silvia Silva da Costa Botelho                     ");
display("  Information Sciences, Elsevier, 2018                                                          ");
display("")
display("#################################################################################################");
display("")

display ("DCMM program stared!")
display ("Begin of Case study 1 (CS1)...");
fflush(stdout); 

%-------------------------------------- Menu for similation selection -------------------------------------- 
display ("Program menu"); 
c = menu("Please select an option", "1) Generate all plots", "2) Generate CS1-A", "3) Generate CS1-B", "4) Generate CS1-C","5) Exit the program"); 

%% For bypassing the menu, uncomment the line below
%c = 1

CS1 = [0 0 0];

switch (c)
  case 1  
    CS1(1:3) = 1;
    display ("Option selected: 1) Generate all plots");
  case 2  
    CS1(1) = 1;
    display ("Option selected: 2) Generate CS1-A"); 
  case 3  
    CS1(2) = 1;
    display ("Option selected: 3) Generate CS1-B"); 
  case 4  
    CS1(3) = 1;
    display ("Option selected: 4) Generate CS1-C");
  case 5
     display ("Option selected: 5) Exiting program");
  otherwise 
    display("Invalid option selected");
endswitch

fflush(stdout); 
if ~exist("plots")&&(c<5)
    
    display("Creating [plots] dir ...");
    mkdir "plots";
    
  endif


display("Configuring shared parameters ...");
fflush(stdout); 

%% Shared simulation parameters

motNeurons = 6;     % Number motivation layer units
nContext = 4;     % Life context simulated (1: Education, 2: Interpersonal relations, 3: Leissure, 4: Null)
nMediators = 3;   % Number of mediator processes (1: Autonomy, 2: Competence, 3: Relatedness)
medNeurons = 11; % Number of mediator layer units
  
sigmaCann = 0.3;  % Sigma CANN gaussian activation
sigmaMed = 0.5;   % Sigma Mediator gaussian activation

%% DCMM parameter setting

setWeightParameters(sigmaCann, sigmaMed, motNeurons, medNeurons);
 
load("parameters/WAuToM.m");
load("parameters/WCoToM.m");
load("parameters/WReToM.m");
load("parameters/WCANN.m");

  
%-------------------------------------- Case study 1-A -------------------------------------- 
if CS1(1) > 0
  display("Simulating CS1-A ... " );
  fflush(stdout); 

  dt = 5.0; # in sec
  
  % Parameters for situational level 
  tS = 100*dt;
  eS = 0.1;
  g1 = g2 = g3 = g4 = 1;
  
  % simulation time
  simTimeMin = 120; 
  timeSec = simTimeMin*60;
  n = timeSec*(1/dt); 

  %% setting Mediator's activation 
   
  step =  (motNeurons-1)/(medNeurons-1);
  preferedUnitValue = [1:step:motNeurons]';  % preferred value encoded by the units
  sigmaNoise = 0.0;                        % gausian Noise for mediator activation
  
  %% --- situational Mediators
  
  % linear mapping from the motivation layer units to the mediation layer units
  step = (motNeurons-1)/(n-1);
  actFun = [1:step:motNeurons];                       % ramp activation function
  actFunForNSteps = ones(medNeurons,n).*actFun;
  
  % firing rate for f mediator neurons (set iqual to a ramp function)
  f_As = f_Cs = f_Rs = generateGaussian(sigmaMed,actFunForNSteps,preferedUnitValue);
  
  %% Setting initial states
  
  % Activation Ms
  h_Ms = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Ms_all = [];
  
  % Activation Mc
  h_Mc = ones(motNeurons,1).*(1/motNeurons); % uniform distribution

  % Applying the inhibition factors
  W_Ms = WCANN .+ eS;
  
  % TNL connections and activation
  WTNLtoAs = WTNLtoCs = WTNLtoRs = zeros(medNeurons,nMediators);
  % unbiased contribution from mediators
  WTNLtoAs(:,1) =  ones(medNeurons,1);
  WTNLtoCs(:,2) =  ones(medNeurons,1);
  WTNLtoRs(:,3) =  ones(medNeurons,1);
  
  f_TNL = ones(nMediators,1); 
  
  %% Simulation begin 
  for i = 1 : n
    
    % Situational level

    h_Ms_all = [h_Ms_all h_Ms];
    f_Ms  = h_Ms;
    
    h_As  = WTNLtoAs*f_TNL.*f_As(:,i);
    h_Cs  = WTNLtoCs*f_TNL.*f_Cs(:,i);
    h_Rs  = WTNLtoRs*f_TNL.*f_Rs(:,i);
       
    h_Ms = (1-(dt/tS))*h_Ms + (dt/tS)*(W_Ms*f_Ms + g1*WAuToM*h_As + g2*WCoToM*h_Cs + g3*WReToM*h_Rs + g4*h_Mc);
    
    h_Ms = h_Ms./sum(h_Ms,1);
     
  endfor
  %% Simulation end
  
  display("Plotting CS1-A ... " );
  fflush(stdout); 
  
  % begin Ploting
  
  % Mediatiors activation plot
  
  figure(1);
  times = (dt.*[0:n-1])/60;
  hold on;
  
  for i = 1 : medNeurons
    plot (times, f_As(i,:), 'b');  
  endfor 
  
  hx = xlabel("Time in min");
  hy = ylabel("Activation");
  ht = title("Time evolution of mediator's firing rate");
  
  print -deps -color plots/CS1a_ActFunMed.eps

  % saving plot data
  plotMat = [times; f_As];
  csvwrite("csv/CS1a_ActFunMed.csv", plotMat);
  
  % Situational motivation plot
  
  figure(2);
  times = (dt.*[0:n-1])/60;
  hold on;

  plot (times, h_Ms_all(1,:), 'r');
  plot (times, h_Ms_all(2,:), 'g');
  plot (times, h_Ms_all(3,:), 'k');
  plot (times, h_Ms_all(4,:), 'c');
  plot (times, h_Ms_all(5,:), 'm');
  plot (times, h_Ms_all(6,:), 'b');
    
  hx = xlabel("Time in min");
  hy = ylabel("Activation");
  ht = title("Time evolution of situational motivation");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  legend (h, "location", "northeastoutside");
  
  print -deps -color plots/CS1a_ActFunMs.eps
  
  % saving plot data
  plotMat = [times; h_Ms_all];
  csvwrite("csv/CS1a_ActFunMs.csv", plotMat);
  
endif 

%-------------------------------------- Case study 1-B -------------------------------------- 
 
if CS1(2) > 0
  display("Simulating CS1-B ... " );
  fflush(stdout); 
  
  dt = 5.0; # in sec
  
  % Parameters for situational level 
  tS = 100*dt;
  eS = 0.1;
  g1 = g2 = g3 = g4 = 1;
  
  % simulation time
  simTimeMin = 120; 
  timeSec = simTimeMin*60;
  n = timeSec*(1/dt); 

  %% setting Mediator's activation 
   
  step =  (motNeurons-1)/(medNeurons-1);
  preferedUnitValue = [1:step:motNeurons]';  % preferred value encoded by the units
  sigmaNoise = 0.0;                        % gausian Noise for mediator activation
  
  %% --- situational Mediators
  
  % linear mapping from the motivation layer units to the mediation layer units
  step = (motNeurons-1)/(n-1);
  actFun = [1:step:motNeurons];                       % ramp activation function
  actFunForNSteps = ones(medNeurons,n).*actFun;
  
  % firing rate for f mediator neurons (set iqual to a ramp function)
  f_As = f_Cs = generateGaussian(sigmaMed,actFunForNSteps,preferedUnitValue);
  actFun = [motNeurons:-step:1];                      % reversed ramp activation function
  actFunForNSteps = ones(medNeurons,n).*actFun;
  f_Rs =  generateGaussian(sigmaMed,actFunForNSteps,preferedUnitValue);
  
  %% Setting initial states
  
  % Activation Ms
  h_Ms_case1 = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Ms_case2 = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Ms_all_case1 = [];
  h_Ms_all_case2 = [];
  
  % Activation Mc
  h_Mc = ones(motNeurons,1).*(1/motNeurons); % uniform distribution

  % Applying the inhibition factors
  W_Ms = WCANN .+ eS;
  
  % TNL connections and activation
  WTNLtoAs = WTNLtoCs = WTNLtoRs = zeros(medNeurons,nMediators);
  % unbiased contribution from mediators
  WTNLtoAs(:,1) =  ones(medNeurons,1);
  WTNLtoCs(:,2) =  ones(medNeurons,1);
  WTNLtoRs(:,3) =  ones(medNeurons,1);
  
  f_TNL_case1 = [0.2 0.3 0.7]';
  f_TNL_case2 = [0.5 0.3 0.3]';
  
  %% Simulation begin 
  for i = 1 : n
    
    % Situational level - case 1

    h_Ms_all_case1 = [h_Ms_all_case1 h_Ms_case1];
    
    h_Ms = h_Ms_case1;
    f_Ms  = h_Ms;
    f_TNL = f_TNL_case1;
    
    h_As  = WTNLtoAs*f_TNL.*f_As(:,i);
    h_Cs  = WTNLtoCs*f_TNL.*f_Cs(:,i);
    h_Rs  = WTNLtoRs*f_TNL.*f_Rs(:,i);
       
    h_Ms = (1-(dt/tS))*h_Ms + (dt/tS)*(W_Ms*f_Ms + g1*WAuToM*h_As + g2*WCoToM*h_Cs + g3*WReToM*h_Rs + g4*h_Mc);
    % normalizing
    h_Ms_case1 = h_Ms./sum(h_Ms,1);

    % Situational level - case 2

    h_Ms_all_case2 = [h_Ms_all_case2 h_Ms_case2];
    
    h_Ms = h_Ms_case2;
    f_Ms  = h_Ms;
    f_TNL = f_TNL_case2;
    
    h_As  = WTNLtoAs*f_TNL.*f_As(:,i);
    h_Cs  = WTNLtoCs*f_TNL.*f_Cs(:,i);
    h_Rs  = WTNLtoRs*f_TNL.*f_Rs(:,i);
       
    h_Ms = (1-(dt/tS))*h_Ms + (dt/tS)*(W_Ms*f_Ms + g1*WAuToM*h_As + g2*WCoToM*h_Cs + g3*WReToM*h_Rs + g4*h_Mc);
    % normalizing
    h_Ms_case2 = h_Ms./sum(h_Ms,1);

     
  endfor
  %% Simulation end

  % begin Ploting
  display("Plotting CS1-B ... " );
  fflush(stdout);
  
  % Mediatiors activation plot
  
  figure(3);
  times = (dt.*[0:n-1])/60;
  hold on;
  
  for i = 1 : medNeurons
    plot (times, f_As(i,:), 'b');  
    plot (times, f_Rs(i,:), 'r');
  endfor 
  
  hx = xlabel("Time in min");
  hy = ylabel("Activation");
  ht = title("Time evolution of mediator's firing rate");
  
  print -deps -color plots/CSb_ActFunMed.eps

  % saving plot data
  plotMat = [times; f_As; f_Cs; f_Rs];
  csvwrite("csv/CSb_ActFunMed.csv", plotMat);
  
  % Situational motivation case 1
  
  figure(4);
  h_Ms_all = h_Ms_all_case1;
  times = (dt.*[0:n-1])/60;
  hold on;

  plot (times, h_Ms_all(1,:), 'r');
  plot (times, h_Ms_all(2,:), 'g');
  plot (times, h_Ms_all(3,:), 'k');
  plot (times, h_Ms_all(4,:), 'c');
  plot (times, h_Ms_all(5,:), 'm');
  plot (times, h_Ms_all(6,:), 'b');
    
  hx = xlabel("Time in min");
  hy = ylabel("Activation");
  ht = title("Time evolution of situational motivation");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  legend (h, "location", "northeastoutside");
  
  print -deps -color plots/CS1b_ActFunMs1.eps
  
  % saving plot data
  plotMat = [times; h_Ms_all];
  csvwrite("csv/CS1b_ActFunMs1.csv", plotMat);

 % Situational motivation case 2
  
  figure(5);
  h_Ms_all = h_Ms_all_case2;
  times = (dt.*[0:n-1])/60;
  hold on;

  plot (times, h_Ms_all(1,:), 'r');
  plot (times, h_Ms_all(2,:), 'g');
  plot (times, h_Ms_all(3,:), 'k');
  plot (times, h_Ms_all(4,:), 'c');
  plot (times, h_Ms_all(5,:), 'm');
  plot (times, h_Ms_all(6,:), 'b');
    
  hx = xlabel("Time in min");
  hy = ylabel("Activation");
  ht = title("Time evolution of situational motivation");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  legend (h, "location", "northeastoutside");
  
  print -deps -color plots/CS1b_ActFunMs1.eps
  
  % saving plot data
  plotMat = [times; h_Ms_all];
  csvwrite("csv/CS1b_ActFunMs1.csv", plotMat);
  
endif 

%-------------------------------------- Case study 1-C -------------------------------------- 

if CS1(3) > 0
  display("Simulating CS1-C ... " );
  fflush(stdout);

  dt = 5.0; # in sec
  
  % Parameters for situational level 
  tS = 100*dt;
  eS_case1 = 0.0;
  eS_case2 = 0.1;
  
  g1 = g2 = g3 = 1;
  g4 = 0;

  % simulation time
  simTimeMin = 120; 
  timeSec = simTimeMin*60;
  n = timeSec*(1/dt); 

  %% setting Mediator's activation 
   
  step =  (motNeurons-1)/(medNeurons-1);
  preferedUnitValue = [1:step:motNeurons]';  % preferred value encoded by the units
  sigmaNoise = 0.0;                        % gausian Noise for mediator activation
  
  %% --- situational Mediators
  
  % linear mapping from the motivation layer units to the mediation layer units
  step = (motNeurons-1)/(n-1);
  actFun = [1:step:motNeurons];                       % ramp activation function
  actFunForNSteps = ones(medNeurons,n).*actFun;
  
  % firing rate for f mediator neurons (set iqual to a ramp function)
  f_As = f_Cs = f_Rs = generateGaussian(sigmaMed,actFunForNSteps,preferedUnitValue);
  
  %% Setting initial states
  
  % Activation Ms
  h_Ms_case1 = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Ms_case2 = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Ms_all_case1 = [];
  h_Ms_all_case2 = [];
  
  % Activation Mc
  h_Mc = ones(motNeurons,1).*(1/motNeurons); % uniform distribution

  % Applying the inhibition factors
  W_Ms_case1 = WCANN .+ eS_case1;
  W_Ms_case2 = WCANN .+ eS_case2;
  
  % TNL connections and activation
  WTNLtoAs = WTNLtoCs = WTNLtoRs = zeros(medNeurons,nMediators);
  % unbiased contribution from mediators
  WTNLtoAs(:,1) =  ones(medNeurons,1);
  WTNLtoCs(:,2) =  ones(medNeurons,1);
  WTNLtoRs(:,3) =  ones(medNeurons,1);
  
  f_TNL = ones(nMediators,1);
  
  middle = round(n/2);
  
  %% Simulation begin 
  for i = 1 : n
    
    if i>middle
       f_TNL = zeros(3,1);
    endif

    % Situational level

    h_As  = WTNLtoAs*f_TNL.*f_As(:,i);
    h_Cs  = WTNLtoCs*f_TNL.*f_Cs(:,i);
    h_Rs  = WTNLtoRs*f_TNL.*f_Rs(:,i);
    
    % case 1 
    
    h_Ms_all_case1 = [h_Ms_all_case1 h_Ms_case1];
    
    h_Ms = h_Ms_case1;
    f_Ms  = h_Ms;
    W_Ms = W_Ms_case1;
  
    h_Ms = (1-(dt/tS))*h_Ms + (dt/tS)*(W_Ms*f_Ms + g1*WAuToM*h_As + g2*WCoToM*h_Cs + g3*WReToM*h_Rs + g4*h_Mc);
    
    h_Ms_case1 = h_Ms./sum(h_Ms,1); % Normalizing 
    
    % case 2
    
    h_Ms_all_case2 = [h_Ms_all_case2 h_Ms_case2];
    
    h_Ms = h_Ms_case2;
    f_Ms  = h_Ms;
    W_Ms = W_Ms_case2;
    
    h_Ms = (1-(dt/tS))*h_Ms + (dt/tS)*(W_Ms*f_Ms + g1*WAuToM*h_As + g2*WCoToM*h_Cs + g3*WReToM*h_Rs + g4*h_Mc);
    h_Ms_case2 = h_Ms./sum(h_Ms,1); % Normalizing
    
  endfor
  
  %% Simulation end

  display("Plotting CS1-C ... " );
  fflush(stdout);

  % begin Ploting
  
  %% Situational motivation - case 1
  figure(6);
  h_Ms_all = h_Ms_all_case1;
  times = (dt.*[0:n-1])/60;
  hold on;

  plot (times, h_Ms_all(1,:), 'r');
  plot (times, h_Ms_all(2,:), 'g');
  plot (times, h_Ms_all(3,:), 'k');
  plot (times, h_Ms_all(4,:), 'c');
  plot (times, h_Ms_all(5,:), 'm');
  plot (times, h_Ms_all(6,:), 'b');

  hx = xlabel("Time in min");
  hy = ylabel("Activation");
  ht = title("Time evolution of situational motivation");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  legend (h, "location", "northeastoutside");
  
  print -deps -color plots/CS1c_ActFunMs1.eps

  % saving plot data
  plotMat = [times; h_Ms_all];
  csvwrite("csv/CS1c_ActFunMs2.csv", plotMat);

  %% Situational motivation - case 2
  
  figure(7);
  h_Ms_all = h_Ms_all_case2;
  times = (dt.*[0:n-1])/60;
  hold on;

  plot (times, h_Ms_all(1,:), 'r');
  plot (times, h_Ms_all(2,:), 'g');
  plot (times, h_Ms_all(3,:), 'k');
  plot (times, h_Ms_all(4,:), 'c');
  plot (times, h_Ms_all(5,:), 'm');
  plot (times, h_Ms_all(6,:), 'b');

  hx = xlabel("Time in min");
  hy = ylabel("Activation");
  ht = title("Time evolution of situational motivation");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  legend (h, "location", "northeastoutside");
  
  print -deps -color plots/CS1c_ActFunMs2.eps

  % saving plot data
  plotMat = [times; h_Ms_all];
  csvwrite("csv/CS1c_ActFunMs2.csv", plotMat);
  
endif 

display ("End of Case study 1 (CS1)");
display ("DCMM program finished")
