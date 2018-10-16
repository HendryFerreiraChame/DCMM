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
#   Description:  this program implenets Case Study 2 (CS2). 
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
display ("Begin of Case study 2 (CS2)...");
fflush(stdout); 

%-------------------------------------- Menu for similation selection -------------------------------------- 
display ("Program menu"); 
 
c = menu("Please select an option", "1) Generate all plots", "2) Generate CS2-A", "3) Generate CS2-B", "4) Generate CS2-C","5) Exit the program"); 

%% For bypassing the menu, uncomment the line below
%c = 1

CS2 = [0 0 0];

switch (c)
  case 1  
    CS2(1:3) = 1;
    display ("Option selected: 1) Generate all plots");
  case 2  
    CS2(1) = 1;
    display ("Option selected: 2) Generate CS2-A"); 
  case 3  
    CS2(2) = 1;
    display ("Option selected: 2) Generate CS2-B");
  case 4  
    CS2(3) = 1;
    display ("Option selected: 2) Generate CS2-C");
  case 5
     display ("Option selected: 5) Exiting program");
  otherwise 
    display("Invalid option selected");
endswitch

fflush(stdout); 

if ~exist("plots")&&(c<5)
    
    display("Creating [plots] dir ...");
    fflush(stdout); 
    mkdir "plots";
    
endif

display("Configuring shared parameters ...");
fflush(stdout); 

%% Shared simulation  parameters

motNeurons = 6;    % Number motivation layer units
nContext = 4;      % Life context simulated (1: Education, 2: Interpersonal relations, 3: Leissure, 4: Null)
nMediators = 3;    % Number of mediator processes (1: Autonomy, 2: Competence, 3: Relatedness)
medNeurons = 11;   % Number of mediator layer units
  
sigmaCann = 0.3;  % Sigma CANN gaussian activation
sigmaMed = 0.5;   % Sigma Mediator gaussian activation

%% DCMM parameter setting 

setWeightParameters(sigmaCann, sigmaMed, motNeurons, medNeurons);
 
load("parameters/WAuToM.m");
load("parameters/WCoToM.m");
load("parameters/WReToM.m");
load("parameters/WCANN.m");


%% Simulation Cases

%-------------------------------------- Case study 2-A -------------------------------------- 
  
if CS2(1) > 0
  display("Simulating CS2-A ... " );
  fflush(stdout); 
  
  dt = 5.0; # in sec
  
  % Parameters for situational level 
  tS = 100*dt;
  eS = 0.1;
  g1 = g2 = g3 = g4 = 1;
  
  % Parameters for contextual level 
  tC = tS*20;  
  eC = eye(nContext).*0.1;
  g5 = g6 = g7 = g8 = g9 = ones(nContext,1);
  
  % Parameters for global level 
  tG = tC*20;
  eG = 0.1;
  g10 = g11 = g12 = g13 = 1;
  
  % simulation time
  simTimeMin = 120;         % time in min 
  timeSec = simTimeMin*60;  % time in sec 
  n = timeSec*(1/dt);       % step number

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

  %% --- Contextual Mediators
  
  f_Ac_k = zeros(medNeurons,n,nContext);
  f_Cc_k = zeros(medNeurons,n,nContext);
  f_Rc_k = zeros(medNeurons,n,nContext); 
  
  motType = 2; %EM1
  actFunForNStepsA = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsC = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsR = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_Ac_k(:,:,1) = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue);
  f_Cc_k(:,:,1) = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);
  f_Rc_k(:,:,1) = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);
  
  %% --- Global Mediators
  
  motType = 5; %EM4
  actFunForNStepsA = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsC = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsR = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_Ag = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue);
  f_Cg = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);
  f_Rg = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);
  
  %% Setting initial states
  
  % Activation Ms
  h_Ms = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Ms_all = [];
  
  % Activation Mc
  motType = 2; %EM1
  h_Mc = [1 7 2 1 0.5 0.2]'; 
  h_Mc = h_Mc ./sum(h_Mc);
  h_Mc_all = [];
  
  % Activation Mc for the k contexts
  h_Mc_K = ones(motNeurons,nContext).*(1/motNeurons); % uniform distribution
  h_Mc_K(:,1) = h_Mc; % only education is considered
  h_Mc_K_all =  [];
  
  % Activation Mg  
  motType = 5; %EM4
  h_Mg = [0.1 0.3 0.5 1.7 6 2.3]'; 
  h_Mg = h_Mg ./sum(h_Mg);
  h_Mg_all = [];
 
  % Applying the inhibition factors
  W_Ms = WCANN .+ eS; 
  W_Mg = WCANN .+ eG; 
 
  % inter-contexts recurrent connexions 
  W_Mc_K = zeros(motNeurons*nContext,motNeurons,nContext);
  index = [0: nContext].*motNeurons;
  
  for i = 2 : nContext+1
    W_Mc_K(index(i-1)+1:index(i),1:motNeurons,i-1) = WCANN;  
  endfor 
  
  % CAL connections and activation 
  WCALtoMc_K = zeros(motNeurons,nContext,nContext);
  % no inter context dynamics is included
  WCALtoMc_K(:,1,1) = ones(motNeurons,1);
  WCALtoMc_K(:,2,2) = ones(motNeurons,1);
  WCALtoMc_K(:,3,3) = ones(motNeurons,1);
  WCALtoMc_K(:,4,4) = ones(motNeurons,1);
  
  % Only education is relevant for this case simulation
  h_CAL = [1 0 0 0]';
  h_CAL = h_CAL./sum(h_CAL);
  
  % TNL connections and activation
  WTNLtoAs = WTNLtoCs = WTNLtoRs = zeros(medNeurons,nMediators);
  % no inter need interaction is set
  WTNLtoAs(:,1) =  ones(medNeurons,1);
  WTNLtoCs(:,2) =  ones(medNeurons,1);
  WTNLtoRs(:,3) =  ones(medNeurons,1);
  
  f_TNL = ones(nMediators,1); 
  
  
  %% Simulation begin 
     
  for i = 1 : n
  
    % Global level
    
    f_Mg  = h_Mg;
       
    h_Mg = (1-(dt/tG))*h_Mg + (dt/tG)*(W_Mg*f_Mg + g10*WAuToM*f_Ag(:,i) + g11*WCoToM*f_Cg(:,i) + g12*WReToM*f_Rg(:,i) + g13*h_Mc);
    
    h_Mg = h_Mg./sum(h_Mg,1);
    h_Mg_all = [h_Mg_all h_Mg];
    
    % Contextual level    

    h_Mc = zeros(motNeurons,1); 

    for k = 1 : nContext
            
      W_Mc_k = reshape (W_Mc_K(:,:,k), motNeurons*nContext,motNeurons);
      f_Mc_k = zeros(motNeurons,1);
      
      for y = 1 : nContext
      
        iL = (y-1)*motNeurons+1;
        rL = iL + motNeurons -1;
        W_Mc_y_k = W_Mc_k(iL:rL,:) .+ eC(k,y);
        f_Mc_y =  h_Mc_K(:,y);
        f_Mc_k  = f_Mc_k + W_Mc_y_k*f_Mc_y;
        
      endfor
  
      f_Ac = reshape (f_Ac_k(:,:,k), medNeurons,n);
      f_Cc = reshape (f_Cc_k(:,:,k), medNeurons,n);
      f_Rc = reshape (f_Rc_k(:,:,k), medNeurons,n);
      
      h_Ac  = (WAuToM*f_Ac(:,i));
      h_Cc  = (WCoToM*f_Cc(:,i));
      h_Rc  = (WReToM*f_Rc(:,i));
      
      h_Mc_k = (1-(dt/tC))*h_Mc_K(:,k) + (dt/tC)*(f_Mc_k + g5(k)*h_Ac + g6(k)*h_Cc + g7(k)*h_Rc + g8(k)*h_CAL(k)*h_Ms + g9(k)*h_CAL(k)*h_Mg);
      
      h_Mc_k = h_Mc_k./(sum(h_Mc_k,1)+1.0e-10);   
      h_Mc_K(:,k) = h_Mc_k;
      
      WCALtoMc_k = reshape (WCALtoMc_K(:,:,k), motNeurons,nContext);     
      h_Mc = h_Mc + WCALtoMc_k*h_CAL.*h_Mc_k;  
      
    endfor
    
    h_Mc_K_all = [h_Mc_K_all h_Mc_K];
    h_Mc  = h_Mc./(sum(h_Mc));
    h_Mc_all = [h_Mc_all h_Mc];
    
    % Situational level

    f_Ms  = h_Ms;
    
    h_As  = WTNLtoAs*f_TNL.*f_As(:,i);
    h_Cs  = WTNLtoCs*f_TNL.*f_Cs(:,i);
    h_Rs  = WTNLtoRs*f_TNL.*f_Rs(:,i);
       
    h_Ms = (1-(dt/tS))*h_Ms + (dt/tS)*(W_Ms*f_Ms + g1*WAuToM*h_As + g2*WCoToM*h_Cs + g3*WReToM*h_Rs + g4*h_Mc);
    
    h_Ms = h_Ms./sum(h_Ms,1);
    h_Ms_all = [h_Ms_all h_Ms];
    
  
  endfor
  % Simulation end
  
  % begin Ploting
  display("Plotting CS2-A ... " );
  fflush(stdout); 
  
  % --- Situational level plot 
  
  figure(1);
  h_Ms_all = h_Ms_all;
  times = (dt.*[0:n-1])/60;
  hold on;

  plot (times, h_Ms_all(1,:), 'r');
  plot (times, h_Ms_all(2,:), 'g');
  plot (times, h_Ms_all(3,:), 'k');
  plot (times, h_Ms_all(4,:), 'c');
  plot (times, h_Ms_all(5,:), 'm');
  plot (times, h_Ms_all(6,:), 'b');
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  hx = xlabel("Time in min");
  hy = ylabel("Normalized activation");
  ht = title("Time evolution of situational motivation activation functions");

  print -deps -color plots/CS2a_ActFunMs.eps

  % saving plot data
  plotMat = [times; h_Ms_all];
  csvwrite("csv/CS2a_ActFunMs.csv", plotMat);
  
  % --- Contextual level plot 
  
  figure(2);
  h_Mc_all = h_Mc_all;
  times = (dt.*[0:n-1])/60;
  hold on;
  plot (times, h_Mc_all(1,:), 'r');
  plot (times, h_Mc_all(2,:), 'g');
  plot (times, h_Mc_all(3,:), 'k');
  plot (times, h_Mc_all(4,:), 'c');
  plot (times, h_Mc_all(5,:), 'm');
  plot (times, h_Mc_all(6,:), 'b');
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  hx = xlabel("Time in min");
  hy = ylabel("Normalized activation");
  ht = title("Time evolution of contextual motivation activation functions");
   
  print -deps -color plots/CS2a_ActFunMc.eps
  
  % saving plot data
  plotMat = [times; h_Mc_all];
  csvwrite("csv/CS2a_ActFunMc.csv", plotMat);

  % --- Global level plot 
  
  figure(3);
  h_Mg_all = h_Mg_all;
  times = (dt.*[0:n-1])/60;
  hold on;
  plot (times, h_Mg_all(1,:), 'r');
  plot (times, h_Mg_all(2,:), 'g');
  plot (times, h_Mg_all(3,:), 'k');
  plot (times, h_Mg_all(4,:), 'c');
  plot (times, h_Mg_all(5,:), 'm');
  plot (times, h_Mg_all(6,:), 'b');
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  hx = xlabel("Time in min");
  hy = ylabel("Normalized activation");
  ht = title("Time evolution of global motivation activation functions");
   
  print -deps -color plots/CS2a_ActFunMg.eps

  % saving plot data
  plotMat = [times; h_Mg_all];
  csvwrite("csv/CS2a_ActFunMg.csv", plotMat);

endif 

%-------------------------------------- Case study 2-B -------------------------------------- 

if CS2(2) > 0
  display("Simulating CS2-B ... " );
  fflush(stdout);  
  
  dt = 5.0; # in sec
  
  % Parameters for situational level 
  tS = 100*dt;
  eS = 0.1;
  g1 = g2 = g3 = g4 = 1;
  
  % Parameters for contextual level 
  tC = tS*20;  
  eC = eye(nContext).*0.1;
  g5 = g6 = g7 = g8 = g9 = ones(nContext,1);
  g5(4) = g6(4) = g7(4) = 0.1;
  g9(4) = 20;
  
  % Parameters for global level 
  tG = tC*20;
  eG = 0.1;
  g10 = g11 = g12 = g13 = 1;
  
  % simulation time
  simTimeMin = 20;          % time in min 
  timeSec = simTimeMin*60;  % time in sec 
  n = timeSec*(1/dt);       % step number

  %% setting Mediator's activation 
   
  step =  (motNeurons-1)/(medNeurons-1);
  preferedUnitValue = [1:step:motNeurons]';  % preferred value encoded by the units
  sigmaNoise = 0.0;                        % gausian Noise for mediator activation
  
  %% --- situational Mediators
  
  % firing rate for f mediator neurons (set iqual to a ramp function)
  motType = 5.8; %IM
  actFunForNStepsA = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsR = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_As = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue); 
  f_Rs = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);
  motType = 1.6; %close to AM
  actFunForNStepsC = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_Cs = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);

  %% --- Contextual Mediators
  
  f_Ac_k = zeros(medNeurons,n,nContext);
  f_Cc_k = zeros(medNeurons,n,nContext);
  f_Rc_k = zeros(medNeurons,n,nContext); 
  
  motType = 2; %EM1
  actFunForNStepsA = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsC = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsR = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_Ac_k(:,:,4) = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue);
  f_Cc_k(:,:,4) = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);
  f_Rc_k(:,:,4) = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);
  
  %% --- Global Mediators
  
  motType = 1; %AM
  actFunForNStepsA = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsC = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsR = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_Ag = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue);
  f_Cg = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);
  f_Rg = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);
  
  %% Setting initial states
  
  % Activation Ms
  h_Ms = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Ms_all = [];
  
  % Activation Mc
  h_Mc = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Mc_all = [];
  
  % Activation Mc for the k contexts
  h_Mc_K = ones(motNeurons,nContext).*(1/motNeurons); % uniform distribution
  h_Mc_K_all =  [];
  
  % Activation Mg  
  motType = 1; %AM
  h_Mg = generateGaussian(0.8,[motType],[1:motNeurons])';  
  h_Mg_all = [];
 
  % Applying the inhibition factors
  W_Ms = WCANN .+ eS; 
  W_Mg = WCANN .+ eG; 
 
  % inter-contexts recurrent connexions 
  W_Mc_K = zeros(motNeurons*nContext,motNeurons,nContext);
  index = [0: nContext].*motNeurons;
  
  for i = 2 : nContext+1
    W_Mc_K(index(i-1)+1:index(i),1:motNeurons,i-1) = WCANN;  
  endfor 
  
  % CAL connections and activation 
  WCALtoMc_K = zeros(motNeurons,nContext,nContext);
  % no inter context dynamics is included
  WCALtoMc_K(:,1,1) = ones(motNeurons,1);
  WCALtoMc_K(:,2,2) = ones(motNeurons,1);
  WCALtoMc_K(:,3,3) = ones(motNeurons,1);
  WCALtoMc_K(:,4,4) = ones(motNeurons,1);
  
  % Only the null context is relevant for this case simulation
  h_CAL = [0 0 0 1]';
  h_CAL = h_CAL./sum(h_CAL);
  
  % TNL connections and activation
  WTNLtoAs = WTNLtoCs = WTNLtoRs = zeros(medNeurons,nMediators);
  % no inter need interaction is set
  WTNLtoAs(:,1) =  ones(medNeurons,1);
  WTNLtoCs(:,2) =  ones(medNeurons,1);
  WTNLtoRs(:,3) =  ones(medNeurons,1);
  
  f_TNL = ones(nMediators,1); 
  
  %% Simulation begin 
     
  for i = 1 : n
    
    % Global level
    
    f_Mg  = h_Mg;
       
    h_Mg = (1-(dt/tG))*h_Mg + (dt/tG)*(W_Mg*f_Mg + g10*WAuToM*f_Ag(:,i) + g11*WCoToM*f_Cg(:,i) + g12*WReToM*f_Rg(:,i) + g13*h_Mc);
    
    h_Mg = h_Mg./sum(h_Mg,1);
    h_Mg_all = [h_Mg_all h_Mg];
    
    % Contextual level
    
    h_Mc = zeros(motNeurons,1); 
    
    for k = 1 : nContext
            
      W_Mc_k = reshape (W_Mc_K(:,:,k), motNeurons*nContext,motNeurons);
      f_Mc_k = zeros(motNeurons,1);
      
      for y = 1 : nContext
      
        iL = (y-1)*motNeurons+1;
        rL = iL + motNeurons -1;
        W_Mc_y_k = W_Mc_k(iL:rL,:) .+ eC(k,y);
        f_Mc_y =  h_Mc_K(:,y);
        f_Mc_k  = f_Mc_k + W_Mc_y_k*f_Mc_y;
        
      endfor
  
      f_Ac = reshape (f_Ac_k(:,:,k), medNeurons,n);
      f_Cc = reshape (f_Cc_k(:,:,k), medNeurons,n);
      f_Rc = reshape (f_Rc_k(:,:,k), medNeurons,n);
      
      h_Ac  = (WAuToM*f_Ac(:,i));
      h_Cc  = (WCoToM*f_Cc(:,i));
      h_Rc  = (WReToM*f_Rc(:,i));
      
      h_Mc_k = (1-(dt/tC))*h_Mc_K(:,k) + (dt/tC)*(f_Mc_k + g5(k)*h_Ac + g6(k)*h_Cc + g7(k)*h_Rc + g8(k)*h_CAL(k)*h_Ms + g9(k)*h_CAL(k)*h_Mg);
      
      h_Mc_k = h_Mc_k./(sum(h_Mc_k,1)+1.0e-10);   
      h_Mc_K(:,k) = h_Mc_k;
      
      WCALtoMc_k = reshape (WCALtoMc_K(:,:,k), motNeurons,nContext);     
      h_Mc = h_Mc + WCALtoMc_k*h_CAL.*h_Mc_k;      
      
    endfor
    
    h_Mc_K_all = [h_Mc_K_all h_Mc_K];
    h_Mc  = h_Mc./(sum(h_Mc));
    h_Mc_all = [h_Mc_all h_Mc];
    
    % Situational level
    
    f_Ms  = h_Ms;
    
    h_As  = WTNLtoAs*f_TNL.*f_As(:,i);
    h_Cs  = WTNLtoCs*f_TNL.*f_Cs(:,i);
    h_Rs  = WTNLtoRs*f_TNL.*f_Rs(:,i);
        
    h_Ms = (1-(dt/tS))*h_Ms + (dt/tS)*(W_Ms*f_Ms + g1*WAuToM*h_As + g2*WCoToM*h_Cs + g3*WReToM*h_Rs + g4*h_Mc);
    
    h_Ms = h_Ms./sum(h_Ms,1);
    h_Ms_all = [h_Ms_all h_Ms];
    
  endfor
  % Simulation end  
  
  % begin Ploting
  display("Plotting CS2-B ... " );
  fflush(stdout); 
  
  % --- Situational level plot 
  
  figure(4);
  h_Ms_all = h_Ms_all;
  times = (dt.*[0:n-1])/60;
  hold on;
  
  plot (times, h_Ms_all(1,:), 'r');
  plot (times, h_Ms_all(2,:), 'g');
  plot (times, h_Ms_all(3,:), 'k');
  plot (times, h_Ms_all(4,:), 'c');
  plot (times, h_Ms_all(5,:), 'm');
  plot (times, h_Ms_all(6,:), 'b');
  
  hx = xlabel("Time in min");
  hy = ylabel("Normalized activation");
  ht = title("Time evolution of situational motivation activation functions");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  legend (h, "location", "northeastoutside");

  print -deps -color plots/CS2b_ActFunMs.eps

  % saving plot data
  plotMat = [times; h_Ms_all];
  csvwrite("csv/CS2b_ActFunMs.csv", plotMat);
  
  % --- Contextual level plot 
  
  figure(5);
  h_Mc_all = h_Mc_all;
  times = (dt.*[0:n-1])/60;
  hold on;
  plot (times, h_Mc_all(1,:), 'r');
  plot (times, h_Mc_all(2,:), 'g');
  plot (times, h_Mc_all(3,:), 'k');
  plot (times, h_Mc_all(4,:), 'c');
  plot (times, h_Mc_all(5,:), 'm');
  plot (times, h_Mc_all(6,:), 'b');
  hx = xlabel("Time in min");
  hy = ylabel("Normalized activation");
  ht = title("Time evolution of contextual motivation activation functions");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  legend (h, "location", "northeastoutside");
  
  print -deps -color plots/CS2b_ActFunMc.eps
  
  % saving plot data
  plotMat = [times; h_Mc_all];
  csvwrite("csv/CS2b_ActFunMc.csv", plotMat);

  % --- Global level plot 
  
  figure(6);
  h_Mg_all = h_Mg_all;
  times = (dt.*[0:n-1])/60;
  hold on;
  plot (times, h_Mg_all(1,:), 'r');
  plot (times, h_Mg_all(2,:), 'g');
  plot (times, h_Mg_all(3,:), 'k');
  plot (times, h_Mg_all(4,:), 'c');
  plot (times, h_Mg_all(5,:), 'm');
  plot (times, h_Mg_all(6,:), 'b');
  hx = xlabel("Time in min");
  hy = ylabel("Normalized activation");
  ht = title("Time evolution of global motivation activation functions");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  legend (h, "location", "northeastoutside");
  
  print -deps -color plots/CS2b_ActFunMg.eps
  
  % saving plot data
  plotMat = [times; h_Mg_all];
  csvwrite("csv/CS2b_ActFunMg.csv", plotMat);

endif

%-------------------------------------- Case study 2-C -------------------------------------- 

if CS2(3) > 0
  display("Simulating CS2-C ... " );
  fflush(stdout); 
  
  dt = 5.0; # in sec
  
  % Parameters for situational level 
  tS = 100*dt;
  eS = 0.1;
  g1 = g2 = g3 = g4 = 1;
  
  % Parameters for contextual level 
  tC = tS*20;  
  eC = eye(nContext).*0.1;
  g5 = g6 = g7 = g8 = g9 = ones(nContext,1);
  g5(3) = g6(3) = g7(3) = g8(3) = g9(3) = 0.07;
  
  % Parameters for global level 
  tG = tC*20;
  eG = 0.1;
  g10 = g11 = g12 = g13 = 1;
  
  % simulation time
  timeDays = 30;             % time in days 
  timeSec = timeDays*24*60;  % time in sec 
  n = timeSec*(1/dt);        % step number

  %% setting Mediator's activation 
   
  step =  (motNeurons-1)/(medNeurons-1);
  preferedUnitValue = [1:step:motNeurons]';  % preferred value encoded by the units
  sigmaNoise = 0.0;                        % gausian Noise for mediator activation
  
  %% --- situational Mediators
  
  % firing rate for f mediator neurons (set iqual to a ramp function)
  motType = 2; %EM1
  actFunForNStepsA = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsC = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsR = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_As = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue);
  f_Cs = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);
  f_Rs = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);

  %% --- Contextual Mediators
  
  f_Ac_k = zeros(medNeurons,n,nContext);
  f_Cc_k = zeros(medNeurons,n,nContext);
  f_Rc_k = zeros(medNeurons,n,nContext); 
  
  % Mediation for Education life context
  motType = 2; %EM1
  actFunForNStepsA = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsC = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsR = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_Ac_k(:,:,1) = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue);
  f_Cc_k(:,:,1) = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);
  f_Rc_k(:,:,1) = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);

  % Mediation for Leissure life context
  motType = 3; %EM2
  actFunForNStepsA = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsC = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsR = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_Ac_k(:,:,3) = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue);
  f_Cc_k(:,:,3) = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);
  f_Rc_k(:,:,3) = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);
  
  %% --- Global Mediators
  
  motType = 5; %EM4
  actFunForNStepsA = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsC = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  actFunForNStepsR = ones(1,n).*motType + normrnd(0, sigmaNoise, 1, n); 
  f_Ag = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue);
  f_Cg = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);
  f_Rg = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);
  
  %% Setting initial states
  
  % Activation Ms
  h_Ms = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Ms_all = [];
  
  % Activation Mc
  h_Mc = ones(motNeurons,1).*(1/motNeurons); % uniform distribution
  h_Mc_all = [];
  
  % Activation Mc for the k contexts
  h_Mc_K = ones(motNeurons,nContext).*(1/motNeurons); % uniform distribution
  
  % set EM4 to education context  
  v = [1 2 5 7 10 6]';
  h_Mc_K(:,1) = v./sum(v);

  % set EM2 to leissure context
  v = [1 4 5 3.5 2 1]';
  h_Mc_K(:,3) = v./sum(v);;

  h_Mc_K_all =  [];

  % Activation Mg  
  h_Mg = [0.1 0.3 0.5 1.7 6 2.3]'; % set EM4
  h_Mg = h_Mg ./sum(h_Mg);  
  h_Mg_all = [];
 
  % Applying the inhibition factors
  W_Ms = WCANN .+ eS; 
  W_Mg = WCANN .+ eG; 
 
  % inter-contexts recurrent connexions 
  W_Mc_K = zeros(motNeurons*nContext,motNeurons,nContext);
  index = [0: nContext].*motNeurons;
  
  for i = 2 : nContext+1
    W_Mc_K(index(i-1)+1:index(i),1:motNeurons,i-1) = WCANN;  
  endfor 
  
  % assign the compensatory effect to leusure (3) from educational (1) 
  % by mirrowing WCANN
  interW = zeros(motNeurons);
  v = [10 8 6 1 0.1 0];
  v = v./sum(v);
  interW(6,:) = v;
  W_Mc_K(1:motNeurons,1:motNeurons,3) = interW;

  % CAL connections and activation 
  WCALtoMc_K = zeros(motNeurons,nContext,nContext);
  % no inter context dynamics is included here
  WCALtoMc_K(:,1,1) = ones(motNeurons,1);
  WCALtoMc_K(:,2,2) = ones(motNeurons,1);
  WCALtoMc_K(:,1,3) = ones(motNeurons,1);
  WCALtoMc_K(:,4,4) = ones(motNeurons,1);
  
  % Only education is relevant for this case simulation
  h_CAL = [1 0 0 0]';
  h_CAL = h_CAL./sum(h_CAL);
  
  % TNL connections and activation
  WTNLtoAs = WTNLtoCs = WTNLtoRs = zeros(medNeurons,nMediators);
  % no inter need interaction is set
  WTNLtoAs(:,1) =  ones(medNeurons,1);
  WTNLtoCs(:,2) =  ones(medNeurons,1);
  WTNLtoRs(:,3) =  ones(medNeurons,1);
  
  f_TNL = ones(nMediators,1); 
  
  %% Simulation begin 
     
  for i = 1 : n
    
    % Global level
    
    f_Mg  = h_Mg;
       
    h_Mg = (1-(dt/tG))*h_Mg + (dt/tG)*(W_Mg*f_Mg + g10*WAuToM*f_Ag(:,i) + g11*WCoToM*f_Cg(:,i) + g12*WReToM*f_Rg(:,i) + g13*h_Mc);
    
    h_Mg = h_Mg./sum(h_Mg,1);
    h_Mg_all = [h_Mg_all h_Mg];
    
    % Contextual level
    h_Mc = zeros(motNeurons,1); 
    
    for k = 1 : nContext
            
      W_Mc_k = reshape (W_Mc_K(:,:,k), motNeurons*nContext,motNeurons);
      f_Mc_k = zeros(motNeurons,1);
      
      for y = 1 : nContext
      
        iL = (y-1)*motNeurons+1;
        rL = iL + motNeurons -1;
        W_Mc_y_k = W_Mc_k(iL:rL,:) .+ eC(k,y);
        f_Mc_y =  h_Mc_K(:,y);
        f_Mc_k  = f_Mc_k + W_Mc_y_k*f_Mc_y;
        
      endfor
  
      f_Ac = reshape (f_Ac_k(:,:,k), medNeurons,n);
      f_Cc = reshape (f_Cc_k(:,:,k), medNeurons,n);
      f_Rc = reshape (f_Rc_k(:,:,k), medNeurons,n);
      
      h_Ac  = (WAuToM*f_Ac(:,i));
      h_Cc  = (WCoToM*f_Cc(:,i));
      h_Rc  = (WReToM*f_Rc(:,i));
      
      h_Mc_k = (1-(dt/tC))*h_Mc_K(:,k) + (dt/tC)*(f_Mc_k + g5(k)*h_Ac + g6(k)*h_Cc + g7(k)*h_Rc + g8(k)*h_CAL(k)*h_Ms + g9(k)*h_CAL(k)*h_Mg);
      
      h_Mc_k = h_Mc_k./(sum(h_Mc_k,1)+1.0e-10);   
      h_Mc_K(:,k) = h_Mc_k;
      
      WCALtoMc_k = reshape (WCALtoMc_K(:,:,k), motNeurons,nContext);     
      
      h_Mc = h_Mc + WCALtoMc_k*h_CAL.*h_Mc_k;        
      
    endfor
    
    h_Mc_K_all = [h_Mc_K_all h_Mc_K];
    h_Mc  = h_Mc./(sum(h_Mc));
    h_Mc_all = [h_Mc_all h_Mc];
    
    % Situational level
    
    f_Ms  = h_Ms;
    
    h_As  = WTNLtoAs*f_TNL.*f_As(:,i);
    h_Cs  = WTNLtoCs*f_TNL.*f_Cs(:,i);
    h_Rs  = WTNLtoRs*f_TNL.*f_Rs(:,i);
       
    h_Ms = (1-(dt/tS))*h_Ms + (dt/tS)*(W_Ms*f_Ms + g1*WAuToM*h_As + g2*WCoToM*h_Cs + g3*WReToM*h_Rs + g4*h_Mc);
    
    h_Ms = h_Ms./sum(h_Ms,1);
    h_Ms_all = [h_Ms_all h_Ms];
    
  endfor
  % Simulation end
  
  % begin Ploting
  display("Plotting CS2-C ... " );
  fflush(stdout); 
  
  % --- Education context plot 
  
   times = (dt.*[0:n-1])/(60*24);
    
  figure(7);
  
  index = [1:nContext:n*nContext];
  h_Mc_k = h_Mc_K_all(:,index);  
  hold on;
  plot (times, h_Mc_k(1,:), 'r');
  plot (times, h_Mc_k(2,:), 'g');
  plot (times, h_Mc_k(3,:), 'k');
  plot (times, h_Mc_k(4,:), 'c');
  plot (times, h_Mc_k(5,:), 'm');
  plot (times, h_Mc_k(6,:), 'b');
  hx = xlabel("Time in days");
  hy = ylabel("Normalized activation");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  ht = title("Time evolution of educational motivation activation functions");
  
  print -deps -color plots/CS2c_ActFunMgEdu.eps
  
  % saving plot data
  plotMat = [times; h_Mc_k];
  csvwrite("csv/CS2c_ActFunMgEdu.csv", plotMat);  
  
  % --- Leussure context plot 
  
  figure(8);
  
  index = [3:nContext:n*nContext];
  h_Mc_k = h_Mc_K_all(:,index);
  hold on;
  plot (times, h_Mc_k(1,:), 'r');
  plot (times, h_Mc_k(2,:), 'g');
  plot (times, h_Mc_k(3,:), 'k');
  plot (times, h_Mc_k(4,:), 'c');
  plot (times, h_Mc_k(5,:), 'm');
  plot (times, h_Mc_k(6,:), 'b');
  hx = xlabel("Time in days");
  hy = ylabel("Normalized activation");
  ht = title("Time evolution of leisure motivation activation functions");
  h = legend("AM ", " EM_1 ", " EM_2 ", " EM_3 ", " EM_4 ", "IM ");
  legend (h, "location", "northeastoutside");
  
  print -deps -color plots/CS2c_ActFunMgLei.eps
  
  % saving plot data
  plotMat = [times; h_Mc_k];
  csvwrite("csv/CS2c_ActFunMgLei.csv", plotMat);

endif 

display ("End of Case study 2 (CS2)");
display ("DCMM program finished")


