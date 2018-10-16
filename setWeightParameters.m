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
#   Description:  This function sets the connections weights of a) the recurrent CANN layers  
#                 and b) the mediators to motivation weights
#
##########################################################################################################


function setWeightParameters (_sigmaMot, _sigmaMed, _motNeurons, _medNeurons)


sigmaCann = _sigmaMot;         % sigma for the CANN recurrent connections 
sigmaMed = _sigmaMed;          % sigma for the mediators connections
motNeurons = _motNeurons;      % Units in the motivation layer 
motPreferred = [1:motNeurons]; % Units preferred varue
medNeurons = _medNeurons;

display ("Setting network parameters ...");
display (strcat("Motivation layers: units[",num2str(motNeurons),"], sigma[", num2str(sigmaCann,"%1.4f"), "]"));
display (strcat("Mediator layers: units[",num2str(medNeurons),"], sigma[", num2str(sigmaMed,"%1.4f"), "]"));
fflush(stdout); 


%% setting the WMedToM

% linear mapping from the motivation layer units to the mediation layer units

step = (motNeurons-1)/(medNeurons-1);
medPreferred = [1:step:motNeurons];                      
   
WAuToM = zeros(motNeurons, medNeurons); 
WCoToM = zeros(motNeurons, medNeurons); 
WReToM = zeros(motNeurons, medNeurons); 

for i=1:medNeurons
  
  WAuToM(:,i) = generateGaussian(sigmaMed,medPreferred(i),motPreferred);
  WCoToM(:,i) = generateGaussian(sigmaMed,medPreferred(i),motPreferred);
  WReToM(:,i) = generateGaussian(sigmaMed,medPreferred(i),motPreferred);

end

%%% CANN recurrent connections

WCANN = zeros(motNeurons);

for i=1:motNeurons 
  WCANN(i,:) = generateGaussian(sigmaCann,motPreferred(i),motPreferred);
end

WCANN = WCANN./sum(WCANN,2); % normalization

% saving data

save "parameters/WCANN.m"  WCANN;
save "parameters/WAuToM.m" WAuToM;
save "parameters/WCoToM.m" WCoToM;
save "parameters/WReToM.m" WReToM;

endfunction