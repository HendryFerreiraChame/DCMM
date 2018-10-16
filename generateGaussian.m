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
#   Description:  This function calculates the Gaussian activation of neurons 
#                 considering the preferred value 'mu', conforming to Eq. (2)
#
##########################################################################################################


function y = generateGaussian (_sigma,_mu,_preferredValue)

   y = (1/(_sigma*sqrt(2*pi))) * exp(-((_preferredValue.-_mu).^ 2) / (2*(_sigma^2)));
     
end