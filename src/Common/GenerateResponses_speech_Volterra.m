function [un,dn,vn] = GenerateResponses_speech_Volterra(Volterra_sys, filename)

% GenerateResponses_speech  Generate Input and Desired Responses with Speech Signal
%
% Arguments:
% b                         Unknown system
% filename                  Filename of speech with variable "un" containing
%                             the speech segment to be used
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd

load(filename,'un');
rand('state',sum(100*clock));
vn = (rand(1,length(un))-0.5)*sqrt(12*1e-6); % white noise with variance -60 dB
a = abs(max(un))/sqrt(2); un = un/a;         % Normalize speech to unit variance
un = un + vn;                                % Add noise to speech

max_len= max(Volterra_sys.M); 

% Fast volterra filtering from {Morhac, M., 1991 - A Fast Algorithm of
%                               Nonlinear Volterra Filtering}

%Volterra_sys.Responses{} are kernels
%Volterra_sys.M is the size of kernels

dn = fastVMcell(un, Volterra_sys.Responses , Volterra_sys.M);
dn = sum(dn,1);

% Adjusting starting index of signals

dn = dn(1:length(un));      
dn = dn(max_len+1:end);
un = un(max_len+1:end);
% 
% % Normalization of u(n) and d(n)
% 
% un = un/std(dn);
% dn = dn/std(dn);

% Normalization of speech signal
  
a = abs(max(dn))/sqrt(2);
un = un/a; dn = dn/a;