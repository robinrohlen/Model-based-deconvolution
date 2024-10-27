function [t,F]=RaikovaForceTwitch5p(fs,Ti,Tlead,Tc,Thr,Ttot,Fmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  
%  fs      Sample rate (Hz)
%  Ti      Time instant of discharge (ms)
%  Tlead   Spike-to-force onset time (ms)
%  Tc      Contraction time (ms)
%  Thr     Half-relaxation time (ms)
%  Ttot    Duration of twitch (ms)
%  Fmax    Maximal amplitude of twitch (a.u.)
%
% Output:
%  t       Time vector (ms)
%  F       Force (a.u.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=linspace(0,Ttot,Ttot*fs*10^(-3));

F=Fmax.*Tc.^(-log(2.0)./(log(Tc)-log(Thr)+Thr./Tc-1.0)).*t.^(log(2.0)./(log(Tc)-log(Thr)+Thr./Tc-1.0)).*exp(log(2.0)./(log(Tc)-log(Thr)+Thr./Tc-1.0)).*exp(-(t.*log(2.0))./(Tc.*(log(Tc)-log(Thr)+Thr./Tc-1.0)));

t=linspace(Ti+Tlead,Ti+Tlead+Ttot,Ttot*fs*10^(-3));