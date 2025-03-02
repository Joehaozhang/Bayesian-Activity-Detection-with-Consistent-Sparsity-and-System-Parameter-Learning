function [PFA_PMD] = PMDPFA(gamma,gamma_hat,N_thr)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
thrMin = 1/1e1;
thrMax = thrMin* 10 ;

threshold = linspace(thrMin,thrMax,N_thr) ;

PFA = zeros(1,length(threshold)) ;
PMD = zeros(1,length(threshold)) ;

for i=1:length(threshold)
    [Pmd,Pfa] = ComputePfaPmdMod(gamma,gamma_hat,threshold(i)) ;
    PMD(i) = Pmd;
    PFA(i) = Pfa;
end


PFA_PMD = [PFA' PMD'] ;
end