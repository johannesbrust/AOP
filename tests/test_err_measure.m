%------------------------- AOP: test err measure -------------------------%
%
% AOP: Algorithms for Oblique Projection Matrices
% J.J. Brust, R.F. Marcia, C.G. Petra
%
% Script to test computing the Frobenius norm of (very) large oblique 
% projection matrices:
%
% W = X inv(Y'X) Y'.
%
% This test is intended to measure the time only.
%
% 10/29/18, J.B.
%
%-------------------------------------------------------------------------%
clc;
clear;

%% Generate data

ns                      = [1000000];
ms                      = [20];

l_ns                    = size(ns,1);
l_ms                    = size(ms,1);

times                   = zeros(l_ns,l_ms);
errs_fro                = zeros(l_ns,l_ms);

%% Nested data
for i = 1:l_ns
    for j = 1:l_ms
       
        n   = ns(i);
        m   = ms(j);
        
        X   = randn(n,m);
        Y   = randn(n,m);

        YX  = Y'*X;
        
        denom   = 1000;
        ratio   = n/denom;
        
        part_rat = 10;
        
        t_err   = tic; 
        
        err_fro = 0;
        for k = 1:part_rat
            
            sidx = (k-1)*denom + 1;
            eidx = k*denom;
            
            M1  = X*linsolve(YX,Y(sidx:eidx,:)');
            M2  = X*linsolve(YX,Y(sidx:eidx,:)');
            
            err_fro = err_fro + norm(M1-M2,'fro'); 
            
        end
        t_err   = toc(t_err);
        
        times(i,j)  = t_err;
        errs_fro    = err_fro;
        
    end 
end

