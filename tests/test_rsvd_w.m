%------------------------- AOP: test rsvd w -------------------------------%
%
% AOP: Algorithms for Oblique Projection Matrices
% J.J. Brust, R.F. Marcia, C.G. Petra
%
% Script to test the randomized svd for the oblique projection,
% W, where
%
% W = X inv(Y'X) Y'.
%
% The rsvd_w function computes 
% 
% U SI V' \approx W,
%
% by exploiting rectangular matrix products (cf. rsvd/rsvd_w.m)
%
% 10/29/18, J.B.
%
%-------------------------------------------------------------------------%
clc;
clear;

addpath('../rsvd_mod');

%% Generate data

t_data  = tic;
mn                      = [10000,20]; % 20,4

X                       = randn(mn);
Y                       = randn(mn);

YX(:,:)                 = Y(:,:)'*X(:,:);

W(:,:)                  = X(:,:)*(YX(:,:)\Y(:,:)');

t_data = toc(t_data);

%% Compute factorization(s)

K                       = mn(2); %mn(2); % mn(1)-mn(2) 

t_wmod                  = tic;
[Uh,SIh,Vh]             = rsvd_w(X,Y,YX,K);
t_wmod                  = toc(t_wmod);

W_mod                   = Uh*(SIh*Vh');

t_wh                    = tic;
[U,SI,V]                = rsvd(W,K);
t_wh                    = toc(t_wh);

W_                      = U*(SI*V');

% Error calculation          
err                     = norm(W_mod-W_,'fro');

