%------------------------- AOP: test rsvd wh -------------------------------%
%
% AOP: Algorithms for Oblique Projection Matrices
% J.J. Brust, R.F. Marcia, C.G. Petra
%
% Script to test the randomized svd for the oblique projection complement,
% \hat{W}, where
%
% \hat{W} = I - X inv(Y'X) Y' = I - W.
%
% The rsvd_wh function computes 
% 
% U SI V' \approx \hat{W},
%
% by exploiting rectangular matrix products (cf. rsvd/rsvd_wh.m)
%
% 10/29/18, J.B.
%
%-------------------------------------------------------------------------%
clc;
clear;

addpath('../rsvd_mod');

%% Generate data

t_data  = tic;
mn                      = [1000,20]; % 20,4

X                       = randn(mn);
Y                       = randn(mn);

YX(:,:)                 = Y(:,:)'*X(:,:);

W(:,:)                  = X(:,:)*(YX(:,:)\Y(:,:)');

diag_idx                = 1:mn(1);

In                      = eye(mn(1));
Wh                      = In-W;

t_data = toc(t_data);

%% Compute factorization(s)

K                       = mn(1)-mn(2); %mn(2); % mn(1)-mn(2) 

t_whmod                 = tic;
[Uwh,SIwh,Vwh]          = rsvd_wh(X,Y,YX,K);
t_whmod                 = toc(t_whmod);

Wh_mod                  = Uwh*(SIwh*Vwh');

t_wh                    = tic;
[Uh,SIh,Vh]             = rsvd(Wh,K);
t_wh                    = toc(t_wh);

Wh_                     = Uh*(SIh*Vh');

err                     = norm(Wh_mod-Wh_,'fro');

%err_d   = norm(diag(d)-X'*X,'fro');