%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AOP: Algorithms for Oblique Projection Matrices
% J.J. Brust, R.F. Marcia, C.G. Petra
%
% test_gm_poisson.m: Script to test the build-in 'poisson' test matrix
%
% Description: Form oblique projections 
%
% W = X inv(Y'X) Y',
%
% where X,Y are constructed from the build in function: "gallery('poisson',...)"
%
% 10/26/18, J.B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

%% 'poisson': Block tri-diagonal matrix

mn      =[20,4]; % 10000,100
In      = eye(mn(1));

A      = gallery('poisson',mn(1));

X(:,:)  = A(:,1:mn(2));
Y(:,:)  = A(:,(mn(2)+1):(2*mn(2)));

YX(:,:) = Y(:,:)'*X(:,:);

W(:,:)  = X(:,:)*(YX(:,:)\Y(:,:)');
Wh(:,:) = In - W;

err_fW  = norm(W*W-W,'fro');
err_fWh = norm(Wh*Wh-Wh,'fro');
