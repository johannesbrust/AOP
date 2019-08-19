%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AOP: Algorithms for Oblique Projection Matrices
% J.J. Brust, R.F. Marcia, C.G. Petra
%
% test_gm_moler.m: Script to test the build-in 'moler' test matrix
%
% Description: Form oblique projections 
%
% W = X inv(Y'X) Y',
%
% where X,Y are constructed from the build in function: "gallery('moler',...)"
%
% 10/26/18, J.B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

%% 'moler': Moler matrix defined as the product of upper-triangular/
% trapezoidal matrices

mn      =[20,4]; % 10000,100
In      = eye(mn(1));

% Extra computation of the matrix AX, AY for the krylov matrix.
A      = gallery('moler',mn(1));

X(:,:)  = A(:,1:mn(2));
Y(:,:)  = A(:,(end-mn(2)+1):end);

YX(:,:) = Y(:,:)'*X(:,:);

W(:,:)  = X(:,:)*(YX(:,:)\Y(:,:)');
Wh(:,:) = In - W;

err_fW  = norm(W*W-W,'fro');
err_fWh = norm(Wh*Wh-Wh,'fro');
