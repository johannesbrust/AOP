%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AOP: Algorithms for Oblique Projection Matrices
% J.J. Brust, R.F. Marcia, C.G. Petra
%
% test_gm_krylov.m: Script to test the build-in 'krylov' test matrix
%
% Description: Form oblique projections 
%
% W = X inv(Y'X) Y',
%
% where X,Y are constructed from the build in function: "gallery('krylov',...)"
%
% 10/26/18, J.B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

%% 'krylov': Krylov subspace matrices B = [x, Ax, A^2x,...]

mn      =[20,4]; % 10000,100
sIn     = speye(mn(1));

% Extra computation of the matrix AX, AY for the krylov matrix.
AX      = gallery('toeppen',mn(1));
AY      = AX + sIn;

X(:,:)  = gallery('krylov',AX,ones(mn(1),1),mn(2));
Y(:,:)  = gallery('krylov',AY,ones(mn(1),1),mn(2));

YX(:,:) = Y(:,:)'*X(:,:);

W(:,:)  = X(:,:)*(YX(:,:)\Y(:,:)');
Wh(:,:) = full(sIn) - W;

err_fW  = norm(W*W-W,'fro');
err_fWh = norm(Wh*Wh-Wh,'fro');
