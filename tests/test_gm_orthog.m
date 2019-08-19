%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AOP: Algorithms for Oblique Projection Matrices
% J.J. Brust, R.F. Marcia, C.G. Petra
%
% test_gm_orthog.m: Script to test the build-in 'orthog' matrix based on 
% a modifed routine.
%
% Description: Form oblique projections 
%
% W = X inv(Y'X) Y',
%
% where X,Y are constructed from the build in function: "gallery('orthog',...)"
%
% 10/26/18, J.B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

%% 'orthog': Orthogonal matrices, or scaled orthogonal matrices

tic;
mn      =[1000,20]; % 20,4
%In      = eye(mn(1));

A       = gallery('orthog',mn(1),-1);

X(:,:)  = A(:,1:mn(2));
d       = sum(X.*X,2);
sd      = spdiags(1./d,0,mn(1),mn(1));

Y(:,:)  = sd*A(:,(mn(2)+1):(2*mn(2)));
%Y(:,:)  = sd*X(:,:);

YX(:,:) = Y(:,:)'*X(:,:);

W(:,:)  = X(:,:)*(YX(:,:)\Y(:,:)');
Wh(:,:) = spdiags(1-diag(W(:,:)),0,-W(:,:));

t_op = toc;

%Wh(:,:) = In - W;

err_fW  = norm(W*W-W,'fro');
err_fWh = norm(Wh*Wh-Wh,'fro');

%err_d   = norm(diag(d)-X'*X,'fro');