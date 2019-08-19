%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AOP: Algorithms for Oblique Projection Matrices
% J.J. Brust, R.F. Marcia, C.G. Petra
%
% test_gm_chebvand.m: Script to test the build-in 'chebvand' test matrix
%
% Description: Form oblique projections 
%
% W = X inv(Y'X) Y',
%
% where X,Y are constructed from the build in function: "gallery('chebvand',...)"
%
% 10/26/18, J.B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

%% 'chebvand': Vandermonde-like matrix. Type 'help private/chebvand' for
% more information

m       = 20;
Ps      = [(1:4)',(5:8)'];

In      = eye(m);
X(:,:)  = gallery('chebvand',m,Ps(:,1));
Y(:,:)  = gallery('chebvand',m,Ps(:,2));

YX(:,:) = Y(:,:)'*X(:,:);

W(:,:)  = X(:,:)*(YX(:,:)\Y(:,:)');
Wh(:,:) = In - W;

err_fW  = norm(W*W-W,'fro');
err_fWh = norm(Wh*Wh-Wh,'fro');
