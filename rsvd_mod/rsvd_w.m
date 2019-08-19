function [U,S,V] = rsvd_w(X,Y,YX,K)
%------------------------AOP: rsvd_w -------------------------------------%
% Modification of random SVD function rsvd (Antoine Liutkus, Inria 2014),
% for oblique projection matrices W. This is used in the article
% "On the Eigendecomposition and Singular Value Decomposition of Oblique Projection Matrices", 
% J.J.Brust, R.F.Marcia, C.G.Petra, 2018.
%
% The oblique projector, W, and its complement, \hat{W}, are defined by
% 
%   W = X(Y'X)^{-1}Y',  \hat{W} = I - W.
% 
%  input:
%  * X  : Rectangular matrix,
%  * Y  : Rectangular matrix,
%  * YX : Square matrix,
%  * K  : Number of components to keep.
%
%  output:
%  * U,S,V : classical output as the builtin svd matlab function
%-------------------------------------------------------------------------%
% J.B., 2018

 
N               = size(X,1);
P               = min(2*K,N);
Omg             = randn(N,P);

Y_              = X*linsolve(YX,Y'*Omg);

W1              = orth(Y_);


opts_l.TRANSA   = true;
B               = (Y*(linsolve(YX,X'*W1,opts_l)))';

%B               = (W1'*X)*linsolve(YX,Y');

[W2,S,V]        = svd(B,'econ');
U               = W1*W2;
K               = min(K,size(U,2));
U               = U(:,1:K);
S               = S(1:K,1:K);
V               = V(:,1:K);
