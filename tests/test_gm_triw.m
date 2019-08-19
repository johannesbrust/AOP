%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AOP: Algorithms for Oblique Projection Matrices
% J.J. Brust, R.F. Marcia, C.G. Petra
%
% test_gm_triw.m: Script to test the build-in triw test matrix
%
% Description: Form oblique projections 
%
% W = X inv(Y'X) Y',
%
% where X,Y are constructed from the build in function: "gallery('triw',...)"
%
% 10/26/18, J.B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 'triw': Upper triangular/trapezoidal matrix
mn      = [20,4];
alphas  = [2,3];

In      = eye(mn(1));
X(:,:)  = gallery('triw',mn,alphas(1));
Y(:,:)  = gallery('triw',mn,alphas(2));

YX(:,:) = Y(:,:)'*X(:,:);

W(:,:)  = X(:,:)*(YX(:,:)\Y(:,:)');
Wh(:,:) = In - W;

err_fW  = norm(W*W-W,'fro');
err_fWh = norm(Wh*Wh-Wh,'fro');
