%-------------------------- AOP EXPERIMENT II  ----------------------------%
% 08/16/18, J.B.,
% Code accompanying the manuscript "On the Eigendecomposition and Singular 
% Value Decomposition of Oblique Projection Matrices", J.J.Brust,
% R.F.Marcia, C.G.Petra, 2018.

% This experiment tests the SVDs of oblique projection matrices
%
%   W - U SI V' = 0,    \hat{W} - \hat{U}\hat{SI}\hat{V}' = 0,
%
% where U, \hat{U}, V, \hat{V} are matrices with orthonormal columns, and 
% SI and \hat{SI} are diagonal matrices. The oblique projector and its
% complement are
%   W = X(Y'X)^{-1}Y',  \hat{W} = I - W.
%
% This experiment tests Algorithm 2, and Algorithm 3 in relation to the
% build-in Matlab command 'svd'.

% Versions:
% 08/16/18, J.B., Setup for large matrices.
% 08/21/18, J.B., Inline testing of U,Uh and V,Vh for orthogonality.
% 08/22/18, J.B., Version to reruns of computations and stores average and
% max. errors.
% 09/09/18, J.B., Preparation for release, and error calculations for small
% or large dimensions.
% 09/14/18, J.B., Tests for larger values of m.
% 09/14/18, J.B., Tests of alternatives for computing the singular vectors of
% \hat{W}.
% 09/17/18, J.B., Modyfing the computation of the singular vectors of both 
% W and \hat{W}. Use multiplications instead of inverses when diagonal
% matrices are used, e.g. use diag(1./d) instead of \diag(d), for a vector
% d = [d1,d2,...,dm]. Preparation of code for release.

%----------------------- Storing of results -----------------------------%
clc;
clear;

fprintf('------- AOP: Algorithms For Oblique Projection Matrices ------- \n');
fprintf('-------      J.J.Brust, R.F.Marica, C.G.Petra, 2018     ------- \n \n');

fprintf('EXPERIMENT II: Singular Vector Computations \n');
fprintf('ALGORITHMS: ALG. 2, ALG.3 and Matlab \n \n');

fprintf('Please enter input A (Times(s)= 1, Error=2). Then hit return. \n');
select_a = input('Input A= ','s');
if isempty(select_a); select_a = '1'; end;

fprintf('Please enter input B (Small=1, Large=2). Then hit return. \n');
select_b = input('Input B= ','s');
if isempty(select_b); select_b = '1'; end;
issmall = 1;

fprintf('\n Starting EXPERIMENT II with: Input A=%s, Input B=%s \n',select_a, select_b);
fprintf(' n \t ALG. 2 \t ALG. 3 \t MTLB.-W \t MTLB.-Wh \n');

if select_b ~= '2'
    ns          = [800;1000;2000;3000;4000;5000];    
else
    ns          = [7500;10000;50000;100000;200000;300000];    
    issmall     = 0;
end

lns         = length(ns);

timeswp     = zeros(lns,1); % Times proposed method for W
timesws     = zeros(lns,1); % Times standard method for W
timeswhp    = zeros(lns,1); % Times proposed method for \hat{W}
timeswhs    = zeros(lns,1); % Times standard method for \hat{W}

errswp      = zeros(lns,1); % Errors proposed method for W, average
errsws      = zeros(lns,1); % Errors standard method for W, average
errswhp     = zeros(lns,1); % Errors proposed method for \hat{W}, average
errswhs     = zeros(lns,1); % Errors standard method for \hat{W}, average

errswpm      = zeros(lns,1); % Errors proposed method for W, max
errswsm      = zeros(lns,1); % Errors standard method for W, max
errswhpm     = zeros(lns,1); % Errors proposed method for \hat{W}, max
errswhsm     = zeros(lns,1); % Errors standard method for \hat{W}, max

ls_opt_trans.TRANSA = true; % Transpose matrix
ls_opt_utri.UT      = true; % Upper triangular matrix
ls_opt_sym.SYM      = true; % Symmetric matrix

% Experiment paramters

nreruns     = 10; % Number of reruns // 10
m           = 4; % Number of columns in X and Y // 4
r           = 8; % Number of columns to compare in the error calculation
del_strt    = 2; % Number of experiments to cut-off from the average 
                 % time/error computations. This is included, because the 
                 % first runs tend to occur start-up times.

if nreruns <= del_strt; del_strt=0; end

om          = ones(m,m);

%---------------------- Computations for varying 'n'---------------------%
for i = 1:lns
    
    n = ns(i);
    
    for ir = 1:nreruns
        
        % Generate data
        X           = randn(n,m);
        Y           = randn(n,m);
        YX          = Y'*X;
        YY          = Y'*Y;
        XX          = X'*X;
        XYm         = [X Y];
                
        % Selection of columns to compare. This is based on problem sizes
        if issmall
            r           = n;
            idx         = 1:r;               
            
            znr         = eye(r);
        else
            idx         = randi(n,r,1);
            
            znr         = zeros(n,r);
            %znr(idx,:)      = eye(r);
            for j = 1:r
                znr(idx(j),j) = 1;
            end
        end
                
        %% Proposed algorithm for W
        % Algorithm 2
        % Orthogonal basis of X = Qp Rp
        
        t_a2    = tic; % Alg.2 timer
        
        Rp      = chol(XX);
        % Temporary buffer Buff
        Buff    = linsolve(YX,Rp',ls_opt_trans);
        
        % Eigendecomposition to compute the singular values, and vectors:
        % Qp' WW' Qp = Vpsi SI Vpsi'.
        [Vpsi, SI2] = eig(Buff'*(YY*Buff));
        
        SI          = sqrt(abs(SI2));
        isi         = 1./diag(SI);
        isis        = spdiags(isi,0,m,m); % Sparse diagonal matrix
        
        % Computing the singular vectors.
           
        U = X*(linsolve(Rp,Vpsi,ls_opt_utri));        
               
        V = Y*(Buff*(Vpsi*isis));        
                        
        % Alg.2 data              
        t_a2            = toc(t_a2);
        if del_strt < ir && del_strt < nreruns
            % Alg.2 time            
            timeswp(i)  = timeswp(i) + t_a2;

            % Alg.2 error
            WIr         = X*linsolve(YX,Y(idx,:)');
            WpcIr       = U*(SI*V(idx,:)');
            ep          = norm(WIr - WpcIr,'fro');        
            errswp(i)   = errswp(i) + ep;
            if errswpm(i) < ep; errswpm(i) = ep; end
        end
        
        %% Proposed algorithm for \hat{W}
        % Algorithm 3
        
        t_a3 = tic; % Alg.3 timer
        
        si      = diag(SI);
        d       = sqrt(abs(1 - 1./(si.^2)));
        d1      = si.*d;
        
        % Singular vectors.
        id1     = 1./d1;
        id1s    = spdiags(id1,0,m,m);
        
        id      = 1./d;
        ids     = spdiags(id,0,m,m);
        
        Uh1     = V*id1s - U*ids;
        Vh1     = V*ids - U*id1s;
        
        %% Projection [X Y]([X Y]'[X Y])^{-1}[X Y]' = Qhpar Qhpar'.
        %QhperP = XYm*linsolve([XX,YX';YX,YY],XYm',ls_opt_sym);        
        
        QhperPs = XYm*linsolve([XX,YX';YX,YY],XYm(idx,:)',ls_opt_sym);
        
        % Alg.3 data
        t_a3            = toc(t_a3);
        if del_strt < ir && del_strt < nreruns
            % Alg.3 time
            timeswhp(i)     = timeswhp(i) + (t_a3+t_a2);

            % Alg.3 error        
            WhIr            = znr - X*linsolve(YX,Y(idx,:)');
            WhpcIr          = Uh1*(SI*Vh1(idx,:)')+(znr-QhperPs(:,:)); % QhperP(:,idx)
            ehp             = norm(WhIr - WhpcIr,'fro');

            errswhp(i)  = errswhp(i) + ehp;
            if errswhpm(i) < ehp; errswhpm(i) = ehp; end
        end
        
        %% Matlab standard computations using the function svd
        if issmall == 1
            
            %% For W
            W               = X*(linsolve(YX,Y'));
            
            tic;            
            [Us,SIs,Vs]     = svd(W);            
            ts              = toc;
            
            % Matlab data
            if del_strt < ir && del_strt < nreruns
                % Matlab time
                timesws(i)      = timesws(i) + ts;
                
                % Matlab error
                WscIr           = Us*(SIs*Vs(idx,:)');
                es              = norm(WIr - WscIr,'fro');

                errsws(i)       = errsws(i) + es;
                if errswsm(i) < es; errswsm(i) = es; end
            end
            
            %% For \hat{W}
            In         = eye(n);
            Wh         = In - X*(linsolve(YX,Y'));
            
            tic;            
            [Uhs,SIhs,Vhs]  = svd(Wh);
            ts              = toc;
            
            % Matlab data
            if del_strt < ir && del_strt < nreruns
                % Matlab time
                timeswhs(i)     = timeswhs(i) + ts;
                
                % Matlab error
                WhscIr          = Uhs*(SIhs*Vhs(idx,:)');
                ehs             = norm(WhIr - WhscIr,'fro');
                
                errswhs(i)      = errswhs(i) + ehs;
                if errswhsm(i) < ehs; errswhsm(i) = ehs; end
            end
            
        end        
    end
    
    denom           = nreruns-del_strt;
    
    timeswp(i)      = timeswp(i)/denom; 
    timesws(i)      = timesws(i)/denom; 
    timeswhp(i)     = timeswhp(i)/denom; 
    timeswhs(i)     = timeswhs(i)/denom; 
    
    errswp(i)       = errswp(i)/denom;
    errsws(i)       = errsws(i)/denom;
    errswhp(i)      = errswhp(i)/denom;
    errswhs(i)      = errswhs(i)/denom;
    
    if issmall == 1
        
        if select_a ~= '2'
        
            fprintf(' %i \t %1.4f \t %1.4f \t %1.4f \t %1.4f \n',n,...
                timeswp(i),timeswhp(i),timesws(i),timeswhs(i));
            
        else
            
            fprintf(' %i \t %1.4f \t %1.4f \t %1.4f \t %1.4f \n',n,...
                errswp(i),errswhp(i),errsws(i),errswhs(i));
            
        end
        
    else
        
       if select_a ~= '2'
        
            fprintf(' %i \t %1.4f \t %1.4f \t %s \t %s \n',n,...
                timeswp(i),timeswhp(i),'------','------');
            
        else
            
            fprintf(' %i \t %1.4f \t %1.4f \t %s \t %s \n',n,...
                errswp(i),errswhp(i),'------','------');
            
        end
        
    end
end

if issmall == 1
    
    times_all   = [timeswp,timesws,timeswhp,timeswhs];
    errs_all    = [errswp,errsws,errswhp,errswhs,...
                   errswpm,errswsm,errswhpm,errswhsm];
    
    save('EXPERIMENT_II_S_SM.mat','times_all','errs_all');
    
else
    
    times_all   = [timeswp,timeswhp];
    errs_all    = [errswp,errswhp,errswpm,errswhpm];
    
    save('EXPERIMENT_II_L_SM.mat','times_all','errs_all');
    
end

