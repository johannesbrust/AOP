%-------------------------- AOP EXPERIMENT I   -------------------------%
% 08/23/18, J.B.,
% Code accompanying the manuscript "On the Eigendecomposition and Singular 
% Value Decomposition of Oblique Projection Matrices", J.J.Brust,
% R.F.Marcia, C.G.Petra, 2018.

% This experiment tests the computation of singular values of oblique 
% projection matrices
%
%   W - U SI V' = 0,    \hat{W} - \hat{U}\hat{SI}\hat{V}' = 0,
%
% where U, \hat{U}, V, \hat{V} are matrices with orthonormal columns 
% (singular vectors), and SI and \hat{SI} are 
% diagonal matrices (singular values). 
% The oblique projector and its complement are
%   W = X(Y'X)^{-1}Y',  \hat{W} = I - W.

% This experiment also computes the singular values by the approach
% described in "Signal Processing Applications of Oblique Projection
% Operators", R.T. Behrens, L.L. Scharf, 1994.

% Versions:
% 08/23/18, J.B., Setup of the script and initial tests.
% 09/09/18, J.B., Storing results of the modified Behrens and Scharf
% approach. Preparation of script for release.
% 09/14/18, J.B., Tests for larger values of m.
% 09/17/18, J.B., Preparation of code for release.
% 09/18/18, J.B., Improvement of the implementation of Algorithm 5, and 
% explicit usage of a basis for the nullspace, Z.

%----------------------- Storing of results -----------------------------%
clc;
clear;

fprintf('------- AOP: Algorithms For Oblique Projection Matrices ------- \n');
fprintf('-------      J.J.Brust, R.F.Marica, C.G.Petra, 2018     ------- \n \n');

fprintf('EXPERIMENT I: Singular Value(s) Computations \n');
fprintf('ALGORITHMS: ALG. 1, ALG.4, ALG. 5 and Matlab \n \n');
fprintf('Please enter input (Times(s)= 1, Error=2). Then hit return. \n');

select = input('Input: ','s');

if isempty(select)
    select = '1';
end
fprintf('\n Starting EXPERIMENT I with input: %s (Times = 1, Errors = 2) \n',select);
fprintf(' n \t ALG. 1 \t ALG. 4 \t ALG. 5 \t MTLB. \n');

ns          = [800;1000;2000;3000;4000;5000];
lns         = length(ns);

times_a1        = zeros(lns,1); % Average times proposed method for 
                                    % singular values of W
times_mt        = zeros(lns,1); % Average times standard method for 
                                    % singular values of W
times_a4        = zeros(lns,1); % Average times B.S. method for 
                                    % singular values of W
times_a5        = zeros(lns,1); % Average times modified B.S. method for 
                                % singular values of W                                   
% Algorithm 1
% Singular value errors
errs_ave_a1_sv  = zeros(lns,1); % Average errors in singular values 
                                    % norm(SI_PROP-SI_STAND) for W
errs_max_a1_sv  = zeros(lns,1); % Max errors in singular values 
                                % norm(SI_PROP-SI_STAND) for W
% Spectral norm errors
errs_ave_a1_sn  = zeros(lns,1); % Average errors in spectral norms 
                                % |max(SI_PROP)-max(SI_STAND)| for W
errs_max_a1_sn  = zeros(lns,1); % Max errors in spectral norms |
                                % max(SI_PROP)-max(SI_STAND)| for W

% Algorithm 4: Behrens and Scharf                                 
% Singular value errors
errs_ave_a4_sv  = zeros(lns,1); % Average errors in singular values 
                                    % norm(SI_PROP-SI_STAND) for \hat{W}, B+S
errs_max_a4_sv = zeros(lns,1); % Max errors in singular values 
                                % norm(SI_PROP-SI_STAND) for \hat{W}, B+S
% Spectral norm errors
errs_ave_a4_sn  = zeros(lns,1); % Average errors in spectral norms 
                                % |max(SI_PROP)-max(SI_STAND)| for W, B+S
errs_max_a4_sn = zeros(lns,1); % Max errors in spectral norms |
                                % max(SI_PROP)-max(SI_STAND)| for W, B+S
% Algorithm 5: Modified Behrens and Scharf                                
% Singular values
errs_ave_a5_sv  = zeros(lns,1); % Average errors in singular values 
                                    % norm(SI_PROP-SI_STAND) for \hat{W}, mod-B+S
errs_max_a5_sv = zeros(lns,1); % Max errors in singular values 
                                % norm(SI_PROP-SI_STAND) for \hat{W}, mod-B+S
% Spectral norm
errs_ave_a5_sn  = zeros(lns,1); % Average errors in spectral norms 
                                % |max(SI_PROP)-max(SI_STAND)| for W, mod-B+S
errs_max_a5_sn = zeros(lns,1); % Max errors in spectral norms |
                                % max(SI_PROP)-max(SI_STAND)| for W, mod-B+S                                

% Experiment parameters

nreruns     = 10; % Number of reruns // 10
m           = 300; % Number of columns in X and Y // 4
del_strt    = 2; % Number of experiments to cut-off from the average 
                 % time/error computations. This is included, because the 
                 % first runs tend to occur start-up times.

if nreruns <= del_strt; del_strt=0; end

ls_opt_trans.TRANSA = true; % Transpose matrix
ls_opt_utri.UT      = true; % Upper triangular matrix
ls_opt_sym.SYM      = true; % Symmetric matrix

%---------------------- Computations for varying 'n'---------------------%
for i = 1:lns
    
    n       = ns(i);
    In      = eye(n); 
    onm2m   = ones(n-2*m,1);
    
    % Reruns            
    for ir = 1:nreruns
    
        % Generate data
        X           = randn(n,m);
        S           = randn(n,m);        
        Sr          = chol(S'*S);
        Qs          = S/Sr;
        
        %Sq          = orth(S);
        Y           = Qs*(Qs'*X);
        
        % Generate Z: A basis of the nullspace of W
        K           = randn(n,(n-m)); 
        Z           = K - S*(linsolve((S'*S),S',ls_opt_sym)*K);
        
        err_sp      = norm(Y'*Z,'fro');
        rk_Z        = rank(Z);
        
        max_rep     = 20;
        iZ          = 1;
        while (1e-5 < err_sp || rk_Z ~= (n-m)) && (iZ < max_rep)
            
            K           = randn(n,(n-m)); 
            Z           = K - S*(linsolve((S'*S),S',ls_opt_sym)*K);
        
            err_sp      = norm(Y'*Z,'fro');
            rk_Z        = rank(Z);
            
            iZ          = iZ + 1;
            
        end
        
        if iZ == max_rep; error(['AOP: EXPERIMENT_I, could not generate null space',...
                                    'matrix Z. Perhaps try rerunning experiment.']); end;
        YX          = Y'*X;

        t_a1 = tic;
        %% Proposed algorithm for W and \hat{W}
        % Algorithm 1

        YY          = Y'*Y;
        XX          = X'*X;
        
        % Orthogonal basis of X = Qp Rp
        Rp = chol(XX);
        
        % Temporary buffer Buff
        Buff = linsolve(YX,Rp',ls_opt_trans);
        
        % Eigendecomposition to compute the singular values only:
        % Qp' WW' Qp = Vpsi SI2 Vpsi'.
        si2         = eig(Buff'*(YY*Buff));
        
        % Singular values of W and \hat{W}
        si          = sqrt(abs(si2));
        
        % Time Alg.1
        t_a1                = toc(t_a1);        
        
        %% Standard computations
        % Matlab
        
        % For W, only computing the non-zero singular values        
        W               = X*(linsolve(YX,Y'));
        
        t_mt = tic;        
        sis             = svds(W,m);
        %[Us,SIs,Vs]     = svd(W);        
        t_mt            = toc(t_mt);        
        
        %% Behrens and Scharf approach
        % Algorithm 4
        t_a4 = tic;
        
        %YY          = Y'*Y;
        XX          = X'*X;
        
        % Orthogonal basis of X = Qp Rp
        Rp = chol(XX);
        Qp = X/Rp;
        
        % Orthogonal basis of perp(Q_Y), where Y = Q_Y R_Y
        %Wy              = In - S*linsolve((S'*S),S',ls_opt_sym);
        
        % Orthonormal basis of Z = Qz Rz
        
        ZZ              = Z'*Z;
        Rz              = chol(ZZ);
        
        Qz              = Z/Rz;
        %Qpery           = orth(Wy);
        
        %[Ly,Dy,py]      = ldl(Wy,'vector');
        %dy              = diag(Dy);
        %idxnzy          = dy > 1e-10;
        %Ry              = (Ly(py(idxnzy),idxnzy)*sqrt(Dy(idxnzy,idxnzy)))';
        %Qpery           = Wy(:,py(idxnzy))/Ry;
        
        % Computing the principal angles and singular values
        %pab             = svds(Qp'*Qpery,m);
        pab             = svds(Qp'*Qz,m);
        
        % Singular values
        sib             = 1./sin(acos(pab));
        
        t_a4 = toc(t_a4);
        
        %% Modified Behrens and Scharf approach
        % Algorithm 5
        t_a5 = tic;

        %XYr             = (Qp'*Y)/Rp;
        
        XX  = X'*X;
        
        % Orthogonal basis of X = Qp Rp
        Rp  = chol(XX);
        
        % Temporary buffer Buff
        Buff = linsolve(YX,Rp',ls_opt_trans);
        
        %sibm2           = eig(XYr);
        
        sibm2           = eig(Rp*Buff);
        %sibm            = 1./sqrt(sibm2);
        sibm            = sqrt(sibm2);
        
        t_a5 = toc(t_a5);
                
        si              = sort(si,1,'descend');
        sis             = sort(sis,1,'descend');
        sib             = sort(sib,1,'descend');        
        sibm            = sort(sibm,1,'descend');

        %% Errors

        e_a1_sv            = norm(si-sis)/n;
        e_a1_sn             = abs(si(1)-sis(1))/n;                
        
        e_a4_sv           = norm(sib-sis)/n;
        e_a4_sn            = abs(sib(1)-sis(1))/n;
        
        e_a5_sv          = norm(sibm-sis)/n;
        e_a5_sn           = abs(sibm(1)-sis(1))/n;
        
        % Data (Time+Errors) of all algorithms        
        if del_strt < ir  
            
            % Times all algorithms
            times_a1(i)      = times_a1(i) + t_a1;            
            times_mt(i)      = times_mt(i) + t_mt;            
            times_a4(i)      = times_a4(i)+t_a4;
            times_a5(i)      = times_a5(i)+t_a5;
        
            % Errors in all singular values 
            errs_ave_a1_sv(i)    = errs_ave_a1_sv(i) + e_a1_sv;        
            if errs_max_a1_sv(i) < e_a1_sv; errs_max_a1_sv(i) = e_a1_sv; end;

            errs_ave_a4_sv(i)   = errs_ave_a4_sv(i) + e_a4_sv;        
            if errs_max_a4_sv(i) < e_a4_sv; errs_max_a4_sv(i) = e_a4_sv; end;

            errs_ave_a5_sv(i)  = errs_ave_a5_sv(i) + e_a5_sv;        
            if errs_max_a5_sv(i) < e_a5_sv; errs_max_a5_sv(i) = e_a5_sv; end;

            % Errors in spectral norms
            errs_ave_a1_sn(i)    = errs_ave_a1_sn(i) + e_a1_sn;        
            if errs_max_a1_sn(i) < e_a1_sn; errs_max_a1_sn(i) = e_a1_sn; end;

            errs_ave_a4_sn(i)   = errs_ave_a4_sn(i) + e_a4_sn;        
            if errs_max_a4_sn(i) < e_a4_sn; errs_max_a4_sn(i) = e_a4_sn; end;

            errs_ave_a5_sn(i)  = errs_ave_a5_sn(i) + e_a5_sn;        
            if errs_max_a5_sn(i) < e_a5_sn; errs_max_a5_sn(i) = e_a5_sn; end;
        end
        
    end
    
    denom               = nreruns-del_strt;
    
    errs_ave_a1_sv(i)        = errs_ave_a1_sv(i)/denom;            
    errs_ave_a1_sn(i)        = errs_ave_a1_sn(i)/denom;
    errs_ave_a4_sv(i)       = errs_ave_a4_sv(i)/denom;
    errs_ave_a4_sn(i)       = errs_ave_a4_sn(i)/denom;
    errs_ave_a5_sv(i)      = errs_ave_a5_sv(i)/denom;
    errs_ave_a5_sn(i)      = errs_ave_a5_sn(i)/denom;
    %errswhsvav(i)       = errswhsvav(i)/nreruns;
    
    times_a1(i)          = times_a1(i)/denom;    
    times_mt(i)          = times_mt(i)/denom;   
    times_a4(i)          = times_a4(i)/denom;
    times_a5(i)          = times_a5(i)/denom;
    %timeswhs(i)         = timeswhs(i)/nreruns;
    
    if select ~= '2'
        fprintf(' %i \t %1.4f \t %1.4f \t %1.4f \t %1.4f \n',n,...
            times_a1(i),times_a4(i),times_a5(i),times_mt(i));
    else
        fprintf(' %i \t %1.4f \t %1.4f \t %1.4f \t %s \n',n,...
            errs_ave_a1_sv(i),errs_ave_a4_sv(i),errs_ave_a5_sv(i),'N/A');
    end
    
end

% Saving data:  Cols errors
%               P-SV|BS-SV|mBS-SV|P-SN|BS-SN|mBS-SN

% Saving data: Cols times
%               P|STD|BS|mBS

errs_ave    = [errs_ave_a1_sv,errs_ave_a4_sv,errs_ave_a5_sv,...
              errs_ave_a1_sn,errs_ave_a4_sn,errs_ave_a5_sn];
          
errs_max    = [errs_max_a1_sv,errs_max_a4_sv,errs_max_a5_sv,...
              errs_max_a1_sn,errs_max_a4_sn,errs_max_a5_sn];          
          
times_all   = [times_a1,times_mt,times_a4,times_a5];

% Preparing data for easier accessing

errs_sv     = [errs_ave_a1_sv,errs_max_a1_sv, errs_ave_a4_sv,errs_max_a4_sv,...
                errs_ave_a5_sv,errs_max_a5_sv];


save('EXPERIMENT_I.mat','errs_ave','errs_max','times_all');

