%-------------------------- AOP EXPERIMENT III  ----------------------------%
% 10/28/18, J.B.,
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
% build-in Matlab command 'svd', and the rsvd (randomized svd) as
% described in "Finding Structure With Randomness: Probabilistic Algorithms
% for Constructing Approximate Matrix Decompositions", Halko, Martinsson, 
% Tropp, 2010. The rsvd implementation is from 
% https://www.mathworks.com/matlabcentral/fileexchange/47835-randomized-singular-value-decomposition, 
% by Antoine Liutkus, 2014.

% Versions:
% 10/28/18, J.B., Initial use of the rsvd function.
% 10/29/18, J.B., Preparation for large-scale tests.

%----------------------- Storing of results -----------------------------%
clc;
clear;
 
fprintf('------- AOP: Algorithms For Oblique Projection Matrices ------- \n');
fprintf('-------      J.J.Brust, R.F.Marica, C.G.Petra, 2018     ------- \n \n');

fprintf('EXPERIMENT III: Singular Vector Computations \n');
fprintf('ALGORITHMS: ALG. 2, ALG.3, Matlab, rSVD \n \n');

% Adding rsvd to the path. Requires the current folder to be AOP/tests.
addpath('rsvd');

fprintf('Please enter input A (Times(s)= 1, Error=2). Then hit return. \n');
select_a = input('Input A= ','s');
if isempty(select_a); select_a = '1'; end

fprintf('Please enter input B (Small=1, Large=2). Then hit return. \n');
select_b = input('Input B= ','s');
if isempty(select_b); select_b = '1'; end
issmall = 1;

fprintf('\n Starting EXPERIMENT III with: Input A=%s, Input B=%s \n',select_a, select_b);
fprintf('n     \t DATA(s) \t ALG. 2 \t ALG. 3 \t MTLB.-W \t MTLB.-Wh \t rSVD.-W \t rSVD.-Wh\n');

if select_b ~= '2'
    ns          = [800;1000;2000;3000;4000;5000];    
else
    ns          = [7500;10000;50000;100000;200000;300000;500000;1000000]; % [7500;10000;50000;100000];    
    issmall     = 0;
end

lns         = length(ns);

timeswp     = zeros(lns,1); % Times proposed method for W
timesws     = zeros(lns,1); % Times standard method for W
timeswr     = zeros(lns,1); % Times randomized method for W
timeswhp    = zeros(lns,1); % Times proposed method for \hat{W}
timeswhs    = zeros(lns,1); % Times standard method for \hat{W}
timeswhr    = zeros(lns,1); % Times randomized method for \hat{W}

errswp      = zeros(lns,1); % Errors proposed method for W, average
errsws      = zeros(lns,1); % Errors standard method for W, average
errswr      = zeros(lns,1); % Errors randomized method for W, average
errswhp     = zeros(lns,1); % Errors proposed method for \hat{W}, average
errswhs     = zeros(lns,1); % Errors standard method for \hat{W}, average
errswhr     = zeros(lns,1); % Errors randomized method for \hat{W}, average

errswpm      = zeros(lns,1); % Errors proposed method for W, max
errswsm      = zeros(lns,1); % Errors standard method for W, max
errswrm      = zeros(lns,1); % Errors randomized method for W, max
errswhpm     = zeros(lns,1); % Errors proposed method for \hat{W}, max
errswhsm     = zeros(lns,1); % Errors standard method for \hat{W}, max
errswhrm     = zeros(lns,1); % Errors randomized method for \hat{W}, max

ls_opt_trans.TRANSA = true; % Transpose matrix
ls_opt_utri.UT      = true; % Upper triangular matrix
ls_opt_sym.SYM      = true; % Symmetric matrix

%% Experiment parameters. This experiment uses deterministic matrices
% from the build-in gallery function for small/medium problems, and
% random data for large problems.

nreruns = 10;
if issmall == 1; nreruns = 1; end
     
m           = 20;%150; % Number of columns in X and Y // 4, 300
r           = 8; % Number of columns to compare in the error calculation
del_strt    = 2; % Number of experiments to cut-off from the average 
                 % time/error computations. This is included, because the 
                 % first runs tend to occur start-up times.

tmp_timers  = zeros(lns,2); % Array to store temporary timing data. For small
                            % problems the storing pattern is consistent.
                            % For large problems it skips all columns,
                            % except the first.
                            % tmp_timers(:,1): Generate data,
                            % tmp_timers(:,2): Form orthogonal matrix A
                            % using gallery,
                            % tmp_timers(:,3): Form explicit W,
                            % tmp_timers(:,4): Form explicit \hat{W}

if nreruns <= del_strt; del_strt=0; end

om          = ones(m,m);

%---------------------- Computations for varying 'n'---------------------%
for i = 1:lns
    
    n = ns(i);
    
    for ir = 1:nreruns
        
        % Generate data. This experiment uses the build-in gallery
        % function. In particular the columns of scaled orthogonal matrices 
        % are used. More information on the orthogonal gallery matrices are 
        % accessible via >>help private/orthog
        
        t_data_all  = tic;

        if issmall == 1

            t_data_g    = tic;
            A           = gallery('orthog',n,-1);        
            t_data_g    = toc(t_data_g);

            X           = A(:,1:m);
            d           = sum(X.*X,2);

            D           = diag(1./d);

            Y           = D*A(:,(m+1):(2*m));
            
            tmp_timers(i,2) = t_data_g;
            
        else          
            
            X           = randn(n,m);
            Y           = randn(n,m);
            
        end
        
        t_data_all      = toc(t_data_all);        
        tmp_timers(i,1) = t_data_all;
        
        YX              = Y'*X;
        YY              = Y'*Y;
        XX              = X'*X;
        XYm             = [X Y];
                
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
        
        %% Randomized algorithms for W. Two type of methods are used:
        % (1) rsvd (for small/medium matrices, because it forms the full
        % projection W). (2) rsvd_w,rsvd_wh (for large matrices, because
        % they exploit rectangular matrix products [cf. rsvd_mod/rsvd_w,rsvd_wh]).
        
        if issmall == 1
        
            %% Computing W
            % Forming the full oblique projection matrix W
            t_wf        = tic; % Timer to form explicit W
            W           = X*(linsolve(YX,Y'));
            t_wf        = toc(t_wf);

            tmp_timers(i,3) = t_wf;

            %% Randomized svd of W. An SVD approximation is computed by fixing
            % the rank to be K = m.

            t_rsvdw     = tic;
            [Ur,Sr,Vr]  = rsvd(W,m);
            t_rsvdw     = toc(t_rsvdw);

            %% Computing \hat{W} 
            % Explicit forming of the complement oblique projection
            % \hat{W} = I - W

            t_whf   = tic;
            Wh      = spdiags(1-diag(W(:,:)),0,-W(:,:));
            t_whf   = toc(t_whf);

            tmp_timers(i,4) = t_whf;

            %In         = eye(n);
            %Wh         = In - X*(linsolve(YX,Y'));

            %% Randomized svd of \hat{W}. An SVD approximation by
            % fixing the rank K = n - m
            K               = n-m;

            t_rsvdwh        = tic;
            [Uhr,Shr,Vhr]   = rsvd(Wh,K);
            t_rsvdwh        = toc(t_rsvdwh);
        
        else
            
            %% Modified randomized svd of W. An SVD approximation is computed by fixing
            % the rank to be K = m.
            K               = m;
             
            t_rsvdw         = tic;
            [Ur,Sr,Vr]      = rsvd_w(X,Y,YX,K);
            t_rsvdw         = toc(t_rsvdw);

            %% Modified randomized svd of \hat{W}. An SVD approximation by
            % fixing the rank K = m. Expect the errors to be larger by
            % this choice of K, because the full singular value spectrum
            % of \hat{W} contains n-m nonzero values.
            K               = m;

            t_rsvdwh        = tic;
            [Uhr,Shr,Vhr]   = rsvd_wh(X,Y,YX,K);
            t_rsvdwh        = toc(t_rsvdwh);
            
        end
        
        % W rsvd data
        if del_strt < ir && del_strt < nreruns
            % rsvd time
            timeswr(i)      = timeswr(i) + t_rsvdw;
                
            % rsvd error
            WscIr           = Ur*(Sr*Vr(idx,:)');
            es              = norm(WIr - WscIr,'fro');

            errswr(i)       = errswr(i) + es;
            if errswrm(i) < es; errswrm(i) = es; end
            
        end
        
        % \hat{W} rsvd data
        if del_strt < ir && del_strt < nreruns
        	% rsvd time
            timeswhr(i)     = timeswhr(i) + t_rsvdwh;
                
            % rsvd error
            WhscIr          = Uhr*(Shr*Vhr(idx,:)');
            ehs             = norm(WhIr - WhscIr,'fro');
                
            errswhr(i)      = errswhr(i) + ehs;
            if errswhrm(i) < ehs; errswhrm(i) = ehs; end
            
        end
        
        %% Matlab standard computations using the function svd
        if issmall == 1
            
            %% For W
            %W               = X*(linsolve(YX,Y'));
            
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
    timeswr(i)      = timeswr(i)/denom;
    timeswhp(i)     = timeswhp(i)/denom; 
    timeswhs(i)     = timeswhs(i)/denom;
    timeswhr(i)     = timeswhr(i)/denom;
    
    errswp(i)       = errswp(i)/denom;
    errsws(i)       = errsws(i)/denom;
    errswr(i)       = errswr(i)/denom;
    errswhp(i)      = errswhp(i)/denom;
    errswhs(i)      = errswhs(i)/denom;
    errswhr(i)      = errswhr(i)/denom;
    
    if issmall == 1
        
        if select_a ~= '2'
        
            fprintf('%i\t %1.4f \t %1.4f \t %1.4f \t %1.4f \t %1.4f \t %1.4f \t %1.4f \n',...
                n,t_data_all,timeswp(i),timeswhp(i),timesws(i),timeswhs(i),timeswr(i),timeswhr(i));
                
        else
            
            fprintf('%i\t %1.4f \t %1.4f \t %1.4f \t %1.4f \t %1.4f \t %1.4f \t %1.4f \n',...
                n,t_data_all,errswp(i),errswhp(i),errsws(i),errswhs(i),errswr(i),errswhr(i));
            
        end
        
    else
        
       if select_a ~= '2'
        
            fprintf('%i\t %1.4f \t %1.4f \t %1.4f \t %s \t %s \t %1.4f \t %1.4f \n',n,...
                t_data_all,timeswp(i),timeswhp(i),'------','------',timeswr(i),timeswhr(i));
            
        else
            
            fprintf('%i\t %1.4f \t %1.4f \t %1.4f \t %s \t %s \t %1.4f \t %1.4f \n',n,...
                t_data_all,errswp(i),errswhp(i),'------','------',errswr(i),errswhr(i));
            
        end
        
    end
end

if issmall == 1
    
    times_all   = [timeswp,timesws,timeswr,timeswhp,timeswhs,timeswhr];
    errs_all    = [errswp,errsws,errswr,errswhp,errswhs,errswhr,...
                   errswpm,errswsm,errswhpm,errswhsm];
               
    data_tbl_w  = [times_all(:,1:3),errs_all(:,1:3)];
                
    data_tbl_wh = [times_all(:,4:6),errs_all(:,4:6)];
    
    save('EXPERIMENT_III_S.mat','times_all','errs_all','data_tbl_w','data_tbl_wh');
    
else
    
    times_all   = [timeswp,timeswr,timeswhp,timeswhr];
    errs_all    = [errswp,errswpm,errswr,errswrm,errswhp,errswhpm,errswhr,errswhrm];
    
    data_tbl    = [times_all,errs_all];
    
    save('EXPERIMENT_III_L.mat','times_all','errs_all','data_tbl');
    
end

