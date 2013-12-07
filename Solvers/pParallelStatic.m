%% Configure algorithm
this    = pMPSOLVER;
maxNIt  = this.MaxNewtonIterations.Value;
nrrTol  = this.NewtonTolerance.Value;
verbose = this.Verbose;

%% Get times points
t    = this.Times;
L    = length(this.Matrices.f(0));
Tmap = map([1 Np], {}, 0:Np-1);
x    = zeros(L, numel(t), Tmap);
tInd = global_ind(x, 2);
Nt   = numel(tInd);
xLoc = local(x);

%% For each time point,
for i = 1:Nt
    j = tInd(i);
    %% Create constant matrices,
    f = this.Matrices.f(t(j));
    K = this.Matrices.K(t(j));
    
    %% Get initial guess,
    if i > 1
        xLoc(:,i) = xLoc(:,i-1);
    end
    
    %% Perform Newton-Raphson iteration,
    nIter  = 1;
    relRes = 1;
    
    while nIter < maxNIt && relRes > nrrTol
        [G, g] = this.Matrices.G(t(j), xLoc(:,i));
        
        r = K * xLoc(:,i) + g - f;
        J = K + G;
        
        if i == 1 && nIter == 1
            [~,~,p,S] = ldl(J, 'vector');
            [~, q]    = sort(p);
            spparms('autoamd', 1);
        end
        
        %Explicit symmetrization ensures MA57 is called for ldl decomposition
        %SYMAMD preording is precalculated since all matrices have the same
        %structure. spparms('autoamd',0) turns of this stage in the sparse solver.
        
        J  = (J + J.');
        r  = 2 * r(p);
        J  = J(p,p);
        dx = J \ r;
        dx = dx(q);
        
        %Equivalent direct call to ldl. Tends to be slower than the explicit
        %symmetrization above. May be faster if (J+J.') is expensive to compute.
        
        %                     r       = r(p);
        %                     J       = J(p,p);
        %                     [L,D,P] = ldl(J);
        %                     dx      = (P' \ (L' \ (D \ (L \ (P \ r)))));
        %                     dx      = dx(q);
        
        xLoc(:,i) = xLoc(:,i) - dx;
        nIter  = nIter + 1;
        relRes = norm(S * r) / norm(S * (g - f));
        
        %                     this.ConvergenceHistory(end+1) = relRes;
    end
end
spparms('autoamd', 1);

x = put_local(x,xLoc);
toc
%% Save Solution
this.X = agg(x);
pMPSOLVER = this;