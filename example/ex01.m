% This M-file shows how to use the MATAMG toolbox.

% Build a matrix A
stencil = [4 -1 -1 -1 -1];
n = 64;
A = stencil2mat(stencil,n);
b = A * ones(length(A),1);

% Generate MATAMG option file
opt = amgset;                               % default option file
opt = amgset(opt,'coarsest',10);            % set the number of levels
opt = amgset(opt,'PreCond','pcg');          % set the Krylov method
opt = amgset(opt,'PrintOnScreen','off');    % turn off the option to print the log on the screen 
opt = amgset(opt,'SaveCsn','on');           % save the set of coarse-grid points
opt = amgset(opt,'CsnType','amg');          % choose the coarsening method

% Initial vector
v = rand(length(A),1);

% Solve a linear system
v = amg(A,v,b,opt);


