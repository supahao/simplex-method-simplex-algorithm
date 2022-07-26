function problem = mpsread(mpsfile,varargin)
%MPSREAD reads LP and MILP optimization data from MPS formatted file.
%
% PROBLEM = mpsread(mpsfile) mpsfile is a string containing a file name.
% The file should contain MPS formatted data. On successful file read,
% PROBLEM is a structure that can be passed directly to either intlinprog
% or linprog functions.
%
% PROBLEM = mpsread(mpsfile, NAME, VALUE) specifies additional
% options for mpsread using one or more name-value pair arguments:
%
% 'ReturnNames'             - Returns variable and constraint names present
%                             in the mpsfile as two string arrays in the
%                             PROBLEM structure.
%
%   See also INTLINPROG, LINPROG.

%   Copyright 2015-2020 The MathWorks, Inc.

defaultReturnNames = false;

p = inputParser;
addRequired(p, 'mpsfile');
addParameter(p,'ReturnNames',defaultReturnNames,@(x) islogical(x) && isscalar(x));

parse(p,mpsfile,varargin{:});

mpsfile = p.Results.mpsfile;
if (isstring(p.Results.mpsfile) && isscalar(p.Results.mpsfile))
    mpsfile = char(p.Results.mpsfile);
end

returnNames = p.Results.ReturnNames;

% Read MPS file
[problem.f,intcon,problem.Aineq,problem.bineq,problem.Aeq,problem.beq, ...
    problem.lb,problem.ub,VariableNames,ConstraintNames] = readMPSfile(mpsfile,returnNames);

% Find indices of variables with integer constraints (non-zero).
problem.intcon = find(intcon ~= 0);
unrestricted_rows = [];
if ~isempty(problem.Aineq)
    % Remove unrestricted rows i.e, constraints with infinite RHS.
    % MPS file may have these rows but they are not needed to solve problems.
    unrestricted_rows = isinf(problem.bineq);
    if nnz(unrestricted_rows) > 0
        problem.Aineq = problem.Aineq(~unrestricted_rows,:);
        problem.bineq = problem.bineq(~unrestricted_rows);
    end
end

% Add solver and options fields depending on problem type.
if isempty(problem.intcon)
    problem.solver = 'linprog';
else
    problem.solver = 'intlinprog';
end
problem.options = optimoptions(problem.solver);

if returnNames
    VariableNames = deblank(string(VariableNames));
    emptyCellInd = cellfun(@isempty, ConstraintNames);
    ConstraintNames(emptyCellInd) = {''};
    ConstraintNames = deblank(string(ConstraintNames));
    problem.variableNames = VariableNames;

    % partition of the constraint names
    problem.constraintNames.eqlin = ConstraintNames(1:size(problem.beq));
    problem.constraintNames.ineqlin = ConstraintNames(size(problem.beq)+1:end);

    greaterThanRowLogicals = endsWith(problem.constraintNames.ineqlin,"_g");
    problem.constraintNames.ineqlin(greaterThanRowLogicals) = ...
        deblank(extractBefore(problem.constraintNames.ineqlin(greaterThanRowLogicals), ...
                strlength(problem.constraintNames.ineqlin(greaterThanRowLogicals)) - 1)) + "_g";

    rangeRowLogicals = endsWith(problem.constraintNames.ineqlin,"_u");
    problem.constraintNames.ineqlin(rangeRowLogicals) = ...
        deblank(extractBefore(problem.constraintNames.ineqlin(rangeRowLogicals), ...
                strlength(problem.constraintNames.ineqlin(rangeRowLogicals)) - 1)) + "_u";
    rangeRowLowerLogicals = circshift(rangeRowLogicals,1);
    problem.constraintNames.ineqlin(rangeRowLowerLogicals) = ...
        extractBefore(problem.constraintNames.ineqlin(rangeRowLogicals), ...
                      strlength(problem.constraintNames.ineqlin(rangeRowLogicals)) - 1) + "_l";

    if nnz(unrestricted_rows) > 0
        problem.constraintNames.ineqlin = problem.constraintNames.ineqlin(~unrestricted_rows,:);
    end
end
