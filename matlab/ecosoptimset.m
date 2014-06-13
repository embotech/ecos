function opts = ecosoptimset(varargin)
% Get options struct for ECOS(QP) solver.
%
% OPTIONS = ECOSOPTIMSET returns a pre-initialized options struct for ECOS.
%
% OPTIONS = ECOSOPTIMSET(PARAMETER1,VALUE1,...,PARAMETERN,VALUEN) returns an
% inialized options struct for ECOS with the parameters PARAMETER set to
% VALUE. Multiple parameters have to be defined in PARAMETER-VALUE pairs.
% PARAMETER must be a string of one of the following keywords:
%
%    feastol, reltol, abstol, maxit and verbose
%
% (c) Alexander Domahid, embotech GmbH, Zurich, Switzerland, 2014.


if( mod(nargin,2) ~= 0)
    error('A parameter list with string-value pairs is expected as argument')
end

% default entries
opts.feastol = 1e-6;
opts.reltol = 1e-7;
opts.abstol = 1e-7;
opts.maxit = 50;
opts.verbose = 1;

% set entries as requested by user
for i = 1:2:nargin
    assert(ischar(varargin{i}),'Parameter %d is not a string');
    opts.(varargin{i}) = varargin{i+1};
end

