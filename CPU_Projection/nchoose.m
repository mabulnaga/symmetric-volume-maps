function W = nchoose(S, I)
% NCHOOSE - all combinations of the elements of a set
%   W = nchoose(S) returns all possible combinations of 0, 1, or more
%   elements of the set S, having N elements.  There are 2^N combinations
%   in total.  W is a cell array and each cell holds one of these
%   combination (as a row vector). 
%   S can be a cell array, and each cell of W will then contain a cell
%   array. W is the powerset of S, as it includes the empty set (0
%   elements) as it first cell.
%
%   For a vector of integers I, W = nchoose(S, I) returns only the sets
%   indicated by the indices I. This might be useful for large sets.
%
%   Examples:
%      nchoose([2 4 6 8])
%        % -> { [] ;
%        %      [2] ; 
%        %      [4] ; 
%        %      [2 4] ; 
%        %      [6] ;
%        %      ...
%        %      [2 6 8] ;
%        %      [4 6 8] ;
%        %      [2 4 6 8]} ; % in total 16 different combinations
% 
%      nchoose([33 22 11], [1 8 4]) 
%        % -> { [] ; [33 22 11] ; [ 33 11]}
%
%   Notes: 
%   - For sets containing more than 18 elements a warning is given, as this
%     can take some time.  Hit Ctrl-C to intterupt calculations.
%   - If S contain non-unique elements (e.g. S = [1 1 2]), NCHOOSE will
%     return non-unique cells. In other words, NCHOOSE treats all elements
%     of S as being unique. One could use NCHOOSE(UNIQUE(S)) to avoid that.
%   - Loosely speaking, NCHOOSE(S) collects all output of multiple calls to
%     NCHOOSEK(S, K) where K is looping from 1 to the number of elements of
%     S. The implementation of NCHOOSE, however, does rely of a different
%     method and is much faster than such a loop.
%   - For more information, see: http://en.wikipedia.org/wiki/Power_set
%   
%   See also NCHOOSEK, PERMS, 
%            PERMN, NCHOOSE2, ALLCOMB on the File Exchange

% Written in Matlab R2018a
% version 3.0 (feb 2018)
% (c) Jos van der Geest
% email: samelinoa@gmail.com

% History
% 1.0, may 2008 - inspired by a post on CSSM
% 2.0, may 2008 - improvemed algorithm
% 2.1, mar 2010 - added a note on "power set", updated ML version
% 2.2, feb 2016 - updated contact info
% 3.0, feb 2018 - changed function to return powerset, implemented the
%                 return of only certain subsets

% Acknowledgements:
% 2.0: Urs Schwarz, for suggesting a significant speed improvement using bitget

narginchk(1,2)
numS = numel(S) ; 
S = reshape(S, 1, numS) ;   % make the set a row vector, for uniform output

Ncomb = 2^numS ;
switch nargin
    case 1
        numW = Ncomb ;
        I = 1:Ncomb ;
    case 2
        if ~isnumeric(I) || ~all(fix(I) == I) || any(I < 1) || any (I > Ncomb)
            error('nchoose:InvalidI', ...
                'Values of K should be integers between 1 and %d.', Ncomb) ;
        end
        numW = numel(I) ;
end

% How many unique combinations of 1 to N elements are there
if numS > 18
    warning('nchoose:LargeOutput', ...
        'There are %d combinations. Please be patient ...', numW) ;
end

% The selection of elements is based on the binary representation of a
% vector of numbers V, which can be calculated using the formula (see
% DEC2BIN):
%    B1 = rem(floor(V(:) * pow2(1-N:0)), 2) 
% (* means matrix multiplication). B(k,:) represents V(k)
% Based on the suggestion by US, this formula can be rewritten as:
%    B2 = bitget(X * (2.^(N-1:-1:0)),N) > 0
% for a specific number X. Although the bits of B2 are in reversed order
% (compared to B1), this is, however, exactly what we want, as we can now
% select the elements of S directly. We do not need to flip the order
% afterwards. 

% Looping over the numbers is faster than creating a whole matrix at once
W = cell(numW, 1) ; % Pre-allocation of output
p2 = 2.^(numS-1:-1:0) ; % This part of the formula can be taken out of the loop
for k = 1:numel(I)
    % calculate the (reversed) binary representation of i
    % select the elements of the set based on this representation
    W{k} = S(bitget((I(k)-1)*p2, numS) > 0) ; 
end

% Two other options, which are much slower and/or efficient
% % 1) looping over nchoosek (slow)
% for i=N:-1:1,
%     W{i} = nchoosek(S,i) ;
% end
% % and now each cell may still have to be split into distinct cells ...
% 
% % 2) in one big step (memory inefficient)
% Q = rem(floor([1:M].'*pow2(1-N:0)),2).' ;
% [R,C] = find(Q) ;
% W = mat2cell(S(R),1,sum(Q)).' ;