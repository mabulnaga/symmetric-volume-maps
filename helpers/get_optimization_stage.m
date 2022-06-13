function [stage] = get_optimization_stage(fixBoundary, count, maxCount, num_times_repair, num_times_uninvert_tets)
%GET_OPTIMIZATION_STAGE: gets which stage of optimization the algorithm is
%at.
%Inputs: 
% fixBoundary: flag, whether the boundary is locked
% count: current optimization count
% maxCount: maximum number of optimization times
% num_times_repair: number of times tetrahedron have been repaired
% num_times_uninvert_tets: number of times to perform the tet repair stage.
%Outputs:
% stage: scalar:
%   stage 0: lock the boundary and optimize the interior
%   stage 1: tetrahedron repair
%   stage 2: optimize all variables.
%   stage 3: post-convergence tetrahedron repair
if(count == 1)
    if(fixBoundary)
        stage = 0;
    else
        stage = 2;
    end
elseif(count > 1)
    if(count == 2 || (mod(count,2) == 0 && num_times_repair < num_times_uninvert_tets))
        %tet inversion stage
        stage = 1;
        %optimization stage
    elseif(count == maxCount)
        stage = 3;
    else
        stage = 2;
    end
end


