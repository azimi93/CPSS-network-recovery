
% INPUT OUTPUT Specifications:
% barA: the system matrix (it must be a square matrix)
% barB: the input matrix of any dimension
%
% Output:
% optVal: the cardinality of an optimal perturbation configuration
% UL_opt_out: the set of leftunmatched state vertices in the maximum
% matching attaining USAN
% UR_opt_out: the set of rightunmatched vertices
% A_opt_perturb: the optimal perturbation configuration where the added
% edges are represented by the non-negative entries
function [optVal, UL_opt_out, UR_opt_out, A_opt_pertub] = OptimalConfiguration(barA, barB)
%It receives the nxn binary matrix barA and nxm binary matrix barB and returns a minimal Feasible
%Perturbation Configuration (FPC)
sparseAbar=sparse(barA);    %Represents the adjacency matrix as edge list
[numberSCCs,SCCsLabels] = graphconncomp(sparseAbar);
dimA = size(barA, 1);       % dimension of the matrix
dimB = size(barB, 2);       % the number of inputs
%To determine the SCC of a digraph representation of barA, i.e., D(barA).
%More precisely, numberSCCs gives the number of different SCCs in D(barA),
%whereas SCCsLabels gives a vector where each entry i states the SCC that a
%state variable i belongs to.
SCCs={};
for i=1:numberSCCs
    SCCs{i}=find(SCCsLabels==i);
end
%SCCs{i} comprises the collection of labels of state variables that belong to
%the same SCC i.
%To determine the non-top linked SCCs we identify which SCCs have no
%path ending in one of its vertices from the vertices in other SCCs (notice that only one vertex from an SCC is 
%enough to perform the test, because in an SCC all vertices are reached from the others)
%% STEP1 decompose intor reachable and unreachable SCCs and determine the source
SourceSCCsIndex=[]; % keeps the indices i of the SCCs{i} that are source SCCs
for i=1:numberSCCs
    minPathToSCCj=inf;
    for j=i+1:numberSCCs
        minPathFromiToj=graphshortestpath(sparseAbar,SCCs{j}(1),SCCs{i}(1));
        %Computes the minimum path from SCCs{j}(1) to SCCs{i}(1), if there
        %is none the minimum path is inf
       if minPathFromiToj< minPathToSCCj 
           minPathToSCCj=minPathFromiToj;
       end
       
    end
    if minPathToSCCj==inf
        %Tests if there exists a finite path from some SCC to SCCs{i}, if
        %not, then SCCs{i} is a non-top linked SCC
        SourceSCCsIndex=[SourceSCCsIndex i];
    end
end
%Re-label of the beta (i.e.,numberNonTopLinkedSCCs) non-top linked SCCs to
%have labels from 1 to beta.
numberSourceSCCs=length(SourceSCCsIndex);
SourceSCCs={};
for i=1:numberSourceSCCs
    SourceSCCs{i}=SCCs{SourceSCCsIndex(i)};
end
ReachableSCCs = [];
inputConnectedStates = find(sum(barB,2) ~= 0);
for i = 1:length(inputConnectedStates)
    if ~ismember(SCCsLabels(inputConnectedStates(i)),ReachableSCCs) %remove redundancy
        ReachableSCCs = [ReachableSCCs SCCsLabels(inputConnectedStates(i))]; % states belongs to these SCCs are input-reachable
    end
end % this calculates the reachable source SCCs only
UnreachableSourceSCCs = setdiff(SourceSCCsIndex,ReachableSCCs);
% the unreachable source SCCs are those satisfying two criteria: 1. they
% are source sccs, 2. they are unreachable, therefore they are those
% SourceSCCs removing the reachable SCCs. 
% UnreachableSourceSCCs specifies the index of SCCs
numUnSourceSCCs = length(UnreachableSourceSCCs) % number of unreachable source SCCs£¬ i.e., the parameter r in the paper.
%% STEP1 + Buiding the cost matrices for MWMM
% building up the whole stream of reachable SCCs
for i=1:length(ReachableSCCs)
    for j=1:numberSCCs
       minPathFromiToj=graphshortestpath(sparseAbar,SCCs{ReachableSCCs(i)}(1),SCCs{j}(1));
       if minPathFromiToj<inf
            %Tests if there exists a finite path from some SCC to SCCs{i}, if
            %not, then SCCs{i} is a non-top linked SCC
            if ~ismember(j,ReachableSCCs)
                ReachableSCCs = [ReachableSCCs j];
            end
        end
    end
end % this contains all the reachable SCCs
numofreachablestates = 0;
Setofreachablestates = [];
for i = 1:numel(ReachableSCCs)
   numofreachablestates = numofreachablestates + numel(SCCs{ReachableSCCs(i)});
   Setofreachablestates = [Setofreachablestates SCCs{ReachableSCCs(i)}];
end
% calculates the number of states reachable from input
%The bipartite graph is described by the concatenated matrix [A I]
%where I comprises the columns with non-zero entries from the variable i to
%xj if xj belongs to the non-top linked SCC i
UnreachAug=inf*ones(numUnSourceSCCs, dimA); % augmentation on the unreachable source SCCs
UnreachSourceStates = []; % finds all unreachable states that are in unreachable state source SCCs
if numUnSourceSCCs % if not zero
    for i=1:numUnSourceSCCs
       UnreachAug(i, SCCs{UnreachableSourceSCCs(i)})=2; % one slack variable pointing to all
       UnreachSourceStates = [UnreachSourceStates SCCs{UnreachableSourceSCCs(i)}];
    end
end
%Cost in barA is inf in their zeros, and unitary cost in there non-zero
%entries. Similarly, the non-zero entries of I have cost 2 and zero entries
%cost inf.
costA=barA;
costA((costA==0))=inf;
costB=barB;
costB((costB==0))=inf;
costB = costB';
costGraph = [costA inf*ones(dimA, dimB); costB inf*ones(dimB, dimB)]; % augment into cost of graph
%wb_rightunmatched has the following structure
%[X->X X->U; U->X U->U; Slack->X Slack->U (which is inf)]
wb_rightunmatched=[costGraph;[UnreachAug inf*ones(numUnSourceSCCs,dimB)]]; % 
%wb_leftunmatched has the following property
%[costA costB(by a vector) n*r matrix for augmentation process]
%% STEP2 & 3: Find UL* and UR* through MWMM
% Assignment rules: states are indexed from 1~n, while input is at n+1, the
% others are slack variables
%Hungarian algorithm to determine the maximum matching incuring in the
%minimum cost is computed using the m-file assignmentoptimal.
%It uses an algorithm suggested by Munkers and implmented by Markus Buehren,
assignmentright = assignmentoptimal(wb_rightunmatched); 
%% Find a maximum matching with the specified left and right-unmatched vertices
UL_opt = find(assignmentright(1:dimA)==0 | assignmentright(1:dimA)>dimA);                           % the set of left unmatched state vertices
UR_opt = union(setdiff(0:dimA,assignmentright(1:dimA+dimB)), assignmentright(dimA+dimB+1:end));     % the set of right-unmatched ~
UR_opt(UR_opt==0) = []; %remove zero
%% if the set of left-unmatched vertices doesn't contain reachable states
if isempty(intersect(UL_opt, Setofreachablestates))
	inputind = inputConnectedStates(1);% the input connnecting to the state, also the state number 
    assignmentright(assignmentright==inputind) = 0; % the state that is connected to inputind in the maximum matching 
    assignmentright(dimA + inputind) = inputind;
end
UL_opt = find(assignmentright(1:dimA)==0 | assignmentright(1:dimA)>dimA);   
% change the maximum matching such that there is at least one reachable
% left-unmatched state vertex
%% Find out the right-unmatched vertex associated with the starting point
NumRightUnmatched = numel(UR_opt) % number of rightunmatched vertices (nr)
AssignabilityIndex = sum((assignmentright(dimA+dimB+1:end)~=0)) %(q)
optVal = NumRightUnmatched + numUnSourceSCCs - AssignabilityIndex; % the minimum number of required perturbations - q in the paper
UL_opt_out = UL_opt; % make a copy for output
UR_opt_out = UR_opt; % make a copy to output
%The output assignment gives a column vector of size n where the row index
%is the state variable and the corresponding entry gives the column that as
%been assigned. Therefore the edges in the maximum matching are given by 
%(assignment(label),label). Therefore, if an assignment(label)==0 it means
%that there is no assignment.
%The entries in assignment(label) that are greater than n (size of state
%space) correspond to the root of an edge in E_{I,X} hence the state
%variable indexed by the line label in assignment is a right-unmatched
%vertex in a non-top linked SCC.
reverseassign = zeros(dimA + dimB, 1);
for original_index = 1:dimA+dimB
    reversed_index = assignmentright(original_index);
    if reversed_index~= 0
        reverseassign(reversed_index) = original_index;
    end
end
%% Find the right unmatched vertex that is the start of the path
UL_reachable = intersect(UL_opt, Setofreachablestates);
startvtx = UL_reachable(1); % select one left-unmatched vertex that is reachable injast k mishe bejaie 1 halate digar ro gozasht
startpt = startvtx; % a pointer used to trace the path
while 1
    EndRightunmatched = reverseassign(startpt);
    if EndRightunmatched == 0
        break;
    else
        startpt = EndRightunmatched; % the startpt is the right unmatched vertex
    end
end
%% Optimal perturbation configuration
% Set of reachable states:  Setofreachablestates
% Set of unreachable states in source SCCs: UnreachSourceStates
A_opt_pertub = zeros(dimA, dimA);
LeftoverUnreach = UnreachableSourceSCCs; % the set of SCCs that are unreachable
if ismember(startpt,UnreachSourceStates) % if it is a right-unmatched source state
    UR_opt(UR_opt == startpt) = [];
    RightUnmatchedSourceUnreachable = intersect(UR_opt, UnreachSourceStates);
    if ~isempty(RightUnmatchedSourceUnreachable)
        leftsel = startvtx;
        rightsel = RightUnmatchedSourceUnreachable(1);
        while true
            A_opt_pertub(leftsel, rightsel) = 1; % assign
            UL_opt(UL_opt == leftsel) = [];                                 % remove added left-unmatched vertex from sets UL_opt.
            UR_opt(UR_opt == rightsel) = [];                                % remove added right-unmatched vertex from sets UR_opt.
            RightUnmatchedSourceUnreachable(RightUnmatchedSourceUnreachable == rightsel)=[];
            LeftoverUnreach(LeftoverUnreach == SCCsLabels(rightsel)) = [];  % unreachable source state SCCs become reachable
            leftsel = TraceLeftUnmatched(assignmentright,rightsel);
            if isempty(RightUnmatchedSourceUnreachable)
                A_opt_pertub(leftsel, startpt) = 1;
                UL_opt(UL_opt == leftsel) = [];   
                LeftoverUnreach(LeftoverUnreach == SCCsLabels(startpt)) = [];
                break; % if there are no leftover unreachable source right unmatched vertex then break
            else
                rightsel = RightUnmatchedSourceUnreachable(1);
            end
        end
    else % there is only one such right unmatched source unreachable 
        A_opt_pertub(startvtx, startpt) = 1;
        UL_opt(UL_opt == startvtx) = [];   
        UR_opt(UR_opt == startpt) = [];   
        LeftoverUnreach(LeftoverUnreach == SCCsLabels(startpt)) = [];
    end
else % if it is not a right-unmatched source state
    leftsel = startvtx;
    if ~isempty(intersect(UR_opt, UnreachSourceStates))
        RightUnmatchedSourceUnreachable = intersect(UR_opt, UnreachSourceStates);
        rightsel = RightUnmatchedSourceUnreachable(1);
        while true
            A_opt_pertub(leftsel, rightsel) = 1; % assign
            UL_opt(UL_opt == leftsel) = [];                                 % remove added left-unmatched vertex from sets UL_opt.
            UR_opt(UR_opt == rightsel) = [];                                % remove added right-unmatched vertex from sets UR_opt.
            RightUnmatchedSourceUnreachable(RightUnmatchedSourceUnreachable == rightsel)=[];
            LeftoverUnreach(LeftoverUnreach == SCCsLabels(rightsel)) = []; % unreachable source state SCCs become reachable
            leftsel = TraceLeftUnmatched(assignmentright,rightsel);
            if isempty(RightUnmatchedSourceUnreachable)
                break; % if there are no leftover unreachable source right unmatched vertex then break
            else
                rightsel = RightUnmatchedSourceUnreachable(1);
            end
        end
    end
end
% connect the remaining right-unmatched vertices
while ~isempty(UR_opt)
    A_opt_pertub(UL_opt(1), UR_opt(1)) = 1; % assign
    UL_opt(1) = []; % remove added left-unmatched vertex from sets UL_opt.
    UR_opt(1) = [];	% remove added right-unmatched vertex from sets UR_opt.
end
% Connect to the unreachable sccs
if ~isempty(LeftoverUnreach)
   for i = 1:numel(LeftoverUnreach) % the set of indices on unreachable source SCCs
      A_opt_pertub(Setofreachablestates(1), SCCs{LeftoverUnreach(i)}(1)) = 1;
   end
end
%% Verification
% to check whether it is valid
assert(sum(sum(A_opt_pertub)) == optVal, 'Error: optval value does not match with theoretical value');
end
