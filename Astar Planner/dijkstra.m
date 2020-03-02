function [path, num_expanded] = dijkstra(map, start, goal, astar)
% DIJKSTRA Find the shortest path from start to goal.
%   PATH = DIJKSTRA(map, start, goal) returns an mx3 matrix, where each row
%   consists of the (x, y, z) coordinates of a point on the path. The first
%   row is start and the last row is goal. If no path is found, PATH is a
%   0x3 matrix. Consecutive points in PATH should not be farther apart than
%   neighboring voxels in the map (e.g. if 5 consecutive points in PATH are
%   co-linear, don't simplify PATH by removing the 3 intermediate points).
%
%   PATH = DIJKSTRA(map, start, goal, astar) finds the path using euclidean
%   distance to goal as a heuristic if astar is true.
%
%   [PATH, NUM_EXPANDED] = DIJKSTRA(...) returns the path as well as
%   the number of nodes that were expanded while performing the search.
%   
% paramaters:
%   map     - the map object to plan in
%   start   - 1x3 vector of the starting coordinates [x,y,z]
%   goal:   - 1x3 vector of the goal coordinates [x,y,z]
%   astar   - boolean use astar or dijkstra


if nargin < 4
    astar = false;
end
%% Initialize
status = zeros(1, size(map.occgrid,1)*size(map.occgrid,2)*size(map.occgrid,3));
parent_inds = zeros(1, size(map.occgrid,1)*size(map.occgrid,2)*size(map.occgrid,3));   % for storing parent ind for each node


startnode_ind = pos2ind(map,start);
startnode_pos = start;  % pos of center of voxel
startnode_G = 0;
startnode_H = compute_H(startnode_pos, goal, astar);
startnode_cost = startnode_G + startnode_H;

Nodes = zeros(7,1000);
Nodes(1,1) = startnode_ind;
Nodes(2,1) = startnode_pos(1);
Nodes(3,1) = startnode_pos(2);
Nodes(4,1) = startnode_pos(3);
Nodes(5,1) = startnode_G;
Nodes(6,1) = startnode_H;
Nodes(7,1) = startnode_cost;

parent_inds(startnode_ind) = startnode_ind;
goal_ind = pos2ind(map,goal);

 % initialize the search node list
% Nodes(1000).ind = []; 
% Nodes(1000).pos = [];
% Nodes(1000).G = [];
% Nodes(1000).H = [];
% Nodes(1000).cost = [];


N = 1;

current= Nodes(:,1);

%% Expanding
t_start = tic; % record running time
while current(1) ~= goal_ind
    if N == 0   % No path found
        path = [];
        num_expanded = sum(status==1);
        return
    end
    [Nodes, N, status, parent_inds]= add_neighbors(map, Nodes, N, current, status, parent_inds, goal, astar);
    [Nodes, min_node, N] = delMin(Nodes, N);
    current = min_node;
end

run_time = toc(t_start) 
%% Recover the path 
path = goal;
current_ind = current(1);
while current_ind ~= startnode_ind
     pos = ind2pos(map,current_ind) + map.res_xyz/2;
     path = [pos; path];
     current_ind = parent_inds(current_ind);
end

path = [start;path];
num_expanded = sum(status==1);

plot_path(map,path);

end

function [Nodes, N, status, parent_inds]= add_neighbors(map, Nodes, N, node, status, parent_inds, goal, astar)

node_pos = [node(2), node(3), node(4)];
sub = pos2sub(map, node_pos);
j = sub(1); i = sub(2); k = sub(3);
for a = -1:1:1
    for b = -1:1:1
        for c = -1:1:1
            if a == 0 && b==0 && c==0   % skip if is itself
                continue
            end
            % check if this neighbor is out of bounds
            if j+a<= 0 || j+a> size(map.occgrid,1) 
                continue
            end
            if i+b<= 0 || i+b> size(map.occgrid,2) 
                continue
            end
            if k+c <= 0 || k+c> size(map.occgrid,3) 
                continue
            end
            
            neighbor = zeros(7,1);
            neighbor_pos = sub2pos(map, [j+a,i+b,k+c])+ map.res_xyz/2;  % pos of center of voxel
            neighbor(1) = pos2ind(map,neighbor_pos);           
            neighbor(2) = neighbor_pos(1);
            neighbor(3) = neighbor_pos(2);
            neighbor(4) = neighbor_pos(3);
            
            if map.occgrid(j+a,i+b,k+c) == 1         % skip if occupied by obstacle
                continue
            end
            %check if the neighbor has been visited before, update the neighbor's G if
            %it's a better path through node
            if status(neighbor(1)) == 1
%                 neighbor.ind).= min(Nodes(neighbor.ind).G, node.G+norm(node.pos-neighbor.pos));
%                 Nodes(neighbor.ind).cost = Nodes(neighbor.ind).G + Nodes(neighbor.ind).F;
                continue     
            end
          
            parent_inds(neighbor(1)) = node(1);
            neighbor(5) = node(5) + norm(node_pos-neighbor_pos);
            neighbor(6) = compute_H(neighbor_pos, goal, astar);
            neighbor(7) = neighbor(5) + neighbor(6);
            status(neighbor(1)) = 1;  % mark as visited
            
            if N == size(Nodes,2)
                Nodes = double_size(Nodes);
            end
            [Nodes, N] = insert(Nodes, N, neighbor);  % add this neighbor to the candidate list
        end
    end
end
end


function result = compute_H(current_pos, goal, astar)
if astar == false
    result = 0;
else
    result = norm(current_pos-goal);
end
end

%% Priority queue implementation
function [array, min, N] = delMin(array, N)
min = array(:,1);
array = exch(array, 1, N);
N = N-1;
array = sink(array, 1, N);
array(:,N+1) = []; 

end

function [array, N] = insert(array, N, node)
N = N+1;
array(:,N) = node;
array = swim(array, N);
end

function array = swim(array, k)
while k>1 && less(array, k, floor(k/2))
    array = exch(array, k, floor(k/2));
end
end

function array = sink(array, k, N)
while 2*k <= N
    j = 2*k;
    if j < N && ~less(array, j, j+1)
        j = j+1;
    end
    if less(array, k, j)
        break
    end
    array = exch(array, k, j);
    k = j;
end
end

function array = exch(array, i, j)
temp = array(:,i);
array(:,i) = array(:,j);
array(:,j) = temp;
end

function bool = less(array, i, j)
if array(7,i)< array(7,j)
    bool = true;
else 
    bool = false;
end
end

function Nodes = double_size(Nodes)
n = size(Nodes,2);
new = zeros(7,n);
Nodes = [Nodes, new];
end

