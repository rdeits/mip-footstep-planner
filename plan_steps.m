lb = [-0.2; 0.18; -0.2; 0];
ub = [0.4; 0.4; 0.2; pi/6];
A_reach = [-diag(ones(4,1)); diag(ones(4,1))];
b_reach = [-lb; ub];
right_foot_lead = true;
nom_forward_step = 0.2;
nom_step_width = 0.26;
nom_step = [0; nom_step_width; 0; 0];
goal_pos = struct('right', [50;-0.15;0.1;0;0;0],...
                  'left',  [50; 0.15;0.1;0;0;0]);
start = [0;0.15;0;0;0;0];
R = [rotmat(-start(6)), zeros(2,2);
     zeros(2,2), eye(2)];
w_goal = 10 * [1;1;1;0;0;0];
w_rel = 1.2 * [1/nom_forward_step^2;1;1;0;0;0];
nsteps = 10;

stones = [0, 0.15, 0;
          0.5, 0.15, 0.01;
          0.5, -0.15, 0.02;
          1, 0.15, 0.03;
          1, -0.15, 0.04;
          1.5, -0.15, 0.05;
          1.5, 0.15, 0.06;
          2, 0, 0.07]';
safe_regions = struct('A', {}, 'b', {}, 'pt', {}, 'normal', {});
% for j = 1:size(stones, 2)
%   [Ai, bi] = poly2lincon(stones(1,j) + [-.15, -.15, .15, .15],...
%                          stones(2,j) + [-.1, .1, .1, -.1]);
%   Ai = [Ai, zeros(size(Ai, 1), 1)];
%   safe_regions(end+1) = struct('A', Ai, 'b', bi, 'pt', [0;0;stones(3,j)], 'normal', [0;0;1]);
% end
[Ai, bi] = poly2lincon([-.2, -.2, 5,5], [-.4, .4, .4, -.4]);
Ai = [Ai, zeros(size(Ai, 1), 1)];
safe_regions(1) = struct('A', Ai, 'b', bi, 'pt', [0,0,0], 'normal', [0,0,1]);
nr = length(safe_regions);

nx = 4 * nsteps;
ns = nsteps * nr;
nvar = nx + ns;

% Normalize the goal weight so that the plans don't stretch out as the goal
% gets farther away
goal_pos.center = mean([goal_pos.right, goal_pos.left],2);
dgoal = norm(goal_pos.center(1:2) - start(1:2));
extra_distance = max(dgoal - (nsteps - 1) * nom_forward_step, 0.01);
w_goal(1:2) = w_goal(1:2) * sqrt(1 / (extra_distance));

if ~right_foot_lead
  r_ndx = 1:2:nsteps;
  l_ndx = 2:2:nsteps;
else
  r_ndx = 2:2:nsteps;
  l_ndx = 1:2:nsteps;
end

x_ndx = reshape(1:nx, 4, nsteps);
s_ndx = reshape(nx + (1:ns), nr, nsteps);

A = [];
b = [];
Aeq = [];
beq = [];
Q = zeros(nvar, nvar);
c = zeros(nvar, 1);
lb = -inf(nvar, 1);
ub = inf(nvar, 1);

for j = 2:nsteps
  Ai = zeros(size(A_reach, 1), nvar);
  if ismember(j, r_ndx)
    rA_reach = A_reach * diag([1,-1,1,-1]) * R;
  else
    rA_reach = A_reach * R;
  end
  Ai(:,x_ndx(:,j)) = rA_reach;
  Ai(:,x_ndx(:,j-1)) = -rA_reach;
  bi = b_reach;
  A = [A; Ai];
  b = [b; bi];
end

w_goal = diag(w_goal([1,2,3,6]));
Q(x_ndx(:,r_ndx(end)), x_ndx(:,r_ndx(end))) = w_goal * w_goal';
xg = reshape(goal_pos.right([1,2,3,6]), [], 1);
c(x_ndx(:,r_ndx(end))) = -2 * xg' * w_goal * w_goal';
Q(x_ndx(:,l_ndx(end)), x_ndx(:,l_ndx(end))) = w_goal * w_goal';
xg = reshape(goal_pos.left([1,2,3,6]), [], 1);
c(x_ndx(:,l_ndx(end))) = -2 * xg' * w_goal * w_goal';

w_rel = diag(w_rel([1,2,3,6]));
% TODO: make this relative to nominal step width
for j = 2:nsteps
  Q(x_ndx(:,j), x_ndx(:,j)) = Q(x_ndx(:,j), x_ndx(:,j)) + R' * w_rel * w_rel' * R;
  Q(x_ndx(:,j-1), x_ndx(:,j)) = Q(x_ndx(:,j-1), x_ndx(:,j)) - R' * w_rel * w_rel' * R;
  Q(x_ndx(:,j), x_ndx(:,j-1)) = Q(x_ndx(:,j), x_ndx(:,j-1)) - R' * w_rel * w_rel' * R;
  Q(x_ndx(:,j-1), x_ndx(:,j-1)) = Q(x_ndx(:,j-1), x_ndx(:,j-1)) + R' * w_rel * w_rel' * R;
  
  if ismember(j, r_ndx)
    nom = diag([1,-1,1,-1]) *nom_step;
  else
    nom = nom_step;
  end
  c(x_ndx(:,j)) = c(x_ndx(:,j)) - (2 * nom' * w_rel * w_rel' * R)';
  c(x_ndx(:,j-1)) = c(x_ndx(:,j-1)) + (2 * nom' * w_rel * w_rel' * R)';
end

for j = 1:nsteps
  Aeqi = zeros(1, nvar);
  Aeqi(1, s_ndx(:,j)) = 1;
  beqi = 1;
  Aeq = [Aeq; Aeqi];
  beq = [beq; beqi];
end

M = 1000;
for j = 1:nsteps
  for r = 1:nr
    A_region = safe_regions(r).A;
    b_region = safe_regions(r).b;
    A_region = [A_region(:,1:2), zeros(size(A_region, 1), 1), A_region(:,3)];
    A_region = [A_region; 
                reshape(safe_regions(r).normal, 1, []), 0;
                -reshape(safe_regions(r).normal, 1, []), 0];
    b_region = [b_region; 
                dot(safe_regions(r).normal, safe_regions(r).pt);
                -dot(safe_regions(r).normal, safe_regions(r).pt)];
    
    Ai = zeros(size(A_region, 1), nvar);
    Ai(:,x_ndx(:,j)) = A_region;
    Ai(:,s_ndx(r,j)) = M;
    bi = b_region + M;
    A = [A; Ai];
    b = [b; bi];
  end
end

lb(x_ndx(:,1)) = start([1,2,3,6]);
lb(x_ndx(4,:)) = start(6);
ub(x_ndx(:,1)) = start([1,2,3,6]);
ub(x_ndx(4,:)) = start(6);

clear model params
model.A = sparse([A; Aeq]);
model.obj = c;
model.sense = [repmat('<', size(A,1), 1); repmat('=', size(Aeq, 1), 1)];
model.rhs = [b; beq];
model.lb = lb;
model.ub = ub;
model.vtype = [repmat('C', nx, 1); repmat('B', ns, 1);];
model.Q = sparse(Q);
params = struct();
params.timelimit = 5;

result = gurobi(model, params)
xstar = result.x;
steps = xstar(x_ndx);
diff(steps, 1, 2)
S = xstar(s_ndx);
figure(1);
clf
% quiver(steps(1,:), steps(2,:), ones(1,size(steps, 2)), zeros(1,size(steps, 2)));
plot(steps(1,r_ndx), steps(2, r_ndx), 'bo')
hold on
plot(steps(1,l_ndx), steps(2,l_ndx), 'ro')
plot(steps(1,:), steps(2,:), 'k:')
for j = 1:length(safe_regions)
  V = iris.thirdParty.polytopes.lcon2vert(safe_regions(j).A(:,1:2), safe_regions(j).b);
  k = convhull(V(:,1), V(:,2));
  patch(V(k,1), V(k,2), 'k', 'FaceAlpha', 0.2);
end
axis equal