lb = [-0.2; 0.18; -0.2; 0];
ub = [0.4; 0.4; 0.2; pi/6];
A_reach = [-diag(ones(4,1)); diag(ones(4,1))];
b_reach = [-lb; ub];
right_foot_lead = true;
nom_forward_step = 0.2;
goal_pos = struct('right', [2;-0.15;0.1;0;0;0],...
                  'left',  [2; 0.15;0.1;0;0;0]); 
start = [0;0.15;0;0;0;0];
w_goal = [1;1;1;0;0;0];
w_rel = 0.1 * [1;1;1;0;0;0];

[A1, b1] = poly2lincon([-1,-1,1,1], [-1,1,1,-1]);
[A2, b2] = poly2lincon([1.3,1.3, 3, 3], [-1,1,1,-1]);
safe_regions(1) = struct('A', A1, 'b', b1, 'pt', [0;0;0], 'normal', [0;0;1]);
safe_regions(2) = struct('A', A2, 'b', b2, 'pt', [2;0;0.1], 'normal', [0;0;1]);
nr = length(safe_regions);

nsteps = 10;
nx = 4 * nsteps;
ns = nsteps * nr;
nvar = nx + ns;

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

for j = 1:length(r_ndx)
  step_ndx = r_ndx(j);
  if step_ndx == 1
    continue
  end
  Ai = zeros(size(A_reach, 1), nvar);
  Ai(:,x_ndx(:,step_ndx)) = A_reach * diag([1,-1,1,-1]);
  Ai(:,x_ndx(:,step_ndx-1)) = -A_reach * diag([1,-1,1,-1]);
  bi = b_reach;
  A = [A; Ai];
  b = [b; bi];
end

for j = 1:length(l_ndx)
  step_ndx = l_ndx(j);
  if step_ndx == 1
    continue
  end
  Ai = zeros(size(A_reach, 1), nvar);
  Ai(:,x_ndx(:,step_ndx)) = A_reach;
  Ai(:,x_ndx(:,step_ndx-1)) = -A_reach;
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

w_rel = diag(w_rel([1,2,3,6])) * diag([1/nom_forward_step, 1, 1, 1]);
% TODO: make this relative to nominal step width
for j = 2:nsteps
  Q(x_ndx(:,j), x_ndx(:,j)) = Q(x_ndx(:,j), x_ndx(:,j)) + w_rel * w_rel';
  Q(x_ndx(:,j-1), x_ndx(:,j)) = Q(x_ndx(:,j-1), x_ndx(:,j)) -w_rel * w_rel';
  Q(x_ndx(:,j), x_ndx(:,j-1)) = Q(x_ndx(:,j), x_ndx(:,j-1)) -w_rel * w_rel';
  Q(x_ndx(:,j-1), x_ndx(:,j-1)) = Q(x_ndx(:,j-1), x_ndx(:,j-1)) + w_rel * w_rel';
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
    Ai = zeros(size(safe_regions(r).A, 1), nvar);
    Ai(:,x_ndx(1:2,j)) = safe_regions(r).A;
    Ai(:,s_ndx(r,j)) = M;
    bi = safe_regions(r).b + M;
    A = [A; Ai];
    b = [b; bi];
  end
end

lb(x_ndx(:,1)) = start([1,2,3,6]);
ub(x_ndx(:,1)) = start([1,2,3,6]);

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

result = gurobi(model, params)
xstar = result.x;
steps = xstar(x_ndx);
S = xstar(s_ndx)
figure(1);
clf
quiver(steps(1,:), steps(2,:), ones(1,size(steps, 2)), zeros(1,size(steps, 2)));
  