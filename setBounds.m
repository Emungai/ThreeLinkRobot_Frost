% %% Set bounds for optimization problem
% 
% % Set Bounds
% model_bounds = model.getLimits();
% model_bounds.states.x.lb = [-10, -1, -pi/2, -pi/2, -pi/2];
% model_bounds.states.x.ub = [10, 2, pi/2, pi/2, pi/2];
% model_bounds.states.dx.lb = -20*ones(1,5);
% model_bounds.states.dx.ub = 20*ones(1,5);
% model_bounds.states.ddx.lb = -100*ones(1,5);
% model_bounds.states.ddx.ub = 100*ones(1,5);
% bounds = struct();
% 
% % Right Stance
% bounds.RightStance = model_bounds;
% 
% bounds.RightStance.time.t0.lb = 0;
% bounds.RightStance.time.t0.ub = 0;
% bounds.RightStance.time.t0.x0 = 0;
% 
% bounds.RightStance.time.tf.lb = 0.25;
% bounds.RightStance.time.tf.ub = 0.75;
% bounds.RightStance.time.tf.x0 = 0.75;
% 
% bounds.RightStance.time.duration.lb = 0.25;
% bounds.RightStance.time.duration.ub = 0.75;
% bounds.RightStance.time.duration.x0 = 0.75;
% 
% bounds.RightStance.inputs.ConstraintWrench.fRightToe.lb = -1000;
% bounds.RightStance.inputs.ConstraintWrench.fRightToe.ub = 1000;
% bounds.RightStance.inputs.ConstraintWrench.fRightToe.x0 = 100;
% 
% bounds.RightStance.inputs.Control.u.lb = -100*ones(2,1);
% bounds.RightStance.inputs.Control.u.ub = 100*ones(2,1);
% bounds.RightStance.inputs.Control.u.x0 = zeros(2,1);
% 
% bounds.RightStance.params.pRightToe.lb = -0*ones(3,1);
% bounds.RightStance.params.pRightToe.ub = 0*ones(3,1);
% bounds.RightStance.params.pRightToe.x0 = zeros(3,1);
% 
% bounds.RightStance.params.aVconstraint.lb = -10*ones(6*2,1); %6*2: 2 represents the number of joints I have; 6 is related to the bezier polynomial degree(5) (the coefficients)
% bounds.RightStance.params.aVconstraint.ub = 10*ones(6*2,1);
% bounds.RightStance.params.aVconstraint.x0 = zeros(6*2,1);
% 
% bounds.RightStance.params.pVconstraint.lb = [-pi/4; 0];
% bounds.RightStance.params.pVconstraint.ub = [0; pi/4];
% bounds.RightStance.params.pVconstraint.x0 = [-pi/8; pi/8];
% 
% bounds.RightStance.time.kp = 100;
% bounds.RightStance.time.kd = 20;
% 
% % Right Impact
% bounds.SwingLegImpact = model_bounds;


%% Set bounds for optimization problem

% Set Bounds 
model_bounds = model.getLimits();
model_bounds.states.x.lb = [ -10, 0.75, -pi/8, 1, 1];
model_bounds.states.x.ub = [10, 10, pi/8, 5, 5];
model_bounds.states.dx.lb = -100*ones(1,5);
model_bounds.states.dx.ub = 100*ones(1,5);
bounds = struct();

% Right Stance
bounds.RightStance = model_bounds;
bounds.RightStance.time.t0.lb = 0;
bounds.RightStance.time.t0.ub = 0;
bounds.RightStance.time.t0.x0 = 0;

bounds.RightStance.time.tf.lb = 0.6;
bounds.RightStance.time.tf.ub = 1;
bounds.RightStance.time.tf.x0 = 1;

bounds.RightStance.time.duration.lb = 0.6;
bounds.RightStance.time.duration.ub = 1;
bounds.RightStance.time.duration.x0 = 1;

bounds.RightStance.inputs.ConstraintWrench.fRightToe.lb = -10000;
bounds.RightStance.inputs.ConstraintWrench.fRightToe.ub = 10000;
bounds.RightStance.inputs.ConstraintWrench.fRightToe.x0 = 100;

bounds.RightStance.inputs.Control.u.lb = -100*ones(2,1);
bounds.RightStance.inputs.Control.u.ub = 100*ones(2,1);
bounds.RightStance.inputs.Control.u.x0 = zeros(2,1);

bounds.RightStance.params.pRightToe.lb = -0*ones(3,1);
bounds.RightStance.params.pRightToe.ub = 0*ones(3,1);
bounds.RightStance.params.pRightToe.x0 = zeros(3,1);

bounds.RightStance.params.atimevc.lb = -10*ones(6*2,1);
bounds.RightStance.params.atimevc.ub = 10*ones(6*2,1);
bounds.RightStance.params.atimevc.x0 = zeros(6*2,1);

bounds.RightStance.params.pVconstraint.lb = [-pi/4; 0];
bounds.RightStance.params.pVconstraint.ub = [0; pi/4];
bounds.RightStance.params.pVconstraint.x0 = [-pi/8; pi/8];
% 
bounds.RightStance.time.kp = 100;
bounds.RightStance.time.kd = 20;

% Right Impact
bounds.SwingLegImpact = model_bounds;

