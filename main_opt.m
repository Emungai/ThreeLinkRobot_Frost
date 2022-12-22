clear; close all;
COMPILE = true;
%% defining needed variables
maxangle=pi/8;

urdf = '3linkedrobot.urdf'; %giving the path to the urdf file

model = ThreeLinkRobot(urdf) 

% model.ExportKinematics('gen/kinematics/')

display({model.Joints(:).Name}.', 'Joints are'); %displays things
display({model.Links(:).Name}.', 'Links are'); %displays things
    

%joint frame is defined by the child link

%% Defining the Dynamics
configureDynamics(model);
% defining the domains(i.e. double support; single support etc)
    %
 domain=copy(model);
 domain.setName('Swing');
 
 q=domain.States.x; %giving you symbolic value of state variable
                    %BasePosX = x
                    %BasePosZ =Z
                    %BaseRotY = joint 1(th1)
                    %Joint3= th2 torso to swing leg
                    %Joint2=th3 stance leg to torso
                    %if you had a 3d robot, you would have 3 pos(BasePos) & 3
                    %rotations(BaseRot)
 dq=domain.States.dx;
 
 stancefoot_sole=ToContactFrame(model.ContactPoints.StanceFoot,'PointContactWithFriction');%in the swing motion, the stance foot is what's contacting the ground
                                                                        %defining
                                                                        %the
                                                                        %type
                                                                        %of
                                                                        %contact
                                                                        %I
                                                                        %have
 fric_coef.mu=0.9;
 domain=addContact(domain,stancefoot_sole,fric_coef);%adding contact to the domain; contact constraint; object is floating but domain has a contact
%  
%add event 
stance_foot_position = domain.getCartesianPosition(stancefoot_sole); stance_foot_position = stance_foot_position([1,3]);
% %manually adding contact
% h = HolonomicConstraint(domain, stance_foot_position, 'stanceFootPosition');
% domain = addHolonomicConstraint(domain, h);

% stance_foot_position = domain.getCartesianPosition(stancefoot_sole);
% h_nsf = UnilateralConstraint(domain,stance_foot_position(3),'leftFootHeight','x');
% domain = addEvent(domain, h_nsf);

 % phase variable: theta1
 theta = q('stance_joint') + q('BaseRotY');
 p = SymVariable('p',[2,1]); %max and min of th1
 tau = (theta - p(2))/(p(1)-p(2)); % theta1
                        %refer to Atlas code it's phasebased
                        %i.e.LeftStance.m
                        %goes between 0 & 1
                        
eventFunc = (1 - tau); %the negative determines the direction that we are looking at; needs to be decreasing when it hits zero
% param = struct;
% param.aVConstraint = zeros(2,6); %linspace(theta_max, theta_min, 6)]; %virtual constraints that depend on tau; makes the bezier polynomial a constant
% param.pVConstraint = p; %defines what the max and min of th1 are
% param.kVConstraint = [100,20];
% param=addParam(domain,p);


%  
 
 %y = [q('joint2')+q('BaseRotY');pi-2*q('BaseRotY')-q('joint2')-q('joint3')]; %th1=q1
                                                                           %th3=q1+q3
                                                                           %th2=pi-q1-q2-q3
                                                                           %q1=BaseRotY; q3=joint2;q2=joint3
y=[q('BaseRotY'); q('swing_joint') + 2*q('BaseRotY')+q('stance_joint')];
%y=[q('BaseRotY');q('stance_joint')+2*q('BaseRotY')+q('swing_joint')+2*pi];
    
y2 = VirtualConstraint(domain,y,'VConstraint','DesiredType','Bezier','PolyDegree',5,'RelativeDegree',2,'PhaseType','StateBased','PhaseVariable',tau,'PhaseParams',p,'Holonomic',true); %relative degree: the amount of time you differentiate q to get input to model(in this case torque); in most cases this will be 2

domain=addVirtualConstraint(domain,y2);

condition = UnilateralConstraint(domain,eventFunc,'tauIsOne',{'x','pVConstraint'}); %defining the condition for switching the foot
domain=addEvent(domain,condition);
 
%% Impact Code
guard = RigidImpact('SwingLegImpact',domain,'tauIsOne');
guard.R = guard.R(:,[1,2,3,5,4]); %swapping th1 with th2 & vice versa

%defining the swing foot as the foot that's 'stuck' to the ground after
%impact
swing_foot_position =  domain.getCartesianPosition(model.ContactPoints.SwingFoot);
h = HolonomicConstraint(domain, swing_foot_position, 'swingFootPosition');
guard.addImpactConstraint(h);


% Example constraint removal
% removeConstraint(nlp.Phase(1),'u_friction_cone_RightToe');

%% Defining Hybrid System
% addpath(genpath(export_path))

system=HybridSystem('System');
% io_control  = IOFeedback('IO');
% param = struct;
% param.aVConstraint = zeros(2,6); %linspace(theta_max, theta_min, 6)]; %virtual constraints that depend on tau; makes the bezier polynomial a constant
% param.pVConstraint = [theta_min, theta_max]; %defines what the max and min of th1 are
% param.kVConstraint = [100,20];

system = addVertex(system,'Swing','Domain',domain);%'Control', io_control) %adding a node
    %'Swing': the names of the vertices to be added
    %'Domain': a character vector specifiy the name of a vertex property for the added vertices.
    %domain: 
%system = setVertexProperties(system,'Swing','Param',param);

srcs = {'Swing'}; %source node
tars = {'Swing'}; %target node
system = addEdge(system, srcs, tars); %sets an edge that starts and ends at the same node
system = setEdgeProperties(system, srcs, tars, ...
    'Guard', {guard}); %sets the name of the vertex or node as guard 
%x0=[0;0;0;-pi/8;pi/8; cos(pi/8);sin(pi/8);0;1;-4];
%x0=[0;0;0;-pi/8;pi/8; zeros(5, 1)];
% logger = system.simulate(0, x0, 1, [],'NumCycle',1);

%% User Costs: 
u=domain.Inputs.Control.u; %saving torque 
u2r=tovector(norm(u).^2); %saving the norm of torque as a mathematica symbolic expression
                          %since optimization expects one equation
u2r_fun=SymFunction(['torque_' domain.Name],u2r,{u});
    %['torque_' domain.Name] = torque_Swing
    %wraps a symbolic expression with additional information
    % for the convenience of compiling and exporting the symbolic
    % expression to a C/C++ source files.
%% Create Optimization Problem
num_grid=10;
%nlp = HybridTrajectoryOptimization('ThreeLink_system',system,num_grid,[],'EqualityConstraintBoundary',1e-4);
nlp = HybridTrajectoryOptimization('ThreeLink_system',system,num_grid,[],'EqualityConstraintBoundary',1e-4);
setBounds;
nlp.configure(bounds);

% Add costs
addRunningCost(nlp.Phase(getPhaseIndex(nlp,'Swing')),u2r_fun,'u');

%Update
nlp.update;

% %% Swing Foot should be 0 at last node
% X = SymVariable('x',[nlp.Phase(1).OptVarTable.x(1).Dimension,1]);
% swingFootHeight = SymFunction('swing_foot_height', nlp.Phase(1).Plant.EventFuncs.leftFootHeight.ConstrExpr, {X});
% addNodeConstraint(nlp.Phase(1), swingFootHeight,{'x'}, 'last', 0, 0, 'Nonlinear');
% 
%  % average velocity
%     velocity_desired =0.3;
%     DOF = 7;
%     T = SymVariable('t',[nlp.Phase(1).OptVarTable.T(1).Dimension,1]);
%     X0 = SymVariable('x0',[nlp.Phase(1).OptVarTable.x(1).Dimension,1]);
%     XF = SymVariable('xF',[nlp.Phase(1).OptVarTable.x(1).Dimension,1]);
%     v_des = SymVariable('vdes');
%     avg_vel = (XF(1) - X0(1))/(T(2)-T(1))-v_des;
%     avg_vel_fun = SymFunction('average_velocity_a', avg_vel, {T,X0,XF},{v_des});
%     
%     avg_vel_cstr = NlpFunction('Name','average_velocity_b',...
%         'Dimension',1,...
%         'lb', 0,...
%         'ub', 0,...
%         'Type','Nonlinear',...
%         'SymFun',avg_vel_fun,...
%         'DepVariables',[nlp.Phase(1).OptVarTable.T(1); nlp.Phase(1).OptVarTable.x(1); nlp.Phase(1).OptVarTable.x(end)],...
%         'AuxData',velocity_desired);
%   
%     addConstraint(nlp.Phase(1), 'average_velocity_c', 'last', avg_vel_cstr);
%% Compile 
export_path = 'gen/';
% domain.compile(export_path);
if COMPILE
    if ~exist([export_path, 'opt\'])
        mkdir([export_path, 'opt\'])
    end
%     model.ExportKinematics([export_path,'kinematics\']);
    compileConstraint(nlp,[],[],[export_path, 'opt\']);
    compileObjective(nlp,[],[],[export_path, 'opt\']);
end

% nlp.Phase(1).removeConstraint('dynamics_equation');
nlp.Phase(1).removeConstraint('u_friction_cone_StanceFootEnd');
% nlp.Phase(1).removeConstraint('u_tauIsOne_Swing');
% % nlp.Phase(1).removeConstraint('u_tauIsOne_Swing');
% % nlp.Phase(1).removeConstraint('u_tauIsOne_Swing');
% % nlp.Phase(1).removeConstraint('u_tauIsOne_Swing');
% 
nlp.Phase(2).removeConstraint('tContDomain'); %time continuously increasing; don't need for periodic
% nlp.Phase(2).removeConstraint('xMinusCont');
% nlp.Phase(2).removeConstraint('xPlusCont');
% nlp.Phase(2).removeConstraint('dxMinusCont');
% nlp.Phase(2).removeConstraint('dxPlusCont');
nlp.Phase(2).removeConstraint('xDiscreteMapSwingLegImpact');

R = guard.R;
x = guard.States.x; %pre-impact
xn = guard.States.xn;
x_diff = R*x-xn;
x_map = SymFunction(['xDiscreteMap' guard.Name],x_diff(3:end),{x,xn});

addNodeConstraint(nlp.Phase(2), x_map, {'x','xn'}, 'first', 0, 0, 'Linear');
% nlp.Phase(2).removeConstraint('dxDiscreteMapSwingLegImpact');
% nlp.Phase(2).removeConstraint('dxDiscreteMapSwingLegImpact');

if COMPILE
    if ~exist([export_path, 'opt\'])
        mkdir([export_path, 'opt\'])
    end
%     model.ExportKinematics([export_path,'kinematics\']);
    compileConstraint(nlp,2,'xDiscreteMapSwingLegImpact',[export_path, 'opt\']);
end
%% Create Ipopt solver
addpath(genpath(export_path));
nlp.update;
solver = IpoptApplication(nlp);

% Run Optimization
tic
% old = load('x0');
% [sol, info] = optimize(solver, old.sol);
[sol, info] = optimize(solver);
toc
[tspan, states, inputs, params] = exportSolution(nlp, sol);

checkConstraints(nlp, sol, 1e-4, 'constraintCheck.txt')
% %% Create initial condition: using fminsearch to minimize the velocity of the stance foot
% velocity_stance_foot = jacobian(stance_foot_position, q)*dq; 
% velocity_stance_foot.export('Vars', {q,dq}, 'File', 'gen/velocity_stance_foot_function')
% 
% q0=[0;0;pi/16;-pi/8;pi/8];
% dq0_first = [cos(pi/8);sin(pi/8);0];
% 
% fun = @(dq0_last) norm(velocity_stance_foot_function(q0,[dq0_first; dq0_last]));
% dq0_last = fminsearch(fun, [0;0]); %minimizing velocity of the stance foot given the initial condition [0;0]
% 
% dq0 = [dq0_first;dq0_last];
% x0 = [q0; dq0];

%% Solve
% logger = system.simulate(0, x0, 3, [],'NumCycle',2);


%% Create display
fig = figure;
disp = frost.Animator.Display(fig, model);
% ax = fig.Children(1);
% view(ax, 0, 0);
% axis(ax, 'manual');
% axis(ax, [-5, 5, -5, 5]);
% hold(ax,'on') 
%% Animate

for  j = 1:length(tspan{1})
    qlog = states{1}.x(:,j)

    disp.update(qlog);

    pause(0.2);
end



% %%
% pos =[];
% for i = 1:length(logger.flow.t)
%     qlog = logger.flow.states.x(:,i);
%     pos(:,i) = p_StanceFootEnd(qlog);
%     
% end
% 
% % plot(logger.flow.t, pos)


