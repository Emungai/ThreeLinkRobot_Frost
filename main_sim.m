clear; close all;
addpath('C:\Users\mungam\Documents\GitHub\frost-dev')
frost_addpath();
%% defining needed variables
maxangle=pi/8;

urdf = '3linkedrobot.urdf'; %giving the path to the urdf file

model = ThreeLinkRobot(urdf) 

model.ExportKinematics('gen/kinematics/')

display({model.Joints(:).Name}.', 'Joints are'); %displays things
display({model.Links(:).Name}.', 'Links are'); %displays things
    

%joint frame is defined by the child link

%% defining the contact frame in coordinate frame terms
% swingjoint_frame=model.Joints(getJointIndices(model,'swing_joint')); %defining the contact frame as the frame attached to the swing foot
%                                                          %the contact point
%                                                          %is the point
%                                                          %where the robot
%                                                          %comes in contact
%                                                          %with the ground;
%                                                          %the stance leg
%                                                          %base joint isn't
%                                                          %a contact point
%                                                          %since we defined
%                                                          %it as floating
% SwingFoot= CoordinateFrame(...
%     'Name','SwingFootEnd',...
%     'Reference',swingjoint_frame,...
%     'Offset',[0,0,-1],... %position of contact point in reference frame
%     'R',[0,0,0]... %z-axis is the normal axis, so there's no need for rotation
%     );
% 
% stancejoint_frame=model.Joints(getJointIndices(model,'stance_joint'));
% StanceFoot= CoordinateFrame(...
%     'Name','StanceFootEnd',...
%     'Reference',stancejoint_frame,...
%     'Offset',[0,0,-1],... %position of contact point in reference frame
%     'R',[0,0,0]... %z-axis is the normal axis, so there's no need for rotation
%     );   
%                                                          
%% defining other frames

% torso_frame=model.Joints(getJointIndices(model,'BaseRotY'));
% Torso=CoordinateFrame(...
%     'Name','TorsoEnd',...
%     'Reference',torso_frame,...
%     'Offset',[0,0,0],...
%     'R',[0,0,0]...
%     );

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
%  fric_coef.mu=0.9;
%  domain=addContact(domain,stancefoot_sole,fric_coef);%adding contact to the domain; contact constraint; object is floating but domain has a contact
 

stance_foot_position = domain.getCartesianPosition(stancefoot_sole); stance_foot_position = stance_foot_position([1,3]);
h = HolonomicConstraint(domain, stance_foot_position, 'stanceFootPosition');
domain = addHolonomicConstraint(domain, h);

 
 % phase variable: theta1
 theta_min = -pi/8;
 theta_max = pi/8;
 theta = q('stance_joint') + q('BaseRotY');
 p = SymVariable('p',[2,1]); %max and min of th1
 tau = (theta - p(1))/(p(2)-p(1)); % theta1
                        %refer to Atlas code it's phasebased
                        %i.e.LeftStance.m
                        %goes between 0 & 1
                        
 % add event
 eventFunc =  1 - (theta - theta_min) / (theta_max - theta_min); %the negative determines the direction that we are looking at; needs to be decreasing when it hits zero
 condition = UnilateralConstraint(domain,eventFunc,'tauIsOne',{'x'}); %defining the condition for switching the foot
 domain=addEvent(domain,condition);
 
 
 %y = [q('joint2')+q('BaseRotY');pi-2*q('BaseRotY')-q('joint2')-q('joint3')]; %th1=q1
                                                                           %th3=q1+q3
                                                                           %th2=pi-q1-q2-q3
                                                                           %q1=BaseRotY; q3=joint2;q2=joint3
y=[q('BaseRotY'); q('swing_joint') + 2*q('BaseRotY')+q('stance_joint')];
%y=[q('BaseRotY');q('stance_joint')+2*q('BaseRotY')+q('swing_joint')+2*pi];
    
y2 = VirtualConstraint(domain,y,'VConstraint','DesiredType','Bezier','PolyDegree',5,'RelativeDegree',2,'PhaseType','StateBased','PhaseVariable',tau,'PhaseParams',p,'Holonomic',true); %relative degree: the amount of time you differentiate q to get input to model(in this case torque); in most cases this will be 2

 domain=addVirtualConstraint(domain,y2);
 
%% Impact Code
guard = RigidImpact('SwingLegImpact',domain,'tauIsOne');
guard.R = guard.R(:,[1,2,3,5,4]); %swapping th1 with th2 & vice versa
guard.addImpactConstraint(struct2array(domain.HolonomicConstraints));

%% Compile 
export_path = 'gen/';
domain.compile(export_path);

%% Defining Hybrid System
addpath(genpath(export_path))

system=HybridSystem('System');
io_control  = IOFeedback('IO');
param = struct;
param.aVConstraint = zeros(2,6); %linspace(theta_max, theta_min, 6)]; %virtual constraints that depend on tau; makes the bezier polynomial a constant
param.pVConstraint = [theta_min, theta_max]; %defines what the max and min of th1 are
param.kVConstraint = [100,20];

system = addVertex(system,'RightStance','Domain',domain,'Control', io_control);
system = setVertexProperties(system,'RightStance','Param',param);

srcs = {'RightStance'};
tars = {'RightStance'};
system = addEdge(system, srcs, tars);
system = setEdgeProperties(system, srcs, tars, ...
    'Guard', {guard});
x0=[0;0;0;-pi/8;pi/8; cos(pi/8);sin(pi/8);0;1;-4];
%x0=[0;0;0;-pi/8;pi/8; zeros(5, 1)];
% logger = system.simulate(0, x0, 1, [],'NumCycle',1);

%% Create initial condition: using fminsearch to minimize the velocity of the stance foot
velocity_stance_foot = jacobian(stance_foot_position, q)*dq; 
velocity_stance_foot.export('Vars', {q,dq}, 'File', 'gen/velocity_stance_foot_function')

q0=[0;0;pi/16;-pi/8;pi/8];
dq0_first = [cos(pi/8);sin(pi/8);0];

fun = @(dq0_last) norm(velocity_stance_foot_function(q0,[dq0_first; dq0_last]));
dq0_last = fminsearch(fun, [0;0]); %minimizing velocity of the stance foot given the initial condition [0;0]

dq0 = [dq0_first;dq0_last];
x0 = [q0; dq0];

%% Solve

logger = system.simulate(0, x0, 3, [],'NumCycle',2);


%% Create display
fig = figure;
disp = frost.Animator.Display(fig, model);
ax = fig.Children(1);
view(ax, 0, 0);
axis(ax, 'manual');
axis(ax, [-1, 2, -1, 2]);

%% Animate
for i = 1:length(logger)
    for  j = 1:length(logger(i).flow.t)
        qlog = logger(i).flow.states.x(:,j);

        disp.update(qlog);

        pause(0.2);
    end
end


%%
pos =[];
for i = 1:length(logger.flow.t)
    qlog = logger.flow.states.x(:,i);
    pos(:,i) = p_StanceFootEnd(qlog);
    
end

plot(logger.flow.t, pos)

