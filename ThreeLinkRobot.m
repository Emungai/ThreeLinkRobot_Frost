classdef ThreeLinkRobot < RobotLinks

    properties
        ContactPoints
        OtherPoints
        DisplayPoints;

    end
    
    methods
        
        function obj = ThreeLinkRobot(urdf)
            
            
            base=get_base_dofs('planar'); %gives 'ground'/'base' joint coordinates; defines it as a floating base model
            %setting the DOF limits; since it's a planar model, it has only 3 DOFs
            limits=[base.Limit];
            [limits.lower]=deal(-5,0,-pi/2); %setting the angle lower limit; angle is defined relative to the initial frame('the straight up coordinate system')
            %defines how much the joint can rotate
            [limits.upper]=deal(5,2,pi/2);%setting the angle upper limit
            [limits.velocity]=deal(15,10,15);%radians/sec
            [limits.effort]=deal(15,10,5); %torque limit, how much force do I need to move the joint
            
            obj = obj@RobotLinks(urdf,base);
            
            
            %% define contact frames
            swingjoint_frame=obj.Joints(getJointIndices(obj,'swing_joint')); %defining the contact frame as the frame attached to the swing foot
                                                                                %the contact point
                                                                                %is the point
                                                                                %where the robot
                                                                                %comes in contact
                                                                                %with the ground;
                                                                                %the stance leg
                                                                                %base joint isn't
                                                                                %a contact point
                                                                                %since we defined
                                                                                %it as floating
            obj.ContactPoints.SwingFoot= CoordinateFrame(...
                'Name','SwingFootEnd',...
                'Reference',swingjoint_frame,...
                'Offset',[0,0,-1],... %position of contact point in reference frame
                'R',[0,0,0]... %z-axis is the normal axis, so there's no need for rotation
                );
            
            stancejoint_frame=obj.Joints(getJointIndices(obj,'stance_joint'));
            obj.ContactPoints.StanceFoot= CoordinateFrame(...
                'Name','StanceFootEnd',...
                'Reference',stancejoint_frame,...
                'Offset',[0,0,-1],... %position of contact point in reference frame
                'R',[0,0,0]... %z-axis is the normal axis, so there's no need for rotation
                );
            
            
            
            %% define other frames
            torso_frame = obj.Joints(getJointIndices(obj,'BaseRotY'));
            obj.OtherPoints.Torso = CoordinateFrame(...
                'Name','TorsoEnd',...
                'Reference',torso_frame,...
                'Offset',[0,0,0.8],...
                'R',[0,0,0]...
                );
            
            %% Define Display Points
            obj.DisplayPoints.SwingFoot = obj.ContactPoints.SwingFoot;
            obj.DisplayPoints.StanceFoot = obj.ContactPoints.StanceFoot;
            obj.DisplayPoints.Torso = obj.OtherPoints.Torso;
            
        end
        
        function ExportKinematics(obj, export_path)
            % Generates code for forward kinematics 
            
            if ~exist(export_path,'dir')
                mkdir(export_path);
                addpath(export_path);
            end
            
            % Compute positions of all joints
            for i = 1:length(obj.Joints)
                position = obj.Joints(i).computeCartesianPosition;
                vars = obj.States.x;
                filename = [export_path, 'p_', obj.Joints(i).Name];
                export(position, 'Vars', vars, 'File', filename);
            end
            
            % Compute positions of contact points
            cp_fields = fields(obj.ContactPoints);
            for i = 1:length(cp_fields)
                position = obj.ContactPoints.(cp_fields{i}).computeCartesianPosition;
                vars = obj.States.x;
                filename = [export_path, 'p_', obj.ContactPoints.(cp_fields{i}).Name];
                export(position, 'Vars', vars, 'File', filename);
            end
            
            % Compute positions of other points
            op_fields = fields(obj.OtherPoints);
            for i = 1:length(op_fields)
                position = obj.OtherPoints.(op_fields{i}).computeCartesianPosition;
                vars = obj.States.x;
                filename = [export_path, 'p_', obj.OtherPoints.(op_fields{i}).Name];
                export(position, 'Vars', vars, 'File', filename);
            end
            
            
        end
    end
    
    
end