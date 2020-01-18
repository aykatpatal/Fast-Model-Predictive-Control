clear all;clc;warning off all;close all;
%% Configure simulation for the below params
% time
t0 = 0;                     % starting time
fps = 10;                   % sampling rate
tf = 10;                    % final time
dt = 1/fps;                 % time step
tspan = t0:dt:sum(tf)-dt;   % time span

%% initial state
q0 = [0;3;-pi/4];
dq0 = [40;-40;0];
ddq0 = [-4;7;0];

q = q0;         % state register for displacement
dq = dq0;       % state register for velocity
ddq = ddq0;     % state register for accleration

%% controller params and limits to take care of while controlling
N = 10;                 % length of window
Q = [1.2;1.2;1];        % Cost matrix on states
R = [0.1;0.1];          % Cost matrix on the control inputs
u_max = [0.47;3.77];    % bound on control input

% calculation of const params
u_min = -u_max;
Qhat = diag(repmat(Q,[N,1]));
Rhat = diag(repmat(R,[N,1]));
d = [repmat(u_max,[N,1]);repmat(u_min,[N,1])];
D = [eye(2*N);-eye(2*N)];

%% Quintic polynomial
poly_order = 5;
numberofCoeffs = 6;
coeff_x = [1;3;1;18;0;0];
coeff_y = [2;2;0.5;20;0;0];

%% Generate trajectory
[p,pd,qd,dpd,dqd,ddpd,ddqd] = trajGeneration(poly_order,coeff_x(1:numberofCoeffs),coeff_y(1:numberofCoeffs));
nbp1 = @(pn,npd)(npd/2*pn-npd/2+1);
nbp2 = @(pn,npd)(npd/2*pn+npd/2);

%% Iteration params
iter = 0;          % iterations
path_num = 1;   
time = 0;          % time counter

% initialized q for and it is used for warm start 
init_q = zeros(20,1);
warm_start = [0;0;0;0;0;0;0;0;0;0];

%% MAIN LOOP
tic %% to measure the time of execution #STARTING TIME
while (iter+N)*dt<=sum(tf) %run loop until all the iteration completes
    % calculations within the window
    u_ref = [];    % create new controller
    
    for i = iter:iter+N-1
        % local time and regulation
        if path_num == 1
            time = i*dt;
        else
            time = i*dt-sum(tf(1:path_num-1));
        end
        
        % if path is piecewise and it's time to change
        if time > tf(path_num)
            path_num = path_num+1;      % next path
            time = time-tf(path_num-1);       % regulate the time
            % create new pd if pn changes
            [p,pd,qd,dpd,dqd,ddpd,ddqd] = Create2DPd(poly_order,...
                coeff_x(nbp1(path_num,numberofCoeffs):nbp2(path_num,numberofCoeffs)),coeff_y(nbp1(path_num,numberofCoeffs):nbp2(path_num,numberofCoeffs)));
        end
        
        % if it is still not the time, but pn+1 in last window
        if  time < 0
            path_num = path_num-1;      % roll back to last path
            time = time+tf(path_num);         % regulate the time
            
            % create new pd if pn changes
            [p,pd,qd,dpd,dqd,ddpd,ddqd] = Create2DPd(poly_order,...
                coeff_x(nbp1(path_num,numberofCoeffs):nbp2(path_num,numberofCoeffs)),coeff_y(nbp1(path_num,numberofCoeffs):nbp2(path_num,numberofCoeffs)));
        end
        
        % ///////////    update reference qr,dqr,vr
        qr(1:2,i+1) = qd(time,tf(path_num));
        dqr(1:2,i+1) = dqd(time,tf(path_num));
        ddqr(1:2,i+1) = ddqd(time,tf(path_num));
        qr(3,i+1) = atan2(dqr(2,i+1),dqr(1,i+1));
        
        if iter == 0
            dqr(3,i+1) = 0;
        else
            dqr(3,i+1) = (qr(3,i+1)-qr(3,iter))/dt;
        end
        vr(i+1) = norm(dqr(1:2,i+1),2);
        % update ur
        u_ref = [u_ref;vr(i+1);dqr(3,i+1)];
    end
    
    % update boundary requirement d
    d = d-[u_ref;-u_ref];
    
    %% we need different A and B at every point as system is not linear
    % update A and B matrix for the window
    %% Note A and B are stored in stack
    for i = 1:N
        A(:,:,i) = [1, 0, -vr(iter+i)*sin(qr(3,iter+i)*dt);
            0, 1, vr(iter+i)*cos(qr(3,iter+i)*dt)
            0, 0 ,1];
        B(:,:,i) = [cos(qr(3,iter+i))*dt, 0;
            sin(qr(3,iter+i))*dt, 0;
            0, dt];
    end
    
    %% A and B are stored in BIGG matrix for optimization 
    % introduce new Ahat(3N,3), Bhat(3N,2N)
    A_hat = repmat(eye(3,3),N,1);
    B_hat = repmat(zeros(3,2),N,N);
    for i = 1:N %% this is how we make big matrix
        B_hat(i*3-2:end,i*2-1:i*2) = repmat(B(:,:,i),N+1-i,1);
        for j = i:N
            A_hat(j*3-2:j*3,:) = A(:,:,i)*A_hat(j*3-2:j*3,:);
            for m = i+1:j
                B_hat(j*3-2:j*3,i*2-1:i*2) = A(:,:,m)*B_hat(j*3-2:j*3,i*2-1:i*2);
            end
        end
    end
    
    %% HERE COMES "optimization"
    % error in states
    e_hat = q(:,iter+1)-qr(:,iter+1);   
    % Hessian matrix: quadratic cost 
    H = 2*(B_hat'*Qhat*B_hat+Rhat);     
    % linear cost
    f = 2*B_hat'*Qhat*A_hat*e_hat;      
    
    % limits and params for optimization 
    x_min = [-25;-25;-3];
    x_max = [25;25;3];
    T = N;
    n = size(x_min,1);
    m = size(u_min,1);
    x = q0;    

    %% Inequality constraint construction
    P = zeros(2*T*(n+m),T*(n+m));
    h = zeros(2*T*(n+m),1);

    for i=1:2*(m+n):size(P,1)-2*(m+n)+1 % THIS IS HOW WE MAKE BIG INEQUALITY CONSTRAINTS
        if i==1
            P(i:i+2*(m+n)-1,i:i+(m+n)-1) = [eye(m) zeros(m,n);(-eye(m)) zeros(m,n);
                zeros(n,m) eye(n); zeros(n,m) -(eye(n))];
        else
            P(i:i+2*(m+n)-1,(i+1)/2:(i+1)/2+(m+n)-1) = [eye(m) zeros(m,n);(-eye(m)) zeros(m,n);...
                zeros(n,m) eye(n); zeros(n,m) -(eye(n))];
        end
    end

    for i=1:2*(m+n):2*T*(m+n)-2*(m+n)+1
       h(i:2*(m+n)+i-1) = [u_max;-u_min;x_max;-x_min]; 
    end

    h_star=h(1:20,1);    
    [row, col] = size(P);
    temp = [];

    for j = 1:5:50
        temp = [temp P(:,j) P(:,j+1)];    
    end
    
    % size of the p is removed by removing the state terms from P
    P = temp;    
    
    
    %% HERE COMES THE SHOW
    
    keq = 1; % barrier parameter initial state
    mu=1/10; % mu factor which will change the Keq by that factor   

    % Objective function with log barrier
    Fun = @(z) (z'*H*z + f'*z)-keq*(-sum(log(h - P*z)));
    
    OPTIONS = optimset('Display','off');
    % Optimizer: Uses Interior point method
    while (keq*length(init_q) >= 10e-3) % iterate until maximum iteration reached
     q_op = fmincon(Fun,init_q,[],[],[],[],[],[],[],OPTIONS);
     init_q=q_op; 
     keq=mu*keq;
    end 
     % now we have best q_op 
     
    %% WARM START
    %update q op for next Horizon 
    init_q = [init_q(3:20); 0 ; 0]; %zeroes appended at the end otherwise ideal control can also be used
    
    % update state q, include q0 so starts from k+2
    q(1,iter+2) = q(1,iter+1)+(u_ref(1)+q_op(1))*cos(q(3,iter+1))*dt;
    q(2,iter+2) = q(2,iter+1)+(u_ref(1)+q_op(1))*sin(q(3,iter+1))*dt;
    q(3,iter+2) = q(3,iter+1)+(u_ref(2)+q_op(2))*dt;
    dq(:,iter+2) =  (q(:,iter+2)-q(:,iter+1))./dt;
    ddq(:,iter+2) =  (dq(:,iter+2)-dq(:,iter+1))./dt;    
    % update counter
    iter=iter+1;
end
toc  %end of time calculation 

%% SHOW plot
tq=0:dt:iter*dt;% time span for q

%           XY Plot
figure(1)
hold on
plot(tspan,qr(1,:),'g-');
plot(tq,q(1,:),'r-');
title('x path compare');
legend('Desired x','Fast MPC}');
xlabel('t');
ylabel('x');

%           Velocity in x Plot
figure(2)
hold on
plot(tspan,dqr(1,:),'g-');
plot(tq,dq(1,:),'r-');
title('Velocity in x');
legend('Desired velocity in x','Fast MPC}');
xlabel('t');
ylabel('x');

%           Acceleration in x
figure(3)
hold on
plot(tspan,ddqr(1,:),'g-');
plot(tq,ddq(1,:),'r-');
title('Acceleration in x');
legend('Desired acceleration in x','Fast LMPC}');
xlabel('t');
ylabel('x');

%           y trajectory
figure(4)
hold on
plot(tspan,qr(2,:),'g-');
plot(tq,q(2,:),'r-');
title('Trajectory in y');
legend('Desired y','Fast MPC}');
xlabel('t');
ylabel('Y');

%           velocity in y
figure(5)
hold on
plot(tspan,dqr(2,:),'g-');
plot(tq,dq(2,:),'r-');
title('Velocity in y');
legend('Desired velocity in y','Fast MPC');
xlabel('t');
ylabel('Y');

%           Acceleration in Y
figure(6)
hold on
plot(tspan,ddqr(2,:),'g-');
plot(tq,ddq(2,:),'r-');
title('Acceleration in Y');
legend('Desired acceleration in Y','Fast MPC');
xlabel('t');
ylabel('Y');

%           Trajectory for \theta
figure(7)
hold on
plot(tspan,qr(3,:),'g-');
plot(tq,q(3,:),'r-');
title('Trajectory for \theta');
legend('Desired \theta', 'Fast MPC');
xlabel('t');
ylabel('\theta');

%           XY trajectory
figure(8)
hold on
plot(qr(1,:),qr(2,:),'g-');
plot(q(1,:),q(2,:),'r-');
%plot(bpx,bpy,'g-x')
title('XY trajectory');
legend('p_d','p_{Fast MPC}','Control Points');
%legend('p_d','p_{LMPC}');
xlabel('X');
ylabel('Y');

%           error in y
figure(9)
hold on
plot(tq,q(2,:)-qr(2,1:size(tq,2)),'g-');
plot(tq,q(1,:)-qr(1,1:size(tq,2)),'r-');
legend('error in y','error in x')
xlabel('t')

%% HELPER FUNCTION TO GENERATE TRAJCTORY
function [polyFunc,pd,qd,varargout]= trajGeneration(poly_order,bpx,bpy)
    if poly_order == 1
        t_span= @(tf)[1,0;1,tf];
        polyFunc = @(tf,bp)(t_span(tf)\bp)';
        pd = @(t,tf,bp)(polyFunc(tf,bp)*[1; t]);
        dpd = @(t,tf,bp)(polyFunc(tf,bp)*[0; 1]);
        ddpd = @(t,tf,bp)(polyFunc(tf,bp)*[0; 0]);
    end
    if poly_order == 2
        t_span= @(tf)[1,0,0,0;0,1,0,0;1,tf,tf^2,tf^3;0,1,2*tf,3*tf^2];
        polyFunc = @(tf,bp)(t_span(tf)\bp)';
        pd = @(t,tf,bp)(polyFunc(tf,bp)*[1; t ;t^2 ;t^3]);
        dpd = @(t,tf,bp)(polyFunc(tf,bp)*[0; 1 ;2*t ;3*t^2]);
        ddpd = @(t,tf,bp)(polyFunc(tf,bp)*[0; 0; 2 ;6*t]);
    end
    if poly_order == 3
        t_span= @(tf)[1,0,0,0;-3/tf,3/tf,0,0;3/tf^2,-6/tf^2,3/tf^2,0;-1/tf^3,3/tf^3,-3/tf^3,1/tf^3];
        polyFunc = @(tf,bp)(t_span(tf)*bp)';
        pd = @(t,tf,bp)(polyFunc(tf,bp)*[1; t ;t^2 ;t^3]);
        dpd = @(t,tf,bp)(polyFunc(tf,bp)*[0; 1 ;2*t ;3*t^2]);
        ddpd = @(t,tf,bp)(polyFunc(tf,bp)*[0; 0; 2 ;6*t]);
        % calculations concerning control points
        bpx(1) = (bpx(1)+bpx(2))/2;
        bpx(4) = (bpx(3)+bpx(4))/2;
        bpy(1) = (bpy(1)+bpy(2))/2;
        bpy(4) = (bpy(3)+bpy(4))/2;
    end
    if poly_order == 5
        t_span= @(tf)[1,0,0,0,0,0;0,1,0,0,0,0;0,0,2,0,0,0;1,tf,tf^2,tf^3,tf^4,tf^5;...
            0,1,2*tf,3*tf^2,4*tf^3,5*tf^4;0,0,2,6*tf,12*tf^2,20*tf^3];
        polyFunc = @(tf,bp)(t_span(tf)\bp)';
        pd = @(t,tf,bp)(polyFunc(tf,bp)*[1; t ;t^2 ;t^3 ;t^4 ;t^5]);
        dpd = @(t,tf,bp)(polyFunc(tf,bp)*[0; 1 ;2*t ;3*t^2 ;4*t^3 ;5*t^4]);
        ddpd = @(t,tf,bp)(polyFunc(tf,bp)*[0; 0; 2 ;6*t ;12*t^2 ;20*t^3]);
    end
    qd = @(t,tf)[pd(t,tf,bpx); pd(t,tf,bpy)];
    dqd=@(t,tf)[dpd(t,tf,bpx); dpd(t,tf,bpy)];
    ddqd=@(t,tf)[ddpd(t,tf,bpx); ddpd(t,tf,bpy)];
    % output
    if nargout > 3
        varargout{1} = dpd;
        varargout{2} = dqd;
        if nargout > 5
            varargout{3} = ddpd;
            varargout{4} = ddqd;
        end
    end
end