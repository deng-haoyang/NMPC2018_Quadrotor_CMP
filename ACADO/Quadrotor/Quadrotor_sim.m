%% PARAMETERS SIMULATION
X0   = [1 0 1 0 1 0 0 0 0];
Xref = [0 0 0 0 0 0 0 0 0];

input.x  = repmat(Xref,N+1,1);
Xref     = repmat(Xref,N,1);
input.od = [];

Uref     = [g 0 0 0];
Uref     = repmat(Uref,N,1);
input.u  = Uref;

input.y  = [Xref(1:N,:) Uref];
input.yN = Xref(N,:);

input.W  = diag([10, 1, 2, 1, 10, 1, 1, 1, 1, 1, 1, 1, 1]);
input.WN = diag([10, 1, 2, 1, 10, 1, 1, 1, 1]);

%% SIMULATION LOOP
% display('------------------------------------------------------------------')
% display('               Simulation Loop'                                    )
% display('------------------------------------------------------------------')
tol        = 5e-3;
maxIterNum = 10;
iter       = 0; 
time       = 0;
KKT_MPC    = [];
INFO_MPC   = [];
state_sim  = X0;
iterNum    = [];
calTime    = [];
timeElapsed  = 0;
controls_MPC = [];

while time(end) < simLength
    % Solve NMPC OCP
    input.x0 = state_sim(end,:);
    iterNumth = 0;
    timeElapsed = 0;
    if  time(end)>=10 && time(end)<=10.04
        input.y(:,[1,3,5])  = (time(end) - 10)/simTs*0.5;
        input.yN(:,[1,3,5]) = (time(end) - 10)/simTs*0.5;
    end
    while iterNumth<maxIterNum
        output = acado_MPCstep(input);
        input.x = output.x;
        input.u = output.u;
        iterNumth = iterNumth + 1;
        timeElapsed = timeElapsed + output.info.cpuTime;
        if output.info.kktValue < tol
            break;
        end
    end
    % Save the MPC step
    KKT_MPC  = [KKT_MPC; output.info.kktValue];
    controls_MPC = [controls_MPC; output.u(1,:)];
    iterNum = [iterNum;iterNumth];
    calTime = [calTime;timeElapsed];
    % Simulate system
    sim_input.x = state_sim(end,:).';
    sim_input.u = output.u(1,:).';
    states = integrate_sys(sim_input);
    state_sim = [state_sim; states.value'];
    iter = iter+1;
    nextTime = iter*simTs; 
    time = [time;nextTime];
end
rec.u = controls_MPC;
rec.x = state_sim;
rec.iter = iterNum;
rec.calTime = calTime;
rec.e = KKT_MPC;
rec.simTime = time;
disp(['All RTI Time: ' num2str(sum(rec.calTime)) 's']);