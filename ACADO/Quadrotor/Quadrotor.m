clc;
clear all;
% close all;

simTs = 0.01;
simLength = 20;
T = 1.0;
N = 48;
EXPORT = 1;

% 1: qpOASES 2: qpDUNES
solver_selection = 2;

DifferentialState x1 x2 x3 x4 x5 x6 x7 x8 x9;
Control u1 u2 u3 u4;

n_XD = length(diffStates);
n_U = length(controls);

%% Differential Equation
g = 9.81;
f = dot([x1;x2;x3;x4;x5;x6;x7;x8;x9]) == [x2;...
    u1*(cos(x7)*sin(x8)*cos(x9) + sin(x7)*sin(x9));...
    x4;...
    u1*(cos(x7)*sin(x8)*sin(x9) - sin(x7)*cos(x9));...
    x6;...
    u1*cos(x7)*cos(x8) - g;...
    (u2*cos(x7) + u3*sin(x7))/cos(x8);...
    -u2*sin(x7) + u3*cos(x7);...
    u2*cos(x7)*tan(x8) + u3*sin(x7)*tan(x8) + u4];

h = [diffStates; controls];
hN = [diffStates];
%% SIMexport
acadoSet('problemname', 'sim');
sim = acado.SIMexport( simTs );
sim.setModel(f);
sim.set( 'INTEGRATOR_TYPE',             'INT_RK4' );
sim.set( 'NUM_INTEGRATOR_STEPS',         1        );

if EXPORT
    sim.exportCode( 'export_SIM' );
    cd export_SIM
    make_acado_integrator('../integrate_sys');
    cd ..
end
%% MPCexport
acadoSet('problemname', 'mpc');

ocp = acado.OCP( 0.0, T, N );
W_mat = eye(n_XD+n_U,n_XD+n_U);
WN_mat = eye(n_XD,n_XD);
W = acado.BMatrix(W_mat);
WN = acado.BMatrix(WN_mat);
ocp.minimizeLSQ( W, h );
ocp.minimizeLSQEndTerm( WN, hN );

ocp.subjectTo( 0.0  <= u1 <= 11.0 );% Bounds
ocp.subjectTo( -1.0 <= u2 <= 1.0 );
ocp.subjectTo( -1.0 <= u3 <= 1.0 );
ocp.subjectTo( -1.0 <= u4 <= 1.0 );
ocp.setModel(f);

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'PRINTLEVEL', 'NONE');
mpc.set( 'INTEGRATOR_TYPE',             'INT_EX_EULER'      );
mpc.set( 'NUM_INTEGRATOR_STEPS',         N                  );
mpc.set( 'LEVENBERG_MARQUARDT', 		 1e-10				);
switch solver_selection
    case 1
        disp('qpOASES selected!');
        % qpOASES dense solver
        mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);
        mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
        mpc.set( 'HOTSTART_QP',                 'YES'             	);
    case 2
        disp('qpDUNES selected!');
        % qpDUNES sparse solver
        mpc.set( 'SPARSE_QP_SOLUTION', 'SPARSE_SOLVER' );
        mpc.set( 'QP_SOLVER', 'QP_QPDUNES' );
    otherwise
        disp('Please select a solver!');
end
if EXPORT
    mpc.exportCode( 'export_MPC' );
    switch solver_selection
        case 1
            disp('qpOASES exported!');
            copyfile('../../../../../../external_packages/qpoases', 'export_MPC/qpoases', 'f')
        case 2
            disp('QP_QPDUNES exported!');
            copyfile('../../../../../../external_packages/qpdunes', 'export_MPC/qpdunes', 'f')
        otherwise
            disp('Please select a solver!');
    end
    cd export_MPC
    make_acado_solver('../acado_MPCstep')
    cd ..
end


