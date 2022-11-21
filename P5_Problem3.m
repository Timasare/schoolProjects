% Computational Methods
% November 2022
% Project 5 - Dynamics
% Name:
% Student number:

clc;clear;close all


% ----------------------------------------------------------------------------------------------------------%
%                                       PRE-PROCESSING
% ----------------------------------------------------------------------------------------------------------%

% Problem 3:
force = -700;
connectivities =dlmread('connectivitiesb.txt');
nodalPositions = 0.1*dlmread('nodalPositionsb.txt');

NeumannBCs = [403 2 force];

% ToDo : Apply fixed boundary conditions on nodes 1,2,3,4,and 5. also
% constrain the movement of nodes 118 and 121 in y-direction.
DirichletBCs = [1 1 0;
                2 1 0;
                3 1 0;
                4 1 0;
                5 1 0;
                1 2 0;
                2 2 0;
                3 2 0;
                4 2 0;
                5 2 0;
                118 2 0;
                121 2 0];

deg = 403; 
% node 403 is the node on which we will apply the external force afterwards. 
% Think about which DoF it will be ;)

nDofsPerNode = 2;

mprop = [ 12e9 0.3 0.01 1.0];
Ey = mprop(1);
v  = mprop(2);
t  = mprop(3);
mass = mprop(4);

% plane strain
D = (Ey/((1+v)*(1-2*v)))*[1-v v 0
    v 1-v 0
    0 0 (1-2*v)/2];
% ----------------------------------------------------------------------------------------------------------%
% You do not need to change anything in the function P5_SetUpFE
%
% Use the variables (detailed below) for the computation in the remainder
% of script.
%
[nNeumannBCs,nDirichletBCs,nNodes,nElements,elementDofs,x1,x2,globalStiffnessMatrix,globalMassMatrix,globalForceVector,globalDisplacementVector]  ...
    = P5_SetUpFE_CST(nDofsPerNode, nodalPositions, connectivities,NeumannBCs, DirichletBCs);

% definition of the variables:
% nNeumannBCs              - number of external NeumannBCs applied
% nDirichletBCs - number of Dirichlet boundary conditions
% nNodes              - number of global nodes
% nElements           - number of elements
% elementDofs         - local-to-global map (connectivities) for the dofs
% x1, x2              - the x1 and x2 coordinates of the element nodes
% global -Force, -Stiffness, -Displacement matrices of the correct size,
% filled with zeros.

% ----------------------------------------------------------------------------------------------------------%
%                                                SOLVING
% ----------------------------------------------------------------------------------------------------------%

for e = 1:nElements
    
    asmNodes = connectivities(e,1:3);  % connectivity in terms of global nodes
    asmDofs  = elementDofs(e,:);       % connectivity in terms of global dofs
    mID      = connectivities(e,3);    % element material id
    
    
    x1e = x1(asmNodes);
    x2e = x2(asmNodes);
    
    ElementStiffnessMatrix = P5_ComputeStiffness_CST(x1e,x2e, D, t);
    ElementMassMatrix      = P5_ComputeLumpedMass_CST(mass);
    
    globalStiffnessMatrix(asmDofs, asmDofs) = globalStiffnessMatrix( asmDofs,asmDofs)...
        + ElementStiffnessMatrix;
    
    globalMassMatrix(asmDofs, asmDofs) = globalMassMatrix(asmDofs,asmDofs)...
        + ElementMassMatrix;
    
end

% ----------------------------------------------------------------------------------------------------------%
%                                        POST-PROCESSING
% ----------------------------------------------------------------------------------------------------------%

%Eigenvalue problem: critical time step:
[Phi,lambda] = eig(full(globalStiffnessMatrix),full(globalMassMatrix));
wmax   = sqrt(max(max(lambda))); %calculate wmax
dtcrit = 2/wmax; % critical timestep for unconstrained problem
%% 

% solve for displacement with time:
% ToDo : fill in the values. Feel free to play around with dt, this is just
% an initial guess. Note that your problem can become unstable!
omega = 10 ;
dt =0.00001;
tend = 3.5;
t = 0:dt:tend;
tstop = 1;


% ToDo: Adapt alpha and beta for part c).
alpha = 0.0001;
beta = 0.1;

    
[K, M, U ] =  ...
    P5_ApplyBCs(globalForceVector,globalDisplacementVector,globalStiffnessMatrix,globalMassMatrix,nDirichletBCs,DirichletBCs,nDofsPerNode);


% ToDo : Initialize C
C = alpha*M +beta*K ;
 
% initialise all U's

Utip = zeros(length(t),1);
Ualpha = U;
UalphaPlus1 = U;
UalphaMinus1 = U;
Fext = zeros(nDofsPerNode*nNodes,1);
tic
for i = 2:length(t)
    % ToDo: Apply the force to the correct DoF!
    dof = 2*deg ; 
        
% ToDo: set Fext(dof) if ... else....
if t(i) < tstop
    Fext(dof) = force*cos(omega*t(i)) ;% calculating Fext at any point in time t(i)
else
Fext(dof) = 0 ;
    
end

% ToDo : Write down the evolution equation for UalphaPlus1
coA  = M + 0.5*dt*C; % create a variable for coefficient of Ualpha+1
coB  = Fext*(dt^2)+(2*M-(dt^2)*K)*Ualpha+(0.5*dt*C-M)*UalphaMinus1; %create a variable for causal elements 
UalphaPlus1 = coA\coB ;
    
% ToDo : Don't forget to update your UalphaMinus1 and Ualpha ! (you
% initialized them before )
UalphaMinus1=Ualpha; %update UalphaMinus1 with current Ualpha value
Ualpha      =UalphaPlus1;%update Ualpha with current UalphaPlus1 value

Utip(i) = UalphaPlus1(deg*2);

end
toc

% Plotting

plot(t,Utip)
grid on
figure;
hold on
pbaspect([20 1 1])
triplot(connectivities(:,1:3),x1,x2,'LineWidth',1)

