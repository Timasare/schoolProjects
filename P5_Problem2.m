% Computational Methods
% November 2022
% Project 5 - Dynamics
% Name:
% Student number:


clc;clear;close all


% ----------------------------------------------------------------------------------------------------------%
%                                       PRE-PROCESSING
% ----------------------------------------------------------------------------------------------------------%
% Input Problem 2:

% h is a parameter for mesh refinement. The lower the parameter, the finer the mesh.
% TODO : adapt h when required in problem c). 
h = 0.5;   
[nodalPositions,connectivities,DirichletBCs,NeumannBCs,mprop] = inputBeam(h); 

nDofsPerNode = 2;
v = 0.3; % Poisson's ratio


% ----------------------------------------------------------------------------------------------------------%
%
% Use the variables (detailed below) for the computation in the remainder
% of script.

% ToDo : Go into P5_SetUpFE_CST and initialize the globalMassMatrix.
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
%% Solving
Ee = mprop(1);
Ae = mprop(3)^2;
m_total = 635; 
t  = 0.01;
% Element mass, used to compute the mass matrices
mass = m_total / nElements;
%z     = Ee/(1-v^2);
%D     = z*[1 v 0; v 1 0; 0 0 (1-v)/2];

for e = 1:nElements
    
    asmNodes = connectivities(e,1:end);  % connectivity in terms of global nodes
    asmDofs  = elementDofs(e,:);         % connectivity in terms of global dofs
    
    x1e = x1(asmNodes);
    x2e = x2(asmNodes);
    
    ElementStiffnessMatrix = P5_ComputeStiffness_CST(x1e,x2e,Ee,Ae);
    
    % TODO : Write and call the ComputeConsistentMass_CST and
    % ComputeLumpedMass_CST files. 
    %ElementMassMatrix       = P5_ComputeConsistentMass_CST(mass);  
    ElementMassMatrix       = P5_ComputeLumpedMass_CST(mass);
    
    globalStiffnessMatrix(asmDofs, asmDofs) = globalStiffnessMatrix(asmDofs,asmDofs)...
        + ElementStiffnessMatrix;
      % TODO : assemble the globalMassMatrix
    globalMassMatrix(asmDofs, asmDofs) = globalMassMatrix(asmDofs,asmDofs)...
        + ElementMassMatrix;
    
  
  
        
end

% TODO : Compute the eigenvalues and eigenvectors. 
[eigenVectors, eigenValues] = eigs(globalStiffnessMatrix, globalMassMatrix, 9, 'smallestabs');
eigenValues=real(eigenValues);
eigenVectors=real(eigenVectors);
% [eigenVectors, eigenValues] = ... ;
% eigenValues=real(eigenValues);
% eigenVectors=real(eigenVectors);

% ----------------------------------------------------------------------------------------------------------%
%                                        POST-PROCESSING
% ----------------------------------------------------------------------------------------------------------%


%% Postprocessing and plotting
figure;
pbaspect([20 1 1])
hold on
for eigenMode = 1:9
    
    subplot(3,3,eigenMode)
    
    axis equal
    hold on
    x_displaced = zeros(nDofsPerNode*nNodes,1);
    y_displaced=zeros(nDofsPerNode*nNodes,1);
    scaleFactor = 10;
    
    
    globalDisplacementVector = eigenVectors(:,eigenMode);
    x = nodalPositions(:,1)';
    y = nodalPositions(:,2)';
    
    for i = 1:nNodes
        x_displaced(i) = x(i) + scaleFactor* globalDisplacementVector(2*i-1);
        y_displaced(i) = y(i) + scaleFactor* globalDisplacementVector(2*i);
    end
    triplot(connectivities(:,1:end),x_displaced,y_displaced,'r','LineWidth',1)
    
    title(sprintf('%1.0f-th eigenmode \n \\omega_%1.0f = %1.2e',...
            eigenMode,eigenMode,sqrt(eigenValues(eigenMode,eigenMode))))
    
end