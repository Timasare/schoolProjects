% Computational Methods
% November 2022
% Project 5 - Dynamics
% Name:
% Student number:

clc;clear;close all;


% ----------------------------------------------------------------------------------------------------------%
%                                       PRE-PROCESSING
% ----------------------------------------------------------------------------------------------------------%

% Read input file Eiffel.m
Eiffel;
nDofsPerNode = 2;

% ----------------------------------------------------------------------------------------------------------%
% You do not need to change anything in the function P5_SetUpFE
%
% Use the variables (detailed below) for the computation in the remainder
% of script.
%
[nNeumannBCs,nDirichletBCs,nNodes,nElements,elementDofs,x1,x2,globalStiffnessMatrix,globalMassMatrix,globalForceVector,globalDisplacementVector]  ...
    = P5_SetUpFE_BAR(nDofsPerNode, nodalPositions, connectivities,NeumannBCs, DirichletBCs);

% definition of the variables:
% nNeumannBCs              - number of external loads applied
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

Ee = mprop(1);
Ae = mprop(2);
mass = 1;

for e = 1:nElements
    
    asmNodes = connectivities(e,1:end);  % connectivity in terms of global nodes
    asmDofs  = elementDofs(e,:);         % connectivity in terms of global dofs
    

    x1e = x1(asmNodes);
    x2e = x2(asmNodes);

    
    ElementStiffnessMatrix = P5_ComputeStiffness_BAR(x1e,x2e,Ee,Ae);
    
    % ToDo : Adapt the script to work with the consistent mass matrix. 
    %ElementMassMatrix      = P5_ComputeConsistentMass_BAR(mass);
    ElementMassMatrix      = P5_ComputeLumpedMass_BAR(mass);

    globalStiffnessMatrix(asmDofs, asmDofs) = globalStiffnessMatrix( asmDofs,asmDofs)...
                                                       + ElementStiffnessMatrix;
                                                   
    globalMassMatrix(asmDofs, asmDofs) = globalMassMatrix( asmDofs,asmDofs)...
                                                       + ElementMassMatrix;
end

globalStiffnessMatrix;
globalMassMatrix;

[eigenVectors, eigenValues] = eigs(globalStiffnessMatrix, globalMassMatrix, 9, 'smallestabs');
eigenValues=real(eigenValues);
eigenVectors=real(eigenVectors);


% ----------------------------------------------------------------------------------------------------------%
%                                        POST-PROCESSING
% ----------------------------------------------------------------------------------------------------------%


%% Postprocessing and plotting  

x1 = nodalPositions(:,1);
x2 = nodalPositions(:,2);

figure;
hold on
for eigenMode = 1:9
   
    subplot(3,3,eigenMode)
    
    axis equal
    hold on
    x_displaced = zeros(nDofsPerNode*nNodes,1);
    y_displaced=zeros(nDofsPerNode*nNodes,1);
    scaleFactor = 3;
 
 globalDisplacementVector = eigenVectors(:,eigenMode);
 x = nodalPositions(:,1)';
 y = nodalPositions(:,2)';
 
for i = 1:nNodes
    x_displaced(i) = x(i) + scaleFactor* globalDisplacementVector(2*i-1);
    y_displaced(i) = y(i) + scaleFactor* globalDisplacementVector(2*i);
end

   for e = 1:nElements
     asmNodes = connectivities(e,1:2);
     x1e = x_displaced(asmNodes);
     x2e = y_displaced(asmNodes);  
     
     x1eoriginal = x1(asmNodes);
     x2eoriginal = x2(asmNodes);  
     
     xtemp = [x1e(1) x1e(2)];
     ytemp = [x2e(1) x2e(2)];
     xtemporiginal = [x1eoriginal(1) x1eoriginal(2)];
     ytemporiginal = [x2eoriginal(1) x2eoriginal(2)];
     xtemp = [x1e(1) x1e(2)];
     ytemp = [x2e(1) x2e(2)];
     
     p1 = plot(xtemp,ytemp,'r');
     p2 = plot(xtemporiginal,ytemporiginal,'k');
     
   end
   title(sprintf('%1.0f-th eigenmode \n \\omega_%1.0f = %1.2e',...
            eigenMode,eigenMode,sqrt(eigenValues(eigenMode,eigenMode))))
   grid on
   


end
