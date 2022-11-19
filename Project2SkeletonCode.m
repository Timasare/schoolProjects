% Computational Method
clc;clear;close all
% default settings 
set(0,'DefaultFigurePosition', [100 500 500 400],...
   'DefaultAxesFontSize',18,'DefaultLegendFontSize',18,...
   'DefaultLineLineWidth', 2);
% ----------------------------------------------------------------------------------------------------------%
%                                       PRE-PROCESSING
% ----------------------------------------------------------------------------------------------------------%

% Input

% Problem 1: Three Bar Truss Example 
% FILL IN matrices/vectors:

%coordinates of the global nodes
% [x1 x2]
nodalPositions = [0.0 0.0; 1 0; 0.5 sqrt(3)/2];                  

% local-to-global-map 
% [node1,node2,material id],
connectivities = [1 2 1; 2 3 1; 1 3 1];         % 1st element
                % last element
                                                 
% element property matrix
% [ E A ] 
E1 = 200e9;
A1 = 0.01;
mprop           = [E1 A1];     	   % first row = material id no. 1
 					         
               
% prescribed loads
% [node no.,dof_x1 -> 1 or  dof_x2 -> 2, force]

Force = 600e3;


NeumannBCs           = [3,2, -Force];               

% boundary conditions
% [node no., dof_x1 -> 1 or  dof_x2 -> 2, disp]
DirichletBCs             = [1, 2, 0; 2, 2, 0 ; 3, 1, 0];

scaleFactor = 1000;  % we are working with small strains - you will need scaling for visualisation 



% Problem 3 use the provided Eiffel input:
%Eiffel; % simply uncomment this 
%(have a look at the provided file for details)

% ----------------------------------------------------------------------------------------------------------%

%SET UP FE problem - all given for you:

nDofsPerNode          = 2;                                                       % number of dofs per Node
nNodes                    = size(nodalPositions,1);                         % number of global nodes
nDofs                       = nNodes * nDofsPerNode;                      % number of global dofs
nElements                = size(connectivities,1);                          % number of elements
nNodesPerElement   = length(connectivities(1,1:2));                % number of nodes per element
nDofsPerElement      = nDofsPerNode*nNodesPerElement;        % number of dofs per element
nNeumannBCs         = size(NeumannBCs,1);                            % number of external loads applied
nDirichletBCs           = size(DirichletBCs,1);                              % number of Dirichlet boundary conditions

% elemental dofs connectivity:
% (local-to-global map) 
elementDofs = zeros(nElements,nDofsPerElement);

for e = 1:nElements
    asmNode1 = connectivities(e,1:nNodesPerElement-1);
    asmNode2 = connectivities(e,nNodesPerElement);
    elementDofs(e,:) = [asmNode1*2-1, asmNode1*2,asmNode2*2-1, asmNode2*2];
end

% extract x1 and x2 for global nodes
x1 = nodalPositions(:,1);
x2 = nodalPositions(:,2);

% let's plot our initial structure:
  hold on
   for e = 1:nElements
     asmNodes = connectivities(e,1:2); % asmNodes = assemble nodes
     x1e = x1(asmNodes);
     x2e = x2(asmNodes);  

     
     xtemp = [x1e(1) x1e(2)];
     ytemp = [x2e(1) x2e(2)];
     
     p1 = plot(xtemp,ytemp,'k');
     
   end

 vxlabels = arrayfun(@(n){sprintf('%d',n)},(1:nNodes)'); %numbering nodes
Hpl = text(x1 +0.02, x2, vxlabels,'FontWeight', 'bold','FontSize',14); %fonts etc on nodes
p3 = plot(x1,x2,'k o', 'LineWidth',1,'MarkerEdgeColor','k',...
  'MarkerFaceColor','k');%markers on nodes

% ----------------------------------------------------------------------------------------------------------%
%                                                SOLVING
% ----------------------------------------------------------------------------------------------------------%

% Initialise global vectors and matrices - TO DO: fill in the sizes

globalStiffnessMatrix           = zeros(nDofs); %SIZE = ?
globalForceVector               = zeros(nDofs,1); %SIZE = ?
globalDisplacementVector   = zeros(nDofs,1); %SIZE = ?


% loop over all elements 
%-> compute elemental stiffness matrix -> assemble global stiffness matrix

for e = 1:nElements
    
    asmNodes = connectivities(e,1:2);   % assemble nodes = connectivity in terms of global nodes
    asmDofs    = elementDofs(e,:);         % assemble dofs  = connectivity in terms of global dofs
    mID           = connectivities(e,3);      % element material id
   % xtemp = [x1e(1) x1e(2)];
    point1 = nodalPositions(asmNodes(1),:);
    point2 = nodalPositions(asmNodes(2),:);
    %ytemp = [x2e(1) x2e(2)];
    l     = norm (point1 - point2);
    theta   = atan2(point2(2)-point1(2), point2(1)-point1(1));
    T     = [cos(theta) sin(theta) 0 0; 0 0 cos(theta) sin(theta)];
    Ke     = (mprop(mID,1)*(mprop(mID,2)/l))*[1 -1;-1 1];
    Ki    = T'*Ke*T;
    globalStiffnessMatrix(asmDofs,asmDofs)= globalStiffnessMatrix(asmDofs,asmDofs)+Ki;
    % TO DO: Compute elemental stiffness matrices and assemble into global system:

                                                   
                                                   
end

% ----------------------------------------------------------------------------------------------------------%

% assemble global force vector: 
    
for loadIndex=1:nNeumannBCs
        dof = 2*(NeumannBCs(loadIndex,1)-1)+ NeumannBCs(loadIndex,2);
        %dof = NeumannBCs(loadIndex,1)*nDofsPerNode;      % TO DO: call position in global force vector to apply external forces
        globalForceVector(dof) = NeumannBCs(loadIndex,3) ;

end

% ----------------------------------------------------------------------------------------------------------%

% make a copy of K to enforce the Dirichlet (essential) BCs 
K = globalStiffnessMatrix;
GlobalForceVector=globalForceVector;

for boundIndex = 1:nDirichletBCs

        % find essential boundary condition dof
        dof = 2*(DirichletBCs(boundIndex,1)-1)+ DirichletBCs(boundIndex,2); % FILL IN
        
        % enforce essential boundary conditions:
        K(dof,:) = zeros(nDofs,1);%K(boundIndex,dof)  ; % FILL IN
        K(dof,dof) = 1;
        GlobalForceVector(dof) =  0 ; % FILL IN
        
end

% solve for displacement:
globalDisplacementVector = K\GlobalForceVector; % FILL IN

% displ_3(f) =globalDisplacementVector(6);
% solve for reaction forces:
globalForceVector = globalStiffnessMatrix*globalDisplacementVector ;       % FILL IN

% ----------------------------------------------------------------------------------------------------------%
%                                        POST-PROCESSING
% ----------------------------------------------------------------------------------------------------------%

% calculate displaced positions
% we introduce a scale factor because the displacements are very small 
% (see input for scale factor)
x1_new = length(x1);
x2_new = length(x2);
for n = 1:nNodes
   x1_new(n) = x1(n) + scaleFactor*globalDisplacementVector(2*n-1);
   x2_new(n) = x2(n) + scaleFactor*globalDisplacementVector(2*n);
end


%plot displaced positions
for l = 1:nElements
     asmNodes = connectivities(l,1:2);
     x1_displaced = x1_new(asmNodes);
     x2_displaced = x2_new(asmNodes);
     p2 = plot(x1_displaced,x2_displaced,'r--');  
end


xlabel('x1')
ylabel('x2')
grid on
legend([p1, p2],{'original structure','displaced structure'})

% ----------------------------------------------------------------------------------------------------------%

% compute strains and stresses :
strainvector = zeros(nElements,1);
stressvector = zeros(nElements,1);

for e = 1:nElements
    mID                    = connectivities(e,3);     % element material id
    asmDofs              = elementDofs(e,:);       % assemble nodes = connectivity in terms of global nodes
    asmNodes = connectivities(e,1:2);
    % TO DO: compute the strain and stress of the element:
    x1_displaced = x1_new(asmNodes)
    x2_displaced = x2_new(asmNodes)
    point1 = nodalPositions(asmNodes(1),:)
    point2 = nodalPositions(asmNodes(2),:)
    strainvector(e)        = norm (x1_displaced - x2_displaced)/norm(point1 - point2)
    % x2_displaced = x2_new(asmNodes)
    stressvector(e)        = strainvector(e)*mprop(mID,1)
    
%    stress(f) = stressvector(1);
%    strain(f) = strainvector(1);
end
% hold off;
% plot(abs(displ_3),F);
%title('D');


% plot(strain, stress);
% grid on
% 
% xlabel("Strain")
% 
% ylabel("Stress")
% title ("Stress against Strain")
% ----------------------------------------------------------------------------------------------------------%

%Reporting

fprintf('%s','K = ')
fprintf('\n')
disp(globalStiffnessMatrix');

fprintf('%s','F = ')
fprintf('\n')
disp(globalForceVector');

fprintf('%s','U = ')
fprintf('\n')
disp(globalDisplacementVector');

fprintf('%s','Stress per element = ')
fprintf('\n')
disp(stressvector');

fprintf('%s','Strain per element = ')
fprintf('\n')
disp(strainvector');




