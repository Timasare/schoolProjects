function [KBc, MBc, UBc] = ...
    P5_ApplyBCs(globalForceVector,globalDisplacementVector,globalStiffnessMatrix,globalMassMatrix,nDirichletBCs,DirichletBCs,nDofsPerNode)


for boundIndex = 1:nDirichletBCs

        % Find essential boundary condition dof
        dof = nDofsPerNode * DirichletBCs(boundIndex,1)- nDofsPerNode + DirichletBCs(boundIndex,2);
        % Enforce essential boundary condition
        globalStiffnessMatrix(dof,:)   = 0;
        globalStiffnessMatrix(dof,dof) = 1.;
        globalMassMatrix(dof,:)   = 0;
        globalMassMatrix(dof,dof) = 1.;
        globalForceVector(dof) =  DirichletBCs(boundIndex,3);
        globalDisplacementVector(dof) = DirichletBCs(boundIndex,3);
        
end

KBc = globalStiffnessMatrix;
MBc = globalMassMatrix;
UBc = globalDisplacementVector;


end
