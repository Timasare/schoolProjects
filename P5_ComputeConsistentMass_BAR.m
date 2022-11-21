function [ElementConsistentMassMatrix] = P5_ComputeConsistentMass_BAR(mass)

% ToDo: Implement this function
ElementConsistentMassMatrix = (mass/6)*[2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2];

end