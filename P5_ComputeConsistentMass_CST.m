function [ElementConsistentMassMatrix] = P5_ComputeConsistentMass_CST(me)

% ToDo : Implement this function. See the lecture notes for help!
ElementConsistentMassMatrix= me/12*[2 0 1 0 1 0;0 2 0 1 0 1; 1 0 2 0 1 0; 0 1 0 2 0 1; 1 0 1 0 2 0; 0 1 0 1 0 2];
end
