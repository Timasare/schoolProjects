function [ElementConsistentMassMatrix] = P5_ComputeLumpedMass_BAR( me)


ElementConsistentMassMatrix = me/2 * eye(4,4);
 

end