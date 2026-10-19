function DOF = getDOF(node, localDOF)
    DOF = 6*(node - 1) + localDOF;
end
