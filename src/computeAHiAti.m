% Computes the action of the inverse Schur complement on a vector:
% AHiAti = @(v)( inv(A*inv(H)*A')*v  )
function AHiAti = computeAHiAti(Li, probData)

    % Assemble the A*inv(H)*A' matrix.
    AHiAt = computeAHiAt(Li, probData);

    [L,D,P] = ldl(AHiAt);
    AHiAti = @(v)( P*(L'\(D\(L\(P'*v)))) );

return
