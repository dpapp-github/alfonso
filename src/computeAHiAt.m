% Computes the Schur complement matrix
% AHiAt = A*inv(H)*A'
function AHiAt = computeAHiAt(Li, probData)

    % Relies on matlab for a single block, but "manually" handles
    % block diagonal Hessians in a for loop.
    if ~iscell(Li)
        LiAt  = Li(probData.At);  % L'\A'
        AHiAt  = LiAt.'*LiAt;
    else
        % This is the readable version of what is happening below. But the code below is more efficient.
        % AHiAt = zeros(size(probData.A,1));
        % idx = 0;
        % for i=1:length(probData.As)
        %     ni = size(probData.As{i},1);
        %     Si = Li{i}(probData.As{i});
        %     AHiAt(probData.rB(:,i),probData.rB(:,i)) = AHiAt(probData.rB(:,i),probData.rB(:,i)) + Si.'*Si;
        %     idx = idx+ni;
        % end
        % AHiAt = sparse(AHiAt);

        % Allocate space for sparse array.
        N   = sum(sum(probData.rB,1).^2);
        ijv = zeros(N,3);

        ijvshift = 0;        % current starting ijv array position (minus 1)
        sidx = 1;            % starting index of the column block
        for i=1:length(probData.Kdims)
            eidx = sidx - 1 + probData.Kdims(i); % ending index of the column block

            % Compute the small Ai*inv(Hi)*Ai' block as Si'*Si, with Si = L'\A(i)'.
            rBi = probData.frB{i};
            Si  = Li{i}(probData.As{i});
            % Distribute the entries in the full A*inv(H)*A' matrix, as in
            % AHiAt(rBi,rBi) = AHiAt(rBi,rBi) + Si.'*Si;
            ni = length(rBi);
            c  = repmat(rBi,ni,1);
            r  = kron(rBi,ones(ni,1));
            M  = Si.'*Si;
            ijv(ijvshift + (1:ni^2),:) = [r c M(:)];

            sidx = eidx + 1;
            ijvshift = ijvshift + ni^2;
        end

        % Assemble the final result.
        % Relies on the property of sparse() that repeated entries are added.
        m = size(probData.A,1);
        AHiAt = sparse(ijv(:,1), ijv(:,2), ijv(:,3), m, m);
    end

return
