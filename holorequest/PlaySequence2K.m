function PlaySequence2K(slm, holograms, order)
%PLAYSEQUENCE2K Play one already-received firing-order sequence on the SLM(s).
%   This is the body of ShootSequencesMsocket2K AFTER the socket read, factored
%   out so a listener can read the sequence itself (and tell a firing-order cell
%   apart from a next-experiment holoRequest struct) before playing. The SLM
%   feed/wait logic is identical to ShootSequencesMsocket2K.

    if isempty(order)
        return
    end
    if numel(slm) ~= numel(holograms)
        disp('Number sequences must equal number SLMs')
        return
    end
    N = numel(slm);

    disp(['received sequence of length ' num2str(length(order{1}))]);
    disp(order)
    if any(cellfun(@max, order) > cellfun(@(x) size(x, 3), holograms))
        disp('ERROR: Sequence error. blanking SLM...')
        blank = zeros(size(holograms, 1), size(holograms, 2));
        for ii = 1:N
            slm(ii).feed(blank);
        end
        return
    end

    timeout = false;
    counter = ones(numel(holograms), 1);

    while ~timeout && any(counter <= (cellfun(@numel, order)'))
        for ii = 1:N
            if size(holograms{ii}, 3) >= counter(ii)
                slm(ii).feed(holograms{ii}(:, :, order{ii}(counter(ii))));
                counter(ii) = counter(ii) + 1;
            end
        end

        outcome = zeros(N, 1);
        for ii = 1:N
            outcome(ii) = slm(ii).wait();
        end

        if all(outcome == -1)
            timeout = true;
        end
    end

    if ~timeout
        disp('completed sequence to the end')
    else
        disp(['timeout while waiting to display hologram order ' num2str(counter - 1)]);
    end
end
