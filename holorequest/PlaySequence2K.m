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

    % One firing order per SLM, checked BEFORE the play loop. A short order used
    % to slip past the range check below (implicit expansion makes the comparison
    % legal) and then throw "index out of bounds" at order{ii} inside the loop --
    % which, called from a listener's serve loop, propagates outside the prime
    % try/catch and kills the listener process while the DAQ keeps gating the
    % laser. A long order was silently truncated at N instead.
    if numel(order) ~= N
        disp(['ERROR: got ' num2str(numel(order)) ' firing order(s) for ' ...
              num2str(N) ' SLM(s). Not playing.'])
        return
    end

    disp(['received sequence of length ' num2str(length(order{1}))]);
    disp(order)
    if any(cellfun(@max, order) > cellfun(@(x) size(x, 3), holograms))
        disp('ERROR: Sequence error. blanking SLM...')
        % Blank must be a FRAME, sized from the hologram stack itself. This read
        % size() of the CELL ARRAY, which is 1xN, so the "blanking" path fed each
        % SLM a 1xN row of zeros rather than a blank image.
        for ii = 1:N
            blank = zeros(size(holograms{ii}, 1), size(holograms{ii}, 2));
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
