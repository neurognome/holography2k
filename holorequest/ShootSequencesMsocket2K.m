function order = ShootSequencesMsocket2K(slm, holograms, control)
%updated 1/19/21 to inlcude output;
% so i think sequences must be a cell array right? 


if numel(slm) ~= numel(holograms)
    disp('Number sequences must equal number SLMs')
    return
end
N = numel(slm);


disp('waiting for socket to send sequence number')
order = control.read(300); 
disp(['received sequence of length ' num2str(length(order))]);

if any(cellfun(@max, order) > cellfun(@(x) size(x,3), holograms))
    disp('ERROR: Sequence error. blanking SLM...')
    blank = zeros(size(holograms,1),size(holograms,2));
    for ii = 1:N
        slm(ii).feed(blank);
    end
    return
end

timeout = false;
counter = 1;

while ~timeout && counter <= max(cellfun(@length, order))
outcome = zeros(N, 1);
    for ii = 1:N
        slm(ii).feed(holograms{ii}(:, :, order{ii}(counter)));
    end

    for ii = 1:N
        outcome(ii) = slm(ii).wait();
    end

    if all(outcome == -1)
        timeout=true;
    end
    counter = counter + 1;
end

if ~timeout
    disp('completed sequence to the end')
else
    disp(['timeout while waiting to display hologram order ' num2str(counter-1)]);
end