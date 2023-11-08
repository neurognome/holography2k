function order = ShootSequencesMsocket2K(slm, sequences, control)
%updated 1/19/21 to inlcude output;

if numel(slm) ~= numel(sequences)
    disp('Number sequences must equal number SLMs')
    return
end
N = numel(slm);

control.io.flush();

sendVar = 'C';
control.io.send(sendVar);


order = [];
disp('waiting for socket to send sequence number')
while isempty(order)
    order = control.io.read(0.5); %msrecv(masterSocket,.5);
end
disp(['received sequence of length ' num2str(length(order))]);

% if any(order>size(sequences,3))
if any(max(order) > cellfun(@(x) size(x,3), sequences))
    disp('ERROR: Sequence error. blanking SLM...')
    blank = zeros(size(sequences,1),size(sequences,2));
    for ii = 1:N
        slm(ii).feed(blank);
    end
    return
end

timeout = false;
counter = 1;

while ~timeout && counter <= length(order)
    for ii = 1:N
        outcome(ii) = slm(ii).feed(sequences{ii}(:, :, order(counter)));
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