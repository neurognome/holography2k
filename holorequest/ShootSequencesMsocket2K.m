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

T=zeros([1 10E5]);
T2=zeros([1 10E5]);

O = zeros([N 10E5]);
% t=tic;

timeout = false;
counter = 1;

useSmallOrder =1;

if useSmallOrder
    for ii = 1:N
        [itemsUsed, ~, smallOrder] = unique(order);
        smallSeq = sequences{ii}(:,:,itemsUsed);

        order = smallOrder;
        sequences{ii} = smallSeq;
    end
end


saveDetails =1;
if saveDetails
    
    % SLM = Setup.SLM;
    while ~timeout && counter<=length(order)
        %disp(['now queuing hologram ' num2str(order(counter))])
        t=tic;
        for ii = 1:N
            outcome(ii) = slm(ii).feed(sequences{ii}(:, :, order(counter)));
        end
        T(counter)=toc(t);
        
        t = tic;
        % outcome = calllib('Blink_C_wrapper', 'ImageWriteComplete', 1, SLM.timeout_ms);
        
        T2(counter)=toc(t);
        O(:, counter) = outcome';
        if all(outcome == -1)
            timeout = true;
        end
        counter = counter+1;
    end
    
else
    while ~timeout && counter<=length(order)
        %disp(['now queuing hologram ' num2str(order(counter))])
        % outcome = Function_Feed_SLM(Setup.SLM, sequences(:,:,order(counter)));
        for ii = 1:N
            outcome(ii) = slm(ii).feed(sequences{ii}(:, :, order(counter)));
        end
        if all(outcome == -1)
            timeout = true;
        end
        counter = counter+1;
    end
end

if ~timeout
    disp('completed sequence to the end')
else
    disp(['timeout while waiting to display hologram order ' num2str(counter-1)]);
end

%t;