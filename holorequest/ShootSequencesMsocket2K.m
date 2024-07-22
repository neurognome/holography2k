function order = ShootSequencesMsocket2K(slm, sequences, control)
%updated 1/19/21 to inlcude output;
% so i think sequences must be a cell array right? 


if numel(slm) ~= numel(sequences)
    disp('Number sequences must equal number SLMs')
    return
end
N = numel(slm);

% for ii = 1:N
%     slm(ii).preload(sequences{ii});
% end
% disp('Preloaded sequences to SLM')

% control.io.flush();
% control.flush()

% sendVar = 'C';
% control.send(sendVar);


disp('waiting for socket to send sequence number')
order = control.read(300); %msrecv(masterSocket,.5); 
% if isempty(order)
%     order = 'done';
%     disp('Donezo')
%     return
% end
disp(['received sequence of length ' num2str(length(order))]);

% counter = 1;
% while counter <= length(order)
%     for ii = 1:N
%         slm(ii).flip(order(counter))
%         disp(['flip' num2str(ii)])
%     end
%     disp(counter)
%     counter = counter + 1;
% end    
% 
% disp('done')
% if any(order>size(sequences,3))
% if ~iscell(order)
%     order = num2cell(order);
%     disp('Manual convert...')
% end
if any(cellfun(@max, order) > cellfun(@(x) size(x,3), sequences))
    disp('ERROR: Sequence error. blanking SLM...')
    blank = zeros(size(sequences,1),size(sequences,2));
    for ii = 1:N
        slm(ii).feed(blank);
    end
    return
end
%
timeout = false;
counter = 1;
% % while counter <= length(order)
% %     for ii = 1:N
% %         slm(ii).flip(order(counter))
% %     end
% %     counter = counter + 1;
% % end
% 
while ~timeout && counter <= max(cellfun(@length, order))
outcome = zeros(N, 1);
    for ii = 1:N
        slm(ii).feed(sequences{ii}(:, :, order{ii}(counter)));
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