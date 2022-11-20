function ys = run_anc(xn,sn,win,hop)
% Inputs:
    % xn: signal with nosie
    % sn: noise source only
    % win: window
    % hop: hop size
    % TODO: I think this will end up taking a 5th input indicated the 
    % noise cancelation method
% Output:
    % ys: denoised xn
    
xn = xn';
sn = sn';
win_len = length(win);
xlen = length(xn);
ovr_len = win_len - hop;
extra = (win_len - ovr_len) - rem( (xlen - ovr_len), (win_len - ovr_len ) );
xn = [ xn, zeros(1,extra) ]; % zero pad input signal for perferct reconstruction
sn = [ sn, zeros(1,extra) ]; % pad as well so no erros
xlen = length(xn);

num_frames = 1 + ( xlen - win_len ) / hop;
ylen = win_len + (num_frames-1)*hop;
ys = zeros(1,ylen);

for frame = 0:num_frames-1
    ind = 1+frame*hop:win_len+frame*hop;
    xn_frame = win .* xn(ind);
    sn_frame = win .* sn(ind); % window?
    
    ys_frame = xn_frame; % output the final error signal 
    
    ys(ind) = ys(ind) + ys_frame; % *window eventually i think
end


end

