function P = conelp_getPerm(K, flag)
% Returns fill-in reducing permutation.

if flag == 0
    P = amd(K);
elseif flag == 1
    P = symamd(K);
else
%     map = metisdice(K, 128);
%     for i = 1:128,
%         figure(1)
%         spy(K(map==i-1,map==i-1))
%         pause(0.1)
%     end
%     save kkt K
%     
%     pause
    P = metis(K);
end

% TODO: add nesdis