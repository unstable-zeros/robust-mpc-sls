function [Delta] = DeltaOperator(n_1, n_2, k, epsilon, ismemory)
%DELTAOPERATOR A random sampling method in the inf-norm ball of radius epsilon
%in the n_1 X n_2*k matrix space. 
% ismemory = 1 if we generate uncertainty operator with memory; = 0 if we
% generate a memoryless operator (only the last n_1 X n_2 block is
% non-empty).
if nargin < 4
    ismemory = 1;
end

dim_1 = n_1; dim_2 = n_2*k;

if ismemory
    if epsilon ~= 0
        Delta = rand(dim_1, dim_2)*2 - 1;
        M = norm(Delta, inf);
        gamma = rand();
        Delta = epsilon*gamma*Delta./M;
    else
        Delta = zeros(dim_1, dim_2);
    end
else
    if epsilon ~= 0
        temp = rand(n_1, n_2)*2 - 1;
        M = norm(temp, inf);
        gamma = rand();
        temp = epsilon*gamma*temp./M;
    else
        temp = zeros(n_1, n_2);
    end
    
    if k > 1
        Delta = [zeros(n_1, n_2*(k-1)) temp];
    else
        Delta = temp;
    end
end

% Delta = [];
% N = 10000;
% for ii = 1:dim_1
%     count = 1;
%     while count < N
%         sample = rand(dim_2, 1)*2*epsilon - epsilon;
%         if norm(sample, 1) <= epsilon
%             break;
%         end
%         count = count + 1;
%     end
%     
%     if count >= N
%         sample = zeros(dim_2, 1);
%     end
%     
%     Delta = [Delta; sample'];
% end


end

