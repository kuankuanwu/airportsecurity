function [x0, seed] = GenInput_DynamicStartification(seed, xmin, xmax)
    global x_set %N*2 matrix (N: solution space)
    global selected_idx %N*1 vector
    global cur_selected_order
    global set_dist %N*1 vector
    
    global Budget
    global Beta1
    global Beta2
    global c1
    global c2
    
    if isempty(x_set) || (length(selected_idx) == length(set_dist))
        %Not initialized yet or number of selected points filled the space
        %entirely
        if isempty(x_set)
            %Not initialized yet ==> Generate the solution space
            for s1 = xmin(1):xmax(1)
                Upper = floor((Budget - Beta1 - c1*s1)/c2);
                Lower = floor(max((Budget - Beta2 - c1*s1)/c2,0));
                for s2 = Lower:Upper
                    x_set = [x_set; [s1 s2]];
                end
            end
        end
        %Generate the first x0
        [seed, u] = u16807d(seed);
        x0(1) = floor(xmin(1) + u*(xmax(1) - xmin(1)) + 0.5);
        Upper = floor((Budget - Beta1 - c1*x0(1))/c2);
        Lower = floor(max((Budget - Beta2 - c1*x0(1))/c2,0));
        [seed, u] = u16807d(seed);
        x0(2) = floor(Lower + u*(Upper - Lower) + 0.5);
        
        s1_dist = abs(x_set(:,1) - x0(1));
        s2_dist = abs(x_set(:,2) - x0(2));
        set_dist = s1_dist + s2_dist;
        
        selected_idx = find(set_dist == 0);
        set_dist(selected_idx) = -1;
        
        cur_selected_order = 1;
    else
        TargetIdx = find(set_dist == max(set_dist));
        if length(TargetIdx) > 1
            [seed, u] = u16807d(seed);
            TmpIdx = floor(1 + u*(length(TargetIdx) - 1) + 0.5);
            TargetIdx = TargetIdx(TmpIdx);
        end
        x0 = x_set(TargetIdx,:);
        
        s1_dist = abs(x_set(:,1) - x0(1));
        s2_dist = abs(x_set(:,2) - x0(2));
        Tmp_dist = s1_dist + s2_dist;
        %Update the new dist
        set_dist = min([set_dist Tmp_dist]')';
        cur_selected_order = cur_selected_order + 1;
        
        selected_idx = [selected_idx TargetIdx];
    end
    
    
end

function [iseed,u16807d]=u16807d(iseed)
%
% A linear congruential pseudorandom number generator                 
% using constant 16807 and modulus (2**31)-1.                       
% iseed = iseed*16807 (mod 2^31 -1)                          
%     INPUT
%        iseed.   Integer.                                            
%                 Chosen from [1,2147483646] on the first call.    
%                 Thereafter, the value returned from the last call
%     OUTPUT
%        iseed.   Integer. To be used in the next call.                     
%        u16807d. Real. A pseudorandom number in (0,1).   

u16807d=0;
while (u16807d<=0 || u16807d>=1)
    iseed = mod (iseed * 16807,2147483647);
    u16807d = iseed / 2147483648;
end
end
