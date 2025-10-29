function [positionLAB,index_interest] = Interphase_fromResults(x_interest,InfoLAB,seed_interest)
% this function computes the position of the interphase at s of interest
% s: parametric coordinate from [1,..., # datum]

dataLAB_x = InfoLAB.LABx;
dataLAB_y = InfoLAB.LABy;
% data length
dat_leng = length(dataLAB_x);
s = (1:1:dat_leng)';                       % indices of data available 1-#Data

x_interest_km = x_interest/1000; % x in km
difX = abs(dataLAB_x-x_interest_km);    % compute the distance of X-data to X-interest

if isempty(seed_interest.area) ~= 1                     % case area of interest
    index_interest = s(difX == min(difX));
    if length(index_interest) > 1
        error('the X of interest has several-Y values')
    end 
    positionLAB = (InfoLAB.maxDepth/1000)+dataLAB_y(index_interest);
elseif isempty(seed_interest.previous_seed) ~= 1                                    % case index of interest
    seed = seed_interest.previous_seed;
    if x_interest == 0
        index_interest = 1;
        positionLAB = (InfoLAB.maxDepth/1000)+dataLAB_y(index_interest);
    else
        index_interest = seed+1;
        % check the value of X is exact
        if round(x_interest_km - dataLAB_x(index_interest),6) ~= 0 
            % check if x_interest_km is in [index_interest; index_interest+1] or 
            % [index_interest; index_interest-1] 
            in_between1 = sign(x_interest_km - dataLAB_x(index_interest)) + sign(x_interest_km - dataLAB_x(index_interest+1));
            if in_between1 == 0
                S_domain = [index_interest index_interest+1];
            else
                S_domain = [index_interest-1 index_interest];
            end 
            X_domain = dataLAB_x(S_domain);
            ds = interp1(X_domain,S_domain,x_interest_km,'linear','extrap')-S_domain(1);
            Y_domain = dataLAB_y(S_domain);
            positionLAB2 = Y_domain(1) + (Y_domain(2) - Y_domain(1))*ds;
            positionLAB = (InfoLAB.maxDepth/1000) + positionLAB2;
            index_interest = max(S_domain);
        else
            positionLAB = (InfoLAB.maxDepth/1000)+dataLAB_y(index_interest);
        end 
        

    end 
    
end 

end 