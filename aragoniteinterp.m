function [k,n] = aragoniteinterp(temp)
%uses data from Burton and Walter, 1987
caln = [0.4, 1.7, 2.4];
calk = [21.8, 40.6, 45.1];
calT = [5, 25, 37];
k = interp1(calT, calk, temp);
n = interp1(calT, caln, temp);

if temp > 37
    k = 45.1;
    n = 2.4;
end

end
