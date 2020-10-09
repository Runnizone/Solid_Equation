function ydata=prebenzenefittingall(n0,xdata,const,const_or_not)

% constract full parameter with n0 and const
a1 = zeros(length(const_or_not),1);
in = 0;   ic = 0;
for ii = 1:length(const_or_not)
    if const_or_not(ii) == 1
        ic = ic + 1;   a1(ii) = const(ic);
    else
        in = in + 1;   a1(ii) = n0(in);
    end
end

% calcuate ydata targetted to zero
ndata = size(xdata,1);
ydata = zeros(ndata,1);
for idata = 1:ndata
    if xdata(idata,4) == 1
        ydata(idata)=(benzenePropsfitV(a1,xdata(idata,1:2))*xdata(idata,4)-xdata(idata,3))^2;
    elseif xdata(idata,4) == 4
        ydata(idata)=(benzenePropsfitalpha(a1,xdata(idata,1:2))*xdata(idata,4)-xdata(idata,3))^2;
    elseif xdata(idata,4) == 5
        ydata(idata)=(benzenePropsfitcp(a1,xdata(idata,1:2))*xdata(idata,4)-xdata(idata,3))^2;
    elseif xdata(idata,4) == 6
        ydata(idata)=(benzenePropsfitcv(a1,xdata(idata,1:2))*xdata(idata,4)-xdata(idata,3))^2;
    elseif xdata(idata,4) == 11
        ydata(idata)=(benzenePropsfitH(a1,xdata(idata,1:2))*xdata(idata,4)-xdata(idata,3))^2;
    else
        error('Flag can only be 1 4 5 6 11')
    end   
end
% 1 means volume (cm3 mol^-1), 
% 4 means isobaric expansivity (K^-1), 
% 5 for isobaric molar heat capacity (J K^-1 mol^-1),
% 6 for isochoric molar heat capacity (J K^-1 mol^-1) 
% 11 for enthalpy.

    

end 