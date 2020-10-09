clear;clc;format long; path(path,[pwd,'\SUB']);
markers = {'s','d','o','^','v'};  
mycolor = [ 0 0 0;   255  0 0;   0 0 255;   95 58 91; 72 113 57;  27 71 116;222 110 38;139 44 42]/256;
fontsize = 12; markersize = 8;  linewidth = 1; 

%%% constants 
nproperty = 12;     Max_Exp_Data = 800;

% read parameters and exp data
[T_K,p_MPa,property,flag,weight] = textread('SUB/expdata.txt','%f%f%f%f%f','headerlines',5); 
[const_or_not,value,~] = textread('SUB/paremeters_set.txt','%f%f%s','headerlines',1); 

% seperate paramters to as fitted n0 and constant const
ndata = length(T_K);
in = 0;   ic = 0;
for ii = 1:length(const_or_not)
    if const_or_not(ii) == 1
        ic = ic + 1;   const(ic) = value(ii);
    else
        in = in + 1;   n0(in) = value(ii);
    end
end
%%% fitting
xdata = [T_K,p_MPa,property,flag,weight];
ydata = zeros(ndata,1);
% opts=optimoptions(@lsqcurvefit,'Algorithm','Levenberg-Marquardt', 'MaxFunctionEvaluations',20000000, 'MaxIterations', 20000000, 'StepTolerance',1e-66,'FunctionTolerance',1e-60,'OptimalityTolerance',1e-60);
opts = optimset('Display','on'); 
n_fit=lsqcurvefit(@(n,xdata) prebenzenefittingall(n,xdata,const,const_or_not),n0,xdata,ydata,[],[],opts);

% constract full parameter with n_fit and const
para_all = value;
in = 0;
for ii = 1:length(const_or_not)
    if const_or_not(ii) == 0
        in = in + 1;   para_all(ii) = n_fit(in);
    end
end

% calcuation using the fitted paramters and seperate different properties
prop_calc = zeros(ndata,1); 
ncount = zeros(nproperty,1);
nindex = zeros(nproperty,Max_Exp_Data);
for idata = 1:ndata
    if xdata(idata,4) == 1
        prop_calc(idata) = benzenePropsfitV(para_all,xdata(idata,1:2));
        ncount(1) = ncount(1) + 1; 
        nindex(1,ncount(1)) = idata;
    elseif xdata(idata,4) == 4
        prop_calc(idata)=benzenePropsfitalpha(para_all,xdata(idata,1:2));   
        ncount(4) = ncount(4) + 1; 
        nindex(4,ncount(4)) = idata;
    elseif xdata(idata,4) == 5
        prop_calc(idata)= benzenePropsfitcp(para_all,xdata(idata,1:2));
        ncount(5) = ncount(5) + 1; 
        nindex(5,ncount(5)) = idata;
    elseif xdata(idata,4) == 6
        prop_calc(idata)= benzenePropsfitcv(para_all,xdata(idata,1:2));
        ncount(6) = ncount(6) + 1; 
        nindex(6,ncount(6)) = idata;
    elseif xdata(idata,4) == 11
        prop_calc(idata)= benzenePropsfitH(para_all,xdata(idata,1:2));
        ncount(11) = ncount(11) + 1; 
        nindex(11,ncount(11)) = idata;
    else
        error('Flag can only be 1 4 5 6 11 now')
    end   
end
prop_exp = property.*weight;
dev = prop_exp - prop_calc;
re_dev = dev ./ prop_calc;

%%% plot the results
iplot = 0; 
for iproperty = 1:nproperty
    if ncount(iproperty) ~= 0
        iplot = iplot + 1;
        subplot(3,2,iplot)
        plot(xdata(nindex(iproperty,1:ncount(iproperty)),2),re_dev(nindex(iproperty,1:ncount(iproperty)))*100,markers{iplot},'markerfacecolor',mycolor(iplot,:),'markeredgecolor',mycolor(iplot,:),'markersize',markersize,'linewidth',linewidth-0.1);
        xlabel('\itp \rm / MPa')
        if iplot == 3 ylabel('(exp - cal)/cal * 100'); end
    end
end

set(gcf,'paperunits','centimeters');
set(gcf,'paperposition',[0 0 9 6]);
print(gcf,'-dtiff','-r600','Figures/devation.tiff');