
function benzenePropsfitH0=benzenePropsfitH0(para) 
%*************************************************************************
% Function to compute properties at specified temperature T and pressure p
% for fitting purposes
% I = +1 for molar volume (cm3 mol^-1)
% abs(I) = 2 for isothermal compressibility (MPa^-1)
% abs(I) = 3 for isentropic compressibility (MPa^-1)
% abs(I) = 4 for isobaric expansivity (K^-1)
% abs(I) = 5 for isobaric molar heat capacity (J K^-1 mol^-1)
% abs(I) = 6 for isochoric molar heat capacity (J K^-1 mol^-1)
% abs(I) = 7 for thermal Grüneisen parameter (-)
% abs(I) = 8 for internal energy (J mol^-1)
% abs(I) = 9 for entropy (J K^-1 mol^-1)
% abs(I) = 10 for Helmholtz energy (J mol^-1)
% abs(I) = 11 for enthalpy (J mol^-1
% abs(I) = 12 for Gibbs energy (J mol^-1)
%*************************************************************************

lim = 9.99999999999999E+307;
Xmax = log(sqrt(lim)) - 1;
%Read model parameters from the cell range passed

    v00 = para(1);
    c1 = para(2);
    c2 = para(3);
    c3 = para(4);
    for j=1:5
        a(j)=para(j+4);
    end
    for j=6:10
        a(j)=1;
    end
    for j=11:15
        a(j)=2;
    end
    for j=16:20
        a(j)=1;
    end
    for j=21:25
        a(j)=2;
    end
    for j=1:25
        Th(j)=para(9+j);
    end
    for j=1:5
        q(j)=para(34+j);
    end
    for j=1:25
        g(j)=para(39+j);
    end
    aa = para(65);
    bb = para(66);
    cc = para(67);
    for j = 6:25
    q(j) = 0;
    end
R = 8.31451;
for j=1:5
    if q(j)<-20||q(j)>20
        benzenePropsfitH0=1000000;
        return
    end
end

%Compute molar volume at given temperature and pressure
 T=278.674;
 p =0.004785;%triple point


    if p > 1 
        eps = 0.0000000001 * p;
    else
        eps = 0.0000000001;
    end
%v = 0.8 * v00                                   %Initial guess for root finder
    v = 70 * exp(-p / 25000);                   %Initial guess for root finder
         
    j = 0;
    perror=100;
    k=1;
    while (perror > eps) && v>40&&k<=100&&v<101&&isnan(v)==0&&isinf(v)==0
       
        %evaluate p and (dp/dv) at given T and v
        

for j = 1:25
    Gamma(j) = g(j) * (v / v00) ^ q(j);
    if abs(q(j)) > 0.0000000001 
        Theta(j) = Th(j) * exp((g(j) / q(j)) * (1 - (v / v00) ^ q(j)));
    else
        Theta(j) = Th(j) * (v00 / v) ^ g(j);
    end
end

z = (v00 / v);%10
lnz = log(z);
pj(1) = z * (c1 * lnz + c2 * lnz ^ 2 + c3 * lnz ^ 3);
pj(2) = a(1) * (R * T * Gamma(1) / v) * Debye3(Theta(1) / T);%the index of pj needs to be j + 1
pcalc = pj(1) + pj(2);

for j = 2:25
    Xi = Theta(j) / T;
    if Xi < Xmax 
        pj(j+1) = a(j) * (R * Theta(j) * Gamma(j) / v) / (exp(Xi) - 1);
    else
        pj(j+1) = 0;
    end
    pcalc = pcalc + pj(j+1);
end

dpj(1) = -(z / v) * (c1 * (1 + lnz) + c2 * lnz * (2 + lnz) + c3 * (lnz ^ 2) * (3 + lnz));
dpj(2) = (pj(2) / v) * (q(1) - 1) - a(1) * R * Theta(1) * ((Gamma(1) / v) ^ 2) * Debye3D(Theta(1) / T);%the index of dpj needs to be j + 1
dpdv = dpj(1) + dpj(2);
Xmax = log(1E+16);
for j = 2:25
    Xi = Theta(j) / T;
    if Xi < Xmax 
        dpj(j+1) = (pj(j+1) / v) * ((q(j) - 1) + Gamma(j) * Xi * exp(Xi) / (exp(Xi) - 1) - Gamma(j));
    else
        dpj(j+1) = (pj(j+1) / v) * ((q(j) - 1) + Gamma(j) * Xi - Gamma(j));
    end 
        dpdv = dpdv + dpj(j+1);
end
z = T / Th(1);
pcalc = pcalc - aa * (cc / v00) * R * Th(1) * (z ^ 4 / (1 + bb * z ^ 2)) * exp(cc * (v - v00) / v00);
dpdv = dpdv - aa * (cc / v00) ^ 2 * R * Th(1) * (z ^ 4 / (1 + bb * z ^ 2)) * exp(cc * (v - v00) / v00);

        j = j + 1;
        perror = p - pcalc;
        vstep = perror / dpdv;
        if abs(vstep) > 0.1 * v00
            vstep = sign(vstep) * 0.2 * v00;      %Limit step size to 10% of v00
        end
        v = v + perror / dpdv;
        perror = abs(perror);
k=k+1;
    end
    if j > 100
        return
    end
if v<=0 || v>100||isnan(v)==1||isinf(v)==1
benzenePropsfitH0=1000000;
    return
end


%Compute anharmonic temperature factor and derivatives, and volume factor
z = T / Th(1);
fz = z ^ 4 / (1 + bb * z ^ 2);
fz1 = 2 * z ^ 3 * (2 + bb * z ^ 2) / (1 + bb * z ^ 2) ^ 2;
fz2 = 2 * z ^ 2 * (6 - 3 * bb * z ^ 2 + bb ^ 2 * z ^ 4) / (1 + bb * z ^ 2) ^ 3;
fv = exp(cc * (v - v00) / v00);

%Compute Cv, (dp/dT) at constant v and Cp
cvj(1) = a(1) * R * (Debye3(Theta(1) / T) - (Theta(1) / T) * Debye3D(Theta(1) / T));
cv = cvj(1);
dpdt = Gamma(1) * cvj(1);
for j =2:25
    cvj(j) = a(j) * R * Einstein(Theta(j) / T);
    cv = cv + cvj(j);
    dpdt = dpdt + Gamma(j) * cvj(j);
end
cv = cv - (aa * R * T / Th(1)) * fz2 * fv;
dpdt = (dpdt / v) - aa * (cc / v00) * R * fz1 * fv;
cp= cv - T * (dpdt ^ 2) / dpdv;

%Compute isothermal compressibility, isentropic compressibility and isobaric expansivity
KappaT = -1 / (v * dpdv);
KappaS = -cv / (cp * v * dpdv);
Alpha = -(dpdt / dpdv) / v;

%Compute U, S, A and G

    U = v00 * (0.5 * c1 * lnz ^ 2 + c2 * (lnz ^ 3) / 3 + 0.25 * c3 * (lnz ^ 4));
    U = U + a(1) * R * T * Debye3(Theta(1) / T);
    for j = 2: 25
        Xi = Theta(j) / T;
        if Xi < Xmax 
            U = U + a(j) * R * Theta(j) / (exp(Xi) - 1);
        end
    end
    U = U + aa * R * Th(1) * (fz - z * fz1) * fv;
    
    Xi = Theta(1) / T;
    if Xi < Xmax 
        S = a(1) * R * ((4/ 3) * Debye3(Xi) - log(1 - exp(-Xi)));
    else
        S = 0;
    end 
    for j = 2:25
        Xi = Theta(j) / T;
        if Xi < Xmax 
            S = S + a(j) * R * (Xi / (exp(Xi) - 1) - log(1 - exp(-Xi)));
        end
    end
    S = S - aa * R * fz1 * fv;
  
H=U + p * v;
benzenePropsfitH0=H;
end