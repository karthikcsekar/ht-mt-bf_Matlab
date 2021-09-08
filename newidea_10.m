% % Calculating Tg, Ts : CONV-DIFF + HT b/w phases + Qr
% % Calculating cCO, cCO2, ri 

clear
clc

rad_oo      = 0.015;    % radius of ore pellet
eps         = 0.45;
shpfct      = 1;
S_term      = 6*(1-eps)/(2*rad_oo*shpfct);

Pressr      = 1; % in atm
Rgasct      = 8.205E-5; % in m3-atm/K/mol
M_CO        = 28.01E-3;
M_CO2       = 44.01E-3;
M_Fe        = 55.845E-3;
M_Fe2O3     = 159.69E-3;

Tg_bot      = 1600; % in Kelvin
Tg_top      = 400;
Ts_bot      = 1600;
Ts_top      = 400;

G1_bot      = 5;%cCObot;
G1_top      = 0.5*G1_bot;
G2_bot      = 0.1*G1_bot;
G2_top      = 0.6*G1_bot;

ri_top      = 0.9999*rad_oo;

% SIMVARS
ncells      = 200+2;
bfradi      = 5;
bfmaxh      = 15; 
var_ri      = 0.9999*rad_oo.*ones(ncells,1);
vol_BF      = pi*(bfradi/2)^2*bfmaxh;           % Volume of entire BF
vol_Fp      = (4/3)*pi*(rad_oo^3);              % Volume of iron ore pellet
numb_F      = eps*vol_BF/vol_Fp;                % Number of ore particles in BF
nFinCV      = numb_F/ncells;

% Initializing variables 
var_Tg      = (Tg_bot:(Tg_top-Tg_bot)/(ncells-1):Tg_top)';
var_Ts      = (Ts_bot:(Ts_top-Ts_bot)/(ncells-1):Ts_top)';
var_G1      = (G1_bot:(G1_top-G1_bot)/(ncells-1):G1_top)';
var_G2      = (G2_bot:(G2_top-G2_bot)/(ncells-1):G2_top)';

% Dependent variables
rFe2O3      = 3950;
rho_Fe      = 7000;
ktherg      = zeros(ncells,1);      % Thermal conductivity
kthers      = 12.5*ones(ncells,1);
rhoofg      = zeros(ncells,1);      % density
rhoofs      = 3950*ones(ncells,1);
muforg      = zeros(ncells,1);      % viscosity
mufors      = zeros(ncells,1);
specpg      = zeros(ncells,1);      % specific heat
specps      = zeros(ncells,1);
mf_CO       = zeros(ncells,1);      % mole fraction
mf_CO2      = zeros(ncells,1);
delHrr      = zeros(ncells,1);      % heat of reaction
reacrF      = zeros(ncells,1);      % reaction rate of Fe2O3 reduction
sfracF      = zeros(ncells,1);      % Fe2O3 reduction degree
htcoef      = zeros(ncells,1);
diffiv      = zeros(ncells,1);

matATg      = zeros(ncells,ncells);
matATs      = zeros(ncells,ncells);
matAri      = zeros(ncells,ncells);
matAG1      = zeros(ncells,ncells);
matAG2      = zeros(ncells,ncells);
rhs_Tg      = zeros(ncells,1);
rhs_Ts      = zeros(ncells,1);
rhs_ri      = zeros(ncells,1);
rhs_G1      = zeros(ncells,1);
rhs_G2      = zeros(ncells,1);
dfracF      = zeros(ncells-1,1);

tol         = 1E-6;
tolS        = 1E-6;
tolQ        = 1E-6;
tolG1       = 1E-6;
itemax      = 100;
itemxS      = 10000;
itemxQ      = 10000;
% 1st Loop
err1        = zeros(itemax,1);
err2        = zeros(itemax,1);
G1err1      = zeros(itemax,1);
G1err2      = zeros(itemxS,1);
% 1st Loop
G2err1      = zeros(itemax,1);
G2err2      = zeros(itemxS,1);
err3        = zeros(itemxS,1);
err4        = zeros(itemxS,1);
% 1st Loop
err5        = zeros(itemxQ,1);
err6        = zeros(itemxQ,1);
G1err3      = zeros(itemxS,1);
G2err3      = zeros(itemax,1);

z                       = linspace(0,bfmaxh,ncells-2)';
zx                      = linspace(0,bfmaxh,ncells)';
dz                      = z(2)-z(1);
tottime                 = 0;

% Iteration begins ...
for varx = 1:1
    zvel_g                  =  0.002;
    zvel_s                  = -0.0000411;
% ------------------------------------
% No Source term and no coupling loop
% ------------------------------------
    tic;
    for ite=1:itemax
        % Calculating all dependent variables from current variable state 
        mf_CO                   = var_G1./(var_G1 + var_G2);
        mf_CO2                  = var_G2./(var_G1 + var_G2);
        sfracF                  = 1 - (var_ri./ri_top).^3;
        [ktherg,muforg]         = thermcond(var_Tg,var_G1,var_G2);
        [rhoofg,rhoofs]         = densityca(var_Tg,var_G1,var_G2,sfracF,Pressr);
        [specpg,specps]         = specifich(var_Tg,var_Ts,var_G1,var_G2,sfracF);
        Mav_s                   = ((1-sfracF)./(1+sfracF)).*M_Fe2O3 + 2.*sfracF.*M_Fe;
        Mav_g                   = mf_CO.*M_CO + mf_CO2.*M_CO2;
%         rhoofg                  = rhoofg./M_CO;         % molar density
%         rhoofs                  = rhoofs./M_Fe2O3;
        
        ktherg(:)               = 0.02;
        rhoofg(:)               = 1;
        rhoofs(:)               = 3950;
        specpg(:)               = 40;
        specps(:)               = 120;
        diffiv(:)               = 2E-4;
        for jj=1:ncells-1
            dfracF(jj) = sfracF(jj+1) - sfracF(jj);
        end

        % Matrix for Ts
        for i=1:ncells
            if i==1
                matATs(i,i)         = 1;
                rhs_Ts(i)           = 2*Ts_bot - var_Ts(i+1);      % Dirichlet
            elseif i==ncells
                matATs(i,i)         = 1;
                rhs_Ts(i)           = 2*Ts_top - var_Ts(i-1);
            else
                vp      = zvel_s;   
                vm      = zvel_s; 
                rm      = 2*rhoofs(i-1)*rhoofs(i)/(rhoofs(i-1) + rhoofs(i));
                cm      = 2*specps(i-1)*specps(i)/(specps(i-1) + specps(i));
                km      = 2*kthers(i-1)*kthers(i)/(kthers(i-1) + kthers(i));
                rp      = 2*rhoofs(i+1)*rhoofs(i)/(rhoofs(i+1) + rhoofs(i));
                cp      = 2*specps(i+1)*specps(i)/(specps(i+1) + specps(i));
                kp      = 2*kthers(i+1)*kthers(i)/(kthers(i+1) + kthers(i));
                
                sP      = 0;
                aP      = -vm*rm*cm + kp/dz + km/dz + sP;
                aE      = -vp*rp*cp + kp/dz;
                aW      = km/dz;
                b       = 0;
                
                matATs(i,i)     = aP;
                matATs(i,i+1)   = -aE;
                matATs(i,i-1)   = -aW;
                rhs_Ts(i)       = b;
            end
        end
        % Inverting matrix for Ts calculation
        Tsold           = var_Ts;
        var_Ts          = sormethod(matATs,var_Ts,rhs_Ts,Ts_bot,Ts_top);
        err2(ite)       = max(abs(var_Ts(2:ncells-1) - Tsold(2:ncells-1)));
        
        % Matrix for Tg
        for i=1:ncells
            if i==1
                matATg(i,i)         = 1;
                rhs_Tg(i)           = 2*Tg_bot - var_Tg(i+1);    % Dirichlet BC
            elseif i==ncells
                matATg(i,i)         = 1;
                rhs_Tg(i)           = 2*Tg_top - var_Tg(i-1);   % Dirichlet BC
            else   
                vp      = zvel_g;   
                vm      = zvel_g; 
                rm      = 2*rhoofg(i-1)*rhoofg(i)/(rhoofg(i-1) + rhoofg(i));
                cm      = 2*specpg(i-1)*specpg(i)/(specpg(i-1) + specpg(i));
                km      = 2*ktherg(i-1)*ktherg(i)/(ktherg(i-1) + ktherg(i));
                rp      = 2*rhoofg(i+1)*rhoofg(i)/(rhoofg(i+1) + rhoofg(i));
                cp      = 2*specpg(i+1)*specpg(i)/(specpg(i+1) + specpg(i));
                kp      = 2*ktherg(i+1)*ktherg(i)/(ktherg(i+1) + ktherg(i));
    
                sP      = 0;%S_term*htcoef(i)*dz;
                aP      = vp*rp*cp + kp/dz + km/dz + sP;
                aE      = kp/dz;
                aW      = vm*rm*cm + km/dz;
                b       = 0;%S_term*var_Ts(i)*dz*htcoef(i);

                matATg(i,i)     = aP;
                matATg(i,i+1)   = -aE;
                matATg(i,i-1)   = -aW;
                rhs_Tg(i)       = b;
            end
        end
        % Inverting matrix for Tg calculation
        Tgold           = var_Tg;
        var_Tg          = sormethod(matATg,var_Tg,rhs_Tg,Tg_bot,Tg_top);
        err1(ite)       = max(abs(var_Tg(2:ncells-1) - Tgold(2:ncells-1)));
        
        % Matrix for cCO
        for i=1:ncells
            if i==1
                matAG1(i,i)     = 1;
                rhs_G1(i)       = 2*G1_bot - var_G1(i+1);
            elseif i==ncells
                matAG1(i,i)     = 1;
                rhs_G1(i)       = var_G1(i-1);
            else
                vp          = zvel_g;
                vm          = zvel_g;
                Dp          = 2*diffiv(i+1)*diffiv(i)/(diffiv(i+1)+diffiv(i));
                Dm          = 2*diffiv(i-1)*diffiv(i)/(diffiv(i-1)+diffiv(i));                
                
                aP          = vp + Dp/dz + Dm/dz;
                aE          = Dp/dz;
                aW          = vm + Dm/dz;
                b           = 0;
                
                matAG1(i,i)     = aP;
                matAG1(i,i+1)   = -aE;
                matAG1(i,i-1)   = -aW;
                rhs_G1(i)       = b;
            end
            
        end
        G1old           = var_G1;
        var_G1          = sormethod_DN(matAG1,var_G1,rhs_G1,G1_bot);
        G1err1(ite)     = max(abs(var_G1(2:ncells-1) - G1old(2:ncells-1)));
        
        % Matrix for cCO2
        for i=1:ncells
            if i==1
                matAG2(i,i)     = 1;
                rhs_G2(i)       = 2*G2_bot - var_G2(i+1);
            elseif i==ncells
                matAG2(i,i)     = 1;
                rhs_G2(i)       = var_G2(i-1);
            else
                vp          = zvel_g;
                vm          = zvel_g;
                Dp          = 2*diffiv(i+1)*diffiv(i)/(diffiv(i+1)+diffiv(i));
                Dm          = 2*diffiv(i-1)*diffiv(i)/(diffiv(i-1)+diffiv(i));                
                
                aP          = vp + Dp/dz + Dm/dz;
                aE          = Dp/dz;
                aW          = vm + Dm/dz;
                b           = 0;
                
                matAG2(i,i)     = aP;
                matAG2(i,i+1)   = -aE;
                matAG2(i,i-1)   = -aW;
                rhs_G2(i)       = b;
            end
            
        end
        G2old           = var_G2;
        var_G2          = sormethod_DN(matAG2,var_G2,rhs_G2,G2_bot);
        G2err1(ite)     = max(abs(var_G2(2:ncells-1) - G2old(2:ncells-1)));

        PG1             = bfmaxh.*zvel_g./diffiv; 
        Peg             = bfmaxh.*zvel_g.*rhoofg.*specpg./ktherg; 
        Pes             = bfmaxh.*zvel_s.*rhoofs.*specps./kthers; 
        
        % CONVERGENCE CHECK
        if max([err1(ite);err2(ite);G1err1(ite);G2err1(ite)])<tol
            disp('Convergence of No-Source No-coupling loop reached @ ')
            disp(ite)
            break
        end
    end
    temp1g = var_Tg;
    temp1s = var_Ts;
    temp11 = var_G1;
    temp12 = var_G2;
    toc;
    tottime = tottime + toc;
% --------------------------------------------------
% Adding Source term - Heat transfer b/w phases &
%                      Reacted fraction addn in cCO
% --------------------------------------------------
    tic;
    for ite=1:itemxS
        % Calculating all dependent variables from current variable state 
        mf_CO                   = var_G1./(var_G1 + var_G2);
        mf_CO2                  = var_G2./(var_G1 + var_G2);
        sfracF                  = 1 - (var_ri./ri_top).^3;
        [ktherg,muforg]         = thermcond(var_Tg,var_G1,var_G2);
        [rhoofg,rhoofs]         = densityca(var_Tg,var_G1,var_G2,sfracF,Pressr);
        Mav_s                   = ((1-sfracF)./(1+sfracF)).*M_Fe2O3 + 2.*sfracF.*M_Fe;
        Mav_g                   = mf_CO.*M_CO + mf_CO2.*M_CO2;
%         rhoofg                  = rhoofg./M_CO;         % molar density
%         rhoofs                  = rhoofs./M_Fe2O3;
        [specpg,specps]         = specifich(var_Tg,var_Ts,var_G1,var_G2,sfracF);
        [htcoef]                = heattrans(muforg,rhoofg,ktherg,zvel_g,rad_oo,...
                                            specpg);
        ktherg(:)               = 0.02;
        rhoofg(:)               = 1;
        rhoofs(:)               = 3950;
        specpg(:)               = 40;
        specps(:)               = 120;
        htcoef(:)               = 8;
        diffiv(:)               = 2E-4;
        for jj=1:ncells-1
            dfracF(jj) = sfracF(jj+1) - sfracF(jj);
        end

        % Matrix for Ts
        for i=1:ncells
            if i==1
                matATs(i,i)         = 1;
                rhs_Ts(i)           = 2*Ts_bot - var_Ts(i+1);      % Dirichlet
            elseif i==ncells
                matATs(i,i)         = 1;
                rhs_Ts(i)           = 2*Ts_top - var_Ts(i-1);
            else
                vp      = zvel_s;   
                vm      = zvel_s; 
                rm      = 2*rhoofs(i-1)*rhoofs(i)/(rhoofs(i-1) + rhoofs(i));
                cm      = 2*specps(i-1)*specps(i)/(specps(i-1) + specps(i));
                km      = 2*kthers(i-1)*kthers(i)/(kthers(i-1) + kthers(i));
                rp      = 2*rhoofs(i+1)*rhoofs(i)/(rhoofs(i+1) + rhoofs(i));
                cp      = 2*specps(i+1)*specps(i)/(specps(i+1) + specps(i));
                kp      = 2*kthers(i+1)*kthers(i)/(kthers(i+1) + kthers(i));
                Qr      = 0;              
                
                sP      = S_term*htcoef(i)*dz;
                aP      = -vm*rm*cm + kp/dz + km/dz + sP;
                aE      = -vp*rp*cp + kp/dz;
                aW      = km/dz;
                b       = S_term*var_Tg(i)*dz*htcoef(i);
                
                matATs(i,i)     = aP;
                matATs(i,i+1)   = -aE;
                matATs(i,i-1)   = -aW;
                rhs_Ts(i)       = b;
            end
        end
        % Inverting matrix for Ts calculation
        Tsold           = var_Ts;
        var_Ts          = sormethod(matATs,var_Ts,rhs_Ts,Ts_bot,Ts_top);
        err3(ite)       = max(abs(var_Ts(2:ncells-1) - Tsold(2:ncells-1)));
        
        % Matrix for Tg
        for i=1:ncells
            if i==1
                matATg(i,i)         = 1;
                rhs_Tg(i)           = 2*Tg_bot - var_Tg(i+1);    % Dirichlet BC
            elseif i==ncells
                matATg(i,i)         = 1;
                rhs_Tg(i)           = 2*Tg_top - var_Tg(i-1);   % Dirichlet BC
            else   
                vp      = zvel_g;   
                vm      = zvel_g; 
                rm      = 2*rhoofg(i-1)*rhoofg(i)/(rhoofg(i-1) + rhoofg(i));
                cm      = 2*specpg(i-1)*specpg(i)/(specpg(i-1) + specpg(i));
                km      = 2*ktherg(i-1)*ktherg(i)/(ktherg(i-1) + ktherg(i));
                rp      = 2*rhoofg(i+1)*rhoofg(i)/(rhoofg(i+1) + rhoofg(i));
                cp      = 2*specpg(i+1)*specpg(i)/(specpg(i+1) + specpg(i));
                kp      = 2*ktherg(i+1)*ktherg(i)/(ktherg(i+1) + ktherg(i));
    
                sP      = S_term*htcoef(i)*dz;
                aP      = vp*rp*cp + kp/dz + km/dz + sP;
                aE      = kp/dz;
                aW      = vm*rm*cm + km/dz;
                b       = S_term*var_Ts(i)*dz*htcoef(i);

                matATg(i,i)     = aP;
                matATg(i,i+1)   = -aE;
                matATg(i,i-1)   = -aW;
                rhs_Tg(i)       = b;
            end
        end
        % Inverting matrix for Tg calculation
        Tgold           = var_Tg;
        var_Tg          = sormethod(matATg,var_Tg,rhs_Tg,Tg_bot,Tg_top);
        err4(ite)       = max(abs(var_Tg(2:ncells-1) - Tgold(2:ncells-1)));
               
        % Matrix for cCO
        for i=1:ncells
            if i==1
                matAG1(i,i)     = 1;
                rhs_G1(i)       = 2*G1_bot - var_G1(i+1);
            elseif i==ncells
                matAG1(i,i)     = 1;
                rhs_G1(i)       = var_G1(i-1);
            else
                vp          = zvel_g;
                vm          = zvel_g;
                Dp          = 2*diffiv(i+1)*diffiv(i)/(diffiv(i+1)+diffiv(i));
                Dm          = 2*diffiv(i-1)*diffiv(i)/(diffiv(i-1)+diffiv(i));   
                mp          = 2*rhoofs(i+1)*rhoofs(i)/(rhoofs(i+1)+rhoofs(i));
                mm          = 2*rhoofs(i-1)*rhoofs(i)/(rhoofs(i-1)+rhoofs(i));
%                 fp          = 2*sfracF(i+1)*sfracF(i)/(sfracF(i+1)+sfracF(i));
%                 fm          = 2*sfracF(i-1)*sfracF(i)/(sfracF(i-1)+sfracF(i));

                aP          = vp + Dp/dz + Dm/dz;
                aE          = Dp/dz;
                aW          = vm + Dm/dz;
                b           = -reacrF(i)*dz;
                
                matAG1(i,i)     = aP;
                matAG1(i,i+1)   = -aE;
                matAG1(i,i-1)   = -aW;
                rhs_G1(i)       = b;
            end
            
        end
        G1old           = var_G1;
        var_G1          = sormethod_DN(matAG1,var_G1,rhs_G1,G1_bot);
        G1err2(ite)     = max(abs(var_G1(2:ncells-1) - G1old(2:ncells-1)));
        
        % Matrix for cCO2
        for i=1:ncells
            if i==1
                matAG2(i,i)     = 1;
                rhs_G2(i)       = 2*G2_bot - var_G2(i+1);
            elseif i==ncells
                matAG2(i,i)     = 1;
                rhs_G2(i)       = var_G2(i-1);
            else
                vp          = zvel_g;
                vm          = zvel_g;
                Dp          = 2*diffiv(i+1)*diffiv(i)/(diffiv(i+1)+diffiv(i));
                Dm          = 2*diffiv(i-1)*diffiv(i)/(diffiv(i-1)+diffiv(i));                
                
                aP          = vp + Dp/dz + Dm/dz;
                aE          = Dp/dz;
                aW          = vm + Dm/dz;
                b           = reacrF(i)*dz;
                
                matAG2(i,i)     = aP;
                matAG2(i,i+1)   = -aE;
                matAG2(i,i-1)   = -aW;
                rhs_G2(i)       = b;
            end
        end
        G2old           = var_G2;
        var_G2          = sormethod_DN(matAG2,var_G2,rhs_G2,G2_bot);
        G2err2(ite)     = max(abs(var_G2(2:ncells-1) - G2old(2:ncells-1)));
        
        % Calculation for var_ri
        for ttt=1:1
            mf_CO       = var_G1./(var_G1 + var_G2);
            mf_CO2      = 1 - mf_CO;
            for i=ncells-1:-1:3
                k1          = Ffunc(var_Tg(i),var_Ts(i),Pressr,mf_CO(i),...
                                        mf_CO2(i),zvel_g,ri_top,var_ri(i));
                k2          = Ffunc(var_Tg(i),var_Ts(i),Pressr,mf_CO(i),...
                                        mf_CO2(i),zvel_g,ri_top,var_ri(i) + k1*dz/2);
                k3          = Ffunc(var_Tg(i),var_Ts(i),Pressr,mf_CO(i),...
                                        mf_CO2(i),zvel_g,ri_top,var_ri(i) + k2*dz/2);
                k4          = Ffunc(var_Tg(i),var_Ts(i),Pressr,mf_CO(i),...
                                        mf_CO2(i),zvel_g,ri_top,var_ri(i) + k3*dz);

                var_ri(i-1) = var_ri(i) + (dz/6)*(k1 + 2*k2 + 2*k3 + k4);
            end
            % For ghost cells
            var_ri(ncells)  = var_ri(ncells-1);
            var_ri(1)       = var_ri(2);
            sfracF          = 1 - (var_ri./ri_top).^3;
        end

        PG1             = bfmaxh.*zvel_g./diffiv; 
        Peg             = bfmaxh.*zvel_g.*rhoofg.*specpg./ktherg; 
        Pes             = bfmaxh.*zvel_s.*rhoofs.*specps./kthers; 

        if max([err3(ite);err4(ite);G1err2(ite);G2err2(ite)])<tolS
            disp('Convergence of Source and 2-Way Coupling loop reached @ ')
            disp(ite)
            break
        end
    end
    temp2g = var_Tg;
    temp2s = var_Ts;
    temp21 = var_G1;
    temp22 = var_G2;
    toc;
    tottime = tottime + toc;

% --------------------------------------------------
% BIG ITERATIVE LOOP
% --------------------------------------------------
    vol_CV      = 1*1*dz;
    tic;
    for ite=1:itemxQ
        % Calculating all dependent variables from current variable state 
        mf_CO                   = var_G1./(var_G1 + var_G2);
        mf_CO2                  = var_G2./(var_G1 + var_G2);
        sfracF                  = 1 - (var_ri./rad_oo).^3;
        [ktherg,muforg]         = thermcond(var_Tg,var_G1,var_G2);
        [rhoofg,rhoofs]         = densityca(var_Tg,var_G1,var_G2,sfracF,Pressr);
        Mav_s                   = ((1-sfracF)./(1+sfracF)).*M_Fe2O3 + 2.*sfracF.*M_Fe;
        Mav_g                   = mf_CO.*M_CO + mf_CO2.*M_CO2;
        rhoofg                  = rhoofg./M_CO;         % molar density
        rhoofs                  = rhoofs./M_Fe2O3;
        [specpg,specps]         = specifich(var_Tg,var_Ts,var_G1,var_G2,sfracF);
        
        [htcoef]                = heattrans(muforg,rhoofg,ktherg,zvel_g,rad_oo,...
                                            specpg);
                                        
        ktherg(:)               = 0.02;
        rhoofg(:)               = 1;
        rhoofs(:)               = 3950;
        specpg(:)               = 40;
        specps(:)               = 120;
        htcoef(:)               = 8;
        diffiv(:)               = 2E-4;
        muforg(:)               = 2E-5;
        [delHrr]                = heatofrxn(var_Ts);
%         [reacrF,Ke]             = reacxnrat(var_Ts,var_Tg,rad_oo,var_G1,var_G2,...
%                                             sfracF,var_ri,Pressr,rhoofg,muforg,...
%                                             zvel_g);
        [reac2F]                = reacxnratold(var_Ts,var_Tg,rad_oo,var_G1,var_G2,sfracF,...
                                        Pressr,rhoofg,muforg,zvel_g,eps); % in kmol/m3-K
        reacrF = nFinCV.*reac2F;
                                        
%         reacrF(:)               = reacrF./vol_CV;
        
        for jj=1:ncells-1
            dfracF(jj) = sfracF(jj+1) - sfracF(jj);
        end
        % Matrix for Ts
        for i=1:ncells
            if i==1
                matATs(i,i)         = 1;
                rhs_Ts(i)           = 2*Ts_bot - var_Ts(i+1);      % Dirichlet
            elseif i==ncells
                matATs(i,i)         = 1;
                rhs_Ts(i)           = var_Ts(i-1);
            else
                vp      = zvel_s;   
                vm      = zvel_s; 
                rm      = 2*rhoofs(i-1)*rhoofs(i)/(rhoofs(i-1) + rhoofs(i));
                cm      = 2*specps(i-1)*specps(i)/(specps(i-1) + specps(i));
                km      = 2*kthers(i-1)*kthers(i)/(kthers(i-1) + kthers(i));
                rp      = 2*rhoofs(i+1)*rhoofs(i)/(rhoofs(i+1) + rhoofs(i));
                cp      = 2*specps(i+1)*specps(i)/(specps(i+1) + specps(i));
                kp      = 2*kthers(i+1)*kthers(i)/(kthers(i+1) + kthers(i));
                Qr      = delHrr(i)*reacrF(i);
                
                sP      = S_term*htcoef(i)*dz;
                aP      = -vm*rm*cm + kp/dz + km/dz + sP;
                aE      = -vp*rp*cp + kp/dz;
                aW      = km/dz;
                b       = S_term*var_Tg(i)*dz*htcoef(i) - 3*Qr*dz;
                
                matATs(i,i)     = aP;
                matATs(i,i+1)   = -aE;
                matATs(i,i-1)   = -aW;
                rhs_Ts(i)       = b;
            end
        end
        % Inverting matrix for Ts calculation
        Tsold           = var_Ts;
        var_Ts          = sormethod(matATs,var_Ts,rhs_Ts,Ts_bot,Ts_top);
        err5(ite)       = max(abs(var_Ts(2:ncells-1) - Tsold(2:ncells-1)));
        
        % Matrix for Tg
        for i=1:ncells
            if i==1
                matATg(i,i)         = 1;
                rhs_Tg(i)           = 2*Tg_bot - var_Tg(i+1);    % Dirichlet BC
            elseif i==ncells
                matATg(i,i)         = 1;
                rhs_Tg(i)           = 2*Tg_top - var_Tg(i-1);   % Dirichlet BC
            else   
                vp      = zvel_g;   
                vm      = zvel_g; 
                rm      = 2*rhoofg(i-1)*rhoofg(i)/(rhoofg(i-1) + rhoofg(i));
                cm      = 2*specpg(i-1)*specpg(i)/(specpg(i-1) + specpg(i));
                km      = 2*ktherg(i-1)*ktherg(i)/(ktherg(i-1) + ktherg(i));
                rp      = 2*rhoofg(i+1)*rhoofg(i)/(rhoofg(i+1) + rhoofg(i));
                cp      = 2*specpg(i+1)*specpg(i)/(specpg(i+1) + specpg(i));
                kp      = 2*ktherg(i+1)*ktherg(i)/(ktherg(i+1) + ktherg(i));
    
                sP      = S_term*htcoef(i)*dz;
                aP      = vp*rp*cp + kp/dz + km/dz + sP;
                aE      = kp/dz;
                aW      = vm*rm*cm + km/dz;
                b       = S_term*var_Ts(i)*dz*htcoef(i);

                matATg(i,i)     = aP;
                matATg(i,i+1)   = -aE;
                matATg(i,i-1)   = -aW;
                rhs_Tg(i)       = b;
            end
        end
        % Inverting matrix for Tg calculation
        Tgold           = var_Tg;
        var_Tg          = sormethod(matATg,var_Tg,rhs_Tg,Tg_bot,Tg_top);
        err6(ite)       = max(abs(var_Tg(2:ncells-1) - Tgold(2:ncells-1)));
        
        % Matrix for cCO
        for i=1:ncells
            if i==1
                matAG1(i,i)     = 1;
                rhs_G1(i)       = 2*G1_bot - var_G1(i+1);
            elseif i==ncells
                matAG1(i,i)     = 1;
                rhs_G1(i)       = var_G1(i-1);
            else
                vp          = zvel_g;
                vm          = zvel_g;
                Dp          = 2*diffiv(i+1)*diffiv(i)/(diffiv(i+1)+diffiv(i));
                Dm          = 2*diffiv(i-1)*diffiv(i)/(diffiv(i-1)+diffiv(i));   
                mp          = 2*rhoofs(i+1)*rhoofs(i)/(rhoofs(i+1)+rhoofs(i));
                mm          = 2*rhoofs(i-1)*rhoofs(i)/(rhoofs(i-1)+rhoofs(i));
                
                aP          = vp + Dp/dz + Dm/dz;
                aE          = Dp/dz;
                aW          = vm + Dm/dz;
                b           = -reacrF(i)*dz;%-(mp*dfracF(i) - mm*dfracF(i-1));
                
                matAG1(i,i)     = aP;
                matAG1(i,i+1)   = -aE;
                matAG1(i,i-1)   = -aW;
                rhs_G1(i)       = b;
            end
            
        end
        G1old           = var_G1;
        var_G1          = sormethod_DN(matAG1,var_G1,rhs_G1,G1_bot);
        G1err3(ite)     = max(abs(var_G1(2:ncells-1) - G1old(2:ncells-1)));
        
        % Matrix for cCO2
        for i=1:ncells
            if i==1
                matAG2(i,i)     = 1;
                rhs_G2(i)       = 2*G2_bot - var_G2(i+1);
            elseif i==ncells
                matAG2(i,i)     = 1;
                rhs_G2(i)       = var_G2(i-1);
            else
                vp          = zvel_g;
                vm          = zvel_g;
                Dp          = 2*diffiv(i+1)*diffiv(i)/(diffiv(i+1)+diffiv(i));
                Dm          = 2*diffiv(i-1)*diffiv(i)/(diffiv(i-1)+diffiv(i));                
                
                aP          = vp + Dp/dz + Dm/dz;
                aE          = Dp/dz;
                aW          = vm + Dm/dz;
                b           = reacrF(i)*dz;%(mp*dfracF(i) - mm*dfracF(i-1));
                
                matAG2(i,i)     = aP;
                matAG2(i,i+1)   = -aE;
                matAG2(i,i-1)   = -aW;
                rhs_G2(i)       = b;
            end
        end
        G2old           = var_G2;
        var_G2          = sormethod_DN(matAG2,var_G2,rhs_G2,G2_bot);
        G2err3(ite)     = max(abs(var_G2(2:ncells-1) - G2old(2:ncells-1)));
        
        % For var_ri
        mf_CO       = var_G1./(var_G1 + var_G2);
        mf_CO2      = 1 - mf_CO;
        for i=ncells-1:-1:3
            k1          = Ffunc(var_Tg(i),var_Ts(i),Pressr,mf_CO(i),...
                                    mf_CO2(i),zvel_g,ri_top,var_ri(i));
            k2          = Ffunc(var_Tg(i),var_Ts(i),Pressr,mf_CO(i),...
                                    mf_CO2(i),zvel_g,ri_top,var_ri(i) + k1*dz/2);
            k3          = Ffunc(var_Tg(i),var_Ts(i),Pressr,mf_CO(i),...
                                    mf_CO2(i),zvel_g,ri_top,var_ri(i) + k2*dz/2);
            k4          = Ffunc(var_Tg(i),var_Ts(i),Pressr,mf_CO(i),...
                                    mf_CO2(i),zvel_g,ri_top,var_ri(i) + k3*dz);

            var_ri(i-1) = var_ri(i) + (dz/6)*(k1 + 2*k2 + 2*k3 + k4);
        end
        % For ghost cells
        var_ri(ncells)  = var_ri(ncells-1);
        var_ri(1)       = var_ri(2);
        sfracF          = 1 - (var_ri./ri_top).^3;
        
        
        % Peclet Numbers calculations
        PG1             = bfmaxh.*zvel_g./diffiv; 
        Peg             = bfmaxh.*zvel_g.*rhoofg.*specpg./ktherg; 
        Pes             = bfmaxh.*zvel_s.*rhoofs.*specps./kthers; 
        
        % Convergence Check
        if max([err5(ite),err6(ite),G1err3(ite),G2err3(ite)])<tolQ
            disp('Convergence of Qr loop reached @ ')
            disp(ite)
            break
        end
    end
    toc;
    tottime = tottime + toc;
end

figure
plot(var_Tg(2:ncells-1),z,var_Ts(2:ncells-1),z,'b+');
title('CONV-DIFF + HT b/w phases + Qr')

figure
plot(temp2g(2:ncells-1),z,temp2s(2:ncells-1),z,'b+');
title('CONV-DIFF + HT b/w phases')

figure
plot(temp1g(2:ncells-1),z,temp1s(2:ncells-1),z,'b+');
title('Pure CONV-DIFF')

figure 
plot(sfracF(2:ncells-1),z);
title('Unreacted Fraction')

figure
plot(var_G1(2:ncells-1),z,temp11(2:ncells-1),z,'+');
title('CO: After and Before Qr')

figure
plot(var_G2(2:ncells-1),z,temp12(2:ncells-1),z,'+');
title('CO2: After and Before Qr')

figure
plot(mf_CO(2:ncells-1),z,mf_CO2(2:ncells-1),z,'+');
title('Mole fraction of CO, CO_2')


figure
t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,reacrF(2:ncells-1),z,'-r')
ax1.XColor = 'r';
ax1.YColor = 'r';
xlabel('Reaction rate (in mol/m^3-s')

ax2 = axes(t);
plot(ax2,delHrr(2:ncells-1),z,'-k')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';

disp('TOTAL TIME: ')
disp(tottime/60)


function [X]             = sormethod(A,X,b,X_bot,X_top)
    tol         = 1E-7;
    n           = length(X);
    err         = 1000;
    k           = 1;
    omega       = 0.75;
    while err>tol
        temp = X;
        for i=1:n
            sig = 0;
            for j=1:n
                if j~=i
                    sig = sig + A(i,j)*X(j);
                end
            end
            X(i) = (1-omega)*X(i) + omega*(b(i) - sig)/A(i,i);
        end
        X(1) = 2*X_bot - X(2);
        X(n) = 2*X_top - X(n-1);
        err  = max(abs(temp(2:n-1)-X(2:n-1)));  
        k    = k + 1;
    end
end
function [X]             = sormethod_DN(A,X,b,X_bot)
    tol         = 1E-7;
    n           = length(X);
    err         = 1000;
    k           = 1;
    omega       = 0.75;
    while err>tol
        temp = X;
        for i=1:n
            sig = 0;
            for j=1:n
                if j~=i
                    sig = sig + A(i,j)*X(j);
                end
            end
            X(i) = (1-omega)*X(i) + omega*(b(i) - sig)/A(i,i);
        end
        X(1) = 2*X_bot - X(2);
        X(n) = X(n-1);
        err  = max(abs(temp(2:n-1)-X(2:n-1)));  
        k    = k + 1;
    end
end
function [X]             = sormethod_ri(A,X,b,X_top)
    tol         = 1E-7;
    n           = length(X);
    err         = 1000;
    k           = 1;
    omega       = 0.75;
    while err>tol
        temp = X;
        for i=1:n
            sig = 0;
            for j=1:n
                if j~=i
                    sig = sig + A(i,j)*X(j);
                end
            end
            X(i) = (1-omega)*X(i) + omega*(b(i) - sig)/A(i,i);
        end
        X(n) = 2*X_top - X(n-1);
        err  = max(abs(temp(2:n-1)-X(2:n-1)));  
        k    = k + 1;
    end

end
function [diffiv]        = diffsivty(var_Tg)
    diffiv      = (1.583E-3).*(var_Tg./1900).^(3/2);
end
function [rhoofg,rhoofs] = densityca(var_Tg,varcCO,vacCO2,sfracF,Pressr)
    rFe2O3          = 3950;
    rho_Fe          = 7000;
    M_CO            = 28.01E-3;
    M_CO2           = 44.01E-3;
    Rgasct          = 8.205E-5;
    n               = length(var_Tg);
    mf_CO           = varcCO./(varcCO + vacCO2);
    mf_CO2          = 1 - mf_CO;
    
    rhoofs          = zeros(n,1);
    rhoofg          = zeros(n,1);
    
    for i=1:n
        rhoofs(i)	= rFe2O3*(1-sfracF(i)) + rho_Fe*sfracF(i);
        rhoofg(i)   = (Pressr/(Rgasct*var_Tg(i)))*(M_CO*mf_CO(i) ...
                                        + M_CO2*mf_CO2(i));
    end
end
function [specpg,specps] = specifich(var_Tg,var_Ts,varcCO,vacCO2,sfracF)
    n               = length(var_Tg);
    
    mf_CO           = varcCO./(varcCO + vacCO2);
    mf_CO2          = vacCO2./(varcCO + vacCO2);
    
    specpg          = zeros(n,1);
    specps          = zeros(n,1);
    
    for i=1:n
        cp_CO       = 4.184*(6.6 + 1.2E-3*var_Tg(i));
        cp_CO2      = 4.184*(10.50 + 2.4E-3*var_Tg(i) - 2E5*var_Tg(i)^(-2));
        cp_Fe2O3    = 4.184*(24.72 + 16.04E-3*var_Ts(i) - 4.23E5*var_Ts(i)^(-2));
        cp_Fe       = 23.97 + 8.36*(var_Ts(i)/1000) + 0.000277*(var_Ts(i)/1000)^2;

        specpg(i)   = (cp_CO*mf_CO(i) + cp_CO2*mf_CO2(i));
        specps(i)   = (cp_Fe2O3*(1-sfracF(i)) +cp_Fe*(sfracF(i)));
    end
end
function [ktherg,muforg] = thermcond(var_Tg,concCO,conCO2)
    M_CO        = 28.01; % in g/mol for mu_i calculations
    M_CO2       = 44.01;
    
    y_CO        = concCO./(concCO+conCO2);
    y_CO2       = 1-y_CO;
    sig_CO      = 3.590;
    sig_CO2     = 3.996;
    
    ebk_CO      = 110;
    ebk_CO2     = 190;
    
    Tgs_CO      = var_Tg./ebk_CO;
    Tgs_CO2     = var_Tg./ebk_CO2;
    
    Omeg_CO     = zeros(length(var_Tg),1);
    Omeg_CO2    = zeros(length(var_Tg),1);
    mu_CO       = zeros(length(var_Tg),1);
    mu_CO2      = zeros(length(var_Tg),1);
    kg_CO       = zeros(length(var_Tg),1);
    kg_CO2      = zeros(length(var_Tg),1);
    z_CO_CO2    = zeros(length(var_Tg),1);
    z_CO2_CO    = zeros(length(var_Tg),1);
    muforg      = zeros(length(var_Tg),1);
    ktherg      = zeros(length(var_Tg),1);
    
    for i=1:length(var_Tg)
        Omeg_CO(i)   = 1.16145/(Tgs_CO(i)^(0.14874)) + ...
                      0.52487/(exp(0.77320*Tgs_CO(i))) + ...
                      2.16178/(exp(2.43787*Tgs_CO(i)));
              
        Omeg_CO2(i)  = 1.16145/(Tgs_CO2(i)^(0.14874)) + ...
                      0.52487/(exp(0.77320*Tgs_CO2(i))) + ...
                      2.16178/(exp(2.43787*Tgs_CO2(i)));
    
        mu_CO(i)     = 2.6693E-6*sqrt(M_CO*var_Tg(i))...
                        /(sig_CO^2*Omeg_CO(i));                
        mu_CO2(i)    = 2.6693E-6*sqrt(M_CO2*var_Tg(i))...
                        /(sig_CO2^2*Omeg_CO2(i));
                    
        kg_CO(i)     = 8.328E-2*sqrt(var_Tg(i)/M_CO)...
                        /(sig_CO^2*Omeg_CO(i));                
        kg_CO2(i)    = 8.328E-2*sqrt(var_Tg(i)/M_CO2)...
                        /(sig_CO2^2*Omeg_CO2(i));
                    
        % Zeta for viscosity
        z_CO_CO2(i) = (1/(sqrt(8)*(1 + M_CO/M_CO2)^0.5))...
                        *(1 + ((mu_CO(i)/mu_CO2(i))^0.5)...
                        *((M_CO2/M_CO)^0.25))^2;
        z_CO2_CO(i) = (1/(sqrt(8)*(1 + M_CO2/M_CO)^0.5))...
                        *(1 + ((mu_CO2(i)/mu_CO(i))^0.5)...
                        *((M_CO/M_CO2)^0.25))^2;
        muforg(i)   = mu_CO(i)*y_CO(i)/(y_CO(i) + z_CO_CO2(i)*y_CO2(i))...
                     +mu_CO2(i)*y_CO2(i)/(y_CO2(i) + z_CO2_CO(i)*y_CO(i));
                 
                 
        % Zeta for thermal conductivity 
        z_CO_CO2(i) = (1/(sqrt(8)*(1 + M_CO/M_CO2)^0.5))...
                        *(1 + ((kg_CO(i)/kg_CO2(i))^0.5)...
                        *((M_CO2/M_CO)^0.25))^2;
        z_CO2_CO(i) = (1/(sqrt(8)*(1 + M_CO2/M_CO)^0.5))...
                        *(1 + ((kg_CO2(i)/kg_CO(i))^0.5)...
                        *((M_CO/M_CO2)^0.25))^2;
        ktherg(i)     = kg_CO(i)*y_CO(i)/(y_CO(i) + z_CO_CO2(i)*y_CO2(i))...
                     +kg_CO2(i)*y_CO2(i)/(y_CO2(i) + z_CO2_CO(i)*y_CO(i)); 
    end    
end
function [delHrr]        = heatofrxn(var_Ts)

    del_H1_std  = -24000;           %in J per mol
    n           = length(var_Ts);
    delHrr      = zeros(n,1);
    for i=1:n
        T_s             = var_Ts(i);
        T_g             = T_s;           % look at Kirchoff' law 
        T_0             = 298;
        cp_CO_dT        = 4.184*(6.6*T_g + 0.6E-3*T_g^2) - ...
                          4.184*(6.6*T_0 + 0.6E-3*T_0^2);

        cp_CO2_dT       = 4.184*(10.50*T_g + 1.2E-3*T_g^2 + 2E5*T_g^(-1)) - ...
                          4.184*(10.50*T_0 + 1.2E-3*T_0^2 + 2E5*T_g^(-1));

        cp_Fe2O3_dT     = 4.184*(24.72*T_s + 8.02E-3*T_s^2 + 4.23E5*T_s^(-1)) - ...
                          4.184*(24.72*T_0 + 8.02E-3*T_0^2 + 4.23E5*T_0^(-1));

        cp_Fe_dT        = 23.97*T_s + 4.18*(T_s^2/1000) + 0.000092*(T_s^3)/(1000)^2 - ...
                          23.97*T_0 + 4.18*(T_0^2/1000) + 0.000092*(T_0^3)/(1000)^2;   

        delHrr(i)       = del_H1_std + 2*cp_Fe_dT + 3*cp_CO2_dT -...
                          cp_Fe2O3_dT - 3*cp_CO_dT;
    end
end
function [htcoef]        = heattrans(muforg,rhoofg,ktherg,zvel_g,rad_oo,specpg)
    n           = length(muforg);
    htcoef      = zeros(n,1);
    for i=1:n
        Re          = 2*rad_oo*zvel_g*rhoofg(i)/muforg(i);
        Pr          = specpg(i)*muforg(i)/ktherg(i);
        Nu          = 2 + 0.6*(Re^(0.5))*(Pr^(1/3));
        htcoef(i)   = ktherg(i)*Nu/(2*rad_oo);
    end
end
function [reacrF,Ke]     = reacxnrat(var_Ts,var_Tg,rad_oo,varcCO,vacCO2,sfracF,...
                              var_ri,Pressr,rhoofg,muforg,zvel_g)

    n           = length(var_Tg);
    mf_CO       = varcCO./(varcCO + vacCO2);
    Ke          = zeros(n,1);
    D           = zeros(n,1);
    kr          = zeros(n,1);
    Deff        = zeros(n,1);
    Re          = zeros(n,1);
    Sc          = zeros(n,1);
    Sh          = zeros(n,1);
    kd          = zeros(n,1);
    ko          = zeros(n,1);
    mf_CO_e     = zeros(n,1);
    dridt       = zeros(n,1);
    reacrF      = zeros(n,1);
    cmtime      = zeros(n,1);
    R           = 82.06E-3*(1/28);
    rho_o       = 900;


    for i=1:n
        if var_Ts(i)<848
            if (sfracF(i)<0.111)
                Ke(i)   = exp(4.91 + 6235/var_Ts(i));
            else
                Ke(i)   = exp(-0.7625 + 543.3/var_Ts(i));
            end
        elseif var_Ts(i)>848
            if (sfracF(i)<0.111)
                Ke(i)   = exp(4.91 + 6235/var_Ts(i));
            elseif (sfracF(i)>0.111 && sfracF(i)<0.333)
                Ke(i)   = exp(2.13 - 2050/var_Ts(i));
            else
                Ke(i)   = exp(-2.642 + 2164/var_Ts(i));
            end
        end
        if var_Tg(i)<848
            D(i)    = 2.592E-6*(var_Tg(i)^1.78)/(Pressr*3600);
        else
            D(i)    = 2.592E-6*(var_Tg(i)^2.00)/(Pressr*3600);
        end
        kr(i)           =(347/3600)*exp(-3460/var_Tg(i));
        Deff(i)         = 0.1*D(i);
        Re(i)           = ((rad_oo/2)*zvel_g*rhoofg(i))/muforg(i);
        Sc(i)           = muforg(i)/(rhoofg(i)*D(i));
        Sh(i)           = 2 + 0.6*(Re(i)^0.5)*(Sc(i)^0.333);
        kd(i)           = 2*Sh(i)*D(i)/rad_oo;      % km term
        ko(i)           = (1/kd(i)) + rad_oo*(rad_oo - var_ri(i))/(var_ri(i)*Deff(i)) ...
                                    + rad_oo^2/(var_ri(i)^2*kr(i)*(1 + 1/Ke(i)));
%         mf_CO_e(i)      = 0.5;%1 + 1/Ke(i)^0.33;  % can make mods here
        mf_CO_e(i)      = 1/(1 + Ke(i)^(-1/3));
        reacrF(i)       = 4*pi*(rad_oo^2)*varcCO(i)*abs(mf_CO(i) - mf_CO_e(i))*ko(i);
        
        % Calculation for rhs of r_i equation
        a           = (R*var_Tg(i)/(kd(i)*(var_ri(i))^2))*(1+1/Ke(i));
        b           = ((rad_oo-var_ri(i))/(rad_oo*var_ri(i)))*...
                        (R*var_Tg(i)/Deff(i))*(1+1/Ke(i));
        c           = R*var_Tg(i)/(kr(i)*var_ri(i)^2);
        SB          = a + b + c;
        Nu          = Pressr*(mf_CO(i) - (1-mf_CO(i))/Ke(i))*...
                        (16/((var_ri(i)^2)*rho_o));
        dridt(i)    = -Nu/SB;
        
        % Calculation of completion time 
        temp        = (16/(rho_o*R*var_Tg(i)*rad_oo))*Pressr*(mf_CO(i) - ...
                      (1-mf_CO(i))/Ke(i));
        cmtime(i)   = (1/kr(i) + (1+1/Ke(i))*(rad_oo/(6*Deff(i)) + ...
                       1/(3*kd(i))))/temp;
    end
end
function [rhs]           = Ffunc(var_Tg,var_Ts,Pressr,mf_CO,mf_CO2,Vg,ro,var_ri)
    M_CO        = 28;
    M_CO2       = 44;                              
    R           = 82.06E-3*(1/M_CO);
    n           = length(var_Tg);
    Ke          = zeros(n,1);
    rho_o       = 900;
    
    F           = 1 - (var_ri/ro)^3;
    if var_Ts<848
        if (F<0.111)
            Ke   = exp(4.91 + 6235/var_Ts);
        else
            Ke   = exp(-0.7625 + 543.3/var_Ts);
        end
    elseif var_Ts>848
        if (F<0.111)
            Ke   = exp(4.91 + 6235/var_Ts);
        elseif (F>0.111 && F<0.333)
            Ke   = exp(2.13 - 2050/var_Ts);
        else
            Ke   = exp(-2.642 + 2164/var_Ts);
        end
    end
    if var_Tg<848
        D    = 2.592E-6*(var_Tg^1.78)/(Pressr*3600);
    else
        D    = 2.592E-6*(var_Tg^2.00)/(Pressr*3600);
    end       

    Deff     = 0.1*D;
    kr       = (347/3600)*exp(-3460/var_Ts);
    rho_g    = (M_CO*mf_CO + M_CO2*mf_CO2)*(273/(var_Tg*22.4));
    mu_g     = 4.960E-3*(var_Tg^1.5)/(var_Tg+103);
    Re      = (2*ro*Vg*rho_g)/mu_g;
    Sc       = mu_g/(rho_g*D);
    Sh       = 2 + 0.6*(Re^0.5)*(Sc^0.333);
    km       = Sh*D/(2*ro);


    A           = R*var_Tg/(km*var_ri^2)*(1 + 1/Ke);
    B           = (1/var_ri - 1/ro)*R*var_Tg*(1 + 1/Ke)/Deff;
    C           = R*var_Tg/(kr*var_ri^2);

    rhs      = (-16/(rho_o*var_ri^2))*Pressr*(mf_CO - ...
                   mf_CO2/Ke)*(A + B + C)^(-1);        
end
function [reacrF]        = reacxnratold(var_Ts,var_Tg,rad_oo,varcCO,vacCO2,sfracF,...
                                        Pressr,rhoofg,muforg,zvel_g,eps)

    n           = length(var_Tg);
    mf_CO       = varcCO./(varcCO + vacCO2);
    Ke          = zeros(n,1);
    D           = zeros(n,1);
    kr          = zeros(n,1);
    Re          = zeros(n,1);
    Sc          = zeros(n,1);
    Sh          = zeros(n,1);
    kd          = zeros(n,1);
    ko          = zeros(n,1);
    mf_CO_e     = zeros(n,1);
    reacrF      = zeros(n,1);

    for i=1:n
        % For equilibrium constant Ke
        if var_Ts(i)<848
            if (sfracF(i)<0.111)
                Ke(i)   = exp(4.91 + 6235/var_Ts(i));
            else
                Ke(i)   = exp(-0.7625 + 543.3/var_Ts(i));
            end
        elseif var_Ts(i)>848
            if (sfracF(i)<0.111)
                Ke(i)   = exp(4.91 + 6235/var_Ts(i));
            elseif (sfracF(i)>0.111 && sfracF(i)<0.333)
                Ke(i)   = exp(2.13 - 2050/var_Ts(i));
            else
                Ke(i)   = exp(-2.642 + 2164/var_Ts(i));
            end
        end
        % For Diffusivity of CO
        if var_Tg(i)<848
            D(i)    = 2.592E-6*(var_Tg(i)^1.78)/(Pressr*3600);
        else
            D(i)    = 2.592E-6*(var_Tg(i)^2.00)/(Pressr*3600);
        end
        Ds1             = D(i)*(0.53 + 0.47*eps)*(0.238*eps + 0.04);
        % For reaction rate constant kr
        kr(i)           =(347/3600)*exp(-3460/var_Tg(i));
        % For mass transfer coefficient kf
        Re(i)           = ((rad_oo/2)*zvel_g*rhoofg(i))/muforg(i);
        Sc(i)           = muforg(i)/(rhoofg(i)*D(i));
        Sh(i)           = 2 + 0.6*(Re(i)^0.5)*(Sc(i)^0.333);
        kd(i)           = 2*Sh(i)*D(i)/rad_oo;      % km term
        
        ko(i)           = (1/kd(i)) + (rad_oo/2)*( (1-sfracF(i))^(1/3) - 1)/Ds1 + ...
                          (kr(i)*((1-sfracF(i))^(2/3))*(1 + 1/Ke(i)))^(-1);
        mf_CO_e(i)      = 1/(1 + Ke(i)^(1/3));
        reacrF(i)       = 4*pi*(rad_oo^2)*273*Pressr*(mf_CO(i) - mf_CO_e(i))...
                          /(22.4*var_Tg(i)*ko(i));
%         if reacrF(i)<0
%             reacrF(i) = 0;
%         end
    end
end
















% function [reacrF,ko] = reacxnnew(var_Ts,var_Tg,rad_oo,varcCO,vacCO2,sfracF,...
%                               var_ri,Pressr,rhoofg,muforg,zvel_g)
% 
%     n       = length(var_Tg);
%     mf_CO       = varcCO./(varcCO + vacCO2);
% 
%     for i=1:n
%         Ke(i)           = exp(7.255 + 3720./var_Ts(i));
%         kc(i)           = exp(-1.445 - 6038./var_Ts(i));
%         De(i)           = exp(4.135 - 14080./var_Ts(i))/Pressr;
%         
%         Re(i)           = ((rad_oo/2)*zvel_g*rhoofg(i))/muforg(i);
%         Sc(i)           = muforg(i)/(rhoofg(i)*D(i));
%         Sh(i)           = 2 + 0.6*(Re(i)^0.5)*(Sc(i)^0.333);
%         kd(i)           = 2*Sh(i)*De(i)/rad_oo;
%         ko(i)           = (1/kd(i)) + rad_oo*(rad_oo - var_ri(i))/(var_ri(i)*De(i)) ...
%                                     + rad_oo^2/(var_ri(i)^2*kc(i)*(1 + 1/Ke(i)));
%         mf_CO_e(i)      = 0.5;
%         reacrF(i)       = 4*pi*(rad_oo^2)*varcCO(i)*abs(mf_CO(i) ...
%                                 - mf_CO_e(i))*ko(i);
%     end
% end

% function [reacrF,dridt,cmtime] = reacxn2(var_Ts,var_Tg,rad_oo,varcCO,vacCO2,sfracF,...
%                               var_ri,Pressr,rhoofg,muforg,zvel_g)
% 
%     n           = length(var_Tg);
%     mf_CO       = varcCO./(varcCO + vacCO2);
%     Ke          = zeros(n,1);
%     D           = zeros(n,1);
%     kr          = zeros(n,1);
%     Deff        = zeros(n,1);
%     Re          = zeros(n,1);
%     Sc          = zeros(n,1);
%     Sh          = zeros(n,1);
%     kd          = zeros(n,1);
%     ko          = zeros(n,1);
%     mf_CO_e     = zeros(n,1);
%     dridt       = zeros(n,1);
%     reacrF      = zeros(n,1);
%     cmtime      = zeros(n,1);
%     R           = 82.06E-3*(1/28);
%     rho_o       = 900;
% 
% 
%     for i=1:n
%         if var_Ts(i)<848
%             if (sfracF(i)<0.111)
%                 Ke(i)   = exp(4.91 + 6235/var_Ts(i));
%             else
%                 Ke(i)   = exp(-0.7625 + 543.3/var_Ts(i));
%             end
%         elseif var_Ts(i)>848
%             if (sfracF(i)<0.111)
%                 Ke(i)   = exp(4.91 + 6235/var_Ts(i));
%             elseif (sfracF(i)>0.111 && sfracF(i)<0.333)
%                 Ke(i)   = exp(2.13 - 2050/var_Ts(i));
%             else
%                 Ke(i)   = exp(-2.642 + 2164/var_Ts(i));
%             end
%         end
%         if var_Tg(i)<848
%             D(i)    = 2.592E-6*(var_Tg(i)^1.78)/(Pressr*3600);
%         else
%             D(i)    = 2.592E-6*(var_Tg(i)^2.00)/(Pressr*3600);
%         end
%         kr(i)           =(347/3600)*exp(-3460/var_Tg(i));
%         Deff(i)         = 0.1*D(i);
%         Re(i)           = ((rad_oo/2)*zvel_g*rhoofg(i))/muforg(i);
%         Sc(i)           = muforg(i)/(rhoofg(i)*D(i));
%         Sh(i)           = 2 + 0.6*(Re(i)^0.5)*(Sc(i)^0.333);
%         kd(i)           = 2*Sh(i)*D(i)/rad_oo;      % km term
%         ko(i)           = (1/kd(i)) + rad_oo*(rad_oo - var_ri(i))/(var_ri(i)*Deff(i)) ...
%                                     + rad_oo^2/(var_ri(i)^2*kr(i)*(1 + 1/Ke(i)));
%         mf_CO_e(i)      = 0.5;%1 + 1/Ke(i)^0.33;  % can make mods here
%         reacrF(i)       = 4*pi*(rad_oo^2)*varcCO(i)*abs(mf_CO(i) - mf_CO_e(i))*ko(i);
% 
%         % Calculation for rhs of r_i equation
%         if i<n
%         a           = (R*var_Tg(i)/(kd(i)*(var_ri(i))^2))*(1+1/Ke(i));
%         b           = ((var_ri(i+1)-var_ri(i))/(var_ri(i+1)*var_ri(i)))*...
%                         (R*var_Tg(i)/Deff(i))*(1+1/Ke(i));
%         c           = R*var_Tg(i)/(kr(i)*var_ri(i)^2);
%         SB          = a + b + c;
%         Nu          = Pressr*(mf_CO(i) - (1-mf_CO(i))/Ke(i))*...
%                         (16/((var_ri(i)^2)*rho_o));
%         dridt(i)    = -Nu/SB;
%         
%         % Calculation of completion time 
%         temp        = (16/(rho_o*R*var_Tg(i)*var_ri(i+1)))*Pressr*(mf_CO(i) - ...
%                       (1-mf_CO(i))/Ke(i));
%         cmtime(i)   = (1/kr(i) + (1+1/Ke(i))*(var_ri(i+1)/(6*Deff(i)) + ...
%                        1/(3*kd(i))))/temp;
%         end
%     end
% end

%     for ite=1:itemri
%     [func_f] = Ffunc(var_Tg,var_Ts,sfracF,Pressr,...
%                                               mf_CO,mf_CO2,zvel_g,...
%                                               rad_oo,var_ri);
%         for i=1:ncells
%             if i==ncells
%                 matAri(i,i)     = 1;
%                 rhs_ri(i)       = rad_oo;
%             else
%                 aP              = -zvel_s;
%                 aE              = -zvel_s;
%                 b               = func_f(i)*dz;
%                 
%                 matAri(i,i)     = aP;
%                 matAri(i,i+1)   = -aE;
%                 rhs_ri(i)       = b;
%             end
%         end
%         riold           = var_ri;
%         var_ri          = sormethod_ri(matAri,var_ri,rhs_ri,rad_oo);
%         rierr1(ite)     = max(abs(var_ri(2:ncells-1) - riold(2:ncells-1)));
%     end
% end





















