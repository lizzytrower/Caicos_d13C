%% Little Ambergris Cay diurnal engine + Suess effect evaluation
% This code was written by Lizzy Trower in Matlab R2021b based on the
% description of the diurnal engine model from Geyman and Maloof (2019). It
% was last updated in February 2024.

%This code requires a few additional functions to run: 
%   *aragoniteinterp (which calculates k and n based on T, using Burton &
%   Walter 1987
%
%   *seawaterdensity (which estimates seawater density as a function of T
%   and salinity using data from Millero & Poisson 1981)
%
%   *phreeqc_diurnalengine_Alk_DIC (which calculates solution chemistry
%   using PHREEQC, to be self-consistent with the calculations for the
%   evaporation model)

%% Part 1: import and plot data from PHREEQC evaporation calcs

evapdata = readtable('PHREEQCevapdata.csv');
evapdata_CO2degas = readtable('PHREEQCevapdata_CO2degas.csv');
Lazardata = readtable('Lazardata.xlsx');
McCaffreydata = readtable('McCaffreydata.xlsx');

figure
tiledlayout(3,2)
nexttile
plot(evapdata.sal,evapdata.Omega_arag_out,'k')
hold on
plot(evapdata_CO2degas.sal,evapdata_CO2degas.Omega_arag_out,'b')
xline(45,'--k')
xlim([35 70])
ylabel('\Omega_{arag}')
xlabel('salinity')

nexttile
plot(evapdata.sal,evapdata.pH_out,'k')
hold on
plot(evapdata_CO2degas.sal,evapdata_CO2degas.pH_out,'b')
scatter(Lazardata.sal,Lazardata.pH,'ok')
scatter(McCaffreydata.sal,McCaffreydata.pH,'dk')
xline(45,'--k')
xlim([35 70])
ylabel('pH')
xlabel('salinity')

nexttile
plot(evapdata.sal,evapdata.DIC_out,'k')
hold on
plot(evapdata_CO2degas.sal,evapdata_CO2degas.DIC_out,'b')
xline(45,'--k')
xlim([35 70])
ylabel('[DIC] (mmol/kg)')
xlabel('salinity')

nexttile
plot(evapdata.sal,evapdata.Alk_out,'k')
hold on
plot(evapdata_CO2degas.sal,evapdata_CO2degas.Alk_out,'b')
scatter(Lazardata.sal,Lazardata.TA,'ok')
xline(45,'--k')
xlim([35 70])
ylabel('Alk (mequiv/kg)')
xlabel('salinity')

nexttile
plot(evapdata.sal,evapdata.Ca_out,'k')
hold on
plot(evapdata_CO2degas.sal,evapdata_CO2degas.Ca_out,'b')
scatter(Lazardata.sal,Lazardata.Ca_mmol_kg,'ok')
scatter(McCaffreydata.sal,McCaffreydata.Ca_mmol_kg,'dk')
xline(45,'--k')
xlim([35 70])
ylabel('[Ca] (mmol/kg)')
xlabel('salinity')

nexttile
plot(evapdata.sal,evapdata.Mg_out,'k')
hold on
plot(evapdata_CO2degas.sal,evapdata_CO2degas.Mg_out,'b')
scatter(Lazardata.sal,Lazardata.Mg_mmol_kg,'ok')
scatter(McCaffreydata.sal,McCaffreydata.Mg_mmol_kg,'dk')
xline(45,'--k')
xlim([35 70])
ylabel('[Mg] (mmol/kg)')
xlabel('salinity')


%% Part 2: use DIC and Alk to drive diurnal engine
%gas exchange calculations
tempC = 26.8; %{°C}
tempK = tempC + 273; %{K}

K0 = exp(-58.0931 + 90.5697*(100/tempK) + 22.294*log(tempK/100) + ...
    35*(0.027766 - 0.025888*(tempK/100) + 0.005078*(tempK/100)^2)); %{mol/L/atm}
K0 = K0*10^3; %{mol/m^3/atm}

Sc = 2116.8 - 136.25*tempC + 4.7353*tempC^2 - 0.092307*tempC^3 + ...
    0.0007555*tempC^4; %{dimensionless}

u = 8; %{m/s}

kCO2 = 0.251*u^2*(Sc/660)^0.5; %{cm/hr}
kCO2 = kCO2/100; %{m/hr}

pCO2_atm = 420; %{uatm}

waterdepth = 1.5; %{m}
waterdensity = seawaterdensity(tempC,35); %{kg/m^3}

%define shape of forcing as a sine function since the equation isn't given
%in Geyman and Maloof (2019)
t_hr = 0:0.1:24;
time_DT = datetime(2023,7,28,0,0,0) + hours(t_hr);
period = 24; %{hr}
offset = 5; %{hr} when the sine curve will cross 0
kappa_p = 140; %{umol/kg/hr}
kappa_p_factor = kappa_p/(period/pi()); %{umol/kg/hr}
photo = kappa_p_factor*sin((2*pi()/period)*(t_hr-offset)); %{umol/kg/hr}

%precipitation kinetics
k_rate = 9*10^-9; %{mol/m^2/s}
k_rate = k_rate*60*60; %{mol/m^2/hr}
[k_BR,n_BR] = aragoniteinterp(tempC);

%carbon isotope fractionations
eps_DIC_ar = 2.7;
eps_g_DIC = -2;
eps_DIC_g = -10.3;
eps_DIC_org = -10;

%pre-allocate space to run diurnal engine for each row of evaporation model
%data
DIC_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
DIC_model(1,:) = evapdata_CO2degas.DIC_out*1000; %{umol/kg}
Alk_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
Alk_model(1,:) = evapdata_CO2degas.Alk_out*1000; %{umol/kg}
pCO2_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
pCO2_model(1,:) = evapdata_CO2degas.pCO2_out; %{uatm}
pH_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
pH_model(1,:) = evapdata_CO2degas.pH_out; %{dimensionless}
Omega_ar_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
Omega_ar_model(1,:) = evapdata_CO2degas.Omega_arag_out; %{dimensionless}
Fcarb_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
Fcarb_model(1,:) = k_rate.*(Omega_ar_model(1,:) - 1).^n_BR*10^6/...
    waterdensity/waterdepth; %{umol/kg/hr}
Fgas_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
Fgas_model(1,:) = kCO2*K0*(pCO2_model(1,:) - pCO2_atm)/waterdensity/...
    waterdepth; %{umol/kg/hr}
d13C_DIC_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
d13C_DIC_model(1,:) = 0.7;
d13C_org_model = -8*ones(length(t_hr),length(evapdata_CO2degas.sal));
d13C_org_mean = -8; %{permil}
d13C_carb_cumulative = zeros(length(t_hr),length(evapdata_CO2degas.sal));
d13C_carb_cumulative_Suess = zeros(length(t_hr),length(evapdata_CO2degas.sal));
d13C_carb_cumulative_Suess2 = zeros(length(t_hr),length(evapdata_CO2degas.sal));
d13C_carb_cumulative(1,:) = d13C_DIC_model(1,:) + eps_DIC_ar;
d13C_carb_cumulative_Suess(1,:) = d13C_carb_cumulative(1,:) + 0.8;
d13C_carb_cumulative_Suess2(1,:) = d13C_carb_cumulative(1,:) + 1.8;
CR_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
CR_model(1,:) = k_BR*(Omega_ar_model(1,:) - 1).^n_BR;

%%

for m = 1:length(evapdata_CO2degas.sal)
%run diurnal engine
for nn = 2:length(t_hr)
    delta_t_model = t_hr(nn) - t_hr(nn-1); %{hr}
    Fcarb_model(nn,m) = k_rate*(Omega_ar_model(nn-1,m) -1)^n_BR*10^6/...
        waterdensity/waterdepth; %{umol/kg/hr}
    Fgas_model(nn,m) = kCO2*K0*(pCO2_model(nn-1,m) - pCO2_atm)/waterdensity/...
        waterdepth; %{umol/kg/hr}
    DIC_model(nn,m) = DIC_model(nn-1,m) - (photo(nn) + ...
        (Fcarb_model(nn,m) + Fgas_model(nn,m)))*delta_t_model;
    Alk_model(nn,m) = Alk_model(nn-1,m) - 2*Fcarb_model(nn,m)*delta_t_model;
    [Alk_out,DIC_out,pH_out,Omega_ar_out,pCO2_out] = ...
        phreeqc_diurnalengine_Alk_DIC(tempC,evapdata_CO2degas.Ca_out(m), ...
        evapdata_CO2degas.Mg_out(m),evapdata_CO2degas.K_out(m), ...
        evapdata_CO2degas.SO4_out(m),evapdata_CO2degas.Na_out(m), ...
        evapdata_CO2degas.Cl_out(m),Alk_model(nn,m)/1000,DIC_model(nn,m)/1000);
    pCO2_model(nn,m) = pCO2_out;
    pH_model(nn,m) = pH_out; %note that this pH is on NBS scale
    Omega_ar_model(nn,m) = Omega_ar_out;
    if Fgas_model(nn,m) >= 0 && photo(nn) >= 0
        d13C_DIC_model(nn,m) = (d13C_DIC_model(nn-1,m)*DIC_model(nn-1,m) - ...
            (d13C_DIC_model(nn-1,m) + eps_DIC_ar)*Fcarb_model(nn,m)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1,m) + eps_DIC_g)*Fgas_model(nn,m)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1,m) + eps_DIC_org)*photo(nn)*delta_t_model)/...
            DIC_model(nn,m);
        d13C_org_model(nn,m) = (d13C_DIC_model(nn-1,m) + eps_DIC_org)*...
            photo(nn)*delta_t_model;
    elseif Fgas_model(nn,m) >= 0 && photo(nn) < 0
        d13C_DIC_model(nn,m) = (d13C_DIC_model(nn-1,m)*DIC_model(nn-1,m) - ...
            (d13C_DIC_model(nn-1,m) + eps_DIC_ar)*Fcarb_model(nn,m)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1,m) + eps_DIC_g)*Fgas_model(nn,m)...
            *delta_t_model - ...
            d13C_org_mean*photo(nn)*delta_t_model)/DIC_model(nn,m);
    elseif Fgas_model(nn,m) < 0 && photo(nn) >= 0
        d13C_DIC_model(nn,m) = (d13C_DIC_model(nn-1,m)*DIC_model(nn-1,m) - ...
            (d13C_DIC_model(nn-1,m) + eps_DIC_ar)*Fcarb_model(nn,m)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1,m) + eps_g_DIC)*Fgas_model(nn,m)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1,m) + eps_DIC_org)*photo(nn)*delta_t_model)/...
            DIC_model(nn,m);
        d13C_org_model(nn,m) = (d13C_DIC_model(nn-1,m) + eps_DIC_org)*...
            photo(nn)*delta_t_model;
    elseif Fgas_model(nn,m) < 0 && photo(nn) < 0
        d13C_DIC_model(nn,m) = (d13C_DIC_model(nn-1,m)*DIC_model(nn-1,m) - ...
            (d13C_DIC_model(nn-1,m) + eps_DIC_ar)*Fcarb_model(nn,m)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1,m) + eps_g_DIC)*Fgas_model(nn,m)...
            *delta_t_model - ...
            d13C_org_mean*photo(nn)*delta_t_model)/DIC_model(nn,m);
    end

    %calculate rates
    %Burton and Walter (1987)
    R_BR = k_BR*(Omega_ar_model(nn,m) - 1).^n_BR; %{umol/m^2/hr}
    CR_model(nn,m) = CR_model(nn-1,m) + R_BR*delta_t_model;

    %calculate instantaneous d13C_carb
    d13C_carb_inst = d13C_DIC_model(nn,m) + eps_DIC_ar;
    d13C_carb_inst_Suess = d13C_carb_inst + 0.8;
    d13C_carb_inst_Suess2 = d13C_carb_inst + 1.8;

%calculate cumulative d13C_carb
    d13C_carb_cumulative(nn,m) = (d13C_carb_cumulative(nn-1,m)*CR_model(nn-1,m) +...
        d13C_carb_inst*R_BR*delta_t_model)/CR_model(nn,m);
    d13C_carb_cumulative_Suess(nn,m) = (d13C_carb_cumulative_Suess(nn-1,m)*CR_model(nn-1,m) +...
        d13C_carb_inst_Suess*R_BR*delta_t_model)/CR_model(nn,m);
    d13C_carb_cumulative_Suess2(nn,m) = (d13C_carb_cumulative_Suess2(nn-1,m)*CR_model(nn-1,m) +...
        d13C_carb_inst_Suess2*R_BR*delta_t_model)/CR_model(nn,m);

end

end

%%
%plot models
figure
tiledlayout(3,2)
colors = parula(6);

nexttile
plot(t_hr,Omega_ar_model(:,1),'Color',colors(1,:))
hold on
plot(t_hr,Omega_ar_model(:,6),'Color',colors(2,:))
plot(t_hr,Omega_ar_model(:,12),'Color',colors(3,:))
plot(t_hr,Omega_ar_model(:,16),'Color',colors(4,:))
plot(t_hr,Omega_ar_model(:,37),'Color',colors(5,:))
xlim([0 24])
xlabel('hour of day')
ylabel('\Omega_{arag}')

nexttile
plot(t_hr,DIC_model(:,1)/1000,'Color',colors(1,:))
hold on
plot(t_hr,DIC_model(:,6)/1000,'Color',colors(2,:))
plot(t_hr,DIC_model(:,12)/1000,'Color',colors(3,:))
plot(t_hr,DIC_model(:,26)/1000,'Color',colors(4,:))
plot(t_hr,DIC_model(:,46)/1000,'Color',colors(5,:))
xlim([0 24])
ylim([1.9 2.7])
xlabel('hour of day')
ylabel('[DIC] (mmol/kg)')

nexttile
plot(t_hr,d13C_DIC_model(:,1),'Color',colors(1,:))
hold on
plot(t_hr,d13C_DIC_model(:,6),'Color',colors(2,:))
plot(t_hr,d13C_DIC_model(:,12),'Color',colors(3,:))
plot(t_hr,d13C_DIC_model(:,16),'Color',colors(4,:))
plot(t_hr,d13C_DIC_model(:,37),'Color',colors(5,:))
xlim([0 24])
xlabel('hour of day')
ylabel('\delta^{13}C_{DIC} (‰)')

nexttile
plot(t_hr,CR_model(:,1),'Color',colors(1,:))
hold on
plot(t_hr,CR_model(:,6),'Color',colors(2,:))
plot(t_hr,CR_model(:,12),'Color',colors(3,:))
plot(t_hr,CR_model(:,16),'Color',colors(4,:))
plot(t_hr,CR_model(:,37),'Color',colors(5,:))
xlim([0 24])
xlabel('hour of day')
ylabel('cumulative CaCO_3 precipitation (\mumol/m^2)')

nexttile
plot(t_hr,d13C_carb_cumulative(:,1),'Color',colors(1,:))
hold on
plot(t_hr,d13C_carb_cumulative(:,6),'Color',colors(2,:))
plot(t_hr,d13C_carb_cumulative(:,12),'Color',colors(3,:))
plot(t_hr,d13C_carb_cumulative(:,16),'Color',colors(4,:))
plot(t_hr,d13C_carb_cumulative(:,37),'Color',colors(5,:))
xlim([0 24])
xlabel('hour of day')
ylabel('cumulative \delta^{13}C_{carb} (‰)')

% nexttile
% histogram(CR_model(end,1:46))
% xlabel('cumulative precipitation')
% ylabel('count')
