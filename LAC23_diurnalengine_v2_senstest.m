%% Little Ambergris Cay diurnal engine sensitivity tests
% This code was written by Lizzy Trower in Matlab R2021b based on the
% description of the diurnal engine model from Geyman and Maloof (2019). It
% was last updated in April 2024.

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


%% Part 2: evaporation sensitivity test
%this sensitivity test uses the outputs from the seawater evaporative model
%to explore the magnitude of diurnal changes when modern seawater is
%concentrated by evaporation

%gas exchange calculations
tempC_evap = 26.8; %{°C}
tempK_evap = tempC_evap + 273; %{K}

K0_evap = exp(-58.0931 + 90.5697*(100/tempK_evap) + 22.294*log(tempK_evap/100) + ...
    35*(0.027766 - 0.025888*(tempK_evap/100) + 0.005078*(tempK_evap/100)^2)); %{mol/L/atm}
K0_evap = K0_evap*10^3; %{mol/m^3/atm}

Sc_evap = 2116.8 - 136.25*tempC_evap + 4.7353*tempC_evap^2 - 0.092307*tempC_evap^3 + ...
    0.0007555*tempC_evap^4; %{dimensionless}

u_evap = 8; %{m/s}

kCO2_evap = 0.251*u_evap^2*(Sc_evap/660)^0.5; %{cm/hr}
kCO2_evap = kCO2_evap/100; %{m/hr}

pCO2_atm = 420; %{uatm}

waterdepth = 1.5; %{m}
waterdensity_evap = seawaterdensity(tempC_evap,evapdata_CO2degas.sal); %{kg/m^3}

%define shape of forcing as a sine function
t_hr = 0:0.1:24;
period = 24; %{hr}
offset = 4; %{hr} when the sine curve will cross 0
kappa_p = 140; %{umol/kg/hr}
kappa_p_factor = kappa_p/(period/pi()); %{umol/kg/hr}
photo = kappa_p_factor*sin((2*pi()/period)*(t_hr-offset)); %{umol/kg/hr}

%precipitation kinetics
k_rate = 9*10^-9; %{mol/m^2/s}
k_rate = k_rate*60*60; %{mol/m^2/hr}
[k_BR,n_BR] = aragoniteinterp(tempC_evap);

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
    waterdensity_evap(1)/waterdepth; %{umol/kg/hr}
Fgas_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
Fgas_model(1,:) = kCO2_evap*K0_evap*(pCO2_model(1,:) - pCO2_atm)/waterdensity_evap(1)/...
    waterdepth; %{umol/kg/hr}
d13C_DIC_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
d13C_DIC_model(1,:) = 0.7;
d13C_org_model = -8*ones(length(t_hr),length(evapdata_CO2degas.sal));
d13C_org_mean = -8; %{permil}
d13C_carb_cumulative = zeros(length(t_hr),length(evapdata_CO2degas.sal));
d13C_carb_cumulative(1,:) = d13C_DIC_model(1,:) + eps_DIC_ar;
CR_model = zeros(length(t_hr),length(evapdata_CO2degas.sal));
CR_model(1,:) = k_BR*(Omega_ar_model(1,:) - 1).^n_BR;

for m = 1:length(evapdata_CO2degas.sal)
%run diurnal engine
for nn = 2:length(t_hr)
    delta_t_model = t_hr(nn) - t_hr(nn-1); %{hr}
    Fcarb_model(nn,m) = k_rate*(Omega_ar_model(nn-1,m) -1)^n_BR*10^6/...
        waterdensity_evap(m)/waterdepth; %{umol/kg/hr}
    Fgas_model(nn,m) = kCO2_evap*K0_evap*(pCO2_model(nn-1,m) - pCO2_atm)/waterdensity_evap(m)/...
        waterdepth; %{umol/kg/hr}
    DIC_model(nn,m) = DIC_model(nn-1,m) - (photo(nn) + ...
        (Fcarb_model(nn,m) + Fgas_model(nn,m)))*delta_t_model;
    Alk_model(nn,m) = Alk_model(nn-1,m) - 2*Fcarb_model(nn,m)*delta_t_model;
    [Alk_out,DIC_out,pH_out,Omega_ar_out,pCO2_out] = ...
        phreeqc_diurnalengine_Alk_DIC(tempC_evap,evapdata_CO2degas.Ca_out(m), ...
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

%calculate cumulative d13C_carb
    d13C_carb_cumulative(nn,m) = (d13C_carb_cumulative(nn-1,m)*CR_model(nn-1,m) +...
        d13C_carb_inst*R_BR*delta_t_model)/CR_model(nn,m);

end

end

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

%% Part 3: photosynthetic forcing sensitivity test
%set up kappa_p vector
kappa_p_vector = 150:50:600;
kappa_p_factor_sens = kappa_p_vector/(period/pi()); %{umol/kg/hr}

tempC_kp = 26.8; %{°C}
tempK_kp = tempC_kp + 273;
sal_kp = 37;
waterdensity_kp = seawaterdensity(tempC_kp,sal_kp); %{kg/m^3}

K0_kp = exp(-58.0931 + 90.5697*(100/tempK_kp) + 22.294*log(tempK_kp/100) + ...
    35*(0.027766 - 0.025888*(tempK_kp/100) + 0.005078*(tempK_kp/100)^2)); %{mol/L/atm}
K0_kp = K0_kp*10^3; %{mol/m^3/atm}

Sc_kp = 2116.8 - 136.25*tempC_kp + 4.7353*tempC_kp^2 - 0.092307*tempC_kp^3 + ...
    0.0007555*tempC_kp^4; %{dimensionless}

u_kp = 8; %{m/s}

kCO2_kp = 0.251*u_kp^2*(Sc_kp/660)^0.5; %{cm/hr}
kCO2_kp = K0_kp/100; %{m/hr}

%pre-allocate space to run diurnal engine for each row of evaporation model
%data
DIC_model_kp = zeros(length(t_hr),length(kappa_p_vector));
DIC_model_kp(1,:) = evapdata_CO2degas.DIC_out(1)*1000; %{umol/kg}
Alk_model_kp = zeros(length(t_hr),length(kappa_p_vector));
Alk_model_kp(1,:) = evapdata_CO2degas.Alk_out(1)*1000; %{umol/kg}
pCO2_model_kp = zeros(length(t_hr),length(kappa_p_vector));
pCO2_model_kp(1,:) = evapdata_CO2degas.pCO2_out(1); %{uatm}
pH_model_kp = zeros(length(t_hr),length(kappa_p_vector));
pH_model_kp(1,:) = evapdata_CO2degas.pH_out(1); %{dimensionless}
Omega_ar_model_kp = zeros(length(t_hr),length(kappa_p_vector));
Omega_ar_model_kp(1,:) = evapdata_CO2degas.Omega_arag_out(1); %{dimensionless}
Fcarb_model_kp = zeros(length(t_hr),length(kappa_p_vector));
Fcarb_model_kp(1,:) = k_rate.*(Omega_ar_model_kp(1,:) - 1).^n_BR*10^6/...
    waterdensity_kp/waterdepth; %{umol/kg/hr}
Fgas_model_kp = zeros(length(t_hr),length(kappa_p_vector));
Fgas_model_kp(1,:) = kCO2_kp*K0_kp*(pCO2_model_kp(1,:) - pCO2_atm)/waterdensity_kp/...
    waterdepth; %{umol/kg/hr}
d13C_DIC_model_kp = zeros(length(t_hr),length(kappa_p_vector));
d13C_DIC_model_kp(1,:) = 0.7;
d13C_org_model_kp = -8*ones(length(t_hr),length(kappa_p_vector));
d13C_org_mean = -8; %{permil}
d13C_carb_cumulative_kp = zeros(length(t_hr),length(kappa_p_vector));
d13C_carb_cumulative_kp(1,:) = d13C_DIC_model_kp(1,:) + eps_DIC_ar;
CR_model_kp = zeros(length(t_hr),length(kappa_p_vector));
CR_model_kp(1,:) = k_BR*(Omega_ar_model_kp(1,:) - 1).^n_BR;

for m = 1:length(kappa_p_vector)
photo_kp = kappa_p_factor_sens(m)*sin((2*pi()/period)*(t_hr-offset)); %{umol/kg/hr}
%run diurnal engine
for nn = 2:length(t_hr)
    delta_t_model_kp = t_hr(nn) - t_hr(nn-1); %{hr}
    Fcarb_model_kp(nn,m) = k_rate*(Omega_ar_model_kp(nn-1,m) -1)^n_BR*10^6/...
        waterdensity_kp/waterdepth; %{umol/kg/hr}
    Fgas_model_kp(nn,m) = kCO2_kp*K0_kp*(pCO2_model_kp(nn-1,m) - pCO2_atm)/waterdensity_kp/...
        waterdepth; %{umol/kg/hr}
    DIC_model_kp(nn,m) = DIC_model_kp(nn-1,m) - (photo_kp(nn) + ...
        (Fcarb_model_kp(nn,m) + Fgas_model_kp(nn,m)))*delta_t_model_kp;
    Alk_model_kp(nn,m) = Alk_model_kp(nn-1,m) - 2*Fcarb_model_kp(nn,m)*delta_t_model_kp;
    [Alk_out,DIC_out,pH_out,Omega_ar_out,pCO2_out] = ...
        phreeqc_diurnalengine_Alk_DIC(tempC_evap,evapdata_CO2degas.Ca_out(1), ...
        evapdata_CO2degas.Mg_out(1),evapdata_CO2degas.K_out(1), ...
        evapdata_CO2degas.SO4_out(1),evapdata_CO2degas.Na_out(1), ...
        evapdata_CO2degas.Cl_out(1),Alk_model_kp(nn,m)/1000,DIC_model_kp(nn,m)/1000);
    pCO2_model_kp(nn,m) = pCO2_out;
    pH_model_kp(nn,m) = pH_out; %note that this pH is on NBS scale
    Omega_ar_model_kp(nn,m) = Omega_ar_out;
    if Fgas_model_kp(nn,m) >= 0 && photo_kp(nn) >= 0
        d13C_DIC_model_kp(nn,m) = (d13C_DIC_model_kp(nn-1,m)*DIC_model_kp(nn-1,m) - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_DIC_ar)*Fcarb_model_kp(nn,m)...
            *delta_t_model_kp - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_DIC_g)*Fgas_model_kp(nn,m)...
            *delta_t_model_kp - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_DIC_org)*photo_kp(nn)*delta_t_model_kp)/...
            DIC_model_kp(nn,m);
        d13C_org_model_kp(nn,m) = (d13C_DIC_model_kp(nn-1,m) + eps_DIC_org)*...
            photo_kp(nn)*delta_t_model_kp;
    elseif Fgas_model_kp(nn,m) >= 0 && photo_kp(nn) < 0
        d13C_DIC_model_kp(nn,m) = (d13C_DIC_model_kp(nn-1,m)*DIC_model_kp(nn-1,m) - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_DIC_ar)*Fcarb_model_kp(nn,m)...
            *delta_t_model_kp - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_DIC_g)*Fgas_model_kp(nn,m)...
            *delta_t_model_kp - ...
            d13C_org_mean*photo_kp(nn)*delta_t_model_kp)/DIC_model_kp(nn,m);
    elseif Fgas_model_kp(nn,m) < 0 && photo_kp(nn) >= 0
        d13C_DIC_model_kp(nn,m) = (d13C_DIC_model_kp(nn-1,m)*DIC_model_kp(nn-1,m) - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_DIC_ar)*Fcarb_model_kp(nn,m)...
            *delta_t_model_kp - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_g_DIC)*Fgas_model_kp(nn,m)...
            *delta_t_model_kp - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_DIC_org)*photo_kp(nn)*delta_t_model_kp)/...
            DIC_model_kp(nn,m);
        d13C_org_model_kp(nn,m) = (d13C_DIC_model_kp(nn-1,m) + eps_DIC_org)*...
            photo_kp(nn)*delta_t_model_kp;
    elseif Fgas_model_kp(nn,m) < 0 && photo_kp(nn) < 0
        d13C_DIC_model_kp(nn,m) = (d13C_DIC_model_kp(nn-1,m)*DIC_model_kp(nn-1,m) - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_DIC_ar)*Fcarb_model_kp(nn,m)...
            *delta_t_model_kp - ...
            (d13C_DIC_model_kp(nn-1,m) + eps_g_DIC)*Fgas_model_kp(nn,m)...
            *delta_t_model_kp - ...
            d13C_org_mean*photo_kp(nn)*delta_t_model_kp)/DIC_model_kp(nn,m);
    end

    %calculate rates
    %Burton and Walter (1987)
    R_BR = k_BR*(Omega_ar_model_kp(nn,m) - 1).^n_BR; %{umol/m^2/hr}
    CR_model_kp(nn,m) = CR_model_kp(nn-1,m) + R_BR*delta_t_model_kp;

    %calculate instantaneous d13C_carb
    d13C_carb_inst = d13C_DIC_model_kp(nn,m) + eps_DIC_ar;

%calculate cumulative d13C_carb
    d13C_carb_cumulative_kp(nn,m) = (d13C_carb_cumulative_kp(nn-1,m)*CR_model_kp(nn-1,m) +...
        d13C_carb_inst*R_BR*delta_t_model_kp)/CR_model_kp(nn,m);
end

end

%plot models
figure
tiledlayout(3,2)
colors = parula(11);

nexttile
plot(t_hr,Omega_ar_model_kp(:,1),'Color',colors(1,:))
hold on
plot(t_hr,Omega_ar_model_kp(:,2),'Color',colors(2,:))
plot(t_hr,Omega_ar_model_kp(:,3),'Color',colors(3,:))
plot(t_hr,Omega_ar_model_kp(:,4),'Color',colors(4,:))
plot(t_hr,Omega_ar_model_kp(:,5),'Color',colors(5,:))
plot(t_hr,Omega_ar_model_kp(:,6),'Color',colors(6,:))
plot(t_hr,Omega_ar_model_kp(:,7),'Color',colors(7,:))
plot(t_hr,Omega_ar_model_kp(:,8),'Color',colors(8,:))
plot(t_hr,Omega_ar_model_kp(:,9),'Color',colors(9,:))
plot(t_hr,Omega_ar_model_kp(:,10),'Color',colors(10,:))
xlim([0 24])
xlabel('hour of day')
ylabel('\Omega_{arag}')

nexttile
plot(t_hr,DIC_model_kp(:,1)/1000,'Color',colors(1,:))
hold on
plot(t_hr,DIC_model_kp(:,2)/1000,'Color',colors(2,:))
plot(t_hr,DIC_model_kp(:,3)/1000,'Color',colors(3,:))
plot(t_hr,DIC_model_kp(:,4)/1000,'Color',colors(4,:))
plot(t_hr,DIC_model_kp(:,5)/1000,'Color',colors(5,:))
plot(t_hr,DIC_model_kp(:,6)/1000,'Color',colors(6,:))
plot(t_hr,DIC_model_kp(:,7)/1000,'Color',colors(7,:))
plot(t_hr,DIC_model_kp(:,8)/1000,'Color',colors(8,:))
plot(t_hr,DIC_model_kp(:,9)/1000,'Color',colors(9,:))
plot(t_hr,DIC_model_kp(:,10)/1000,'Color',colors(10,:))
xlim([0 24])
xlabel('hour of day')
ylabel('[DIC] (mmol/kg)')

nexttile
plot(t_hr,d13C_DIC_model_kp(:,1),'Color',colors(1,:))
hold on
plot(t_hr,d13C_DIC_model_kp(:,2),'Color',colors(2,:))
plot(t_hr,d13C_DIC_model_kp(:,3),'Color',colors(3,:))
plot(t_hr,d13C_DIC_model_kp(:,4),'Color',colors(4,:))
plot(t_hr,d13C_DIC_model_kp(:,5),'Color',colors(5,:))
plot(t_hr,d13C_DIC_model_kp(:,6),'Color',colors(6,:))
plot(t_hr,d13C_DIC_model_kp(:,7),'Color',colors(7,:))
plot(t_hr,d13C_DIC_model_kp(:,8),'Color',colors(8,:))
plot(t_hr,d13C_DIC_model_kp(:,9),'Color',colors(9,:))
plot(t_hr,d13C_DIC_model_kp(:,10),'Color',colors(10,:))
xlim([0 24])
xlabel('hour of day')
ylabel('\delta^{13}C_{DIC} (‰)')

nexttile
plot(t_hr,CR_model_kp(:,1),'Color',colors(1,:))
hold on
plot(t_hr,CR_model_kp(:,2),'Color',colors(2,:))
plot(t_hr,CR_model_kp(:,3),'Color',colors(3,:))
plot(t_hr,CR_model_kp(:,4),'Color',colors(4,:))
plot(t_hr,CR_model_kp(:,5),'Color',colors(5,:))
plot(t_hr,CR_model_kp(:,6),'Color',colors(6,:))
plot(t_hr,CR_model_kp(:,7),'Color',colors(7,:))
plot(t_hr,CR_model_kp(:,8),'Color',colors(8,:))
plot(t_hr,CR_model_kp(:,9),'Color',colors(9,:))
plot(t_hr,CR_model_kp(:,10),'Color',colors(10,:))
xlim([0 24])
xlabel('hour of day')
ylabel('cumulative CaCO_3 precipitation (\mumol/m^2)')

nexttile
plot(t_hr,d13C_carb_cumulative_kp(:,1),'Color',colors(1,:))
hold on
plot(t_hr,d13C_carb_cumulative_kp(:,2),'Color',colors(2,:))
plot(t_hr,d13C_carb_cumulative_kp(:,3),'Color',colors(3,:))
plot(t_hr,d13C_carb_cumulative_kp(:,4),'Color',colors(4,:))
plot(t_hr,d13C_carb_cumulative_kp(:,5),'Color',colors(5,:))
plot(t_hr,d13C_carb_cumulative_kp(:,6),'Color',colors(6,:))
plot(t_hr,d13C_carb_cumulative_kp(:,7),'Color',colors(7,:))
plot(t_hr,d13C_carb_cumulative_kp(:,8),'Color',colors(8,:))
plot(t_hr,d13C_carb_cumulative_kp(:,9),'Color',colors(9,:))
plot(t_hr,d13C_carb_cumulative_kp(:,10),'Color',colors(10,:))
xlim([0 24])
xlabel('hour of day')
ylabel('cumulative \delta^{13}C_{carb} (‰)')

%% Part 4: temperature sensitivity test

tempC_vector = 24:2:38;
tempK_T = tempC_vector + 273; %{K}
sal_T = 37;

waterdensity_T = seawaterdensity(tempC_vector,sal_T); %{kg/m^3}

K0_T = exp(-58.0931 + 90.5697*(100./tempK_T) + 22.294*log(tempK_T./100) + ...
    35*(0.027766 - 0.025888*(tempK_T./100) + 0.005078*(tempK_T./100).^2)); %{mol/L/atm}
K0_T = K0_T*10^3; %{mol/m^3/atm}

Sc_T = 2116.8 - 136.25.*tempC_vector + 4.7353.*tempC_vector.^2 - 0.092307.*tempC_vector.^3 + ...
    0.0007555.*tempC_vector.^4; %{dimensionless}

u_T = 8; %{m/s}

kCO2_T = 0.251*u_T^2.*(Sc_T./660).^0.5; %{cm/hr}
kCO2_T = kCO2_T/100; %{m/hr}

[k_T,n_T] = aragoniteinterp(tempC_vector);

%pre-allocate space to run diurnal engine for each row of evaporation model
%data
DIC_model_T = zeros(length(t_hr),length(tempC_vector));
DIC_model_T(1,:) = evapdata_CO2degas.DIC_out(1)*1000; %{umol/kg}
Alk_model_T = zeros(length(t_hr),length(tempC_vector));
Alk_model_T(1,:) = evapdata_CO2degas.Alk_out(1)*1000; %{umol/kg}
pCO2_model_T = zeros(length(t_hr),length(tempC_vector));
pCO2_model_T(1,:) = evapdata_CO2degas.pCO2_out(1); %{uatm}
pH_model_T = zeros(length(t_hr),length(tempC_vector));
pH_model_T(1,:) = evapdata_CO2degas.pH_out(1); %{dimensionless}
Omega_ar_model_T = zeros(length(t_hr),length(tempC_vector));
Omega_ar_model_T(1,:) = evapdata_CO2degas.Omega_arag_out(1); %{dimensionless}
Fcarb_model_T = zeros(length(t_hr),length(tempC_vector));
Fcarb_model_T(1,:) = k_rate.*(Omega_ar_model_T(1,:) - 1).^n_T(1)*10^6/...
    waterdensity_T(1)/waterdepth; %{umol/kg/hr}
Fgas_model_T = zeros(length(t_hr),length(tempC_vector));
Fgas_model_T(1,:) = kCO2_T(1)*tempK_T(1)*(pCO2_model_T(1,:) - pCO2_atm)/waterdensity_T(1)/...
    waterdepth; %{umol/kg/hr}
d13C_DIC_model_T = zeros(length(t_hr),length(tempC_vector));
d13C_DIC_model_T(1,:) = 0.7;
d13C_org_model_T = -8*ones(length(t_hr),length(tempC_vector));
d13C_org_mean = -8; %{permil}
d13C_carb_cumulative_T = zeros(length(t_hr),length(tempC_vector));
d13C_carb_cumulative_T(1,:) = d13C_DIC_model_T(1,:) + eps_DIC_ar;
CR_model_T = zeros(length(t_hr),length(tempC_vector));
CR_model_T(1,:) = k_T(1)*(Omega_ar_model_T(1,:) - 1).^n_T(1);

photo = kappa_p_factor*sin((2*pi()/period)*(t_hr-offset)); %{umol/kg/hr}

for m = 1:length(tempC_vector)

%run diurnal engine
for nn = 2:length(t_hr)
    delta_t_model_T = t_hr(nn) - t_hr(nn-1); %{hr}
    Fcarb_model_T(nn,m) = k_rate*(Omega_ar_model_T(nn-1,m) -1)^n_T(m)*10^6/...
        waterdensity_T(m)/waterdepth; %{umol/kg/hr}
    Fgas_model_T(nn,m) = kCO2_T(m)*tempK_T(m)*(pCO2_model_T(nn-1,m) - pCO2_atm)/waterdensity_T(m)/...
        waterdepth; %{umol/kg/hr}
    DIC_model_T(nn,m) = DIC_model_T(nn-1,m) - (photo(nn) + ...
        (Fcarb_model_T(nn,m) + Fgas_model_T(nn,m)))*delta_t_model_T;
    Alk_model_T(nn,m) = Alk_model_T(nn-1,m) - 2*Fcarb_model_T(nn,m)*delta_t_model_T;
    [Alk_out,DIC_out,pH_out,Omega_ar_out,pCO2_out] = ...
        phreeqc_diurnalengine_Alk_DIC(tempC_vector(m),evapdata_CO2degas.Ca_out(1), ...
        evapdata_CO2degas.Mg_out(1),evapdata_CO2degas.K_out(1), ...
        evapdata_CO2degas.SO4_out(1),evapdata_CO2degas.Na_out(1), ...
        evapdata_CO2degas.Cl_out(1),Alk_model_T(nn,m)/1000,DIC_model_T(nn,m)/1000);
    pCO2_model_T(nn,m) = pCO2_out;
    pH_model_T(nn,m) = pH_out; %note that this pH is on NBS scale
    Omega_ar_model_T(nn,m) = Omega_ar_out;
    if Fgas_model_T(nn,m) >= 0 && photo(nn) >= 0
        d13C_DIC_model_T(nn,m) = (d13C_DIC_model_T(nn-1,m)*DIC_model_T(nn-1,m) - ...
            (d13C_DIC_model_T(nn-1,m) + eps_DIC_ar)*Fcarb_model_T(nn,m)...
            *delta_t_model_T - ...
            (d13C_DIC_model_T(nn-1,m) + eps_DIC_g)*Fgas_model_T(nn,m)...
            *delta_t_model_T - ...
            (d13C_DIC_model_T(nn-1,m) + eps_DIC_org)*photo(nn)*delta_t_model_T)/...
            DIC_model_T(nn,m);
        d13C_org_model_T(nn,m) = (d13C_DIC_model_T(nn-1,m) + eps_DIC_org)*...
            photo(nn)*delta_t_model_T;
    elseif Fgas_model_T(nn,m) >= 0 && photo(nn) < 0
        d13C_DIC_model_T(nn,m) = (d13C_DIC_model_T(nn-1,m)*DIC_model_T(nn-1,m) - ...
            (d13C_DIC_model_T(nn-1,m) + eps_DIC_ar)*Fcarb_model_T(nn,m)...
            *delta_t_model_T - ...
            (d13C_DIC_model_T(nn-1,m) + eps_DIC_g)*Fgas_model_T(nn,m)...
            *delta_t_model_T - ...
            d13C_org_mean*photo(nn)*delta_t_model_T)/DIC_model_T(nn,m);
    elseif Fgas_model_T(nn,m) < 0 && photo(nn) >= 0
        d13C_DIC_model_T(nn,m) = (d13C_DIC_model_T(nn-1,m)*DIC_model_T(nn-1,m) - ...
            (d13C_DIC_model_T(nn-1,m) + eps_DIC_ar)*Fcarb_model_T(nn,m)...
            *delta_t_model_T - ...
            (d13C_DIC_model_T(nn-1,m) + eps_g_DIC)*Fgas_model_T(nn,m)...
            *delta_t_model_T - ...
            (d13C_DIC_model_T(nn-1,m) + eps_DIC_org)*photo(nn)*delta_t_model_T)/...
            DIC_model_T(nn,m);
        d13C_org_model_T(nn,m) = (d13C_DIC_model_T(nn-1,m) + eps_DIC_org)*...
            photo(nn)*delta_t_model_T;
    elseif Fgas_model_T(nn,m) < 0 && photo(nn) < 0
        d13C_DIC_model_T(nn,m) = (d13C_DIC_model_T(nn-1,m)*DIC_model_T(nn-1,m) - ...
            (d13C_DIC_model_T(nn-1,m) + eps_DIC_ar)*Fcarb_model_T(nn,m)...
            *delta_t_model_T - ...
            (d13C_DIC_model_T(nn-1,m) + eps_g_DIC)*Fgas_model_T(nn,m)...
            *delta_t_model_T - ...
            d13C_org_mean*photo(nn)*delta_t_model_T)/DIC_model_T(nn,m);
    end

    %calculate rates
    %Burton and Walter (1987)
    R_BR = k_T(m)*(Omega_ar_model_T(nn,m) - 1).^n_T(m); %{umol/m^2/hr}
    CR_model_T(nn,m) = CR_model_T(nn-1,m) + R_BR*delta_t_model_T;

    %calculate instantaneous d13C_carb
    d13C_carb_inst = d13C_DIC_model_T(nn,m) + eps_DIC_ar;

%calculate cumulative d13C_carb
    d13C_carb_cumulative_T(nn,m) = (d13C_carb_cumulative_T(nn-1,m)*CR_model_T(nn-1,m) +...
        d13C_carb_inst*R_BR*delta_t_model_T)/CR_model_T(nn,m);

end

end
%%
%plot models
figure
tiledlayout(3,2)
colors = parula(9);

nexttile
plot(t_hr,Omega_ar_model_T(:,1),'Color',colors(1,:))
hold on
plot(t_hr,Omega_ar_model_T(:,2),'Color',colors(2,:))
plot(t_hr,Omega_ar_model_T(:,3),'Color',colors(3,:))
plot(t_hr,Omega_ar_model_T(:,4),'Color',colors(4,:))
plot(t_hr,Omega_ar_model_T(:,5),'Color',colors(5,:))
plot(t_hr,Omega_ar_model_T(:,6),'Color',colors(6,:))
plot(t_hr,Omega_ar_model_T(:,7),'Color',colors(7,:))
plot(t_hr,Omega_ar_model_T(:,8),'Color',colors(8,:))
xlim([0 24])
ylim([3 7])
yline(1)
xlabel('hour of day')
ylabel('\Omega_{arag}')

nexttile
plot(t_hr,DIC_model_T(:,1)/1000,'Color',colors(1,:))
hold on
plot(t_hr,DIC_model_T(:,2)/1000,'Color',colors(2,:))
plot(t_hr,DIC_model_T(:,3)/1000,'Color',colors(3,:))
plot(t_hr,DIC_model_T(:,4)/1000,'Color',colors(4,:))
plot(t_hr,DIC_model_T(:,5)/1000,'Color',colors(5,:))
plot(t_hr,DIC_model_T(:,6)/1000,'Color',colors(6,:))
plot(t_hr,DIC_model_T(:,7)/1000,'Color',colors(7,:))
plot(t_hr,DIC_model_T(:,8)/1000,'Color',colors(8,:))
xlim([0 24])
ylim([1.9 2.2])
xlabel('hour of day')
ylabel('[DIC] (mmol/kg)')

nexttile
plot(t_hr,d13C_DIC_model_T(:,1),'Color',colors(1,:))
hold on
plot(t_hr,d13C_DIC_model_T(:,2),'Color',colors(2,:))
plot(t_hr,d13C_DIC_model_T(:,3),'Color',colors(3,:))
plot(t_hr,d13C_DIC_model_T(:,4),'Color',colors(4,:))
plot(t_hr,d13C_DIC_model_T(:,5),'Color',colors(5,:))
plot(t_hr,d13C_DIC_model_T(:,6),'Color',colors(6,:))
plot(t_hr,d13C_DIC_model_T(:,7),'Color',colors(7,:))
plot(t_hr,d13C_DIC_model_T(:,8),'Color',colors(8,:))
xlim([0 24])
ylim([0.4 1.2])
xlabel('hour of day')
ylabel('\delta^{13}C_{DIC} (‰)')

nexttile
plot(t_hr,CR_model_T(:,1),'Color',colors(1,:))
hold on
plot(t_hr,CR_model_T(:,2),'Color',colors(2,:))
plot(t_hr,CR_model_T(:,3),'Color',colors(3,:))
plot(t_hr,CR_model_T(:,4),'Color',colors(4,:))
plot(t_hr,CR_model_T(:,5),'Color',colors(5,:))
plot(t_hr,CR_model_T(:,6),'Color',colors(6,:))
plot(t_hr,CR_model_T(:,7),'Color',colors(7,:))
plot(t_hr,CR_model_T(:,8),'Color',colors(8,:))
xlim([0 24])
xlabel('hour of day')
ylabel('cumulative CaCO_3 precipitation (\mumol/m^2)')

nexttile
plot(t_hr,d13C_carb_cumulative_T(:,1),'Color',colors(1,:))
hold on
plot(t_hr,d13C_carb_cumulative_T(:,2),'Color',colors(2,:))
plot(t_hr,d13C_carb_cumulative_T(:,3),'Color',colors(3,:))
plot(t_hr,d13C_carb_cumulative_T(:,4),'Color',colors(4,:))
plot(t_hr,d13C_carb_cumulative_T(:,5),'Color',colors(5,:))
plot(t_hr,d13C_carb_cumulative_T(:,6),'Color',colors(6,:))
plot(t_hr,d13C_carb_cumulative_T(:,7),'Color',colors(7,:))
plot(t_hr,d13C_carb_cumulative_T(:,8),'Color',colors(8,:))
xlim([0 24])
ylim([3.2 3.7])
xlabel('hour of day')
ylabel('cumulative \delta^{13}C_{carb} (‰)')

%% Part 5: integrating one rare day with many normal days
%in this test, we will run 15 normal days, 1 rare day (k_p = 600, T = 35),
%then 15 more normal days and compare with 31 normal days

t_month = 0:0.1:24*31;

%first we'll do 31 normal days
%pre-allocate space to run diurnal engine for each row of evaporation model
%data
tempC_norm = 26.8; %{°C}
tempK_norm = tempC_norm + 273; %{K}
sal_norm = 37;

K0_norm = exp(-58.0931 + 90.5697*(100/tempK_norm) + 22.294*log(tempK_norm/100) + ...
    35*(0.027766 - 0.025888*(tempK_norm/100) + 0.005078*(tempK_norm/100)^2)); %{mol/L/atm}
K0_norm = K0_norm*10^3; %{mol/m^3/atm}

Sc_norm = 2116.8 - 136.25*tempC_norm + 4.7353*tempC_norm^2 - 0.092307*tempC_norm^3 + ...
    0.0007555*tempC_norm^4; %{dimensionless}

u_norm = 8; %{m/s}

kCO2_norm = 0.251*u_norm^2*(Sc_norm/660)^0.5; %{cm/hr}
kCO2_norm = kCO2_norm/100; %{m/hr}

pCO2_atm = 420; %{uatm}

waterdepth = 1.5; %{m}
waterdensity_norm = seawaterdensity(tempC_norm,sal_norm); %{kg/m^3}

photo_month = kappa_p_factor*sin((2*pi()/period)*(t_month-offset)); %{umol/kg/hr}

DIC_model_month_norm = zeros(length(t_month),1);
DIC_model_month_norm(1) = evapdata_CO2degas.DIC_out(1)*1000; %{umol/kg}
Alk_model_month_norm = zeros(length(t_month),1);
Alk_model_month_norm(1) = evapdata_CO2degas.Alk_out(1)*1000; %{umol/kg}
pCO2_model_month_norm = zeros(length(t_month),1);
pCO2_model_month_norm(1) = evapdata_CO2degas.pCO2_out(1); %{uatm}
pH_model_month_norm = zeros(length(t_month),1);
pH_model_month_norm(1) = evapdata_CO2degas.pH_out(1); %{dimensionless}
Omega_ar_model_month_norm = zeros(length(t_month),1);
Omega_ar_model_month_norm(1) = evapdata_CO2degas.Omega_arag_out(1); %{dimensionless}
Fcarb_model_month_norm = zeros(length(t_month),1);
Fcarb_model_month_norm(1) = k_rate.*(Omega_ar_model_month_norm(1) - 1).^n_BR*10^6/...
    waterdensity_norm/waterdepth; %{umol/kg/hr}
Fgas_model_month_norm = zeros(length(t_month),1);
Fgas_model_month_norm(1) = kCO2_norm*K0_norm*(pCO2_model_month_norm(1) - pCO2_atm)/waterdensity_norm/...
    waterdepth; %{umol/kg/hr}
d13C_DIC_model_month_norm = zeros(length(t_month),1);
d13C_DIC_model_month_norm(1) = 0.7;
d13C_org_model_month_norm = -8*ones(length(t_month),1);
d13C_org_mean = -8; %{permil}
d13C_carb_cumulative_month = zeros(length(t_month),1);
CR_model_month_norm = zeros(length(t_month),1);
CR_model_month_norm(1) = k_BR*(Omega_ar_model_month_norm(1) - 1).^n_BR;

%run diurnal engine
for nn = 2:length(t_month)
    delta_t_model_month_norm = t_month(nn) - t_month(nn-1); %{hr}
    Fcarb_model_month_norm(nn) = k_rate*(Omega_ar_model_month_norm(nn-1) -1)^n_BR*10^6/...
        waterdensity_norm/waterdepth; %{umol/kg/hr}
    Fgas_model_month_norm(nn) = kCO2_norm*K0_norm*(pCO2_model_month_norm(nn-1) - pCO2_atm)/waterdensity_norm/...
        waterdepth; %{umol/kg/hr}
    DIC_model_month_norm(nn) = DIC_model_month_norm(nn-1) - (photo_month(nn) + ...
        (Fcarb_model_month_norm(nn) + Fgas_model_month_norm(nn)))*delta_t_model_month_norm;
    Alk_model_month_norm(nn) = Alk_model_month_norm(nn-1) - 2*Fcarb_model_month_norm(nn)*delta_t_model_month_norm;
    [Alk_out,DIC_out,pH_out,Omega_ar_out,pCO2_out] = ...
        phreeqc_diurnalengine_Alk_DIC(tempC_norm,evapdata_CO2degas.Ca_out(1), ...
        evapdata_CO2degas.Mg_out(1),evapdata_CO2degas.K_out(1), ...
        evapdata_CO2degas.SO4_out(1),evapdata_CO2degas.Na_out(1), ...
        evapdata_CO2degas.Cl_out(1),Alk_model_month_norm(nn)/1000,DIC_model_month_norm(nn)/1000);
    pCO2_model_month_norm(nn) = pCO2_out;
    pH_model_month_norm(nn) = pH_out; %note that this pH is on NBS scale
    Omega_ar_model_month_norm(nn) = Omega_ar_out;
    if Fgas_model_month_norm(nn) >= 0 && photo_month(nn) >= 0
        d13C_DIC_model_month_norm(nn) = (d13C_DIC_model_month_norm(nn-1)*DIC_model_month_norm(nn-1) - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_DIC_ar)*Fcarb_model_month_norm(nn)...
            *delta_t_model_month_norm - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_DIC_g)*Fgas_model_month_norm(nn)...
            *delta_t_model_month_norm - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_DIC_org)*photo_month(nn)*delta_t_model_month_norm)/...
            DIC_model_month_norm(nn);
        d13C_org_model_month_norm(nn) = (d13C_DIC_model_month_norm(nn-1) + eps_DIC_org)*...
            photo_month(nn)*delta_t_model_month_norm;
    elseif Fgas_model_month_norm(nn) >= 0 && photo_month(nn) < 0
        d13C_DIC_model_month_norm(nn) = (d13C_DIC_model_month_norm(nn-1)*DIC_model_month_norm(nn-1) - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_DIC_ar)*Fcarb_model_month_norm(nn)...
            *delta_t_model_month_norm - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_DIC_g)*Fgas_model_month_norm(nn)...
            *delta_t_model_month_norm - ...
            d13C_org_mean*photo_month(nn)*delta_t_model_month_norm)/DIC_model_month_norm(nn);
    elseif Fgas_model_month_norm(nn) < 0 && photo_month(nn) >= 0
        d13C_DIC_model_month_norm(nn) = (d13C_DIC_model_month_norm(nn-1)*DIC_model_month_norm(nn-1) - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_DIC_ar)*Fcarb_model_month_norm(nn)...
            *delta_t_model_month_norm - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_g_DIC)*Fgas_model_month_norm(nn)...
            *delta_t_model_month_norm - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_DIC_org)*photo_month(nn)*delta_t_model_month_norm)/...
            DIC_model_month_norm(nn);
        d13C_org_model_month_norm(nn) = (d13C_DIC_model_month_norm(nn-1) + eps_DIC_org)*...
            photo_month(nn)*delta_t_model_month_norm;
    elseif Fgas_model_month_norm(nn) < 0 && photo_month(nn) < 0
        d13C_DIC_model_month_norm(nn) = (d13C_DIC_model_month_norm(nn-1)*DIC_model_month_norm(nn-1) - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_DIC_ar)*Fcarb_model_month_norm(nn)...
            *delta_t_model_month_norm - ...
            (d13C_DIC_model_month_norm(nn-1) + eps_g_DIC)*Fgas_model_month_norm(nn)...
            *delta_t_model_month_norm - ...
            d13C_org_mean*photo_month(nn)*delta_t_model_month_norm)/DIC_model_month_norm(nn);
    end

    %calculate rates
    %Burton and Walter (1987)
    R_BR = k_BR*(Omega_ar_model_month_norm(nn) - 1).^n_BR; %{umol/m^2/hr}
    CR_model_month_norm(nn) = CR_model_month_norm(nn-1) + R_BR*delta_t_model_month_norm;

    %calculate instantaneous d13C_carb
    d13C_carb_inst = d13C_DIC_model_month_norm(nn) + eps_DIC_ar;

%calculate cumulative d13C_carb
    d13C_carb_cumulative_month(nn) = (d13C_carb_cumulative_month(nn-1)*CR_model_month_norm(nn-1) +...
        d13C_carb_inst*R_BR*delta_t_model_month_norm)/CR_model_month_norm(nn);

end

%now we'll do 15 normal days, one warm high productivity day, then 15
%normal days

tempC_month_rare = zeros(length(t_month),1);
tempC_month_rare(1:24*10*15) = 26.8;
tempC_month_rare(24*10*15+1:24*10*16) = 35;
tempC_month_rare(24*10*16+1:end) = 26.8;

sal_month_rare = 37;

[k_rare,n_rare] = aragoniteinterp(tempC_month_rare);

tempK_month_rare = tempC_month_rare + 273; %{K}

K0_month_rare = exp(-58.0931 + 90.5697*(100./tempK_month_rare) + 22.294*log(tempK_month_rare./100) + ...
    35*(0.027766 - 0.025888*(tempK_month_rare./100) + 0.005078*(tempK_month_rare./100).^2)); %{mol/L/atm}
K0_month_rare = K0_month_rare*10^3; %{mol/m^3/atm}

Sc_month_rare = 2116.8 - 136.25*tempC_month_rare + 4.7353*tempC_month_rare.^2 - 0.092307*tempC_month_rare.^3 + ...
    0.0007555*tempC_month_rare.^4; %{dimensionless}

u_month_rare = 8; %{m/s}

kCO2_month_rare = 0.251*u_month_rare.^2.*(Sc_month_rare./660).^0.5; %{cm/hr}
kCO2_month_rare = kCO2_month_rare./100; %{m/hr}

pCO2_atm = 420; %{uatm}

waterdepth = 1.5; %{m}
waterdensity_month_rare = seawaterdensity(tempC_month_rare,sal_month_rare); %{kg/m^3}

kappa_p_norm = 140; %{umol/kg/hr}
kappa_p_high = 600; %{umol/kg/hr}
kappa_p_factor_rare = zeros(length(t_month),1);
kappa_p_factor_rare (1:24*10*15) = kappa_p_norm/(period/pi()); %{umol/kg/hr}
kappa_p_factor_rare (24*10*15+1:24*10*16) = kappa_p_high/(period/pi()); %{umol/kg/hr}
kappa_p_factor_rare (24*10*16+1:end) = kappa_p_norm/(period/pi()); %{umol/kg/hr}
photo_month_rare = kappa_p_factor_rare'.*sin((2*pi()/period).*(t_month-offset)); %{umol/kg/hr}

DIC_model_month_rare = zeros(length(t_month),1);
DIC_model_month_rare(1) = evapdata_CO2degas.DIC_out(1)*1000; %{umol/kg}
Alk_model_month_rare = zeros(length(t_month),1);
Alk_model_month_rare(1) = evapdata_CO2degas.Alk_out(1)*1000; %{umol/kg}
pCO2_model_month_rare = zeros(length(t_month),1);
pCO2_model_month_rare(1) = evapdata_CO2degas.pCO2_out(1); %{uatm}
pH_model_month_rare = zeros(length(t_month),1);
pH_model_month_rare(1) = evapdata_CO2degas.pH_out(1); %{dimensionless}
Omega_ar_model_month_rare = zeros(length(t_month),1);
Omega_ar_model_month_rare(1) = evapdata_CO2degas.Omega_arag_out(1); %{dimensionless}
Fcarb_model_month_rare = zeros(length(t_month),1);
Fcarb_model_month_rare(1) = k_rate.*(Omega_ar_model_month_rare(1) - 1).^n_rare(1)*10^6/...
    waterdensity_month_rare(1)/waterdepth; %{umol/kg/hr}
Fgas_model_month_rare = zeros(length(t_month),1);
Fgas_model_month_rare(1) = kCO2_month_rare(1)*K0_month_rare(1)*(pCO2_model_month_rare(1) - pCO2_atm)/waterdensity_month_rare(1)/...
    waterdepth; %{umol/kg/hr}
d13C_DIC_model_month_rare = zeros(length(t_month),1);
d13C_DIC_model_month_rare(1) = 0.7;
d13C_org_model_month_rare = -8*ones(length(t_month),1);
d13C_org_mean = -8; %{permil}
d13C_carb_cumulative_month_rare = zeros(length(t_month),1);
CR_model_month_rare = zeros(length(t_month),1);
CR_model_month_rare(1) = k_BR*(Omega_ar_model_month_rare(1) - 1).^n_rare(1);

%dissolution kinetics - Walter and Morse 1985 
n_diss = 2.5; 
k_diss = 10^3; %{umol/m^2/hr)

%set mean d13C_carb that's used for dissolving carbonate
d13C_carb_mean = 4; %{permil}

%run diurnal engine
for nn = 2:length(t_month)
    delta_t_model_month_rare = t_month(nn) - t_month(nn-1); %{hr}
    if Omega_ar_model_month_rare(nn-1) == 1
        Fcarb_model_month_rare(nn) = 0;
    elseif Omega_ar_model_month_rare(nn-1) > 1
        Fcarb_model_month_rare(nn) = k_rate*(Omega_ar_model_month_rare(nn-1) -1)^n_rare(nn)*10^6/...
            waterdensity_month_rare(nn)/waterdepth; %{umol/kg/hr}
    elseif Omega_ar_model_month_rare(nn-1) < 1
        Fcarb_model_month_rare(nn) = -k_diss*(1 - Omega_ar_model_month_rare(nn-1))^n_diss/...
            waterdensity_month_rare(nn)/waterdepth; %{umol/kg/hr}
    end
    Fgas_model_month_rare(nn) = kCO2_month_rare(nn)*K0_month_rare(nn)*(pCO2_model_month_rare(nn-1) - pCO2_atm)/waterdensity_month_rare(nn)/...
        waterdepth; %{umol/kg/hr}
    DIC_model_month_rare(nn) = DIC_model_month_rare(nn-1) - (photo_month_rare(nn) + ...
        (Fcarb_model_month_rare(nn) + Fgas_model_month_rare(nn)))*delta_t_model_month_rare;
    Alk_model_month_rare(nn) = Alk_model_month_rare(nn-1) - 2*Fcarb_model_month_rare(nn)*delta_t_model_month_rare;
    [Alk_out,DIC_out,pH_out,Omega_ar_out,pCO2_out] = ...
        phreeqc_diurnalengine_Alk_DIC(tempC_month_rare(nn),evapdata_CO2degas.Ca_out(1), ...
        evapdata_CO2degas.Mg_out(1),evapdata_CO2degas.K_out(1), ...
        evapdata_CO2degas.SO4_out(1),evapdata_CO2degas.Na_out(1), ...
        evapdata_CO2degas.Cl_out(1),Alk_model_month_rare(nn)/1000,DIC_model_month_rare(nn)/1000);
    pCO2_model_month_rare(nn) = pCO2_out;
    pH_model_month_rare(nn) = pH_out; %note that this pH is on NBS scale
    Omega_ar_model_month_rare(nn) = Omega_ar_out;
    if Fcarb_model_month_rare(nn) >= 0
        if Fgas_model_month_rare(nn) >= 0 && photo_month_rare(nn) >= 0
            d13C_DIC_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1)*DIC_model_month_rare(nn-1) - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_ar)*Fcarb_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_g)*Fgas_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_org)*photo_month_rare(nn)*delta_t_model_month_rare)/...
                DIC_model_month_rare(nn);
            d13C_org_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1) + eps_DIC_org)*...
                photo_month_rare(nn)*delta_t_model_month_rare;
        elseif Fgas_model_month_rare(nn) >= 0 && photo_month_rare(nn) < 0
            d13C_DIC_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1)*DIC_model_month_rare(nn-1) - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_ar)*Fcarb_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_g)*Fgas_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                d13C_org_mean*photo_month_rare(nn)*delta_t_model_month_rare)/DIC_model_month_rare(nn);
        elseif Fgas_model_month_rare(nn) < 0 && photo_month_rare(nn) >= 0
            d13C_DIC_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1)*DIC_model_month_rare(nn-1) - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_ar)*Fcarb_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_g_DIC)*Fgas_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_org)*photo_month_rare(nn)*delta_t_model_month_rare)/...
                DIC_model_month_rare(nn);
            d13C_org_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1) + eps_DIC_org)*...
                photo_month_rare(nn)*delta_t_model_month_rare;
        elseif Fgas_model_month_rare(nn) < 0 && photo_month_rare(nn) < 0
            d13C_DIC_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1)*DIC_model_month_rare(nn-1) - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_ar)*Fcarb_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_g_DIC)*Fgas_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                d13C_org_mean*photo_month_rare(nn)*delta_t_model_month_rare)/DIC_model_month_rare(nn);
        end
     elseif Fcarb_model_month_rare(nn) < 0
        if Fgas_model_month_rare(nn) >= 0 && photo_month_rare(nn) >= 0
            d13C_DIC_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1)*DIC_model_month_rare(nn-1) - ...
                d13C_carb_mean*Fcarb_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_g)*Fgas_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_org)*photo_month_rare(nn)*delta_t_model_month_rare)/...
                DIC_model_month_rare(nn);
            d13C_org_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1) + eps_DIC_org)*...
                photo_month_rare(nn)*delta_t_model_month_rare;
        elseif Fgas_model_month_rare(nn) >= 0 && photo_month_rare(nn) < 0
            d13C_DIC_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1)*DIC_model_month_rare(nn-1) - ...
                d13C_carb_mean*Fcarb_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_g)*Fgas_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                d13C_org_mean*photo_month_rare(nn)*delta_t_model_month_rare)/DIC_model_month_rare(nn);
        elseif Fgas_model_month_rare(nn) < 0 && photo_month_rare(nn) >= 0
            d13C_DIC_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1)*DIC_model_month_rare(nn-1) - ...
                d13C_carb_mean*Fcarb_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_g_DIC)*Fgas_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_DIC_org)*photo_month_rare(nn)*delta_t_model_month_rare)/...
                DIC_model_month_rare(nn);
            d13C_org_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1) + eps_DIC_org)*...
                photo_month_rare(nn)*delta_t_model_month_rare;
        elseif Fgas_model_month_rare(nn) < 0 && photo_month_rare(nn) < 0
            d13C_DIC_model_month_rare(nn) = (d13C_DIC_model_month_rare(nn-1)*DIC_model_month_rare(nn-1) - ...
                d13C_carb_mean*Fcarb_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                (d13C_DIC_model_month_rare(nn-1) + eps_g_DIC)*Fgas_model_month_rare(nn)...
                *delta_t_model_month_rare - ...
                d13C_org_mean*photo_month_rare(nn)*delta_t_model_month_rare)/DIC_model_month_rare(nn);
        end
    end

    %calculate rates
    %Burton and Walter (1987)
    if Omega_ar_model_month_rare(nn) >= 1
        R_BR = k_rare(nn)*(Omega_ar_model_month_rare(nn) - 1).^n_rare(nn); %{umol/m^2/hr}
    elseif Omega_ar_model_month_rare(nn) < 1
        R_BR = -k_diss*(1 - Omega_ar_model_month_rare(nn))^n_diss; %{umol/m^2/hr}
    end
    CR_model_month_rare(nn) = CR_model_month_rare(nn-1) + R_BR*delta_t_model_month_rare;

    %calculate instantaneous d13C_carb
    d13C_carb_inst = d13C_DIC_model_month_rare(nn) + eps_DIC_ar;
%calculate cumulative d13C_carb
    d13C_carb_cumulative_month_rare(nn) = (d13C_carb_cumulative_month_rare(nn-1)*CR_model_month_rare(nn-1) +...
        d13C_carb_inst*R_BR*delta_t_model_month_rare)/CR_model_month_rare(nn);

end

%plot models
figure
tiledlayout(3,2)
colors = parula(3);

nexttile
plot(t_month/24,photo_month,'Color',colors(1,:))
hold on
plot(t_month/24,photo_month_rare,'Color',colors(2,:))
xlabel('day')
ylabel('F_{photo} (\mumol/kg/hr)')

nexttile
plot(t_month/24,(Omega_ar_model_month_rare-Omega_ar_model_month_norm),'k')
ylim([-2 5])
xlabel('day')
ylabel('\Delta \Omega_{arag}')

nexttile
plot(t_month/24,(DIC_model_month_rare/1000-DIC_model_month_norm/1000),'k')
xlabel('day')
ylabel('\Delta [DIC] (mmol/kg)')

nexttile
plot(t_month/24,(d13C_DIC_model_month_rare-d13C_DIC_model_month_norm),'k')
xlabel('day')
ylabel('\Delta\delta^{13}C_{DIC} (‰)')

nexttile
plot(t_month/24,(CR_model_month_rare-CR_model_month_norm),'k')
xlabel('day')
ylabel('\Delta precipitation (\mumol/m^2)')

nexttile
plot(t_month/24,(d13C_carb_cumulative_month_rare-d13C_carb_cumulative_month),'k')
ylim([-0.1 0.5])
xlabel('day')
ylabel('\Delta cumulative \delta^{13}C_{carb} (‰)')

