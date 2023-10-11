%% Part 1: load in data
clear

%load data from dataset spreadsheets
LAC23_P_data = readtable('DatasetS1_platformdata.xlsx');
LAC23_M_data = readtable('DatasetS2_matdata.xlsx');

%create time of day duration vectors
time = timeofday(LAC23_P_data.time);
timeM = timeofday(LAC23_M_data.time);

%% Part 2: plot raw Alk vs DIC data from platform
figure
plot(LAC23_P_data.DIC_mmol_kg_,LAC23_P_data.Alk_mequiv_kg_,'k')
hold on
scatter(LAC23_P_data.DIC_mmol_kg_,LAC23_P_data.Alk_mequiv_kg_,[],...
    hours(time),'filled')
box on
axis equal
ylabel('TA (mequiv/kg)')
xlabel('[DIC] (mmol/kg)')
caxis([0 24])
colorbar

%% Part 3: calculate max and min DIC and gas exchange
%find maximum and minimum DIC values and associated timepoints
[DIC_max,i_max] = max(LAC23_P_data.DIC_mmol_kg_);
t_max = time(i_max);
[DIC_min,i_min] = min(LAC23_P_data.DIC_mmol_kg_);
t_min = time(i_min);
deltaDIC = DIC_max - DIC_min; %{mmol/kg}

%unit transformation for deltaDIC
deltaDIC = deltaDIC*1000; %{umol/kg)

%find pCO2 values at max/min DIC timepoints
pCO2_DICmax = LAC23_P_data.pCO2_uatm_(i_max);
pCO2_DICmin = LAC23_P_data.pCO2_uatm_(i_min);
pCO2_atm = 420; %{uatm}
delta_pCO2 = zeros(length(LAC23_P_data.DIC_mmol_kg_),1);
for n = 1:length(LAC23_P_data.DIC_mmol_kg_)
    delta_pCO2(n) = LAC23_P_data.pCO2_uatm_(n) - pCO2_atm; %{uatm}
end

%unit transformation for delta_pCO2
delta_pCO2 = delta_pCO2*10^-6; %{atm}

%gas exchange calculations
tempC = 25; %{°C}
tempK = tempC + 273; %{K}

K0 = exp(-58.0931 + 90.5697*(100/tempK) + 22.294*log(tempK/100) + ...
    35*(0.027766 - 0.025888*(tempK/100) + 0.005078*(tempK/100)^2)); %{mol/L/atm}
K0 = K0*10^3; %{mol/m^3/atm}

Sc = 2116.8 - 136.25*tempC + 4.7353*tempC^2 - 0.092307*tempC^3 + ...
    0.0007555*tempC^4; %{dimensionless}

u = 8; %{m/s}

kCO2 = 0.251*u^2*(Sc/660)^0.5; %{cm/hr}
kCO2 = kCO2/100; %{m/hr}

Fgas = kCO2*K0*delta_pCO2; %{mol/m^2/hr}

waterdepth = 1.5; %{m}
waterdensity = seawaterdensity(tempC,35); %{kg/m^3}

Fgas = Fgas/waterdensity/waterdepth; %{mol/kg/hr}
Fgas = Fgas*10^6; %{umol/kg/hr}

%calculate gas exchange for each time step, accounting for the longer
%duration of some steps where there is missing data
deltaCO2_vector = zeros(length(Fgas),1);
for p = 1:length(LAC23_P_data.pCO2_uatm_)
    if p<length(LAC23_P_data.pCO2_uatm_)
        delta_t = hours(time(p+1)-time(p));
    else
        delta_t = hours(time(1)-time(p));
    end
    if delta_t < 0
        delta_t = delta_t + 24;
    end
    deltaCO2_vector(p) = Fgas(p)*delta_t;
end

%calculate net gas exchange
deltaCO2_net = sum(deltaCO2_vector(1:i_min)) +...
    sum(deltaCO2_vector(i_max:end));

%calculate alkalinity change
Alk_DICmax = LAC23_P_data.Alk_mequiv_kg_(i_max);
Alk_DICmin = LAC23_P_data.Alk_mequiv_kg_(i_min);
delta_Alk = Alk_DICmax - Alk_DICmin; %{mmol/kg}
delta_Alk = delta_Alk*1000; %{umol/kg}

%calculate kappa_p; this definition is slightly different than what's used
%in Geyman and Maloof (2019) because we assume no alkalinity change
%associated with photosynthesis
kappa_p = deltaDIC + deltaCO2_net - 0.5*delta_Alk; %{umol/kg}

%% Part 4: set up and run engine to compare with data
%define shape of forcing as a sine function sinece the equation isn't given
%in Geyman and Maloof (2019)
t_hr = 0:0.1:24;
time_DT = datetime(2023,7,28,0,0,0) + hours(t_hr);
period = 24; %{hr}
offset = 5; %{hr} when the sine curve will cross 0
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

%set up carbonate chemistry variables
DIC_model = zeros(length(t_hr),1);
DIC_model(1) = mean(LAC23_P_data.DIC_mmol_kg_)*1000; %{umol/kg}
Alk_model = zeros(length(t_hr),1);
Alk_model(1) = mean(LAC23_P_data.Alk_mequiv_kg_)*1000; %{umol/kg}
pCO2_model = zeros(length(t_hr),1);
pCO2_model(1) = mean(LAC23_P_data.pCO2_uatm_); %{uatm}
pH_model = zeros(length(t_hr),1);
pH_model(1) = mean(LAC23_P_data.pH); %{dimensionless}
Omega_ar_model = zeros(length(t_hr),1);
Omega_ar_model(1) = mean(LAC23_P_data.Omega_ar_CO2SYS_);
Fcarb_model = zeros(length(t_hr),1);
Fcarb_model(1) = k_rate*(Omega_ar_model(1) - 1)^n_BR*10^6/waterdensity/waterdepth; %{umol/kg/hr}
Fgas_model = zeros(length(t_hr),1);
Fgas_model(1) = kCO2*K0*(pCO2_model(1) - pCO2_atm)/waterdensity/waterdepth; %{umol/kg/hr}
d13C_DIC_model = zeros(length(t_hr),1);
d13C_DIC_model(1) = mean(LAC23_P_data.d13C_DIC_permil_);
d13C_org_model = -8*ones(length(t_hr),1);
d13C_org_mean = -8; %{permil}

%run diurnal engine
for nn = 2:length(t_hr)
    delta_t_model = t_hr(nn) - t_hr(nn-1); %{hr}
    Fcarb_model(nn) = k_rate*(Omega_ar_model(nn-1) -1)^n_BR*10^6/waterdensity/waterdepth; %{umol/kg/hr}
    Fgas_model(nn) = kCO2*K0*(pCO2_model(nn-1) - pCO2_atm)/waterdensity/waterdepth; %{umol/kg/hr}
    DIC_model(nn) = DIC_model(nn-1) - (photo(nn) + ...
        (Fcarb_model(nn) + Fgas_model(nn)))*delta_t_model;
    Alk_model(nn) = Alk_model(nn-1) - 2*Fcarb_model(nn)*delta_t_model;
    recalcCO2SYS = CO2SYS(DIC_model(nn),Alk_model(nn),2,1,35,tempC,nan,0,nan,...
    0,0,0,0,1,4,1,1,1);
    pCO2_model(nn) = recalcCO2SYS(4);
    pH_model(nn) = recalcCO2SYS(3);
    Omega_ar_model(nn) = recalcCO2SYS(18);
    if Fgas_model(nn) >= 0 && photo(nn) >= 0
        d13C_DIC_model(nn) = (d13C_DIC_model(nn-1)*DIC_model(nn-1) - ...
            (d13C_DIC_model(nn-1) + eps_DIC_ar)*Fcarb_model(nn)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1) + eps_DIC_g)*Fgas_model(nn)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1) + eps_DIC_org)*photo(nn)*delta_t_model)/...
            DIC_model(nn);
        d13C_org_model(nn) = (d13C_DIC_model(nn-1) + eps_DIC_org)*...
            photo(nn)*delta_t_model;
    elseif Fgas_model(nn) >= 0 && photo(nn) < 0
        d13C_DIC_model(nn) = (d13C_DIC_model(nn-1)*DIC_model(nn-1) - ...
            (d13C_DIC_model(nn-1) + eps_DIC_ar)*Fcarb_model(nn)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1) + eps_DIC_g)*Fgas_model(nn)...
            *delta_t_model - ...
            d13C_org_mean*photo(nn)*delta_t_model)/DIC_model(nn);
    elseif Fgas_model(nn) < 0 && photo(nn) >= 0
        d13C_DIC_model(nn) = (d13C_DIC_model(nn-1)*DIC_model(nn-1) - ...
            (d13C_DIC_model(nn-1) + eps_DIC_ar)*Fcarb_model(nn)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1) + eps_g_DIC)*Fgas_model(nn)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1) + eps_DIC_org)*photo(nn)*delta_t_model)/...
            DIC_model(nn);
        d13C_org_model(nn) = (d13C_DIC_model(nn-1) + eps_DIC_org)*...
            photo(nn)*delta_t_model;
    elseif Fgas_model(nn) < 0 && photo(nn) < 0
        d13C_DIC_model(nn) = (d13C_DIC_model(nn-1)*DIC_model(nn-1) - ...
            (d13C_DIC_model(nn-1) + eps_DIC_ar)*Fcarb_model(nn)...
            *delta_t_model - ...
            (d13C_DIC_model(nn-1) + eps_g_DIC)*Fgas_model(nn)...
            *delta_t_model - ...
            d13C_org_mean*photo(nn)*delta_t_model)/DIC_model(nn);
    end
end

%% Part 5: Make comparison plots

%load in Trower et al. (2018) data to plot comparisons
Trower2018data = readtable('Trower 2018 data.xlsx');

%add tide times
lowtide1 = datetime('2023-07-28 10:09:00','InputFormat', ...
    'yyyy-MM-dd HH:mm:ss');
hightide1 = datetime('2023-07-28 17:22:00','InputFormat', ...
    'yyyy-MM-dd HH:mm:ss');
lowtide2 = datetime('2023-07-28 23:35:00','InputFormat', ...
    'yyyy-MM-dd HH:mm:ss');
hightide2 = datetime('2023-07-29 04:50:00','InputFormat', ...
    'yyyy-MM-dd HH:mm:ss');
lowtide3 = datetime('2023-07-29 11:03:00','InputFormat', ...
    'yyyy-MM-dd HH:mm:ss');

%add civil dawn and dusk times
dawn = datetime('2023-07-29 05:52:00','InputFormat', ...
    'yyyy-MM-dd HH:mm:ss');
dusk = datetime('2023-07-28 19:53:00','InputFormat', ...
    'yyyy-MM-dd HH:mm:ss');

%make figure with data and model fits
figure
tiledlayout(3,2)
nexttile
plot(t_hr,DIC_model/1000,'--k')
hold on
scatter(hours(time),LAC23_P_data.DIC_mmol_kg_,'.k')
scatter(hours(timeM),LAC23_M_data.DIC_mmol_kg_,'.r')
yline(mean(LAC23_P_data.DIC_mmol_kg_),'k')
yline(mean(Trower2018data.DIC_mmol_kg_),'b')
xline(hours(timeofday(hightide1)),'-.k')
xline(hours(timeofday(hightide2)),'-.k')
xline(hours(timeofday(dawn)),'--g')
xline(hours(timeofday(dusk)),'--g')
ylim([1.9 2.2])
xlabel('hour of day')
ylabel('[DIC] (mmol/kg)')
xlim([0 24])

nexttile
plot(t_hr,Alk_model/1000,'--k')
hold on
scatter(hours(time),LAC23_P_data.Alk_mequiv_kg_,'.k')
scatter(hours(timeM),LAC23_M_data.Alk_mequiv_kg_,'.r')
yline(mean(LAC23_P_data.Alk_mequiv_kg_),'k')
yline(mean(Trower2018data.Alk_mmol_kg_),'b')
xlim([0 24])
ylim([2.4 2.55])
xlabel('hour of day')
ylabel('Alk (mmol/kg)')
xline(hours(timeofday(lowtide1)),'--k')
xline(hours(timeofday(hightide1)),'-.k')
xline(hours(timeofday(lowtide2)),'--k')
xline(hours(timeofday(hightide2)),'-.k')
xline(hours(timeofday(lowtide3)),'--k')
xline(hours(timeofday(dawn)),'--g')
xline(hours(timeofday(dusk)),'--g')

nexttile
plot(t_hr,pH_model,'--k')
hold on
scatter(hours(time),LAC23_P_data.pH,'.k')
scatter(hours(timeM),LAC23_M_data.pH,'.r')
yline(mean(LAC23_P_data.pH),'k')
yline(mean(Trower2018data.pH),'b')
xline(hours(timeofday(hightide1)),'-.k')
xline(hours(timeofday(hightide2)),'-.k')
xline(hours(timeofday(dawn)),'--g')
xline(hours(timeofday(dusk)),'--g')
xlim([0 24])
ylim([7.9 8.3])
xlabel('hour of day')
ylabel('pH')

nexttile
plot(t_hr,pCO2_model,'--k')
hold on
scatter(hours(time),LAC23_P_data.pCO2_uatm_,'.k')
scatter(hours(timeM),LAC23_M_data.pCO2_uatm_,'.r')
yline(mean(LAC23_P_data.pCO2_uatm_),'k')
[Trower2018_co2sys,Trower2018_co2sys_headers] = CO2SYS(Trower2018data.pH,...
    Trower2018data.Alk_mmol_kg_*1000,3,1,Trower2018data.salinity_ppt_,25,nan,0,nan,...
    0,0,0,0,1,4,1,1,1);
yline(mean(Trower2018_co2sys(:,4)),'b')
xline(hours(timeofday(hightide1)),'-.k')
xline(hours(timeofday(hightide2)),'-.k')
xline(hours(timeofday(dawn)),'--g')
xline(hours(timeofday(dusk)),'--g')
xlim([0 24])
ylim([200 550])
xlabel('hour of day')
ylabel('pCO_2 (\muatm)')

nexttile
plot(t_hr,Omega_ar_model,'--k')
hold on
scatter(hours(time),LAC23_P_data.Omega_ar_CO2SYS_,'.k')
scatter(hours(timeM),LAC23_M_data.Omega_ar_CO2SYS_,'.r')
yline(mean(LAC23_P_data.Omega_ar_CO2SYS_),'--')
yline(mean(Trower2018_co2sys(:,18)),'b')
xline(hours(timeofday(hightide1)),'-.k')
xline(hours(timeofday(hightide2)),'-.k')
xline(hours(timeofday(dawn)),'--g')
xline(hours(timeofday(dusk)),'--g')
xlim([0 24])
ylim([3 5.5])
xlabel('hour of day')
ylabel('\Omega_{ar}')

nexttile
plot(t_hr,d13C_DIC_model,'--k')
hold on
scatter(hours(time),LAC23_P_data.d13C_DIC_permil_,'.k')
scatter(hours(timeM),LAC23_M_data.d13C_DIC_permil_,'.r')
yline(mean(LAC23_P_data.d13C_DIC_permil_),'--')
yline(mean(Trower2018data.d13C_DIC),'k')
xline(hours(timeofday(hightide1)),'-.k')
xline(hours(timeofday(hightide2)),'-.k')
xline(hours(timeofday(dawn)),'--g')
xline(hours(timeofday(dusk)),'--g')
xlim([0 24])
ylim([0.1 1.3])
xlabel('hour of day')
ylabel('\delta^{13}C_{DIC}')

%plot model drivers
figure
tiledlayout(3,1)
nexttile
plot(t_hr,photo)
box on
ylabel('F_{photo} (\mumol/kg/hr)')
xlabel('hour of day')
yline(0)
xlim([0 24])

nexttile
plot(t_hr,Fcarb_model)
box on
ylabel('F_{carb} (\mumol/kg/hr)')
xlabel('hour of day')
yline(0)
xlim([0 24])
ylim([0 0.25])

nexttile
plot(t_hr,Fgas_model)
box on
ylabel('F_{gas} (\mumol/kg/hr)')
xlabel('hour of day')
yline(0)
hold on
xlim([0 24])


%% Part 6: estimate d13C_carb

%sort data so it's in order of time of day
[time_sort,sort_order] = sort(time);
LAC23_P_data_sort = LAC23_P_data(sort_order,:);

%calculate rates
%Burton and Walter (1987)
R_BR = k_BR*(LAC23_P_data_sort.Omega_ar_CO2SYS_ - 1).^n_BR; %{umol/m^2/hr}
%Geyman and Maloof (2019) rate constant
R_GM = k_rate*(LAC23_P_data_sort.Omega_ar_CO2SYS_ - 1).^n_BR; %{mol/m^2/hr}
R_GM = R_GM*10^6; %{umol/m^2/hr}

%calculate cumulative precipitation for both rate laws
CR_BR = zeros(length(R_BR),1);
CR_BR(1) = R_BR(1);
CR_GM = zeros(length(R_GM),1);
CR_GM(1) = R_GM(1);
for n = 2:length(R_BR)
    delta_t_precip = hours(time_sort(n)-time_sort(n-1));
    CR_BR(n) = CR_BR(n-1) + R_BR(n)*delta_t_precip;
    CR_GM(n) = CR_GM(n-1) + R_GM(n)*delta_t_precip;
end

%calculate instantaneous d13C_carb
d13C_carb_inst = LAC23_P_data_sort.d13C_DIC_permil_ + eps_DIC_ar;
d13C_carb_inst_Suess = d13C_carb_inst + 0.8;
d13C_carb_inst_Suess2 = d13C_carb_inst + 1.8;

%calculate cumulative d13C_carb
d13C_cumulative_BR = zeros(length(CR_BR),1);
d13C_cumulative_BR(1) = d13C_carb_inst(1);
d13C_cumulative_GM = zeros(length(CR_BR),1);
d13C_cumulative_GM(1) = d13C_carb_inst(1);
d13C_cumulative_Suess = zeros(length(CR_BR),1);
d13C_cumulative_Suess(1) = d13C_carb_inst_Suess(1);
d13C_cumulative_Suess2 = zeros(length(CR_BR),1);
d13C_cumulative_Suess2(1) = d13C_carb_inst_Suess2(1);
for m = 2:length(CR_BR)
    delta_t_d13C = hours(time_sort(m)-time_sort(m-1));
    d13C_cumulative_BR(m) = (d13C_cumulative_BR(m-1)*CR_BR(m-1) + ...
        d13C_carb_inst(m)*R_BR(m)*delta_t_d13C)/CR_BR(m);
    d13C_cumulative_GM(m) = (d13C_cumulative_GM(m-1)*CR_GM(m-1) + ...
        d13C_carb_inst(m)*R_GM(m)*delta_t_d13C)/CR_GM(m);
    d13C_cumulative_Suess(m) = (d13C_cumulative_Suess(m-1)*CR_BR(m-1) + ...
        d13C_carb_inst_Suess(m)*R_BR(m)*delta_t_d13C)/CR_BR(m);
    d13C_cumulative_Suess2(m) = (d13C_cumulative_Suess2(m-1)*CR_BR(m-1) + ...
        d13C_carb_inst_Suess2(m)*R_BR(m)*delta_t_d13C)/CR_BR(m);
end

%make figure comparing cumulative d13C_carb with different rate laws; this
%figure shows that it doesn't matter which rate constant we choose for the
%cumulative calculation
figure
tiledlayout(2,1)
nexttile
scatter(time_sort,d13C_cumulative_BR,'.k')
box on
xlabel('time of day')
ylabel('cumulative \delta^{13}C (‰)')
title('Burton and Walter (1987) rate constant')
yline(3.46)

nexttile
scatter(time_sort,d13C_cumulative_GM,'.k')
box on
xlabel('time of day')
ylabel('cumulative \delta^{13}C (‰)')
title('Geyman and Maloof (2019) rate constant')
yline(3.46)

%make figure that compares d13C_carb estimates with and without the Suess
%effect
figure
tiledlayout(3,1)
nexttile
yyaxis left
scatter(time_sort,R_BR,'.')
box on
ylim([100 350])
ylabel('R_a_r (\mumol/m^2/hr)')
yyaxis right
scatter(time_sort,CR_BR,'.')
ylabel('cumulative precip (\mumol/m^2)')
xticks(timeofday(datetime(zeros(6,1),zeros(6,1),zeros(6,1),(0:4:20)', ...
    zeros(6,1),zeros(6,1))))

nexttile
scatter(time_sort,LAC23_P_data_sort.d13C_DIC_permil_,'.k')
hold on
scatter(time_sort,LAC23_P_data_sort.d13C_DIC_permil_ + 0.8, ...
    '.','MarkerEdgeColor',[0.3010 0.7450 0.9330])
d13C_Suess2 = LAC23_P_data_sort.d13C_DIC_permil_ + 1.8;
scatter(time_sort,d13C_Suess2, '.','MarkerEdgeColor',[0 0.4470 0.7410])
box on
ylabel('\delta^{13}C_{DIC}')
ylim([0 3.5])
legend('measured values','Suess 1',"Suess 2",'location','southoutside')
xticks(timeofday(datetime(zeros(6,1),zeros(6,1),zeros(6,1),(0:4:20)', ...
    zeros(6,1),zeros(6,1))))

nexttile
scatter(time_sort,d13C_cumulative_BR,'.k')
hold on
scatter(time_sort,d13C_cumulative_Suess,'.','MarkerEdgeColor', ...
    [0.3010 0.7450 0.9330])
scatter(time_sort,d13C_cumulative_Suess2,'.','MarkerEdgeColor', ...
    [0 0.4470 0.7410])
box on
ylabel('cumulative \delta^{13}C_{carb}')
ylim([3 5.5])
xticks(timeofday(datetime(zeros(6,1),zeros(6,1),zeros(6,1),(0:4:20)', ...
    zeros(6,1),zeros(6,1))))
legend('measured values','Suess 1','Suess 2','location','southoutside')

%% Part 7: load and plot Suess effect data

BATS = readtable('BATS.csv','Range','A36:J508');

d13C_num = BATS.Var8;
d13C_num(isnan(d13C_num)) = [];
date_d13C = BATS.Date;
date_d13C(isnan(BATS.Var8)) = [];
d13C_mm = movmean(d13C_num,5);

% Load Caribbean sclerosponge data
Bahamas_Swart2002 = readtable('bahamas-sclerosponge-2002.txt',...
    NumHeaderLines=79);
Jamaica_Bohm2002_1 = readtable('jamaica_ce95-2-1.txt',...
    NumHeaderLines=49);
Jamaica_Bohm2002_2 = readtable('jamaica_ce95-2-2.txt',...
    NumHeaderLines=49);
Jamaica_Bohm2002_3 = readtable('jamaica_ce96-1.txt',...
    NumHeaderLines=49);
Bahamas_Lazareth2000 = readtable('lazareth2000-d13c_noaa.txt',...
    NumHeaderLines=90);

% Plot data
figure
tiledlayout(3,1)
nexttile
colors = parula(6);
scatter(Bahamas_Swart2002.Year,Bahamas_Swart2002.d13C,[],colors(1,:),'.')
hold on
scatter(Jamaica_Bohm2002_1.YearAD,Jamaica_Bohm2002_1.d13C,[],colors(2,:),'.')
scatter(Jamaica_Bohm2002_2.YearAD,Jamaica_Bohm2002_2.d13C,[],colors(3,:),'.')
scatter(Jamaica_Bohm2002_3.YearAD,Jamaica_Bohm2002_3.d13C,[],colors(4,:),'.')
scatter(Bahamas_Lazareth2000.age,Bahamas_Lazareth2000.d13C,[],colors(5,:),'.')
all_date = cat(1,Jamaica_Bohm2002_1.YearAD,...
    Jamaica_Bohm2002_2.YearAD,Jamaica_Bohm2002_3.YearAD,...
    Bahamas_Lazareth2000.age);
[all_date_sort,sortorder] = sort(all_date);
all_d13C = cat(1,Jamaica_Bohm2002_1.d13C,...
    Jamaica_Bohm2002_2.d13C,Jamaica_Bohm2002_3.d13C,...
    Bahamas_Lazareth2000.d13C);
all_d13C_sort = all_d13C(sortorder);
all_mm = movmean(all_d13C_sort,10);
plot(all_date_sort,all_mm,'r')
box on
ylim([3.5 5])
legend('Bahamas Swart 2002',...
    'Jamaica Bohm 2002 1','Jamaica Bohm 2002 2','Jamaic Bohm 2002 3', ...
    'Jamaica Lazareth 2000','location','southwest')
xlim([1800 2000])
ylabel('\delta^{13}C (‰)')

nexttile
scatter(date_d13C,d13C_num,'.k')
hold on
plot(date_d13C,d13C_mm,'k')
BATS_mdl = fitlm(datenum(date_d13C),d13C_num);
test_dates = datetime(1989:2023,7*ones(1,35),28*ones(1,35));
BATS_mdl_y = predict(BATS_mdl,datenum(test_dates'));
plot(test_dates,BATS_mdl_y,'r')
box on
ylabel('\delta^{13}C_{DIC} (‰)')
ylim([0.5 1.6])
xlim([min(date_d13C) max(date_d13C)])

nexttile
plot(test_dates(test_dates<max(date_d13C)), ...
    (BATS_mdl_y(test_dates<max(date_d13C)) - min(BATS_mdl_y)))
hold on
plot(test_dates(test_dates>max(date_d13C)), ...
    (BATS_mdl_y(test_dates>max(date_d13C)) - min(BATS_mdl_y)),'--r')
plot(datetime(all_date_sort,ones(length(all_date_sort),1), ...
    ones(length(all_date_sort),1)),all_mm - 3.1,'r')
box on
xlim([datetime(1800,1,1) datetime(2023,7,28)])
ylabel('\delta^{13}C_{DIC} (‰) relative to modern')
