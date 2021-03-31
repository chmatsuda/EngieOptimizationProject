%% Load data and plot

clear all;
load('ExcelImported2.mat');
%%
t = [1:8760]';
ElecDemand = table2array(ElecDemand);%kW
RenewableEnergy = table2array(RenewableEnergy);
%% Can we use a battery for a 2MW + 4.5 MW max PV distribution center if the peak is 3.9MW

PVInstalled = 4500; %How many kW of PV are installed, adjust this accordingly
WindInstalled = 0; %3.6MW of Wind installed = 3600kW 

GridProduction = ones(8760,1)*2000; %kW
PV = RenewableEnergy(:,1)*PVInstalled; %Kw
Wind = RenewableEnergy(:,2)*WindInstalled; %Kw
%TotalRen = PV+Wind; %kW
PVnGrid = PV+GridProduction;
%% Plotting the production

figure(1)
plot(t,ElecDemand,'b',t,PV,'g',t,GridProduction,'r');
legend("Demand","PV","Grid");
ylabel('Power(kW)')

figure(2)
plot(t,ElecDemand,'b',t,PVnGrid,'g');
legend("Demand","PV + Grid");
ylabel('Power(kW)')
%% Compute the energy that is consumed over the time period and the energy produced:
% E = P * t

%Assuming that each time step is exactly 60 seconds:
%Each time step is 1 hour = 3600 s
energy_consumed = sum(ElecDemand)*3600/10^6;%% kJ converted to GJ
PV_energy_prod = sum(PV) *3600/10^6; %%GJ
Grid_energy_prod = sum(GridProduction) *3600/10^6; %GJ
total_energy_prod = PV_energy_prod + Grid_energy_prod;

%The energy consumed is 80,209 GJ and energy produced with PV + 2MW grid is
%80,804 GJ so hypothetically could work... for the yearly production/needs.
%Doesn't account for instantaneous needs.
%% What percentage of time is it self sufficient?

%If we don't have a battery
sstime = 0;
PowerStillNeeded = zeros(8760,1);
for i = 1:8760
    if ElecDemand(i) < (PV(i)+GridProduction(i))
       sstime = sstime+1; 
    else
        PowerStillNeeded(i) = ElecDemand(i)- (PV(i)+GridProduction(i));
    end
    
end
percentage = sstime/8760.;
X = sprintf('The system is self-sufficient for %d percent of the time',percentage);
disp(X)
% 51.31% of the time, it is self sufficient
%% Max Power of Battery in Discharge Mode in MW

%Above, we found the difference in demand and production in the variable
%PowerStillNeeded.


%Finding the maximum difference between production and demand,
MaxPower= max(PowerStillNeeded)/1000; %MW
X = sprintf('Maximum battery capacity needed at a given moment: %d MW',MaxPower);
disp(X)
%Gives 1.5725 MW as the max
%% Max Size of Battery in MWh

for k = 1:15000
    Capacities = linspace(1,15000,15000);
    sstime2=0; %Time With Power
    MaxCap = Capacities(k); %MWh
    BatteryE = MaxCap;
    count=0;
    PercentWanted = 1.0;
    for j = 1:8760
        PowerDiff = ((PV(j)+GridProduction(j)) - ElecDemand(j))/1000; %MW
          if PowerDiff >= 0 % production > demand
            sstime2 = sstime2 + 1;
          elseif PowerDiff <0
            EnergyNeeded = abs(PowerDiff); %%MWh
            if EnergyNeeded <= BatteryE
                sstime2 = sstime2+ 1;
            elseif BatteryE < EnergyNeeded
                if BatteryE >0
                    sstime2 = sstime2 + BatteryE/EnergyNeeded;
                end
            end
        end

        BatteryE = ChargeDrawBat(BatteryE,PowerDiff,MaxCap);
        bat(j) = BatteryE;
    end
    
    percentage2 = sstime2/8760;
    suff(k) = percentage2;
    caps(k) = Capacities(k);
    percentage2;
    Capacities(k);
    if percentage2 >= PercentWanted
        X = sprintf('Percentage of time in which demand is covered by supply: %d percent', percentage2*100);
        disp(X)
        CapacityNeeded = Capacities(k);
        X = sprintf('Total battery capacity needed: %d MWh', CapacityNeeded);
        disp(X)
        
        break;
    end
    
    
    
end
%%
figure(3)
plot(suff,caps)
ylabel('Needed Capacity (MWh)')
xlabel('Percent self sufficient')
%% Cost
% Battery CAPEX


LiIonCAPEX = CapacityNeeded*300e3; % 300 euros for per kWh of capacity
fprintf('Battery CAPEX: %d euros',LiIonCAPEX)
%% 
% Cost of purchasing from grid

CommodityPrices = table2array(CommodityPrices);%kW
MarketPrice = CommodityPrices(:,1); %euros/kWh
PurchasedElectricty = zeros(8760);
for i = 1:8760
    CostOn = MarketPrice(i)*2e3 ; %assuming 2MW is continuously supplied
    PurchasedElectricity(i) = CostOn;
end
PurchasedElectricity = sum(PurchasedElectricity);    
X = sprintf('Cost of electricity from grid: %d euros',PurchasedElectricity);
disp(X)
%% 
% Cost of charging battery and profit from selling excess PV production

BatteryChargeHeld = zeros(8760,1);
MaxCap = CapacityNeeded*10^3; % Max capacity of battery set to capacity needed in kWh
BatteryChargeHeld(1) = MaxCap; % Assume we start with a full battery.
BatteryChargeCost = zeros(8760,1);
ExcessSale = zeros(8760,1);
for j = 1:8760
    BatteryChargeHeld(j+1) = BatteryChargeHeld(j)+PV(j)+GridProduction(j)-ElecDemand(j);
    if BatteryChargeHeld(j+1) <= MaxCap
        BatteryChargeCost(j) = PV(j)*20.5; % 20.5 euros per kWh charge 
        ExcessSale(j) = 0; % if battery charge held < max capacity, then no excess PV to be sold
    else
        BatteryChargeCost(j) = (MaxCap - BatteryChargeHeld(j))*20.5; % 20.5 euros per kWh charge 
        ExcessSale(j) = (BatteryChargeHeld(j)+PV(j) - MaxCap)*0.3*MarketPrice(j); % selling at 30% of commodity price
        BatteryChargeHeld(j+1) = MaxCap;
    end
end
BatteryChargeCost = sum(BatteryChargeCost);
X = sprintf('Annual cost of charging battery: %d euros',BatteryChargeCost);
disp(X)
ExcessSale = sum(ExcessSale);
Y = sprintf('Annual money made from selling excess energy produced at 30 percent market price: %d euros', ExcessSale);
disp(Y)
ChargeatYearEnd = BatteryChargeHeld(8760)
%% 
% Total Costs for 2021-2024

TotalCost = LiIonCAPEX + PurchasedElectricity + BatteryChargeCost - ExcessSale
%%
function BatteryEnergy = ChargeDrawBat(BatteryEnergy, PowerDifference, MaxCap)
    EnergyDiff = PowerDifference;
    BatteryEnergy = BatteryEnergy + EnergyDiff; % MWh
        if BatteryEnergy > MaxCap
            BatteryEnergy = MaxCap;
        end
        if BatteryEnergy <0
            BatteryEnergy = 0;
        end
    
end