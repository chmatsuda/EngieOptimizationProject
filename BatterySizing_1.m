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
energy_consumed = sum(ElecDemand)* 3600/1000/1000%% kW converted to G
PV_energy_prod = sum(PV) *3600/1000/1000 %%GJ
Grid_energy_prod = sum(GridProduction) *3600/1000/1000 %GJ
total_energy_prod = PV_energy_prod + Grid_energy_prod;

%The energy consumed is 80,209 GJ and energy produced with PV + 2MW grid is
%83,020 GJ so hypothetically could work...

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
percentage = sstime/8760
% 51.31 of the time, it is self sufficient

%% Max Power of Battery in Discharge Mode in MW
%Above, we found the difference in demand and production in the variable
%PowerStillNeeded.


%Finding the maximum difference between production and demand,
MaxPower= max(PowerStillNeeded)/1000 %MW
%Gives 1.573 MW as the max

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
            if EnergyNeeded < BatteryE
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
    percentage2
    Capacities(k)
    if percentage2 >= PercentWanted
        disp(percentage2)
        disp(Capacities(k))
        break;
    end
    
    
    
end

%%
figure(3)
plot(suff,caps)
ylabel('Needed Capacity (MWh)')
xlabel('Percent self sufficient')


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


