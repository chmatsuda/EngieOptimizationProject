%% Load data and plot
clear all;
load('ExcelImported2.mat');
%%
t = [1:8760]';
ElecDemand = table2array(ElecDemand);%kW


%% Can we use a battery for a 2MW + 4.5 MW max PV distribution center if the peak is 3.9MW

RenewableEnergy = table2array(RenewableEnergy);
PVInstalled = 4000; %How many kW of PV are installed, adjust this accordingly
WindInstalled = 0; %3.6MW of Wind installed = 3600kW 

GridProduction = ones(8760,1)*2000; %kW
PV = RenewableEnergy(:,1)*PVInstalled; %Kw
Wind = RenewableEnergy(:,2)*WindInstalled; %Kw
%TotalRen = PV+Wind; %kW

%% Plotting the production
figure(1)
plot(t,ElecDemand,'b',t,PV,'g',t,GridProduction,'r');
legend("Demand","PV","Grid");
ylabel('Power(kW)')


%% Compute the energy that is consumed over the time period and the energy produced:
% E = P * t
%Assuming that each time step is exactly 60 seconds:
energy_consumed = ElecDemand* 3600;%% kW converted to kJ
total_energy_cons = sum(energy_consumed)%% kJ
PV_energy_prod = PV *3600; %%kJ
Grid_energy_prod = GridProduction *3600;
total_ren_prod = sum(PV_energy_prod) + sum(Grid_energy_prod)%%kJ
%The energy consumed is 8 E10 KJ and the energy produced by wind and 2MW grid is 8.08E10 so hypothetically yes 


% %% What is the percentage of time during which the community is self-sufficient in terms of electricity?
% % At each time step, check if the demand is less than the production. Then
% % it is self sufficient
% sstime = 0;
% PowerStillNeeded = zeros(size(timesteps));
% for i = 1:305280
%     if Demand(i) < Wind(i)
%        sstime = sstime+timesteps(i); 
%     else
%         PowerStillNeeded(i) = Demand(i)- Wind(i);
%     end
%     
% end
% 
% totaltime = sum(timesteps);
% 
% percentage = sstime/totaltime
% % = 54.38% of the time, it is self sufficient, 45.62% of the time it is NOT
% % sufficient
% 
% %% What would the maximum power needed for a battery in the discharge mode? What is the unit?
% 
% %Above, we found the difference in demand and production in the variable
% %PowerStillNeeded.
% 
% %Finding the maximum difference between production and demand,
% MaxDiff = max(PowerStillNeeded)
% %Gives 26.65 MW as the max
% %% Size of Battery in MWh (energy)
% % At each time step, check if the demand is less than the production. Then
% % it is self sufficient
% 
% 
% for k = 1:2000
%     Capacities = linspace(1,2000,2000);
%     sstime2=0; %Time With Power
%     MaxCap = Capacities(k); %MWh
%     BatteryE = MaxCap;
%     count=0;
%     PercentWanted = .85;
%     for j = 1:1:size(timesteps,1) 
%         PowerDiff = Wind(j) - Demand(j);
%         time =timesteps(j);
%           if PowerDiff >= 0
%             sstime2 = sstime2 + timesteps(j);
%         elseif PowerDiff <0
%             EnergyNeeded = abs(PowerDiff) * timesteps(j)/3600; %%MWh
%             if EnergyNeeded < BatteryE
%                 sstime2 = sstime2+ timesteps(j);
%             elseif BatteryE < EnergyNeeded
%                 if BatteryE >0
%                     sstime2 = sstime2 + BatteryE/EnergyNeeded *timesteps(j);
%                 end
%             end
%         end
% 
%         BatteryE = ChargeDrawBat(BatteryE,PowerDiff,MaxCap,time);
%         bat(j) = BatteryE;
%     end
%     
%     percentage2 = sstime2/totaltime;
%     suff(k) = percentage2;
%     caps(k) = Capacities(k);
%     percentage2
%     Capacities(k)
%     if percentage2 >= PercentWanted
%         disp(percentage2)
%         disp(Capacities(k))
%         break;
%     end
%     
% end
% %%
% figure(2)
% plot(suff,caps)
% ylabel('Needed Capacity (MWh)')
% xlabel('Percent self sufficient')
% %% Taking a cost of P400 per kWh, what would be the price to be 100% self sufficient?
% cost = caps*400*1000;
% plot(suff,cost)
% ylabel('Total Cost (Pounds)')
% xlabel('Percent self sufficient')
% %632,000,000 pounds
