clear all
close all
%% Load all the data
load('T1_Data_Community.mat'); %We load the data: t = time, Demand = demand from community, Wind = Wind Production from community
%Demand and Wind are assumed to be in MW.
TimeResolution = 1; %(time resolution will be 1 minute interval after cleaning the data)
Tcoeff = TimeResolution/60; %to take into account the energy calculation
%1W for 1/2hour represents 0.5 Wh
Timediff = seconds(diff(t)); % computes the time interval between successive recordings
Timediff = [Timediff;60]; %to take into account the last value that we consider equal to 60 seconds
count = 0; %index used to populate the missing data
SizeVectors = size(t,1);
%% Clean the data
Demand(Demand < 0) = 0; %remove data less than 0
Wind(Wind < 0) = 0; %remove data less than 0

%Remove all data missing data for longer 2 hours
for i= 2:size(t,1)
    
   if Timediff(i-1)> 2*60*60 %if a time interval between 2 consecutive data is above 2 hours, 
       %we discard data within this time interval as the wind could change consequently and change the analysis we are doing
       t =  [t(1:i-1+count,1);t(i-1+count,1)+minutes(1);t(i+count,1)-minutes(1);t(i+count:SizeVectors+count,1)]; %we create a new time vector 
       Demand =  [Demand(1:i-1+count,1);0;0;Demand(i+count:SizeVectors+count,1)];%we create a new Demand vector with 0 at the boundaries of the time interval where we want to discard data
       Wind =  [Wind(1:i-1+count,1);0;0;Wind(i+count:SizeVectors+count,1)]; %we create a new Wind vector with 0 at the boundaries of the time interval where we want to discard data
        count = count+2;
   end
    
end

%Clean the timesteps


NewTime = t(1,1):minutes(1):t(size(t,1),1);
NewTime=NewTime';
NewDemand = interp1(t,Demand,NewTime); % interpolation for all data for which time step is different from 1 minute.
NewWind = interp1(t,Wind,NewTime);
Wind = NewWind;
Demand = NewDemand;
t = NewTime;
Timediff = seconds(diff(t));
Timediff = [Timediff;60]; %to take into account the last value that we consider equal to 60 seconds

%% Compute the current percentage of time without need for grid (Case without batteries)
BCap = zeros(size(Demand));
BCap(1,1)=0; %if the capacity is 0 MWh, what is the downtime in percentage
Pbat = zeros(size(Demand)); %Battery Power
    for i = 2:size(Demand,1)
          BCap(i) = BCap(i-1) +  (Wind(i) - Demand(i))*Tcoeff ;
          if BCap(i)> BCap(1,1)
             BCap(i) = BCap(1,1); 
          end
          if BCap(i)<0
             BCap(i)=0; 
          
          end
          Pbat(i) = (BCap(i) - BCap(i-1))/(Tcoeff);
    end
A = Wind - Pbat - Demand; 
PercentageTimewithoutPower= sum(A<-0.000001)/size(Demand,1)*100;
X = sprintf('If there is no battery, %d %% of the time consumption will be supplied by renewables',(100.-PercentageTimewithoutPower));
disp(X)
% figure('Name','Power Flows without battery');
% plot(t,Wind,t,Demand,t,BCap/BCap(1,1)*max(Wind),t,Pbat,t,A)


%% Compute the battery capacity to get 100% of demand met by Renewables

BCap = 0; %Capacity of the battery
maxBCap = 0; 
maxPower = 0;
signmaxPower = 1;

%We run the data backward (no need to do that, but that is a personal preference) to get the size of the battery so it covers 100%
%of the consumption
for i = 0:size(Demand,1)-1
      BCap = BCap + (Demand(size(Demand,1)-i) - Wind(size(Demand,1)-i))*Tcoeff;
    if BCap <0
        BCap = 0;
    end
    if BCap > maxBCap
        maxBCap = BCap;
    end
    if abs((Demand(size(Demand,1)-i) - Wind(size(Demand,1)-i))) > maxPower
        maxPower = abs((Demand(size(Demand,1)-i) - Wind(size(Demand,1)-i)));
        signmaxPower = ((Demand(size(Demand,1)-i) - Wind(size(Demand,1)-i)))/abs((Demand(size(Demand,1)-i) - Wind(size(Demand,1)-i)));
    end
end


X = sprintf('maximum battery capacity needed to be off-grid 100 %% of the time : %d MWh, with a maximum Power: %d MW ',maxBCap, signmaxPower*maxPower );
disp(X)
% maxBCap %.1fs the maximum battery capacity needed
% maxPower %.1fs the maximum power put in/from the battery without any smart charging (battery takes all the difference between Wind and Demand)
% signmaxPower; %+1 if the max power is obtained for discharge, -1 if it is obtained for discharge

TimeDuration = seconds(t(size(t),1)-t(1,1));
answer = questdlg('Do you want to compute the required battery capacity for a particular sel-consumption percentage, or compute a curve Self-consumption percentage vs battery capacity (takes several minutes)?', ...
	'choice of study',...
    'one battery capacity', 'curve','cancel');

switch answer
    case 'one battery capacity'
%% Compute the required battery capacity to reach the wanted  % of time with self consumption
%         PercentSelfConsumptionTime = input('What percentage of self-consumption time do you want to achieve (enter a number between 55 (for 55%) and 100 (for 100%))?');
        prompt = {'What percentage of self-consumption time do you want to achieve (enter a number between 55 (for 55%) and 100 (for 100%))?'};
        dlgtitle = 'Input User';
        dims = [1 80];
        definput = {'55'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        PercentSelfConsumptionTime = str2num(answer{1,1});
        %         PercentSelfConsumptionTime = 85; %percentage of self consumption Time that one want to achieve
        %Now, we run forward with various capacities for the battery until the
        %self-consumption time is above the target given by PercentSelfConsumptionTime
        BCap = zeros(size(Demand));
        Pbat = zeros(size(Demand));
        PercentageTimewithoutPower = 0;
        

        k = 0;
            while PercentageTimewithoutPower <= (1-PercentSelfConsumptionTime/100)*100
                BCap(1,1)=round(maxBCap+0.5)-k; %We decrease the battery capacity until the demand is met by renewable+battery the right percentage of time
                TimeWithoutPower = 0;
                for i = 2:size(Demand,1)
                      BCap(i) = BCap(i-1) +  (Wind(i) - Demand(i))*Tcoeff ;
                      if BCap(i)> BCap(1,1)
                         BCap(i) = BCap(1,1); 
                      end
                      if BCap(i)<0
                         BCap(i)=0; 
                      end
                      Pbat(i) = (BCap(i) - BCap(i-1))/(Tcoeff);
                      A(i) =Wind(i) - Pbat(i) - Demand(i);   %we alsodiscard the times where we do not have any data
                       
                        TimeWithoutPower = TimeWithoutPower + (A(i)<-0.00001)*Timediff(i-1);
                      
                end
            PercentageTimewithoutPower= TimeWithoutPower/TimeDuration(1,1)*100;
            k = k+1;
            end

        BatteryCapacityNeeded = BCap(1,1);  %Battery Capacity in MWh


        X = sprintf('maximum battery capacity needed to be off-grid %.1f %% of the time : %.1f MWh',100-PercentageTimewithoutPower, BatteryCapacityNeeded);
        disp(X)

        figure('Name','Power Flows with battery');
        hold on
        yyaxis left
        plot(t,Wind,'-','color','#0072BD')
        plot(t,Demand,'-','color','#D95319')
        plot(t,Pbat,'-','color','#EDB120')
        plot(t,A,'-','color','#77AC30')
        yyaxis right 
        plot(t,BCap)
        legend('Wind Production', 'Demand', 'P battery', 'Wind-battery-Demand','SoC')
        hold off

    

end


