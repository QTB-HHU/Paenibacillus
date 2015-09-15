%%To create all plots, use these scripts iteratively.
%First perform the analysis using ScrumPy and the Plotting.py script, which will create the necessary Data Files
%First initialize the Experimental data

Colors = [0 0 0; %Biomass
255 192 203; % FHL
180 205 20; % Acetone
139 69 19; %Acetate
255 140 0; % Ethanol
0 0 255; %Butanediol
255 0 255; %Nitrate
255 0 0; %Nitrite
148 0 211]/255; % Ammonia

%Identifiers
Biomass = 1;
FHL = 2;
ACETONE = 3;
Acetate = 4;
Ethanol = 5;
Butanediol = 6;
Nitrite = 7;
Nitrate = 8;
Ammonia = 9;


%% Experimental Data

%Experimental Values in mmol / mmol Glucose
AmmoniaExp = [   0.3033628736	0.0896131711; %Biomass
    0 0;
    0.0527975364	0.0096205995; %Acetone    
    0.0853943177	0.004006848;    %Acetate
0.586059998	0.0445454352; %Ethanol
0.6091781782	0.0136487511]; %Butanediol

NitrateExp = [
    0.2196482761	0.0541129334; %Biomass
    0 0;
    0.0592858547	0.0153267301; %Acetone
    0.4771383713	0.0536121659; %Acetate
    0.0263384881	0.0028124995; %Ethanol
    0.5976356476	0.0337292533]; %Butanediol


%% Figure 3
figure
[data,~,colheaders] = importdata('FinalAmmonium.csv')

for i=1:numel(data.data(:,1))
if data.data(i,2) < AmmoniaExp(Biomass,1)
Maintenance = (data.data(i,1) + data.data(i-1,1))/2;
break
end
end
hold off
p = plot(data.data(:,1),data.data(:,2:10),'LineWidth',2)
set(gca,'YLim',[0,2.5])
hold on

p(Biomass).Color = [0 0 0]/255;
p(FHL).Color = [255 192 203]/255;
p(ACETONE).Color = [180 205 20]/255;
p(Acetate).Color = [139 69 19]/255;
p(Ethanol).Color = [255 140 0]/255;
p(Butanediol).Color  = [0 0 255]/255;
p(Nitrite).Color = [255 0 255]/255;
p(Nitrate).Color = [255 0 0]/255;
p(Ammonia).Color = [148 0 211]/255;
p(ACETONE).LineStyle = '--';
p(Butanediol).LineStyle = '--';
p(Nitrite).LineStyle = '-.';
p(Nitrate).LineStyle = '--';
children = get(gca,'children')
delete(children(10-[ Nitrate Nitrite ]));
legend(data.colheaders([[2:7 10]]))
hold on
errorbar(Maintenance,AmmoniaExp(Biomass,1),AmmoniaExp(Biomass,2),'Color',Colors(Biomass,:),'Marker','o','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(Ethanol,1),AmmoniaExp(Ethanol,2),'Color',Colors(Ethanol,:),'Marker','*','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(Butanediol,1),AmmoniaExp(Butanediol,2),'Color',Colors(Butanediol,:),'Marker','s','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(Acetate,1),AmmoniaExp(Acetate,2),'Color',Colors(Acetate,:),'Marker','d','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(ACETONE,1),AmmoniaExp(ACETONE,2),'Color',Colors(ACETONE,:),'Marker','x','MarkerSize',4);

xlabel('Maintenance ATP per biomass carbon','FontSize',14)
ylabel('Flux activity in mmol/mmol glucose','FontSize',14)

%% Supplementary Figure 6
[data,~,colheaders] = importdata('PREDEHOGFix.csv')
figure

for i=1:numel(data.data(:,1))
if data.data(i,2) < AmmoniaExp(Biomass,1)
Maintenance = (data.data(i,1) + data.data(i-1,1))/2;
break
end
end
hold off
p = plot(data.data(:,1),data.data(:,2:10),'LineWidth',2)
hold on
%Plot the Experimental Values
p(Biomass).Color = [0 0 0]/255;
p(FHL).Color = [255 192 203]/255;
p(ACETONE).Color = [180 205 20]/255;
p(Acetate).Color = [139 69 19]/255;
p(Ethanol).Color = [255 140 0]/255;
p(Butanediol).Color  = [0 0 255]/255;
p(Nitrite).Color = [255 0 255]/255;
p(Nitrate).Color = [255 0 0]/255;
p(Ammonia).Color = [148 0 211]/255;
p(ACETONE).LineStyle = '--';
p(Butanediol).LineStyle = '--';
p(Nitrite).LineStyle = '-.';
p(Nitrate).LineStyle = '--';
children = get(gca,'children')
delete(children(10-[ Nitrate Nitrite ]));
legend(data.colheaders([[2:7 10]]))
hold on
errorbar(Maintenance,AmmoniaExp(Biomass,1),AmmoniaExp(Biomass,2),'Color',Colors(Biomass,:),'Marker','o','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(Ethanol,1),AmmoniaExp(Ethanol,2),'Color',Colors(Ethanol,:),'Marker','*','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(Butanediol,1),AmmoniaExp(Butanediol,2),'Color',Colors(Butanediol,:),'Marker','s','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(Acetate,1),AmmoniaExp(Acetate,2),'Color',Colors(Acetate,:),'Marker','d','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(ACETONE,1),AmmoniaExp(ACETONE,2),'Color',Colors(ACETONE,:),'Marker','x','MarkerSize',4);

xlabel('Maintenance ATP per biomass carbon','FontSize',14)
ylabel('Flux activity in mmol/mmol glucose','FontSize',14)


%% Supplementary Figure 7
[data,~,colheaders] = importdata('CarbonScan.csv')
figure
Bestfit = 100;
BestFitVal = 0;

for i=1:numel(data.data(:,1))
    currentfit = (AmmoniaExp(Butanediol,1) - data.data(i,Butanediol+1))^2 + (AmmoniaExp(Ethanol,1) - data.data(i,Ethanol+1))^2 + (AmmoniaExp(Biomass,1) - data.data(i,Biomass+1))^2
    if currentfit < Bestfit
        BestFitVal = data.data(i,1);
        Bestfit = currentfit;
    end
end

hold off
p = plot(data.data(:,1),data.data(:,2:10),'LineWidth',2)
hold on

p(Biomass).Color = [0 0 0]/255;
p(FHL).Color = [255 192 203]/255;
p(ACETONE).Color = [180 205 20]/255;
p(Acetate).Color = [139 69 19]/255;
p(Ethanol).Color = [255 140 0]/255;
p(Butanediol).Color  = [0 0 255]/255;
p(Nitrite).Color = [255 0 255]/255;
p(Nitrate).Color = [255 0 0]/255;
p(Ammonia).Color = [148 0 211]/255;
p(ACETONE).LineStyle = '--';
p(Butanediol).LineStyle = '--';
p(Nitrite).LineStyle = '-.';
p(Nitrate).LineStyle = '--';
children = get(gca,'children')
delete(children(10-[ Nitrate Nitrite ]));
legend(data.colheaders([[2:7 10]]))
hold on
errorbar(BestFitVal,AmmoniaExp(Biomass,1),AmmoniaExp(Biomass,2),'Color',Colors(Biomass,:),'Marker','o','MarkerSize',4);
errorbar(BestFitVal,AmmoniaExp(Ethanol,1),AmmoniaExp(Ethanol,2),'Color',Colors(Ethanol,:),'Marker','*','MarkerSize',4);
errorbar(BestFitVal,AmmoniaExp(Butanediol,1),AmmoniaExp(Butanediol,2),'Color',Colors(Butanediol,:),'Marker','s','MarkerSize',4);
errorbar(BestFitVal,AmmoniaExp(Acetate,1),AmmoniaExp(Acetate,2),'Color',Colors(Acetate,:),'Marker','d','MarkerSize',4);
errorbar(BestFitVal,AmmoniaExp(ACETONE,1),AmmoniaExp(ACETONE,2),'Color',Colors(ACETONE,:),'Marker','x','MarkerSize',4);

xlabel('CO2 output','FontSize',14)
ylabel('Flux activity in mmol/mmol glucose','FontSize',14)


%% Figure 4b
figure
[data,~,colheaders] = importdata('FinalAmmoniumNoFHLDEHOG.csv')

for i=1:numel(data.data(:,1))
if data.data(i,2) < AmmoniaExp(Biomass,1)
Maintenance = (data.data(i,1) + data.data(i-1,1))/2;
break
end
end
hold off
p = plot(data.data(:,1),data.data(:,2:10),'LineWidth',2)
set(gca,'YLim',[0,2.5])

hold on

p(Biomass).Color = [0 0 0]/255;
p(FHL).Color = [255 192 203]/255;
p(ACETONE).Color = [180 205 20]/255;
p(Acetate).Color = [139 69 19]/255;
p(Ethanol).Color = [255 140 0]/255;
p(Butanediol).Color  = [0 0 255]/255;
p(Nitrite).Color = [255 0 255]/255;
p(Nitrate).Color = [255 0 0]/255;
p(Ammonia).Color = [148 0 211]/255;
p(ACETONE).LineStyle = '--';
p(Butanediol).LineStyle = '--';
p(Nitrite).LineStyle = '-.';
p(Nitrate).LineStyle = '--';
children = get(gca,'children');
delete(children(10-[ Nitrate Nitrite ]));
legend(data.colheaders([[2:7 10]]))
hold on
errorbar(Maintenance,AmmoniaExp(Biomass,1),AmmoniaExp(Biomass,2),'Color',Colors(Biomass,:),'Marker','o','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(Ethanol,1),AmmoniaExp(Ethanol,2),'Color',Colors(Ethanol,:),'Marker','*','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(Butanediol,1),AmmoniaExp(Butanediol,2),'Color',Colors(Butanediol,:),'Marker','s','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(Acetate,1),AmmoniaExp(Acetate,2),'Color',Colors(Acetate,:),'Marker','d','MarkerSize',4);
errorbar(Maintenance,AmmoniaExp(ACETONE,1),AmmoniaExp(ACETONE,2),'Color',Colors(ACETONE,:),'Marker','x','MarkerSize',4);

xlabel('Maintenance ATP per biomass carbon','FontSize',14)
ylabel('Flux activity in mmol/mmol glucose','FontSize',14)


%% Figure 4a
figure 
[data,~,colheaders] = importdata('RedoxScan.csv')


hold off
p = plot(data.data(:,1),data.data(:,2:10),'LineWidth',2)
hold on

p(Biomass).Color = [0 0 0]/255;
p(FHL).Color = [255 192 203]/255;
p(ACETONE).Color = [180 205 20]/255;
p(Acetate).Color = [139 69 19]/255;
p(Ethanol).Color = [255 140 0]/255;
p(Butanediol).Color  = [0 0 255]/255;
p(Nitrite).Color = [255 0 255]/255;
p(Nitrate).Color = [255 0 0]/255;
p(Ammonia).Color = [148 0 211]/255;
p(ACETONE).LineStyle = '--';
p(Butanediol).LineStyle = '--';
p(Nitrite).LineStyle = '-.';
p(Nitrate).LineStyle = '--';
children = get(gca,'children');
delete(children([10-ACETONE 10-Nitrate 10-Nitrite 10-Ammonia]));


Bestfit = 100;
BestFitVal = 0;

for i=1:numel(data.data(:,1))
    currentfit = (AmmoniaExp(Butanediol,1) - data.data(i,Butanediol+1))^2 + (AmmoniaExp(Ethanol,1) - data.data(i,Ethanol+1))^2
    if currentfit < Bestfit
        BestFitVal = data.data(i,1);
        Bestfit = currentfit;
    end
end


legend(data.colheaders(setdiff(1:9, [ACETONE Nitrate Nitrite Ammonia] )+1));
hold on
errorbar(BestFitVal,AmmoniaExp(Biomass,1),AmmoniaExp(Biomass,2),'Color',Colors(Biomass,:),'Marker','o','MarkerSize',4);
errorbar(BestFitVal,AmmoniaExp(Ethanol,1),AmmoniaExp(Ethanol,2),'Color',Colors(Ethanol,:),'Marker','*','MarkerSize',4);
errorbar(BestFitVal,AmmoniaExp(Butanediol,1),AmmoniaExp(Butanediol,2),'Color',Colors(Butanediol,:),'Marker','s','MarkerSize',4);
errorbar(BestFitVal,AmmoniaExp(Acetate,1),AmmoniaExp(Acetate,2),'Color',Colors(Acetate,:),'Marker','d','MarkerSize',4);

xlabel('Additional reductant requirement per carbon','FontSize',14)
ylabel('Flux activity in mmol/mmol glucose','FontSize',14)
set(gca,'XLim',[-0.7,1]);


%% Figure 5
figure 
[data,~,colheaders] = importdata('FinalNitrateECTCADEHOG.csv')


hold off
p = plot(data.data(:,1),data.data(:,2:10),'LineWidth',2)
set(gca,'YLim',[-0.5,2.5])

hold on
%Plot the Experimental Values
p(Biomass).Color = [0 0 0]/255;
p(FHL).Color = [255 192 203]/255;
p(ACETONE).Color = [180 205 20]/255;
p(Acetate).Color = [139 69 19]/255;
p(Ethanol).Color = [255 140 0]/255;
p(Butanediol).Color  = [0 0 255]/255;
p(Nitrite).Color = [255 0 255]/255;
p(Nitrate).Color = [255 0 0]/255;
p(Ammonia).Color = [148 0 211]/255;
p(ACETONE).LineStyle = '--';
p(Butanediol).LineStyle = '--';
p(Nitrite).LineStyle = '-.';
p(Nitrate).LineStyle = '--';
children = get(gca,'children');
delete(children([10-FHL]));


xlabel('Maintenance ATP per biomass carbon','FontSize',14)
ylabel('Flux activity in mmol/mmol glucose','FontSize',14)
legend(data.colheaders(setdiff(1:9, [FHL] )+1));

for i=1:numel(data.data(:,1))
if data.data(i,2) < NitrateExp(Biomass,1)
Maintenance = (data.data(i,1) + data.data(i-1,1))/2;
break
end
end

hold on

errorbar(Maintenance,NitrateExp(Biomass,1),NitrateExp(Biomass,2),'Color',Colors(Biomass,:),'Marker','o','MarkerSize',4);
errorbar(Maintenance,NitrateExp(Ethanol,1),NitrateExp(Ethanol,2),'Color',Colors(Ethanol,:),'Marker','*','MarkerSize',4);
errorbar(Maintenance,NitrateExp(Butanediol,1),NitrateExp(Butanediol,2),'Color',Colors(Butanediol,:),'Marker','s','MarkerSize',4);
errorbar(Maintenance,NitrateExp(Acetate,1),NitrateExp(Acetate,2),'Color',Colors(Acetate,:),'Marker','d','MarkerSize',4);
errorbar(Maintenance,NitrateExp(ACETONE,1),NitrateExp(ACETONE,2),'Color',Colors(ACETONE,:),'Marker','x','MarkerSize',4);
