clear;
close all;

cd ~/Downloads/
ncdisp('b.e11.B1850C5CN.f09_g16.005.pop.h.SST.040001-049912.nc');

cd /Users/lucaswashburn/github/El-Nino-Project

T = ncread('b.e11.B1850C5CN.f09_g16.005.pop.h.SST.040001-049912.nc','SST');
time = ncread('b.e11.B1850C5CN.f09_g16.005.pop.h.SST.040001-049912.nc','time');
TAREA = ncread('b.e11.B1850C5CN.f09_g16.005.pop.h.SST.040001-049912.nc','TAREA');
TLAT = ncread('b.e11.B1850C5CN.f09_g16.005.pop.h.SST.040001-049912.nc','TLAT');
TLONG =  ncread('b.e11.B1850C5CN.f09_g16.005.pop.h.SST.040001-049912.nc','TLONG');
time_matrix = datevec(time);
size(T);

Nino_12_avgs = [];
Nino_3_avgs = [];
Nino_34_avgs = [];
Nino_4_avgs = [];

Monthly_Data_Matrix = [];

%permute, in order for x and y coordinates to be the right side up
areas_permuted = permute(TAREA, [2,1]);
lat_permuted = permute(TLAT, [2,1]);
long_permuted = permute(TLONG, [2,1]);



%x-axis: 1° = 0.888888 data points
%y-axis: 1° = 2.133333 data points

%Niño 1+2 (0-10S, 90W-80W):



lat0 = -10
lat1 = 0
long0 = Convert_to_360_longitude(-90)
long1 = Convert_to_360_longitude(-80)

%bottom left coordinates
[m,long_coord] = min(abs(long_permuted(1,:)-long0));
[m,lat_coord] = min(abs(lat_permuted(:,long_coord)-lat0));
coord00 = [long_coord, lat_coord]

%top left coordinates
[m,long_coord] = min(abs(long_permuted(1,:)-long0));
[m,lat_coord] = min(abs(lat_permuted(:,long_coord)-lat1));
coord01 = [long_coord, lat_coord]

%top right coordinates
[m,long_coord] = min(abs(long_permuted(1,:)-long1));
[m,lat_coord] = min(abs(lat_permuted(:,long_coord)-lat1));
coord11 = [long_coord, lat_coord]

%bottom right coordinates
[m,long_coord] = min(abs(long_permuted(1,:)-long1));
[m,lat_coord] = min(abs(lat_permuted(:,long_coord)-lat0));
coord10 = [long_coord, lat_coord]





for counter = 1:size(T,4)
    
    surface = T(:,:,1,counter);
    temperatures_permuted = permute(surface,[2 1]);               
    %pcolor(temperatures_permuted(:,:));
    
    

    %shading flat;

    Nino_12_temperatures = temperatures_permuted([coord00(2):coord11(2)],[coord00(1):coord11(1)]);
    Nino_12_areas = areas_permuted([coord00(2):coord11(2)],[coord00(1):coord11(1)]);
%     Nino_12_temperatures = temperatures_permuted([coord00(1):coord11(1)],[coord00(2):coord11(2)]);
%     Nino_12_areas = areas_permuted([coord00(1):coord11(1)],[coord00(2):coord11(2)]);
    Nino_12_avgs(counter,1) = time_matrix(counter, 1);
    %find weighted mean
    Nino_12_avgs(counter,2) = sum(Nino_12_temperatures.*Nino_12_areas,'all','omitnan')/sum(Nino_12_temperatures./Nino_12_temperatures.*Nino_12_areas,'all','omitnan');
    Nino_12_avgs(counter,3) = counter;
    
    %get data by month
    Monthly_Data_Matrix(time_matrix(counter, 2),time_matrix(counter, 1) + 1 - time_matrix(1,1),1) = Nino_12_avgs(counter,2);

end

%calculate monthly averages
for month_number = 1:12
    monthly_averages(month_number) = sum(Monthly_Data_Matrix(month_number,:))/length(nonzeros(Monthly_Data_Matrix(month_number,:)));
end

%calculate anomaly
for month_number = 1:12
    for year_number = 1:time_matrix(size(time),1) + 1 - time_matrix(1,1);
        if Monthly_Data_Matrix(month_number,year_number,1) ~= 0
            Monthly_Data_Matrix(month_number,year_number,2) = Monthly_Data_Matrix(month_number,year_number,1) - monthly_averages(month_number);
        end
    end
end

for counter = 1:size(T,4)
    month_number = time_matrix(counter,2)
    %calculate anomaly
    Nino_12_avgs(counter,4) = Nino_12_avgs(counter,2) - monthly_averages(month_number);
    
end

for counter = 1
    Nino_12_avgs(counter,5) = mean(nonzeros([Nino_12_avgs(counter,4) Nino_12_avgs(counter+1,4) Nino_12_avgs(counter+2,4)]),'all','omitnan');
end

for counter = 2
    Nino_12_avgs(counter,5) = mean(nonzeros([Nino_12_avgs(counter-1,4) Nino_12_avgs(counter,4) Nino_12_avgs(counter+1,4) Nino_12_avgs(counter+2,4)]),'all','omitnan');
end

for counter = 3:size(T,4)-2
    Nino_12_avgs(counter,5) = mean(nonzeros([Nino_12_avgs(counter-2,4) Nino_12_avgs(counter-1,4) Nino_12_avgs(counter,4) Nino_12_avgs(counter+1,4) Nino_12_avgs(counter+2,5)]),'all','omitnan');
end

for counter = size(T,4)-1
    Nino_12_avgs(counter,5) = mean(nonzeros([Nino_12_avgs(counter-2,4) Nino_12_avgs(counter-1,4) Nino_12_avgs(counter,4) Nino_12_avgs(counter+1,4)]),'all','omitnan');
end

for counter = size(T,4)
    Nino_12_avgs(counter,5) = mean(nonzeros([Nino_12_avgs(counter-2,4) Nino_12_avgs(counter-1,4) Nino_12_avgs(counter,4)]),'all','omitnan');
end

Nino_12_smoothed_STDEV = std(Nino_12_avgs(:,5));

%normalized for standard deviation
for counter = 1:size(T,4)
    Nino_12_avgs(counter,6) = Nino_12_avgs(counter,5)/Nino_12_smoothed_STDEV;
end

% %Niño 3 (5N-5S, 150W-90W): 
% lat5N = 183;
% lat5S = 161;
% long90W = 275;
% long150W = 222;
% 
% Nino_3 = surface_permuted([lat5S:lat5N],[long150W:long90W]);
% Nino_3_avgs(counter,1) = counter;
% Nino_3_avgs(counter,2) = mean(Nino_3,'all','omitnan');
% 
% %Niño 3.4 (5N-5S, 170W-120W):
% lat5N = 183;
% lat5S = 161;
% long120W = 248;
% long170W = 204;
% 
% Nino_34 = surface_permuted([lat5S:lat5N],[long170W:long120W]);
% Nino_34_avgs(counter,1) = counter;
% Nino_34_avgs(counter,2) = mean(Nino_34,'all','omitnan');
% 
% %Niño 4 (5N-5S, 160E-150W): 
% lat5N = 183;
% lat5S = 161;
% long161E = 177;
% long150W = 222;
% 
% Nino_4 = surface_permuted([lat5S:lat5N],[long161E:long150W]);
% Nino_4_avgs(counter,1) = counter;
% Nino_4_avgs(counter,2) = mean(Nino_4,'all','omitnan');





figure;
plot(Nino_12_avgs(:,1),Nino_12_avgs(:,2))
title('Nino 1+2 Index SST vs. Time')
ylabel('Sea Surface Temperature (°C)')

figure;
plot(Nino_12_avgs(:,1),Nino_12_avgs(:,4))
title('Nino 1+2 Index SST vs. Time')
ylabel('Sea Surface Temperature (°C)')

figure;
plot(Nino_12_avgs(:,1),Nino_12_avgs(:,5))
title('Nino 1+2 Index SST vs. Time')
ylabel('Sea Surface Temperature (°C)')

figure;
hold on;
plot(Nino_12_avgs(:,1),Nino_12_avgs(:,6))
fplot(0, [time_matrix(1,1),time_matrix(length(time),1)]);
title('Nino 1+2 Index SST vs. Time')
ylabel('Standard Deviations from Monthly Climatological SST')
xlabel('year')


% figure;
% plot(Nino_3_avgs(:,1),Nino_3_avgs(:,2))
% title('Nino 3 Index SST vs. Time')
% ylabel('Sea Surface Temperature (°C)')
% figure;
% plot(Nino_34_avgs(:,1),Nino_34_avgs(:,2))
% title('Nino 3.4 Index SST vs. Time')
% ylabel('Sea Surface Temperature (°C)')
% xlabel('time in months')
% figure;
% plot(Nino_4_avgs(:,1),Nino_4_avgs(:,2))
% title('Nino 4 Index SST vs. Time')
% ylabel('Sea Surface Temperature (°C)')
% 
% 
