clear;
close all;

cd ~/Downloads/;
ncdisp('b.e11.B1850C5CN.f09_g16.005.pop.h.SST.040001-049912.nc');

T = ncread('b.e11.B1850C5CN.f09_g16.005.pop.h.SST.040001-049912.nc','SST');
time = ncread('b.e11.B1850C5CN.f09_g16.005.pop.h.SST.040001-049912.nc','time');
time_matrix = datevec(time);
size(T);

Nino_12_avgs = [];
Nino_3_avgs = [];
Nino_34_avgs = [];
Nino_4_avgs = [];

steez = [];




%x-axis: 1° = 0.888888 data points
%y-axis: 1° = 2.133333 data points

%Niño 1+2 (0-10S, 90W-80W):
lat0 = 172;
lat10S = 151;
long90W = 275;
long80W = 284;

for counter = 1:size(T,4)
    
    surface = T(:,:,1,counter);
    surface_permuted = permute(surface,[2 1]);               
    pcolor(surface_permuted(:,:));

    shading flat;

    Nino_12 = surface_permuted([lat10S:lat0],[long90W:long80W]);
    Nino_12_avgs(counter,1) = time_matrix(counter, 1);
    Nino_12_avgs(counter,2) = mean(Nino_12,'all','omitnan');
    Nino_12_avgs(counter,3) = counter;
    %get data by month
    steez(time_matrix(counter, 2),time_matrix(counter, 1) + 1 - time_matrix(1,1),1) = Nino_12_avgs(counter,2);

end

%calculate monthly averages
for month_number = 1:12
    monthly_averages(month_number) = sum(steez(month_number,:))/length(nonzeros(steez(month_number,:)));
end

%calculate anomaly
for month_number = 1:12
    for year_number = 1:time_matrix(size(time),1) + 1 - time_matrix(1,1);
        if steez(month_number,year_number,1) ~= 0
            steez(month_number,year_number,2) = steez(month_number,year_number,1) - monthly_averages(month_number);
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




steez = [];

%Niño 3
lat5N = 183;
lat5S = 161;
long90W = 275;
long150W = 222;

for counter = 1:size(T,4)
    
    surface = T(:,:,1,counter);
    surface_permuted = permute(surface,[2 1]);               
    pcolor(surface_permuted(:,:));

    shading flat;

    Nino_3 = surface_permuted([lat5S:lat5N],[long150W:long90W]);
    Nino_3_avgs(counter,1) = time_matrix(counter, 1);
    Nino_3_avgs(counter,2) = mean(Nino_3,'all','omitnan');
    Nino_3_avgs(counter,3) = counter;
    %get data by month
    steez(time_matrix(counter, 2),time_matrix(counter, 1) + 1 - time_matrix(1,1),1) = Nino_3_avgs(counter,2);

end

%calculate monthly averages
for month_number = 1:12
    monthly_averages(month_number) = sum(steez(month_number,:))/length(nonzeros(steez(month_number,:)));
end

%calculate anomaly
for month_number = 1:12
    for year_number = 1:time_matrix(size(time),1) + 1 - time_matrix(1,1);
        if steez(month_number,year_number,1) ~= 0
            steez(month_number,year_number,2) = steez(month_number,year_number,1) - monthly_averages(month_number);
        end
    end
end

for counter = 1:size(T,4)
    month_number = time_matrix(counter,2)
    %calculate anomaly
    Nino_3_avgs(counter,4) = Nino_3_avgs(counter,2) - monthly_averages(month_number);
    
end

for counter = 1
    Nino_3_avgs(counter,5) = mean(nonzeros([Nino_3_avgs(counter,4) Nino_3_avgs(counter+1,4) Nino_3_avgs(counter+2,4)]),'all','omitnan');
end

for counter = 2
    Nino_3_avgs(counter,5) = mean(nonzeros([Nino_3_avgs(counter-1,4) Nino_3_avgs(counter,4) Nino_3_avgs(counter+1,4) Nino_3_avgs(counter+2,4)]),'all','omitnan');
end

for counter = 3:size(T,4)-2
    Nino_3_avgs(counter,5) = mean(nonzeros([Nino_3_avgs(counter-2,4) Nino_3_avgs(counter-1,4) Nino_3_avgs(counter,4) Nino_3_avgs(counter+1,4) Nino_3_avgs(counter+2,5)]),'all','omitnan');
end

for counter = size(T,4)-1
    Nino_3_avgs(counter,5) = mean(nonzeros([Nino_3_avgs(counter-2,4) Nino_3_avgs(counter-1,4) Nino_3_avgs(counter,4) Nino_3_avgs(counter+1,4)]),'all','omitnan');
end

for counter = size(T,4)
    Nino_3_avgs(counter,5) = mean(nonzeros([Nino_3_avgs(counter-2,4) Nino_3_avgs(counter-1,4) Nino_3_avgs(counter,4)]),'all','omitnan');
end

Nino_12_smoothed_STDEV = std(Nino_3_avgs(:,5));

%normalized for standard deviation
for counter = 1:size(T,4)
    Nino_3_avgs(counter,6) = Nino_3_avgs(counter,5)/Nino_12_smoothed_STDEV;
end


Nino_34_avgs = [];



steez = [];






%Niño 3.4
lat5N = 183;
lat5S = 161;
long120W = 248;
long170W = 204;

for counter = 1:size(T,4)
    
    surface = T(:,:,1,counter);
    surface_permuted = permute(surface,[2 1]);               
    pcolor(surface_permuted(:,:));

    shading flat;

    Nino_34 = surface_permuted([lat5S:lat5N],[long170W:long120W]);
    Nino_34_avgs(counter,1) = time_matrix(counter, 1);
    Nino_34_avgs(counter,2) = mean(Nino_34,'all','omitnan');
    Nino_34_avgs(counter,3) = counter;
    %get data by month
    steez(time_matrix(counter, 2),time_matrix(counter, 1) + 1 - time_matrix(1,1),1) = Nino_34_avgs(counter,2);

end

%calculate monthly averages
for month_number = 1:12
    monthly_averages(month_number) = sum(steez(month_number,:))/length(nonzeros(steez(month_number,:)));
end

%calculate anomaly
for month_number = 1:12
    for year_number = 1:time_matrix(size(time),1) + 1 - time_matrix(1,1);
        if steez(month_number,year_number,1) ~= 0
            steez(month_number,year_number,2) = steez(month_number,year_number,1) - monthly_averages(month_number);
        end
    end
end

for counter = 1:size(T,4)
    month_number = time_matrix(counter,2)
    %calculate anomaly
    Nino_34_avgs(counter,4) = Nino_34_avgs(counter,2) - monthly_averages(month_number);
    
end

for counter = 1
    Nino_34_avgs(counter,5) = mean(nonzeros([Nino_34_avgs(counter,4) Nino_34_avgs(counter+1,4) Nino_34_avgs(counter+2,4)]),'all','omitnan');
end

for counter = 2
    Nino_34_avgs(counter,5) = mean(nonzeros([Nino_34_avgs(counter-1,4) Nino_34_avgs(counter,4) Nino_34_avgs(counter+1,4) Nino_34_avgs(counter+2,4)]),'all','omitnan');
end

for counter = 3:size(T,4)-2
    Nino_34_avgs(counter,5) = mean(nonzeros([Nino_34_avgs(counter-2,4) Nino_34_avgs(counter-1,4) Nino_34_avgs(counter,4) Nino_34_avgs(counter+1,4) Nino_34_avgs(counter+2,5)]),'all','omitnan');
end

for counter = size(T,4)-1
    Nino_34_avgs(counter,5) = mean(nonzeros([Nino_34_avgs(counter-2,4) Nino_34_avgs(counter-1,4) Nino_34_avgs(counter,4) Nino_34_avgs(counter+1,4)]),'all','omitnan');
end

for counter = size(T,4)
    Nino_34_avgs(counter,5) = mean(nonzeros([Nino_34_avgs(counter-2,4) Nino_34_avgs(counter-1,4) Nino_34_avgs(counter,4)]),'all','omitnan');
end

Nino_12_smoothed_STDEV = std(Nino_34_avgs(:,5));

%normalized for standard deviation
for counter = 1:size(T,4)
    Nino_34_avgs(counter,6) = Nino_34_avgs(counter,5)/Nino_12_smoothed_STDEV;
end

Nino_4_avgs = [];


steez = [];


%Niño 3.4
lat5N = 183;
lat5S = 161;
long161E = 177;
long150W = 222;

for counter = 1:size(T,4)
    
    surface = T(:,:,1,counter);
    surface_permuted = permute(surface,[2 1]);               
    pcolor(surface_permuted(:,:));

    shading flat;

    Nino_4 = surface_permuted([lat5S:lat5N],[long161E:long150W]);
    Nino_4_avgs(counter,1) = time_matrix(counter, 1);
    Nino_4_avgs(counter,2) = mean(Nino_4,'all','omitnan');
    Nino_4_avgs(counter,3) = counter;
    %get data by month
    steez(time_matrix(counter, 2),time_matrix(counter, 1) + 1 - time_matrix(1,1),1) = Nino_4_avgs(counter,2);

end

%calculate monthly averages
for month_number = 1:12
    monthly_averages(month_number) = sum(steez(month_number,:))/length(nonzeros(steez(month_number,:)));
end

%calculate anomaly
for month_number = 1:12
    for year_number = 1:time_matrix(size(time),1) + 1 - time_matrix(1,1);
        if steez(month_number,year_number,1) ~= 0
            steez(month_number,year_number,2) = steez(month_number,year_number,1) - monthly_averages(month_number);
        end
    end
end

for counter = 1:size(T,4)
    month_number = time_matrix(counter,2)
    %calculate anomaly
    Nino_4_avgs(counter,4) = Nino_4_avgs(counter,2) - monthly_averages(month_number);
    
end

for counter = 1
    Nino_4_avgs(counter,5) = mean(nonzeros([Nino_4_avgs(counter,4) Nino_4_avgs(counter+1,4) Nino_4_avgs(counter+2,4)]),'all','omitnan');
end

for counter = 2
    Nino_4_avgs(counter,5) = mean(nonzeros([Nino_4_avgs(counter-1,4) Nino_4_avgs(counter,4) Nino_4_avgs(counter+1,4) Nino_4_avgs(counter+2,4)]),'all','omitnan');
end

for counter = 3:size(T,4)-2
    Nino_4_avgs(counter,5) = mean(nonzeros([Nino_4_avgs(counter-2,4) Nino_4_avgs(counter-1,4) Nino_4_avgs(counter,4) Nino_4_avgs(counter+1,4) Nino_4_avgs(counter+2,5)]),'all','omitnan');
end

for counter = size(T,4)-1
    Nino_4_avgs(counter,5) = mean(nonzeros([Nino_4_avgs(counter-2,4) Nino_4_avgs(counter-1,4) Nino_4_avgs(counter,4) Nino_4_avgs(counter+1,4)]),'all','omitnan');
end

for counter = size(T,4)
    Nino_4_avgs(counter,5) = mean(nonzeros([Nino_4_avgs(counter-2,4) Nino_4_avgs(counter-1,4) Nino_4_avgs(counter,4)]),'all','omitnan');
end

Nino_12_smoothed_STDEV = std(Nino_4_avgs(:,5));

%normalized for standard deviation
for counter = 1:size(T,4)
    Nino_4_avgs(counter,6) = Nino_4_avgs(counter,5)/Nino_12_smoothed_STDEV;
end



figure;
hold on;
plot(Nino_12_avgs(:,1),Nino_12_avgs(:,6))
plot(Nino_3_avgs(:,1),Nino_3_avgs(:,6))
plot(Nino_34_avgs(:,1),Nino_34_avgs(:,6))
plot(Nino_4_avgs(:,1),Nino_4_avgs(:,6))
fplot(0, [time_matrix(1,1),time_matrix(length(time),1)]);
title('Nino 1+2 Index SST vs. Time')
ylabel('Standard Deviations from Monthly Climatological SST')
xlabel('year')