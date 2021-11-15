function [long360] = Convert_to_360_longitude(long)

if long >= 0
    long360 = long;
end

if long < 0
    long360 = 360+long;
end

