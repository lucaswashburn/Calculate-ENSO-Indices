function [long360] = Convert_to_360_long(long)

if long >= 0
    long360 = long;
end

if long < 0
    long360 = 360+long;
end

