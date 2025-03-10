% clear all;close all;clc;
% ======================================================================
% Create initial time table, 
% daily and 3-Hourly for picking files from 2011 10/01-12-31.
% ======================================================================
% Daily
% Arrange Month
TTD(1:(31),1) = 10;
TTD((31+1):(31+30),1) = 11;
TTD((31+30+1):(31+30+31),1) = 12;
% Arrange Date
TTD(1:(31),2) = [1:31];
TTD((31+1):(31+30),2) = [1:30];
TTD((31+30+1):(31+30+31),2) = [1:31];
% ======================================================================
% 3-Hourly
% Arrange Month
TT(1:(31*8),1) = 10;
TT((31*8+1):(31*8+30*8),1) = 11;
TT((31*8+30*8+1):(31*8+30*8+31*8),1) = 12;
% Arrange Date
for i = 1:31;
    TT((8*(i-1))+1:(8*(i-1))+8,2) = i;
end
for i = 32:61;
    TT((8*(i-1))+1:(8*(i-1))+8,2) = i-31;
end
for i = 62:92;
    TT((8*(i-1))+1:(8*(i-1))+8,2) = i-61;
end
% Arrange Hour
for i = 1:92;
    TT((8*(i-1))+1:(8*(i-1))+8,3) = (0:3:21);
end
% ======================================================================
clear i;