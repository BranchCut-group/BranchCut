function [hmsStr, h, m, s] = sec2hms(sec)
%
% Convert seconds of day to hour, minute, seconds string
%
% Creator John Merryman - sarmap
% Date: 22 Mar 2017
%
% Usage:
%   dateStr = sec2hms(secs)
%
% secs     : (input)   Seconds elapsed since midnight.
% hmsStr   : (output)  'HH:MM:SS.SSSSSSS'.
%
% Examples:
%    sec2hms(49565.821180)
%    ans =  '13:46:05.8211800'
%

secondsPerHour = 3600;
minutesPerHour = 60;
secondsPerMinute = 60;

decimalHours = sec / secondsPerHour;
integerHours = floor(decimalHours);
h = integerHours;

decimalMinutes = (decimalHours - integerHours)*minutesPerHour;
integerMinutes = floor(decimalMinutes);
m = integerMinutes;

decimalSeconds = (decimalMinutes - integerMinutes)*secondsPerMinute;
s = decimalSeconds;

hmsStr = sprintf('%02d:%02d:%010.7f', h, m, s);
