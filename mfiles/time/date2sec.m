function s = date2sec(dateStr, refDateStr)
%
% Date string to seconds from reference date
%
% Creator John Merryman - sarmap
% Date: 22 Mar 2017
%
% Usage:
%   s = date2sec(dateStr, [refDateStr])
%
% dateStr  : (input)  Date string in SARscape format
% refDate  : (input)  Reference date (default = '01-Jan-1970').
% s        : (output) Seconds elapsed since reference date.
%
% Examples:
%     s1 = date2sec('08-OCT-2015 08:32:42.063759');
%     s2 = date2sec('08-OCT-2015');
%
if nargin < 2
    refDateStr = '01-Jan-1990';
end

% Get year, month, day, h, m ,s
[Y1, M1, D1, H1, MN1, S1] = datevec(dateStr);
% Compute seconds properly if they were specified:
C = strsplit(dateStr,':');
if size(C,2) == 3
    S1 = str2num(C{1,size(C,2)});
end

dateStr1 = datestr([Y1 M1 D1 0 0 0]);
dateStr2 = refDateStr;

d1 = datenum(dateStr1);
d2 = datenum(dateStr2);

s = 86400 * (d1 - d2) + 3600 * H1 + 60 * MN1 + S1;
