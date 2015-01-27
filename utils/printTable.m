function [] = printTable(colVal, colPrintType, colWidth)
% colVal and colPrintType are nrCol size cells containing info per column
% the type is either string or double, if colPrintType(col) is not 'string' then
% it indicates the print type (%4.0d or %2.1e etc.)
% if bottomsep then print horizontal separator under the line, etc



nrCol = size(colWidth,2);
horPad = floor((64 - sum(colWidth + 1) + 1)/2); % starting position (horizontal shift of table)
vertSep = '|';
horSep = '-';

sprintTable = sprintf([repmat(' ', 1, horPad), vertSep]);
for col = 1:nrCol
	if strcmp(colPrintType{col}, 'string')
		temp = [' ', colVal{col}];

	else
		temp = [' ', sprintf(colPrintType{col}, colVal{col})];
	end
	sprintTable = [sprintTable, temp, repmat(' ', 1,...
		colWidth(col) - size(temp, 2)), vertSep];
end

announce(sprintTable, 0, 0)

end