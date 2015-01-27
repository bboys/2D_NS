function [] = printTableSep(colWidth)
% print +-----+------------+ corresponding to colWidth


nrCol = size(colWidth,2);
horPad = floor((64 - sum(colWidth + 1) + 1)/2); % starting position (horizontal shift of table)
horSep = '-';
corChar = '+';

sprintTable = [repmat(' ', 1, horPad), corChar];
for col = 1:nrCol
	sprintTable = [sprintTable, sprintf([repmat(horSep, 1,...
	 colWidth(col)), corChar])];
end

announce(sprintTable, 0, 0)


end