function [] = announce(inputText, printTop, printBottom, txtColor)
% prints inputText in following style
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%									%
% 	blablallba .................	%
% 	............................	%
%									%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
	printBottom = 0;
end
if nargin < 2
	printTop = 0;
end

horPadding = 2; % horizontal padding
vertPadding = 0; % vertical lines
textWidth = 70 - 2 - 2*horPadding; % nr of characters per line

horChar = '-';
vertChar = '.';
cornerChar = 'o';

if nargin < 4
	txtColor = [0 0 0];
end


nrChar = size(inputText, 2);
lastChar = 0; % last character printed

if printTop == 1
	fprintf([cornerChar, repmat(horChar,1 , 2*horPadding + textWidth),cornerChar,...
	 '\n'])
end
fprintf(repmat([vertChar, repmat(' ',1 , 2*horPadding + textWidth),vertChar,...
 '\n'], vertPadding, 1))
while lastChar < nrChar
	lineText = inputText(lastChar + 1:min(nrChar, lastChar + textWidth));

	if strcmp(lineText(1), ' ')
		lastChar = lastChar + 1;
		lineText = inputText(lastChar + 1:min(nrChar, lastChar + textWidth));
	end
	if ~(lastChar + textWidth + 1 > nrChar)
		if ~strcmp(inputText(lastChar + textWidth + 1), ' ')
			if strcmp(inputText(lastChar + textWidth), ' ') 
				lineText = [lineText(1:end - 1), ' '];
			else
				lineText = [lineText(1:end - 1), '-'];
			end
			lastChar = lastChar - 1;	
		end
	else
		lineText = [lineText, repmat(' ', 1, textWidth - size(lineText, 2))];
	end


	lastChar = lastChar + textWidth;
	% fprintf([vertChar, repmat(' ',1 , horPadding), lineText,...
	% 	repmat(' ',1 , horPadding), vertChar, '\n'])
	fprintf([vertChar, repmat(' ',1 , horPadding)])
	cprintf(txtColor, lineText)
	fprintf(repmat(' ',1 , horPadding))
	cprintf([0 0 0], [vertChar, '\n'])

end
fprintf(repmat([vertChar, repmat(' ',1 , 2*horPadding + textWidth),vertChar,...
 '\n'], vertPadding, 1))
if printBottom == 1
	fprintf([cornerChar, repmat(horChar,1 , 2*horPadding + textWidth),cornerChar,...
	 '\n'])
end


end