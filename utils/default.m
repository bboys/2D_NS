function output = default(query, defInput, comments, outputType)
% if query is class char
% outputs user input or default value 'defInput' (for parameter input)

% else query is a cell array containing the usual text explaining which choise
% to make, and multiple options to choose from 
% note the output is the corresponding integer for outputType = 'int'
% and the string is given for outputType = 'string'

horChar = '.';
horPad = 2;
textWidth = 70 - 2 - 2*horPad;
hlColor = [0.5 0 0]; % color of highlighted text

% allow output to be user string
if nargin < 4
	outputType = 'int';
end

fprintf(['.', repmat(' ', 1, horPad)])
if strcmp(class(query), 'char')
	printedChar = 0;
	temp = sprintf('%d', defInput);

	if size(temp,2) > 5
		temp = sprintf('%2.1e', defInput);
	end

	% string to display default value
	temp = [' (default value is ',temp,'): '];
	fprintf([query, temp])
	printedChar = printedChar + size(temp,2);

	% request input
	userInp = input('', 's'); % request input string

	if isempty(userInp) 
		output = defInput;
	elseif strcmp(outputType, 'int') == 1
		printedChar = printedChar + size(userInp, 2);
		output = str2num(userInp);  % convert to number/double
	else
		printedChar = printedChar + size(userInp, 2);
		output = userInp;		
	end

	if strcmp(outputType, 'int') == 1
		% prevent large precision numbers to be printed
		temp = sprintf('%d', output);
		if size(temp,2) > 5
			temp = sprintf('%2.1e', output);
		end
	end

	inputSize = size(temp, 2);
	fprintf([repmat('\b', 1, 1 + printedChar), ': ',...
		repmat(horChar, 1, textWidth - size(query,2) - inputSize - 3),...
		 ' '])% ,temp])
	cprintf(hlColor, [temp, '0']) % weird bug
	fprintf('\b')

elseif strcmp(class(query), 'cell')
	printedChar = 0; % keep track of nr of printed chars such that they can be removed
	nrOptions = size(query,2) - 1;
	if nargin < 3
		% no comments were passed
		comments = repmat({''}, 1, nrOptions);
	end

	temp = ', the options are:';
	fprintf([query{1}, temp, '\n'])
	printedChar = printedChar + size(temp,2);

	for opt = 1:nrOptions

		temp = [' %3d: ', query{opt + 1}];
		% fprintf(temp, opt)
		cprintf(hlColor, sprintf(temp, opt))
		printedChar = printedChar + size(temp,2);
		if ~strcmp(comments{opt}, '')
			temp = [' (',comments{opt},')'];
			fprintf(temp)
			printedChar = printedChar + size(temp,2);

		end
		if opt == defInput

			temp = ' (default choise)';
			fprintf(temp)
			printedChar = printedChar + size(temp,2);

		end
		fprintf('\n')
		
	end

	% request input
	userInp = input('');
	if isempty(userInp) | ~strcmp(class(userInp), 'double')
		output = defInput;
	else
		% check if userInp is in options
		if isempty(find([1:nrOptions] == userInp(1)))
			userInp = defInput;
		end

		output = userInp;
	end
	inputSize = size(sprintf(query{output + 1}), 2);

	% remove userinput from screen
	fprintf(repmat('\b',1 , 1 + size(userInp,2) + printedChar + nrOptions + 1))
	fprintf([': ', repmat(horChar, 1, textWidth - size(query{1} ,2) - inputSize - 3),...
		' ']) %,query{output + 1}])
	cprintf(hlColor, query{output + 1})

end
fprintf(repmat(' ', 1, horPad))
cprintf([0 0 0], '.\n')

end

