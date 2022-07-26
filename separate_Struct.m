function varargout = separate_Struct(myStruct)

probFieldNames = fieldnames(myStruct);
% Extract values of all the fields
for i = 1:length(probFieldNames) - 1 % Last field is the 'solver' field
    varargout{i} = myStruct.(probFieldNames{i});
end
% varargout{end+1} = options; % Stuff options as the last field







end