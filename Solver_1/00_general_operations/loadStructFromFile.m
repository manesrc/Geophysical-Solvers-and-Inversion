function varargout = loadStructFromFile(fileName,environmentName)
%LOADSTRUCTFROMFILE Load a structure from a .mat file and return its fields.
%   VARARGOUT = LOADSTRUCTFROMFILE(FILENAME, ENVIRONMENTNAME) loads the
%   structure stored in the .mat file specified by FILENAME under the name
%   ENVIRONMENTNAME. It returns the fields of the structure as separate
%   output arguments.

  varargout = struct2cell(load(fileName,environmentName));
end