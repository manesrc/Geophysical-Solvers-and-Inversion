function varargout = loadStructFromFile(fileName,environmentName)
  varargout = struct2cell(load(fileName,environmentName));
end