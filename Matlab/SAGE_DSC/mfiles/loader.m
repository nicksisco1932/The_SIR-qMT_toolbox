function [vol] = loader(path)
    info = niftiinfo(path);
    vol = double(niftiread(info));
end