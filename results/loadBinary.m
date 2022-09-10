function [ out ] = loadBinary(fname, type, sz)
    fileid = fopen(fname, 'r');
    out = fread(fileid, [prod(sz), 1], type);
    fclose(fileid);
    out = reshape(out, sz);
end
