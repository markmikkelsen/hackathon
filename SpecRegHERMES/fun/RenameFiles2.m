function RenameFiles2(exp, subj)

flist.sdat = dir('*.sdat');
flist.spar = dir('*.spar');
if isempty(flist.sdat)
    flist.sdat = dir('*.SDAT');
    flist.spar = dir('*.SPAR');
end

fnames(1) = cellstr(flist.sdat(1).name);
fnames(2) = cellstr(flist.spar(1).name);
fnames(3) = cellstr(flist.sdat(2).name);
fnames(4) = cellstr(flist.spar(2).name);

new_fnames = {[exp '_' subj '_HERMES_act.sdat'], [exp '_' subj '_HERMES_act.spar'], ...
              [exp '_' subj '_HERMES_ref.sdat'], [exp '_' subj '_HERMES_ref.spar']};

for ii = 1:length(fnames)
    java.io.File(fnames{ii}).renameTo(java.io.File(new_fnames{ii}));
end