% Recoleccion de datos desde xls obtenidos del Registro Civil, agregados y
% guardados en un archivo .mat
%
% 2020 Vicente Parot
% Wellman Center for Photomedicine
% Harvard Medical School
%
% Instituto de Ingeniería Biológica y Médica
% Pontificia Universidad Católica de Chile
%
%% read tables
dirlist = dir('*.xlsx');
[~, idx] = sort({dirlist.name});
dirlist = dirlist(idx);
defs = [];
for it = 1:numel(dirlist)
    disp(it)
    defs = [defs;
        readtable(dirlist(it).name)
        ];
end
disp done
fecha = datestr(datetime(dirlist(end).name(end-14:end-5),'inputformat','yyyy-MM-dd'),'DD/mm');
save defunciones defs fecha
