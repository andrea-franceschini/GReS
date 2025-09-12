function classes = findSubClasses(superclassName, folder)
    path = fullfile(getappdata(0, 'gres_root'), 'Code/discretizer/', folder, '*.m');
    files = dir(path);
    classes = {};

    for k = 1:numel(files)
        name = erase(files(k).name, '.m');
        mc = meta.class.fromName(name);
        if isempty(mc)
            continue;
        end
        if inheritsFrom(mc, superclassName)
            classes{end+1} = name; 
        end
    end
end


function result = inheritsFrom(mc, targetName)
    % Recursively check if mc inherits from targetName
    result = false;
    if isempty(mc.SuperclassList)
        return;
    end
    for s = mc.SuperclassList
        if strcmp(s.Name, targetName)
            result = true;
            return;
        else
            result = inheritsFrom(s, targetName);
            if result
                return;
            end
        end
    end
end
