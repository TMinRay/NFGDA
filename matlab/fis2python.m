function [incoef,outcoef,rulelogic] = fis2python(sugfis,fn)
    fis = sugfis.EvalfisData;
    if ~strcmp(fis.type, 'sugeno') || ~strcmp(fis.andMethod, 'prod') || ~strcmp(fis.orMethod, 'probor')
        error("fis.type = %s:: should be 'sugeno'\n fis.andMethod:: = %s should be 'prod'\n fis.orMethod = %s:: should be 'probor'",fis.type,fis.andMethod,fis.orMethod);
    end
    incoef = nan(length(fis.rule),length(fis.input(1).mf(1).params),length(fis.input));
    outcoef = nan(length(fis.rule),length(fis.output(1).mf(1).params),length(fis.output));
    rulelogic = nan(length(fis.rule),2);
    for ir = 1:length(fis.rule)
        dest= fis.rule(ir).consequent;
        if ~min(fis.rule(1).antecedent==fis.rule(1).consequent,[],'all')
            error('fis:InvalidValue', 'fis.rule(%d).consequent ~= fis.rule(%d).antecedent',ir,ir);
            error('%s ~= %s',num2str(fis.rule(1).consequent),num2str(fis.rule(1).antecedent));
        end
        rulelogic(ir,1) = fis.rule(ir).weight;
        rulelogic(ir,2) = fis.rule(ir).connection;
        for iv = 1:size(incoef,3)
            incoef(dest,:,iv) = fis.input(iv).mf(ir).params;
        end
        for iv = 1:size(outcoef,3)
            outcoef(dest,:,iv) = fis.output(iv).mf(ir).params;
        end
    end
    save(fn,"incoef","outcoef","rulelogic");
end