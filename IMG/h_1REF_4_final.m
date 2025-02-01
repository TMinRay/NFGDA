run ../main/NF00_header

for cindex=1:numel(ttable(:,1));
    
    PUTDAT=[ ttable(cindex,:) ];
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
        a=zeros(401,401);
        moriSCRout=[matPATH '/IMG/Z/oriz' PUTDAT num2str(m,'%02i') '.mat'];    
        load(moriSCRout, 'originscore');
        mLINEout=[matPATH '/LINE/Z/linez' PUTDAT num2str(m,'%02i') '.mat'];    
        linez=max(originscore,[],3);
        % for ii=1:401
        %     for jj=1:401
        %         a(ii,jj)=max(originscore(ii,jj,:));
        %     end
        % end
        linez=a;
        save(mLINEout, 'linez');
        clear originscore 
     end
end





   
