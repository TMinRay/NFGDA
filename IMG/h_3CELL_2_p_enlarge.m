run ../main/NF00_header

load('tmpCELLdatax.mat','datacx');
load('tmpCELLdatay.mat','datacy');
load('tmpCELLdatax2.mat','datacx2');
load('tmpCELLdatay2.mat','datacy2');
indexnum=numel(datacx);
indexnum2=numel(datacx2);
crsize=5;
numINT =crsize+4;

for cindex=1:numel(ttable(:,1));
    PUTDAT=[ ttable(cindex,:) ];
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
        mLINEout=[matPATH '/LINE/CELLZ/cellz' PUTDAT num2str(m,'%02i') '.mat'];
        load(mLINEout, 'CELLline');
        mLINEout2=[matPATH '/LINE/CELLZ/widecellz' PUTDAT num2str(m,'%02i') '.mat'];
        a2 = CELLline;
        [row_indices, col_indices] = find(a2>cellcsrthresh);
        inINT = (row_indices>numINT) & (row_indices<=401-numINT)...
            &(col_indices>numINT) & (col_indices<=401-numINT);
        row_indices = row_indices(inINT);
        col_indices = col_indices(inINT);
        cridx = row_indices.' + datacy;
        ccidx = col_indices.' + datacx;
        c_indices = sub2ind(size(a2), cridx, ccidx);
        cbox=double(a2(c_indices)>cellcsrthresh);
        numtcb=sum(cbox,1);
        cbr = numtcb/indexnum;
        row_indices = row_indices(cbr<1);
        col_indices = col_indices(cbr<1);
        cridx = row_indices.' + datacy2;
        ccidx = col_indices.' + datacx2;
        c_indices = sub2ind(size(a2), cridx, ccidx);
        a2(c_indices)=1;

        % cellscr2=CELLline;
        % tnum=1;
        % indcell=(cellscr2>cellcsrthresh).*tnum;
        
        % for ii=1+numINT:401-numINT
        %     for jj=1+numINT:401-numINT
        %         if (indcell(ii,jj)==tnum)
        %         cbox=[];
        %         for kind=1:indexnum
        %             cadd=[indcell(ii+datacy(kind),jj+datacx(kind))];
        %             cbox=[cbox; cadd;];
        %             clear cadd;
        %         end
        %         indcb=((cbox)==tnum);
        %         numtcb=sum(indcb(:));
        %         cbr=numtcb/indexnum;
        %         if (cbr<1)
        %             for kind=1:indexnum2
        %                 cellscr2(ii+datacy2(kind),jj+datacx2(kind))=1;
        %             end

        %         end
        %         end
        %     end
        % end
        widecellz=double(a2>0.5);
        save(mLINEout2, 'widecellz');
    end
end
