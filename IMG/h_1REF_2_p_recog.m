run ../main/NF00_header


load('tmpREFcdatax.mat','datacx');
load('tmpREFcdatay.mat','datacy');
load('tmpREFsdatax.mat','datasx');
load('tmpREFsdatay.mat','datasy');
indexcnum=numel(datacx);
indexsnum=numel(datasx);
refcentersize=8;
refapart=5;
 
for cindex=1:numel(ttable(:,1));
    
    PUTDAT=[ ttable(cindex,:) ];
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
        mROTout=[matPATH '/IMG/Z/rotz' PUTDAT num2str(m,'%02i') '.mat'];
        load(mROTout, 'rotitp');
        mSCRout=[matPATH '/IMG/Z/scrz' PUTDAT num2str(m,'%02i') '.mat'];
        llscore=[];
        ssscore=[];
        totscore2=zeros(401,401);
        totscore=zeros(401,401,rotnum); 
        for i=1:rotnum
            indi=rotdegree*(i-1);
            a2=rotitp(:,:,i);
            clscore=zeros(401,401);
            sdscore=zeros(401,401);
            inda2=(rotitp(:,:,i)>thrREF);
            numINT =refcentersize+2;
            nunapart=refapart;  
            % totscore(:,1:numINT,i)=nan;
            % totscore(:,401-numINT+1:401,i)=nan;
            % totscore(1:numINT,:,i)=nan;
            % totscore(401-numINT+1:401,:,i)=nan; 
            for ii=1+numINT:401-numINT       
                for jj=1+numINT:401-numINT
                    % cbox=[];
                    % sbox=[];        
                    if (inda2(ii,jj)==1)
                        row_indices = ii + datacy; 
                        col_indices = jj + datacx;
                        lin_indices = sub2ind(size(a2), row_indices, col_indices);
                        cbox=a2(lin_indices);

                        row_indices = ii + datasy; 
                        col_indices = jj + datasx;
                        lin_indices = sub2ind(size(a2), row_indices, col_indices);
                        % for kind=1:indexcnum
                        %     cadd=[a2(ii+datacy(kind),jj+datacx(kind))];
                        %     cbox=[cbox; cadd;];
                        %     clear cadd;
                        % end
                        % for kind=1:indexsnum
                        %     sadd=[a2(ii+datasy(kind),jj+datasx(kind))];
                        %     sbox=[sbox; sadd;];
                        %     clear sadd;
                        % end 
                        indcb=((cbox)>thrREF);   
                        % indsb=(~isnan(sbox));  
                        numtcb=sum(indcb(:));
                        % numtsb=sum(indsb(:));
                        cbr=numtcb/indexcnum;
                        % sbr=numtsb/indexsnum;
                        if (cbr>=0.5 )
                        kkk=0;   
                        h_1REF_2_SCRfunc;
                        else
                        sdscore(ii,jj)=0; clscore(ii,jj)=0;  
                        end %% for cbr>=0.5  
                    else
                        sdscore(ii,jj)=0; clscore(ii,jj)=0; 
                    end %%%%%%%%%%%%%%%% for inda2
                    % clear cbox sbox;
                end %end Azimuth
            end
            pretotscore=(sdscore+clscore);
            totscore2=pretotscore./(3*17+1*18);
            totscore2(totscore2<0)=0;
            % for ki=1:401
            %     for kj=1:401
            %         if totscore2(ki,kj)<0
            %         totscore2(ki,kj)=0;
            %         end
            %     end
            % end
            totscore(:,:,i)=totscore2(:,:);
        end
        save(mSCRout, 'totscore');
    end
end
