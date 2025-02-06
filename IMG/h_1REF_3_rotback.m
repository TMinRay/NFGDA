run ../main/NF00_header
[gcol_indices,grow_indices]=meshgrid(1:401);
grow_indices = grow_indices(:);
gcol_indices = gcol_indices(:);

for cindex=1:numel(ttable(:,1));
    
    PUTDAT=[ ttable(cindex,:) ];
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
        moriSCRout=[matPATH '/IMG/Z/oriz' PUTDAT num2str(m,'%02i') '.mat'];
        mSCRout=[matPATH '/IMG/Z/scrz' PUTDAT num2str(m,'%02i') '.mat'];
        load(mSCRout, 'totscore');
        originscore=zeros(401,401,rotnum);
        for i=1:rotnum
            a2=totscore(:,:,i);
            % notnana2=a2>0;
            % for ii=1:401
            %     for jj=1:401
            %         if notnana2(ii,jj)==1
            %             ycord=-100+0.5*(ii-1);
            %             xcord=-100+0.5*(jj-1);
            %             if (xcord>-100 & xcord<100 & ycord<100 & ycord>-100)
            %                     origindeg=1*rotbackrad*(i-1);
            %                     mat11=cos(origindeg);
            %                     mat12=sin(origindeg);
            %                     mat21=-1*mat12;
            %                     mat22=mat11;
            %                     backprocess=[mat11 mat12; mat21 mat22;];
            %                     rotcord=[xcord ycord]';
            %                     oldcord= backprocess*rotcord;
            %                     xnew=oldcord(1);
            %                     ynew=oldcord(2);
            %                     oldi=round((ynew+100)/0.5+1);
            %                     oldj=round((xnew+100)/0.5+1);
            %                 if (oldi>1 & oldi<400 & oldj<400 & oldj>1)
            %                     originscore(oldi,oldj,i)=a2(ii,jj);
            %                 end
            %             end
            %         clear oldj oldi xnew ynew ycord xcord ;  
            %         end
            %     end
            % end

            row_indices = grow_indices(a2(:)>0);
            col_indices = gcol_indices(a2(:)>0);
            ycord=-100+0.5*(row_indices-1);
            xcord=-100+0.5*(col_indices-1);
            origindeg=1*rotbackrad*(i-1);
            mat11 = cos(origindeg);
            mat12 = sin(origindeg);
            mat21 = -mat12;
            mat22 = mat11;
            backprocess = [mat11 mat12; mat21 mat22;];
            rotcord = [xcord ycord]';
            oldcord = backprocess*rotcord;
            xnew=oldcord(1,:);
            ynew=oldcord(2,:);
            oldi=round((ynew.'+100)/0.5+1);
            oldj=round((xnew.'+100)/0.5+1);
            mappxl = (xcord>-100 & xcord<100 & ycord<100 & ycord>-100) ...
                & (oldi>1 & oldi<400 & oldj<400 & oldj>1);
            row_indices = row_indices(mappxl);
            col_indices = col_indices(mappxl);
            oldi = oldi(mappxl);
            oldj = oldj(mappxl);
            lin_indices = sub2ind(size(a2), row_indices, col_indices);
            lin_old_indices = sub2ind(size(a2), oldi, oldj);
            buf = zeros(401,401);
            buf(lin_old_indices)=a2(lin_indices);
            originscore(:,:,i) = buf;
        end
        save(moriSCRout, 'originscore');
        clear originscore;
   end
 end
 



   