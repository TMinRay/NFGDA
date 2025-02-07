run ../main/NF00_header
[gcol_indices,grow_indices]=meshgrid(1:401);
grow_indices = grow_indices(:);
gcol_indices = gcol_indices(:);

for cindex=1:numel(ttable(:,1));
    
     PUTDAT=[ ttable(cindex,:) ];
     startm=startt(cindex)+1;
     endm=endt(cindex);
  
    for m=startm:endm
    
        mSCRout=[matPATH '/IMG/DELZ/scrdelz' PUTDAT num2str(m,'%02i') '.mat'];
        load(mSCRout, 'totscore');
        moriSCRout=[matPATH '/IMG/DELZ/oridelz' PUTDAT num2str(m,'%02i') '.mat'];
        originscore=zeros(401,401,rotnum);
        for i=1:rotnum
            a2=totscore(:,:,i);
            % notnana2=a2>0;
            origindeg=1*rotbackrad*(i-1);
            originscore(:,:,i) = rot_score_back(a2,origindeg,grow_indices,gcol_indices);
            % for ii=1:401
            %     for jj=1:401
            %         if notnana2(ii,jj)==1
            %         ycord=-100+0.5*(ii-1);
            %         xcord=-100+0.5*(jj-1);
            %         if (xcord>-100 & xcord<100 & ycord<100 & ycord>-100)
            %         origindeg=1*rotbackrad*(i-1);
            %         mat11=cos(origindeg);
            %         mat12=sin(origindeg);
            %         mat21=-1*mat12;
            %         mat22=mat11;
            %         backprocess=[mat11 mat12; mat21 mat22;];
            %         rotcord=[xcord ycord]';
            %         oldcord= backprocess*rotcord;
            %         xnew=oldcord(1);
            %         ynew=oldcord(2);
            %         oldi=round((ynew+100)/0.5+1);
            %         oldj=round((xnew+100)/0.5+1);
            %         if (oldi>1 & oldi<400 & oldj<400 & oldj>1)
            %         originscore(oldi,oldj,i)=a2(ii,jj);
            %         end
            %         end
            %         clear oldj oldi xnew ynew ycord xcord ;
            %         end 
            %     end
            % end

        end
        save(moriSCRout, 'originscore');
        clear originscore
    end
 end
 



   