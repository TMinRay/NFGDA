NF00_header;

totscore=zeros(1,7);
% totscore(threshi,1)=threshH(threshi);
 lscrfile2=[matPATH '/EVAL/neweval_exmigfa.mat'];
load(lscrfile2,'lscr')   
% % % % 
% % % % cindex m heat miss false dtdL fdtdL totL
% % % %    1      2     3    4      5       6     7    8

[totnum varinum]=size(lscr);

% % % % exclude time steps used in making training datasets!!
% % % % exclude time steps used in making training datasets!!
% % % % in header.m those will be defined as
% % % % trncasest=1;
% % % % trncasend=2;
% % % % trnstartt=[13 7];
% % % % trnendt=[13 7];
% % % % % % means making training dataset based on case #1 timestep 13
% % % % % % and case #2 timestep 2
% % % % % %  the excluding should be
% % % for i=1:totnum
% % %     if lscr(totnum,1)==13 || lscr(totnum,2)==7 --> check this part
% % %         lscr(totnum,3:9)=0;   --> no need to change --> excluding
% % % scores for timestep used in making training dataset
% % %     end
% % % end
% % % % % % % % % % % % % % % % % % % 
% % % % % lscr(totnum,1) --> 1st case
% % % % % lscr(totnum,1) == 13 --> 13th timestep
% % % % % 
% % % % % lscr(totnum,2) --> 2nd case
% % % % % lscr(totnum,2) == 7 --> 7th timestep



% for i=1:totnum
%     if lscr(totnum,1)==1 || lscr(totnum,2)==13
%         lscr(totnum,3:9)=0;        
%     end
% end
for i=1:totnum
    if lscr(i,1)==1 && lscr(i,2)==2
        lscr(i,3:9)=0;        
    end
end

hit=lscr(:,3);
miss=lscr(:,4);
fa=lscr(:,5);
cn=lscr(:,6);
dtdL=lscr(:,7);
falseL=lscr(:,8);
totL=lscr(:,9);

shit=sum(hit(:));
smiss=sum(miss(:));
sfa=sum(fa(:));
scn=sum(cn(:));
sdtdL=sum(dtdL(:));
sfalseL=sum(falseL(:));
stotL=sum(totL(:));

S=sprintf('%s', ['MIGFA POD:' num2str((shit/(shit+smiss))*100) ...
',PFD:' num2str((sfa/(shit+sfa))*100) ...
',PC:' num2str(((shit+scn)/(shit+smiss+sfa+scn))*100) ...
',PLD:' num2str((sdtdL/stotL)*100) ...
',PFD:' num2str((sfalseL/(sfalseL+sdtdL))*100)]);
disp(S)
% % 
