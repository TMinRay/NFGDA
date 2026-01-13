function [buf] = rot_score_back(a2,origindeg,grow_indices,gcol_indices)
    mat11 = cos(origindeg);
    mat12 = sin(origindeg);
    mat21 = -mat12;
    mat22 = mat11;
    backprocess = [mat11 mat12; mat21 mat22;];

    row_indices = grow_indices(a2(:)>0);
    col_indices = gcol_indices(a2(:)>0);
    ycord=-100+0.5*(row_indices-1);
    xcord=-100+0.5*(col_indices-1);
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
end