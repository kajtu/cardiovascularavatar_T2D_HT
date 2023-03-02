function title = groupTitle(groups,p,patNums,realPatNums,number)
% get a title for which group the subject belongs to based on home blood pressure measurement and T2d or not
pat = num2str(p);
if groups.T2D_HT_home(p)
    if realPatNums
        title = [patNums{p} ' T2D+HT'];
    elseif nargin >4
        title = [number ' T2D+HT'];
    else
        title = [pat ' T2D+HT'];
    end
elseif groups.T2D_NT_home(p)
    if realPatNums
        title = [patNums{p} ' T2D'];
    elseif nargin >4
        title = [number ' T2D'];
    else
        title = [pat ' T2D'];
    end
elseif groups.C_HT_home(p)
    if realPatNums
        title = [patNums{p} ' HT'];
    elseif nargin >4
        title = [number ' HT'];
    else
        title = [pat ' HT'];
    end
else
    if realPatNums
        title = [patNums{p} ' Control'];
    elseif nargin >4
        title = [number ' Control'];
    else
        title = [pat ' Control'];
    end
end

end