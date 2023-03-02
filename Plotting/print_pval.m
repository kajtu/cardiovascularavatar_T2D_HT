function print_pval(t,xpos,ypos,xlinepos,pval,indSign,reject,rejectCorrected, pvalsize,ymax,ytext)

plot(t,xlinepos, [1 1]*ypos*ymax, '-k')

if isempty(rejectCorrected)
    if ~isnan(reject(indSign)) && reject(indSign)
        if pval < 0.001
            text(t,xpos,ymax*(ypos+ytext),'p<0.001*','FontSize',pvalsize,'FontName','Calibri')
        else
            text(t,xpos,ymax*(ypos+ytext),sprintf('p=%0.3f*',pval),'FontSize',pvalsize,'FontName','Calibri')
        end
    elseif ~isnan(pval)
        text(t,xpos,ymax*(ypos+ytext),sprintf('p=%0.3f',pval),'FontSize',pvalsize,'FontName','Calibri') 
    end
else
    if ~isnan(rejectCorrected(indSign)) && rejectCorrected(indSign) %significant also after correction
        if pval < 0.001
            text(t,xpos,ymax*(ypos+ytext),'p<0.001*','FontSize',pvalsize,'FontName','Calibri')
        else
            text(t,xpos,ymax*(ypos+ytext),sprintf('p=%0.3f*',pval),'FontSize',pvalsize,'FontName','Calibri')
        end
    elseif ~isnan(reject(indSign)) && reject(indSign) %significant only before correction
        text(t,xpos,ymax*(ypos+ytext),sprintf('p=%0.3f',pval),'FontSize',pvalsize,'FontName','Calibri')
    elseif ~isnan(pval) %not significant
        text(t,xpos,ymax*(ypos+ytext),sprintf('p=%0.3f',pval),'FontSize',pvalsize,'FontName','Calibri')
    end
end

end