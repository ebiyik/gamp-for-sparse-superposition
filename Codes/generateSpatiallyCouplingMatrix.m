function J = generateSpatiallyCouplingMatrix(wLeft, wRight, Lr, Lc, Jright)
    % Generates the variance matrix for spatially coupling

    %% initialization
    J = zeros(Lr, Lc);
    
    %% first row
    J(1,1) = 1;
    J(1,2:wRight+1) = Jright;
    nextRowNo = 2;
    
    %% next rows till the usual rows where there are (wLeft + 1) blocks with variance 1 and wRight block with variance Jright
    for rowNo=1:wLeft-1
        J(nextRowNo,1:rowNo+1) = 1;
        J(nextRowNo,rowNo+2:rowNo+1+wRight) = Jright;
        nextRowNo = nextRowNo + 1;
    end
    
    %% the usual rows where there are (wLeft + 1) blocks with variance 1 and wRight block with variance Jright
    for rowNo = 1:Lc-wLeft-wRight+1
        J(nextRowNo,rowNo:rowNo+wLeft) = 1;
        J(nextRowNo,rowNo+wLeft+1:rowNo+wLeft+wRight) = Jright;
        nextRowNo = nextRowNo + 1;
    end
    
    %% last rows
    for rowNo = 1:wRight
        J(nextRowNo,Lr-wLeft-wRight+rowNo:Lr-wRight+rowNo) = 1;
        J(nextRowNo,Lr-wRight+rowNo+1:end) = Jright;
        nextRowNo = nextRowNo + 1;
    end
end