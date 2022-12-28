function [centers] = centersFromGrid(sep,maxGrid,minGrid)

            spacing = sep;
            x1Range = minGrid(1):spacing:maxGrid(1);
            x2Range = minGrid(2):spacing:maxGrid(2);

            [X1,X2] = meshgrid(x1Range,x2Range);
            sz = size(X1);

            centers = zeros(sz(1)*sz(2),2);
            idx = 1;
            for ii = 1:sz(1)
                for jj = 1:sz(2)
                    centers(idx,:) = [X1(ii,jj),X2(ii,jj)];
                    idx = idx+1;
                end
            end

end