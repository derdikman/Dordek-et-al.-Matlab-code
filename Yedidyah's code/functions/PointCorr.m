function out_mat=PointCorr(Rxx,RxxR)
nan_mat=Rxx .* RxxR;                                             
        notnan_inds = find(~isnan(nan_mat));  %normalized to the number of nontnan components (average)
        n=length(notnan_inds);
        
        if n < 2
            out_mat = NaN;
            
        end
  
        sigma_x_y =sum(nan_mat(notnan_inds));
        sigma_x =      sum(Rxx(notnan_inds));
        sigma_y =      sum(RxxR(notnan_inds));
        sigma_x2 = sum(Rxx(notnan_inds).^2);
        sigma_y2 = sum(RxxR(notnan_inds).^2);
        
        out_mat = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
                                    sqrt(n*sigma_x2-sigma_x.^2) ./ ...
                                    sqrt(n*sigma_y2-sigma_y.^2);
                   
end