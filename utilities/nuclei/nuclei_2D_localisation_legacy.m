function labels = nuclei_2D_localisation_legacy(map2D, dims, spotsize)
%% nuclei_2D_localisation_legacy(map2D, dims, spotsize)
%
% Locate and segment nuclei in 2D from map2D image using gLog filters
% Remove objects on edges and small chunks of pixels
% Relabel everything after cleaning
%
% Inputs:
%   map2D - image to segment
%   dims - dimensions of the image (in pixels)
%   spotsize - largest possible width of a spot, in pixels

        imfilt = zeros(dims(1:2));

        sigmax = resXY/0.064*spotsize:spotsize:resXY/0.064*spotsize^2;
        sigmay = resXY/0.064*spotsize:spotsize:resXY/0.064*spotsize^2;
        theta = 0:pi/4:pi/4;
        
        reps=min(10, length(sigmax));
        fprintf(['\t[',repmat(' ',1,reps),']']);
        did = 0;

        for i_sigx = 1:length(sigmax)
            for i_sigy = 1:length(sigmay)
               for i_theta = 1:length(theta)
                % Generate 3D gaussian kernel
                    cellsize = max(max(sigmax), max(sigmay));

                    wxyz = ceil(cellsize*7);
                    cxyz = wxyz/2+0.5;
                    [xx, yy] = meshgrid(1:wxyz);

                    a = cos(theta(i_theta))^2/(2*sigmax(i_sigx)^2) + sin(theta(i_theta))^2/(2*sigmay(i_sigy)^2);
                    b = -sin(2*theta(i_theta))/(4*sigmax(i_sigx)^2) + sin(2*theta(i_theta))/(4*sigmay(i_sigy)^2);
                    c = sin(theta(i_theta))^2/(2*sigmax(i_sigx)^2) + cos(theta(i_theta))^2/(2*sigmay(i_sigy)^2);

                    filter = exp(-(a.*(xx - cxyz).^2 + 2*b.*(xx - cxyz).*(yy - cxyz) + c.*(yy - cxyz).^2));
        %             filter = exp(-( ((xx-cxyz).^2)/(2*sigmax(i_sigx)^2) + ((yy-cxyz).^2)/(2*sigmay(i_sigy)^2)));
                    filter = filter./sum(filter(:));
                    filter = del2(filter);


                    imfilt = imfilt + imfilter(double(map2D), filter, 'symmetric').*(1+log(sigmax(i_sigx))).*(1+log(sigmay(i_sigy)));
                end
            end

            if sum(i_sigx == round(1:(length(sigmax)/reps):length(sigmax))) == 1
                fprintf(repmat('\b',1,reps+1-did)) %send the cursor back to the start
                did = did +1;
                fprintf('-');
                fprintf(repmat(' ',1,reps-did));
                fprintf(']');
            end
        end

        seg = 1-mat2gray(imfilt);
        seg = imtophat(seg, strel('disk', sigmax(end), 0));
        
        locs = imregionalmax(seg); %medfilt2(seg, [cellsize cellsize]));
        if Sig_qual
            segseg = FastFCMeans(im2uint8(seg), 3, 1.1);
            locs(segseg < 2) = 0;
            segseg = segseg > 1;
        else
            segseg = FastFCMeans(im2uint8(seg), 3, 1.1);
            locs(segseg < 3) = 0;
            segseg = segseg > 1;
        end
        
        % First separation step
        ws = watershed(-imgaussfilt(seg, spotsizedef));
        segseg(ws == 0) = 0;
        
        % Check if noisy segmentation (big connected objects)
        CC = bwconncomp(segseg);
        LL = regionprops(CC, 'Area', 'Solidity');
        if corr(cat(1,LL.Area), cat(1,LL.Solidity)) < -0.5
            % Possibly bad segmentation - try to separate
            fm = FastFCMeans(im2uint8(seg), 3, 1.1);
            ws = watershed(3-fm);
            segseg(ws == 0) = false;
            
            CC = bwconncomp(segseg);
            LL = regionprops(CC, 'Area', 'Solidity');
            if corr(cat(1,LL.Area), cat(1,LL.Solidity)) < -0.5
                % Gnaan
                locs = imregionalmax(seg);
                segseg = FastFCMeans(im2uint8(seg), 3, 1.1);
                locs(segseg < 3) = 0;
                segseg = segseg == 3;
            end
        end
        clear CC LL
        
        % Check if noisy segmentation (background kept)
        if sum(segseg(:))/double(dims(1))/double(dims(2)) > .8 % 80% of image (in 2D) segmented
            locs = imregionalmax(seg);
            segseg = FastFCMeans(im2uint8(seg), 3, 1.1);
            locs(segseg < 3) = 0;
            segseg = segseg == 3;
        end

        [px, py] = find(locs);
        candidates = [px, py];
        T = clusterdata(candidates, 'criterion', 'distance', 'cutoff', spotsize);
        nbLocs = max(T);
        newlocs = zeros(nbLocs, 2);
        for i = 1:nbLocs
            newlocs(i,:) = round(mean(candidates(T == i, :), 1));
        end

        locs = zeros(dims(1), dims(2));
        locs(sub2ind([dims(1), dims(2)], newlocs(:,1), newlocs(:,2))) = 1;
        
        ridges = watershed(~segseg);
        map = double(ridges>0);
        map(locs(:)>0) = 2;
        
        labels = watershed(-map);
end
