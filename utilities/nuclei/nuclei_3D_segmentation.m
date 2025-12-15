function signal = nuclei_3D_segmentation(img3D, labels, dims, nbObject, maxradius, rad_dilate)

    labels3Dtubes = repmat(labels, 1, 1, dims(3));
    LL = regionprops(labels, 'BoundingBox');
    padding = 15;
    padse = strel('disk', 5, 0);

    if nbObject < 256
        map3D = zeros(dims, 'uint8');
    else
        map3D = zeros(dims, 'uint16');
    end
    
    for ilab = 1:size(LL,1)
        BB = floor(LL(ilab).BoundingBox);
        BBrow = max(1, BB(2)+1-padding):min(BB(2)+BB(4)+padding, dims(1));
        BBcol = max(1, BB(1)+1-padding):min(BB(1)+BB(3)+padding, dims(2));

        output4seg = nan(dims);
        locmap = imdilate(labels3Dtubes == ilab, padse);
        output4seg(locmap) = img3D(locmap);
        sub3D = output4seg(BBrow, BBcol, :);

        class3D = FastFCMeans(im2uint8(sub3D), 3);

        if nbObject < 256
            map3D(BBrow, BBcol, :) = map3D(BBrow, BBcol, :) + uint8(class3D == 3).*ilab;
        else
            map3D(BBrow, BBcol, :) = map3D(BBrow, BBcol, :) + uint16(class3D == 3).*ilab;
        end
    end

    conts = zeros(dims);
    for i_z = 1:dims(3)
        conts(:,:,i_z) = edge(map3D(:,:,i_z), 0);
    end
    conts = conts.*labels3Dtubes;

    fprintf('\t Nuclei located and labeled: %d found !\n', nbObject);

    if nbObject > 255
        signal = zeros(dims, 'uint16');
    else
        signal = zeros(dims, 'uint8');
    end

    for i_nuc = 1:nbObject
        pointMatrix = find(conts == i_nuc);
        if length(pointMatrix) < 4
            continue
        end
        [px, py, pz] = ind2sub(dims, pointMatrix);
        dt = delaunayTriangulation([px, py, pz]);  %# Create a Delaunay triangulation
        if isempty(dt.ConnectivityList)
            continue
        end

        minx = min(px); maxx = max(px);
        miny = min(py); maxy = max(py);
        minz = min(pz); maxz = max(pz);
        [YY,XX,ZZ] = meshgrid(miny:maxy, minx:maxx, minz:maxz);   %# Create a mesh of coordinates for your volume
        simplexIndex = pointLocation(dt,XX(:),YY(:),ZZ(:));       %# Find index of simplex that each point is inside
        
        maskloc = ~isnan(simplexIndex);                           %# Points outside the convex hull have a simplex index of NaN
        maskloc = double(reshape(maskloc,size(XX)));

        locs = signal(minx:maxx, miny:maxy,minz:maxz) == 0;
        maskloc(locs) = maskloc(locs).*i_nuc;
        
        locs_sig = signal(minx:maxx, miny:maxy,minz:maxz);
        maskloc(~locs) = locs_sig(~locs);

        signal(minx:maxx, miny:maxy,minz:maxz) = maskloc;
    end
    
    signal = imclearborder(signal, 8);

    % Remove nuclei masks that are to big
    nucs2kill = repmat(bwareaopen(max(signal, [], 3), ceil((maxradius/2)^2*pi)), 1, 1, dims(3));
    signal(nucs2kill>0) = 0;
    
    fprintf('\t Nuclei fully segmented\n');
    
    % Check for empty image
    valid = zeros(nbObject, 1);
    bg = median(img3D(signal == 0));
    for i = 1:nbObject
        valid(i) = median(img3D(signal == i));
    end
    
    if sum(valid(~isnan(valid)) > bg*1.05) == 0
        fprintf('\tThe image does not contain any nuclei\n');
        signal = zeros(dims);
    end
    
    % Dilate nuclei masks if required by user
    if rad_dilate > 0
        signal = imdilate(signal, strel('sphere', rad_dilate));
    end
    
    %  Relabel if needed
    gnaan = unique(signal(signal>0));
    if length(gnaan) < nbObject
        nbObject = length(gnaan);
        for i = 1:nbObject
            signal(signal == gnaan(i)) = i;
        end
        fprintf('Nuclei re-labeled, correct number of nuclei segmented: %d\n', nbObject);
    end
end